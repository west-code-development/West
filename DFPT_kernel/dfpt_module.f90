!
! Copyright (C) 2015-2021 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Marco Govoni
!
!-----------------------------------------------------------------------
MODULE dfpt_module
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-------------------------------------------------------------------------
    SUBROUTINE dfpt (m,dvg,dng,tr2,iq)
      !-----------------------------------------------------------------------
      !
      USE kinds,                 ONLY : DP
      USE io_global,             ONLY : stdout
      USE wvfct,                 ONLY : nbnd,et
      USE fft_base,              ONLY : dffts
      USE gvect,                 ONLY : gstart
      USE wavefunctions,         ONLY : evc,psic
      USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
      USE mp_global,             ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm
      USE mp_world,              ONLY : world_comm
      USE buffers,               ONLY : get_buffer
      USE noncollin_module,      ONLY : noncolin,npol
      USE pwcom,                 ONLY : current_spin,isk,npw,npwx,lsda,current_k,ngk,igk_k
      USE cell_base,             ONLY : omega
      USE control_flags,         ONLY : gamma_only
      USE uspp,                  ONLY : nkb,vkb
      USE uspp_init,             ONLY : init_us_2
      USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
      USE fft_at_gamma,          ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
      USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
      USE io_push,               ONLY : io_push_title
      USE types_bz_grid,         ONLY : k_grid,q_grid,compute_phase
      USE westcom,               ONLY : iuwfc,lrwfc,npwqx,npwq,igq_q,fftdriver,&
                                      & l_frac_occ,nbnd_occ,nbnd_occ_full,occupation,docc_thr
      USE distribution_center,   ONLY : band_group,kpt_pool
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN), OPTIONAL :: iq
      INTEGER, INTENT(IN) :: m
      COMPLEX(DP), INTENT(IN) :: dvg(npwqx,m)
      COMPLEX(DP), INTENT(OUT) :: dng(npwqx,m)
      REAL(DP),INTENT(IN) :: tr2
      !
      ! Workspace
      !
      INTEGER :: ipert, ig, ir, ibnd, jbnd, ibnd2, lbnd, iks, ikqs, ikq, ik, is, iks_g
      INTEGER :: nbndval, nbndval_full, ierr
      INTEGER :: npwkq
      !
      REAL(DP) :: g0(3)
      REAL(DP) :: docc
      REAL(DP), ALLOCATABLE :: eprec(:)
      REAL(DP), ALLOCATABLE :: eprec_loc(:)
      REAL(DP), ALLOCATABLE :: et_loc(:)
      !
      COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:),dpsi(:,:)
      COMPLEX(DP), ALLOCATABLE :: aux_r(:),aux_g(:),dpsi_frac(:),dvr(:)
      COMPLEX(DP), ALLOCATABLE :: dpsic(:)
      !
      COMPLEX(DP), ALLOCATABLE :: evckmq(:,:)
      COMPLEX(DP), ALLOCATABLE :: phase(:)
      !
      TYPE(bar_type) :: barra
      !
      LOGICAL :: l_dost
      !
      CHARACTER(LEN=512) :: title
      !
      CALL mp_barrier( world_comm )
      !
      CALL report_dynamical_memory( )
      !
      l_dost = ( tr2 >= 0._DP )
      !
      IF( l_dost ) THEN
         WRITE(title,'(a,es14.6)') "Sternheimer eq. solver... with threshold = ", tr2
      ELSE
         WRITE(title,'(a,es14.6)') "Sternheimer eq. solver... with lite-solver"
      ENDIF
      CALL io_push_title(TRIM(ADJUSTL(title)))
      !
      dng=0.0_DP
      !
      CALL start_bar_type( barra, 'dfpt', MAX(m,1) * kpt_pool%nloc )
      !
      IF (.NOT.gamma_only) THEN
         ALLOCATE( evckmq(npwx*npol,nbnd) )
         ALLOCATE( phase(dffts%nnr) )
      ENDIF
      !
      IF (l_frac_occ .and. .not. gamma_only) THEN
         CALL errore("dfpt", "fraction occupation only implemented for gamma-only case", 1)
      ENDIF
      !
      DO iks = 1, kpt_pool%nloc ! KPOINT-SPIN LOOP
         !
         iks_g = kpt_pool%l2g(iks)
         ik = k_grid%ip(iks_g)
         is = k_grid%is(iks_g)
         !
         ! ... Set k-point, spin, kinetic energy, needed by Hpsi
         !
         current_k = iks
         IF ( lsda ) current_spin = isk(iks)
         CALL g2_kin(iks)
         !
         ! ... More stuff needed by the hamiltonian: nonlocal projectors
         !
         IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), k_grid%p_cart(1,ik), vkb )
         !
         nbndval = nbnd_occ(iks)
         IF (l_frac_occ) nbndval_full = nbnd_occ_full(iks)
         !
         CALL band_group%init( nbndval, 'b', 'band_group', .FALSE. )
         !
         ! ... Number of G vectors for PW expansion of wfs at k
         !
         npw = ngk(iks)
         !
         ! ... Read wavefuctions at k in G space, for all bands, and store them in evc
         !
         IF(kpt_pool%nloc > 1) THEN
            IF ( my_image_id == 0 ) CALL get_buffer( evc, lrwfc, iuwfc, iks )
            CALL mp_bcast( evc, 0, inter_image_comm )
         ENDIF
         !
         IF (gamma_only) THEN
            !
            ikqs = iks
            g0 = 0.0_DP
            !
         ELSE
            !
            ! ... Find G0 and compute phase
            !
            !CALL k_grid%find( k_grid%p_cart(:,ik) - q_grid%p_cart(:,iq), is, 'cart', ikqs, g0 )   !M
            CALL k_grid%find( k_grid%p_cart(:,ik) - q_grid%p_cart(:,iq), 'cart', ikq, g0 )
            ikqs = k_grid%ipis2ips(ikq,is)
            !
            CALL compute_phase( g0, 'cart', phase )
            !
            ! ... Number of G vectors for PW expansion of wfs at [k-q]
            !
            npwkq = ngk(ikqs)
            !
            ! ... Set wavefunctions at [k-q] in G space, for all bands, and store them in evckmq
            !
            IF ( my_image_id == 0 ) CALL get_buffer( evckmq, lrwfc, iuwfc, ikqs )
            CALL mp_bcast( evckmq, 0, inter_image_comm )
            !
         ENDIF
         !
         !
         ALLOCATE(eprec(nbndval))
         ALLOCATE(eprec_loc(band_group%nloc))
         ALLOCATE(et_loc(band_group%nloc))
         CALL set_eprec(nbndval,evc(1,1),eprec)
         !
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)
            eprec_loc(lbnd) = eprec(ibnd)
            et_loc(lbnd) = et(ibnd,ikqs)
         ENDDO
         !
         ALLOCATE(dvpsi(npwx*npol,band_group%nloc))
         ALLOCATE(dpsi(npwx*npol,band_group%nloc))
         !
         DO ipert = 1, m
            !
            ALLOCATE( aux_g(npwqx) ); aux_g = 0._DP
            ALLOCATE( aux_r(dffts%nnr) ); aux_r = 0._DP
            !
            DO CONCURRENT (ig = 1:npwq)
               aux_g(ig) = dvg(ig,ipert)
            ENDDO
            !
            ! ... inverse Fourier transform of the perturbation: (q+)G ---> R
            !
            IF (gamma_only) THEN
               CALL single_invfft_gamma(dffts,npwq,npwqx,aux_g,aux_r,TRIM(fftdriver))
            ELSE
               CALL single_invfft_k(dffts,npwq,npwqx,aux_g,aux_r,'Wave',igq_q(1,iq))
            ENDIF
            !
            ! The perturbation is in aux_r
            !
            dvpsi = 0._DP
            dpsi  = 0._DP
            !
            IF(gamma_only) THEN
               !
               ! double bands @ gamma
               DO lbnd = 1,band_group%nloc-MOD(band_group%nloc,2),2
                  !
                  ibnd = band_group%l2g(lbnd)
                  ibnd2 = band_group%l2g(lbnd+1)
                  !
                  CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,ibnd2),psic,'Wave')
                  DO CONCURRENT (ir=1:dffts%nnr)
                     psic(ir) = psic(ir) * REAL(aux_r(ir),KIND=DP)
                  ENDDO
                  CALL double_fwfft_gamma(dffts,npw,npwx,psic,dvpsi(1,lbnd),dvpsi(1,lbnd+1),'Wave')
                  !
               ENDDO
               !
               ! single band @ gamma
               IF( MOD(band_group%nloc,2) == 1 ) THEN
                  !
                  lbnd = band_group%nloc
                  ibnd = band_group%l2g(lbnd)
                  !
                  CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),psic,'Wave')
                  DO CONCURRENT (ir=1:dffts%nnr)
                     psic(ir) = CMPLX( REAL(psic(ir),KIND=DP) * REAL(aux_r(ir),KIND=DP), 0._DP, KIND=DP)
                  ENDDO
                  CALL single_fwfft_gamma(dffts,npw,npwx,psic,dvpsi(1,lbnd),'Wave')
                  !
               ENDIF
               !
            ELSE
               !
               DO lbnd = 1,band_group%nloc
                  !
                  ibnd = band_group%l2g(lbnd)
                  !
                  ! ... inverse Fourier transform of wfs at [k-q]: (k-q+)G ---> R
                  !
                  CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1,ibnd),psic,'Wave',igk_k(1,ikqs))
                  !
                  ! ... construct right-hand-side term of Sternheimer equation:
                  ! ... product of wavefunction at [k-q], phase and perturbation in real space
                  !
                  DO CONCURRENT (ir = 1:dffts%nnr)
                     psic(ir) = psic(ir) * phase(ir) * aux_r(ir)
                  ENDDO
                  !
                  ! Fourier transform product of wf at [k-q], phase and
                  ! perturbation of wavevector q: R ---> (k+)G
                  !
                  CALL single_fwfft_k(dffts,npw,npwx,psic,dvpsi(1,lbnd),'Wave',igk_k(1,iks))
                  !
                  ! dv|psi> is in dvpsi
                  !
               ENDDO
               !
               IF (noncolin) THEN
                  DO lbnd = 1,band_group%nloc
                     !
                     ibnd = band_group%l2g(lbnd)
                     !
                     CALL single_invfft_k(dffts,npwkq,npwx,evckmq(npwx+1,ibnd),psic,'Wave',igk_k(1,ikqs))
                     !
                     DO CONCURRENT (ir = 1:dffts%nnr)
                        psic(ir) = psic(ir) * phase(ir) * aux_r(ir)
                     ENDDO
                     !
                     CALL single_fwfft_k(dffts,npw,npwx,psic,dvpsi(npwx+1,lbnd),'Wave',igk_k(1,iks))
                     !
                  ENDDO
               ENDIF
               !
            ENDIF
            !
            ! keep dV in R space for the calculation of an extra term in frac occ case
            IF ( l_frac_occ ) THEN
               ALLOCATE( dvr(dffts%nnr) )
               dvr = aux_r
            ENDIF
            !
            DEALLOCATE( aux_g )
            DEALLOCATE( aux_r )
            !
            ! - P_c | dvpsi >
            !
            CALL apply_alpha_pc_to_m_wfcs( nbndval, band_group%nloc, dvpsi, (-1._DP,0._DP) )
            !
            CALL precondition_m_wfcts( band_group%nloc, dvpsi, dpsi, eprec_loc )
            !
            IF( l_dost) THEN
               !
               ! The Sternheimer operator is (H_k - E_(k-q) + alpha * P_v)
               ! The Hamiltonian is evaluated at the k-point current_k in h_psi
               ! (see also PHonon/PH/cch_psi_all.f90, where H_(k+q) is evaluated)
               !
               CALL linsolve_sternheimer_m_wfcts (nbndval, band_group%nloc, dvpsi, dpsi, et_loc, eprec_loc, tr2, ierr )
               !
               IF(ierr/=0) THEN
                  WRITE(stdout, '(7X,"** WARNING : PERT ",i8," iks ",I8," not converged, ierr = ",i8)') ipert,iks,ierr
               ENDIF
               !
            ENDIF
            !
            ALLOCATE( aux_r(dffts%nnr) )
            !
            aux_r=0._DP
            !
            IF(gamma_only) THEN
               !
               ! double band @ gamma
               DO lbnd = 1,band_group%nloc
                  !
                  ibnd = band_group%l2g(lbnd)
                  !
                  CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),dpsi(1,lbnd),psic,'Wave')
                  !
                  IF (l_frac_occ) THEN
                     !
                     ALLOCATE( dpsi_frac(dffts%nnr) )
                     !
                     DO jbnd = MAX(ibnd, nbndval_full), nbndval
                        !
                        docc = occupation(ibnd, iks) - occupation(jbnd, iks)
                        IF ( ABS(docc) < docc_thr ) CYCLE
                        !
                        CALL compute_pt1_dpsi(ibnd, jbnd, iks, dvr, dpsi_frac)
                        !
                        DO ir=1,dffts%nnr
                           psic(ir) = psic(ir) + (docc / occupation(ibnd,iks)) * CMPLX(0._DP, REAL(dpsi_frac(ir),KIND=DP), KIND=DP)
                        ENDDO
                        !
                     ENDDO
                     !
                     DEALLOCATE( dpsi_frac )
                     !
                  ENDIF
                  !
                  DO CONCURRENT (ir=1:dffts%nnr)
                     aux_r(ir) = aux_r(ir) + CMPLX( occupation(ibnd,iks) * REAL( psic(ir),KIND=DP) * AIMAG( psic(ir)) , 0.0_DP, KIND=DP)
                  ENDDO
                  !
               ENDDO
               !
            ELSE
               !
               ALLOCATE( dpsic(dffts%nnr) )
               !
               DO lbnd = 1,band_group%nloc
                  !
                  ibnd = band_group%l2g(lbnd)
                  !
                  ! inverse Fourier transform of wavefunction at [k-q]: (k-q+)G ---> R
                  !
                  CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1,ibnd),psic,'Wave',igk_k(1,ikqs))
                  !
                  ! inverse Fourier transform of perturbed wavefunction: (k+)G ---> R
                  !
                  CALL single_invfft_k(dffts,npw,npwx,dpsi(1,lbnd),dpsic,'Wave',igk_k(1,iks))
                  !
                  DO CONCURRENT (ir = 1: dffts%nnr)
                     aux_r(ir) = aux_r(ir) + CONJG( psic(ir) * phase(ir) ) * dpsic(ir)
                  ENDDO
                  !
               ENDDO
               !
               IF (noncolin) THEN
                  DO lbnd = 1,band_group%nloc
                     !
                     ibnd = band_group%l2g(lbnd)
                     !
                     CALL single_invfft_k(dffts,npwkq,npwx,evckmq(npwx+1,ibnd),psic,'Wave',igk_k(1,ikqs))
                     !
                     CALL single_invfft_k(dffts,npw,npwx,dpsi(npwx+1,lbnd),dpsic,'Wave',igk_k(1,iks))
                     !
                     DO CONCURRENT (ir = 1: dffts%nnr)
                        aux_r(ir) = aux_r(ir) + CONJG( psic(ir) * phase(ir) ) * dpsic(ir)
                     ENDDO
                     !
                  ENDDO
               ENDIF
               !
               DEALLOCATE( dpsic )
               !
            ENDIF
            !
            ! Sum up aux_r from band groups
            !
            CALL mp_sum(aux_r,inter_bgrp_comm)
            !
            ! The perturbation is in aux_r
            !
            ALLOCATE( aux_g(npwqx) )
            !
            IF(gamma_only) THEN
               CALL single_fwfft_gamma(dffts,npwq,npwqx,aux_r,aux_g,TRIM(fftdriver))
            ELSE
               CALL single_fwfft_k(dffts,npwq,npwqx,aux_r,aux_g,'Wave',igq_q(1,iq))
            ENDIF
            !
            DO CONCURRENT( ig = 1: npwq) ! pert acts only on body
               dng(ig,ipert) = dng(ig,ipert) + 2._DP * k_grid%weight(iks_g) * aux_g(ig) / omega
            ENDDO
            !
            DEALLOCATE( aux_g )
            DEALLOCATE( aux_r )
            !
            IF ( l_frac_occ ) DEALLOCATE( dvr )
            !
            CALL update_bar_type( barra, 'dfpt', 1 )
            !
         ENDDO ! ipert
         !
         IF( m == 0 ) CALL update_bar_type( barra, 'dfpt', 1 )
         !
         DEALLOCATE( eprec )
         DEALLOCATE( eprec_loc )
         DEALLOCATE( et_loc )
         DEALLOCATE( dpsi )
         DEALLOCATE( dvpsi )
         !
      ENDDO ! K-POINT and SPIN
      !
      IF ( gamma_only ) THEN
         IF ( gstart == 2 ) dng(1,1:m) = CMPLX( 0._DP, 0._DP, KIND=DP )
      ELSE
         IF ( gstart == 2 .AND. q_grid%l_pIsGamma(iq) ) dng(1,1:m) = CMPLX( 0._DP, 0._DP, KIND=DP )
         DEALLOCATE( evckmq )
         DEALLOCATE( phase )
      ENDIF
      !
      CALL mp_sum(dng,inter_pool_comm)
      !
      CALL mp_barrier( world_comm )
      !
      CALL stop_bar_type( barra, 'dfpt' )
      !
    END SUBROUTINE
    !
    !
    SUBROUTINE compute_pt1_dpsi(ibnd, jbnd, iks, dvr, dpsi)
      !
      ! Given perturbation potential dvr (in R space), compute 1st order
      ! perturbation to ith KS orbital projected into jth orbital
      ! |dpsi_i> = <psi_j| dV | psi_i> / (ei - ej) |dpsi_j>
      !
      ! If ibnd and jbnd are degenerate, raise an error
      !
      ! The total perturbed KS orbital should contain a summation over j 
      ! with appropriate prefactors
      !
      USE kinds,                 ONLY : DP
      USE io_global,             ONLY : stdout
      USE wvfct,                 ONLY : et
      USE fft_base,              ONLY : dffts
      USE wavefunctions,  ONLY : evc,psic
      USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
      USE mp_global,             ONLY : inter_image_comm,inter_pool_comm,my_image_id,intra_bgrp_comm
      USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
      USE fft_at_gamma,          ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
      USE fft_interfaces,        ONLY : fwfft, invfft
      USE control_flags,         ONLY : gamma_only
      USE westcom,               ONLY : l_frac_occ,nbnd_occ,nbnd_occ_full,occupation,de_thr
      USE pwcom,                 ONLY : npw,npwx
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN)  :: ibnd, jbnd, iks
      COMPLEX(DP),INTENT(IN)   :: dvr(dffts%nnr)
      COMPLEX(DP),INTENT(OUT)  :: dpsi(dffts%nnr)
      !
      INTEGER  :: ir
      REAL(DP) :: vji, de
      COMPLEX(DP), ALLOCATABLE :: aux_r(:)
      !
      IF ( .not. gamma_only .or. .not. l_frac_occ ) THEN
         CALL errore("compute_pt1_dpsi", "this routine should only be called in gamma_only and frac occ case", 1)
      ENDIF
      !
      dpsi = 0._DP
      !
      de = et(ibnd,iks) - et(jbnd,iks)
      !
      IF ( ABS(de) < de_thr ) CALL errore("compute_pt1_dpsi", "degenerate orbitals", 1)
      !
      ALLOCATE( aux_r(dffts%nnr) )
      aux_r = 0._DP
      CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,jbnd),aux_r,'Wave')
      !
      vji = 0._DP
      DO ir=1,dffts%nnr
         vji = vji + DIMAG(aux_r(ir)) * REAL(dvr(ir),KIND=DP) * REAL(aux_r(ir),KIND=DP)
      ENDDO
      vji = vji / (dffts%nr1 * dffts%nr2 * dffts%nr3)
      CALL mp_sum( vji, intra_bgrp_comm )
      !
      !PRINT*, ibnd, jbnd, de, vji
      !
      DO ir=1,dffts%nnr
         dpsi(ir) = (vji / de) * CMPLX(DIMAG(aux_r(ir)), 0._DP, KIND=DP)
      ENDDO
      !
      DEALLOCATE( aux_r )
      !
    END SUBROUTINE
    !
END MODULE
