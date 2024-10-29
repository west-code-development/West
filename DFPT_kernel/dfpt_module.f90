!
! Copyright (C) 2015-2024 M. Govoni
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
      USE fft_base,              ONLY : dffts
      USE gvect,                 ONLY : gstart
      USE mp,                    ONLY : mp_sum,mp_bcast
      USE mp_global,             ONLY : inter_image_comm,my_image_id,inter_pool_comm,nbgrp,my_bgrp_id,&
                                      & inter_bgrp_comm,intra_bgrp_comm
      USE buffers,               ONLY : get_buffer
      USE noncollin_module,      ONLY : noncolin,npol
      USE pwcom,                 ONLY : current_spin,isk,npw,npwx,lsda,current_k,ngk,igk_k,nbnd
      USE cell_base,             ONLY : omega
      USE control_flags,         ONLY : gamma_only
      USE uspp,                  ONLY : nkb,vkb
      USE uspp_init,             ONLY : init_us_2
      USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
      USE fft_at_gamma,          ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                      & double_invfft_gamma
      USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
      USE io_push,               ONLY : io_push_title
      USE types_bz_grid,         ONLY : k_grid,q_grid,compute_phase
      USE westcom,               ONLY : iuwfc,lrwfc,npwqx,npwq,igq_q,fftdriver,l_frac_occ,nbnd_occ,&
                                      & nbnd_occ_full,occupation,docc_thr,de_thr
      USE distribution_center,   ONLY : kpt_pool,band_group
      USE wavefunctions,         ONLY : evc,psic
      USE wvfct,                 ONLY : et
#if defined(__CUDA)
      USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu,allocate_linsolve_gpu,&
                                      & deallocate_linsolve_gpu,reallocate_ps_gpu
      USE cublas
#endif
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
      INTEGER :: ipert, ig, ir, ibnd, jbnd, lbnd, iks, ikqs, ikq, ik, is, iks_g
      INTEGER :: nbndval, nbndval_full, nbndval_frac, ierr
      INTEGER :: npwkq
      INTEGER :: dffts_nnr,band_group_nloc
      !
      REAL(DP) :: g0(3)
      REAL(DP) :: docc, de, this_occ
      REAL(DP), ALLOCATABLE :: eprec(:)
      REAL(DP), ALLOCATABLE :: eprec_loc(:)
      REAL(DP), ALLOCATABLE :: et_loc(:)
      REAL(DP), ALLOCATABLE :: psi_dvpsi(:,:)
      !
      COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:),dpsi(:,:)
      COMPLEX(DP), ALLOCATABLE :: aux_r(:),aux_g(:)
      COMPLEX(DP), ALLOCATABLE :: dpsic(:)
      !
      COMPLEX(DP), ALLOCATABLE :: evckmq(:,:)
#if defined(__CUDA)
      ATTRIBUTES(PINNED) :: evckmq
#endif
      COMPLEX(DP), ALLOCATABLE :: phase(:)
      !
      TYPE(bar_type) :: barra
      !
      LOGICAL :: l_dost
      !
      CHARACTER(LEN=512) :: title
      !
      COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
      !
      IF (l_frac_occ .AND. .NOT. gamma_only) THEN
         CALL errore('dfpt', 'fraction occupation only implemented for gamma-only case', 1)
      ENDIF
      !
      ! Allocation
      !
      nbndval = MAXVAL(nbnd_occ)
      dffts_nnr = dffts%nnr
      !
      CALL band_group%init(nbndval,'b','band_group',.FALSE.)
      !
      ALLOCATE(eprec(nbndval))
      ALLOCATE(eprec_loc(band_group%nloc))
      ALLOCATE(et_loc(band_group%nloc))
      ALLOCATE(dvpsi(npwx*npol,band_group%nloc))
      ALLOCATE(dpsi(npwx*npol,band_group%nloc))
      ALLOCATE(aux_r(dffts%nnr))
      ALLOCATE(aux_g(npwqx))
      !$acc enter data create(eprec,eprec_loc,et_loc,dvpsi,dpsi,aux_r,aux_g)
      !
      IF (.NOT. gamma_only) THEN
         ALLOCATE(dpsic(dffts%nnr))
         ALLOCATE(phase(dffts%nnr))
         ALLOCATE(evckmq(npwx*npol,nbnd))
         !$acc enter data create(dpsic,phase,evckmq)
      ENDIF
      !
#if defined(__CUDA)
      CALL allocate_gpu()
      CALL allocate_linsolve_gpu(band_group%nloc)
#endif
      !
      l_dost = ( tr2 >= 0._DP )
      !
      IF( l_dost ) THEN
         WRITE(title,'(A,ES14.6)') 'Sternheimer eq. solver... with threshold = ', tr2
      ELSE
         WRITE(title,'(A,ES14.6)') 'Sternheimer eq. solver... with lite-solver'
      ENDIF
      CALL io_push_title(TRIM(ADJUSTL(title)))
      !
      dng(:,:) = zero
      !
      CALL start_bar_type( barra, 'dfpt', MAX(m,1) * kpt_pool%nloc )
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
#if defined(__CUDA)
         IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), k_grid%p_cart(1,ik), vkb, .TRUE. )
#else
         IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), k_grid%p_cart(1,ik), vkb, .FALSE. )
#endif
         !
         nbndval = nbnd_occ(iks)
         IF (l_frac_occ) THEN
            nbndval_full = nbnd_occ_full(iks)
            nbndval_frac = nbndval - nbndval_full
            ALLOCATE(psi_dvpsi(nbndval_frac,band_group%nloc))
            !$acc enter data create(psi_dvpsi)
         ENDIF
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
            !$acc update device(evc)
         ENDIF
         !
         IF (gamma_only) THEN
            !
            ikqs = iks
            g0 = 0._DP
            !
         ELSE
            !
            ! ... Find G0 and compute phase
            !
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
            !$acc update device(evckmq,phase)
         ENDIF
         !
         CALL set_eprec(nbndval,evc,eprec)
         !
         band_group_nloc = band_group%nloc
         !
         !$acc parallel loop present(eprec_loc,eprec,et_loc,et)
         DO lbnd = 1,band_group_nloc
            !
            ! ibnd = band_group%l2g(lbnd)
            !
            ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
            eprec_loc(lbnd) = eprec(ibnd)
            et_loc(lbnd) = et(ibnd,ikqs)
            !
         ENDDO
         !$acc end parallel
         !
         DO ipert = 1, m
            !
            aux_g(1:npwq) = dvg(1:npwq,ipert)
            !$acc update device(aux_g)
            !
            ! ... inverse Fourier transform of the perturbation: (q+)G ---> R
            !
            IF (gamma_only) THEN
               CALL single_invfft_gamma(dffts,npwq,npwqx,aux_g,aux_r,TRIM(fftdriver))
            ELSE
               CALL single_invfft_k(dffts,npwq,npwqx,aux_g,aux_r,'Wave',igq_q(:,iq))
            ENDIF
            !
            ! The perturbation is in aux_r
            !
            IF(gamma_only) THEN
               !
               ! double bands @ gamma
               DO lbnd = 1,band_group%nloc-MOD(band_group%nloc,2),2
                  !
                  ibnd = band_group%l2g(lbnd)
                  jbnd = band_group%l2g(lbnd+1)
                  !
                  CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
                  !
                  !$acc parallel loop present(aux_r)
                  DO ir = 1,dffts_nnr
                     psic(ir) = psic(ir)*REAL(aux_r(ir),KIND=DP)
                  ENDDO
                  !$acc end parallel
                  !
                  CALL double_fwfft_gamma(dffts,npw,npwx,psic,dvpsi(:,lbnd),dvpsi(:,lbnd+1),'Wave')
                  !
               ENDDO
               !
               ! single band @ gamma
               IF( MOD(band_group%nloc,2) == 1 ) THEN
                  !
                  lbnd = band_group%nloc
                  ibnd = band_group%l2g(lbnd)
                  !
                  CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
                  !
                  !$acc parallel loop present(aux_r)
                  DO ir = 1,dffts_nnr
                     psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(aux_r(ir),KIND=DP),KIND=DP)
                  ENDDO
                  !$acc end parallel
                  !
                  CALL single_fwfft_gamma(dffts,npw,npwx,psic,dvpsi(:,lbnd),'Wave')
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
                  CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1:npwx,ibnd),psic,'Wave',igk_k(:,ikqs))
                  !
                  ! ... construct right-hand-side term of Sternheimer equation:
                  ! ... product of wavefunction at [k-q], phase and perturbation in real space
                  !
                  !$acc parallel loop present(phase,aux_r)
                  DO ir = 1,dffts_nnr
                     psic(ir) = psic(ir)*phase(ir)*aux_r(ir)
                  ENDDO
                  !$acc end parallel
                  !
                  ! Fourier transform product of wf at [k-q], phase and
                  ! perturbation of wavevector q: R ---> (k+)G
                  !
                  CALL single_fwfft_k(dffts,npw,npwx,psic,dvpsi(1:npwx,lbnd),'Wave',igk_k(:,iks))
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
                     CALL single_invfft_k(dffts,npwkq,npwx,evckmq(npwx+1:npwx*2,ibnd),psic,'Wave',igk_k(:,ikqs))
                     !
                     !$acc parallel loop present(phase,aux_r)
                     DO ir = 1,dffts_nnr
                        psic(ir) = psic(ir)*phase(ir)*aux_r(ir)
                     ENDDO
                     !$acc end parallel
                     !
                     CALL single_fwfft_k(dffts,npw,npwx,psic,dvpsi(npwx+1:npwx*2,lbnd),'Wave',igk_k(:,iks))
                     !
                  ENDDO
               ENDIF
               !
            ENDIF
            !
            IF(l_frac_occ) THEN
               !
               ! Compute <psi_j| dV |psi_i>
               !
               CALL glbrak_gamma(evc(1,nbndval_full+1),dvpsi,psi_dvpsi,npw,npwx,nbndval_frac,&
               & band_group%nloc,nbndval_frac,npol)
               !$acc update host(psi_dvpsi)
               !
               CALL mp_sum(psi_dvpsi,intra_bgrp_comm)
               !
            ENDIF
            !
            ! - P_c | dvpsi >
            !
#if defined(__CUDA)
            CALL reallocate_ps_gpu(nbndval,band_group%nloc)
#endif
            !
            CALL apply_alpha_pc_to_m_wfcs( nbndval, band_group%nloc, dvpsi, (-1._DP,0._DP) )
            !
            CALL precondition_m_wfcts( band_group%nloc, dvpsi, dpsi, eprec_loc )
            !
            IF (l_dost) THEN
               !
               ! The Sternheimer operator is (H_k - E_(k-q) + alpha * P_v)
               ! The Hamiltonian is evaluated at the k-point current_k in h_psi
               ! (see also PHonon/PH/cch_psi_all.f90, where H_(k+q) is evaluated)
               !
               CALL linsolve_sternheimer_m_wfcts( nbndval, band_group%nloc, dvpsi, dpsi, et_loc, eprec_loc, tr2, ierr )
               !
               IF(ierr /= 0) &
                  WRITE(stdout, '(7X,"** WARNING : PERT ",I8," iks ",I8," not converged, ierr = ",I8)') ipert,iks,ierr
               !
            ENDIF
            !
            IF(l_frac_occ) THEN
               !
               ! Add to dpsi: \sum_j <psi_j| dV | psi_i> / (e_i - e_j) |psi_j>
               !
               DO lbnd = 1,band_group%nloc
                  !
                  ibnd = band_group%l2g(lbnd)
                  !
                  DO jbnd = nbndval_full+1,nbndval
                     !
                     IF(jbnd <= ibnd) THEN
                        psi_dvpsi(jbnd-nbndval_full,lbnd) = 0._DP
                        CYCLE
                     ENDIF
                     !
                     this_occ = occupation(ibnd,iks)
                     docc = occupation(ibnd,iks) - occupation(jbnd,iks)
                     !
                     IF(ABS(docc) < docc_thr) THEN
                        psi_dvpsi(jbnd-nbndval_full,lbnd) = 0._DP
                        CYCLE
                     ENDIF
                     !
                     de = et(ibnd,iks) - et(jbnd,iks)
                     IF(ABS(de) < de_thr) CALL errore('dfpt','fractional occupation degenerate orbitals',1)
                     !
                     psi_dvpsi(jbnd-nbndval_full,lbnd) = psi_dvpsi(jbnd-nbndval_full,lbnd) * docc / this_occ / de
                     !
                  ENDDO
                  !
               ENDDO
               !
               !$acc update device(psi_dvpsi)
               !$acc host_data use_device(evc,psi_dvpsi,dpsi)
               CALL DGEMM('N','N',2*npw,band_group%nloc,nbndval_frac,1._DP,evc(1,nbndval_full+1),2*npwx,&
               & psi_dvpsi,nbndval_frac,1._DP,dpsi,2*npwx)
               !$acc end host_data
               !
            ENDIF
            !
            !$acc kernels present(aux_r)
            aux_r(:) = zero
            !$acc end kernels
            !
            IF(gamma_only) THEN
               !
               ! double band @ gamma
               DO lbnd = 1,band_group%nloc
                  !
                  ibnd = band_group%l2g(lbnd)
                  !
                  CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),dpsi(:,lbnd),psic,'Wave')
                  !
                  this_occ = occupation(ibnd,iks)
                  !$acc parallel loop present(aux_r)
                  DO ir = 1,dffts_nnr
                     aux_r(ir) = aux_r(ir)+CMPLX(this_occ*REAL(psic(ir),KIND=DP)*AIMAG(psic(ir)),KIND=DP)
                  ENDDO
                  !$acc end parallel
                  !
               ENDDO
               !
            ELSE
               !
               DO lbnd = 1,band_group%nloc
                  !
                  ibnd = band_group%l2g(lbnd)
                  !
                  ! inverse Fourier transform of wavefunction at [k-q]: (k-q+)G ---> R
                  !
                  CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1:npwx,ibnd),psic,'Wave',igk_k(:,ikqs))
                  !
                  ! inverse Fourier transform of perturbed wavefunction: (k+)G ---> R
                  !
                  CALL single_invfft_k(dffts,npw,npwx,dpsi(1:npwx,lbnd),dpsic,'Wave',igk_k(:,iks))
                  !
                  !$acc parallel loop present(aux_r,phase,dpsic)
                  DO ir = 1,dffts_nnr
                     aux_r(ir) = aux_r(ir)+CONJG(psic(ir)*phase(ir))*dpsic(ir)
                  ENDDO
                  !$acc end parallel
                  !
               ENDDO
               !
               IF (noncolin) THEN
                  DO lbnd = 1,band_group%nloc
                     !
                     ibnd = band_group%l2g(lbnd)
                     !
                     CALL single_invfft_k(dffts,npwkq,npwx,evckmq(npwx+1:npwx*2,ibnd),psic,'Wave',igk_k(:,ikqs))
                     !
                     CALL single_invfft_k(dffts,npw,npwx,dpsi(npwx+1:npwx*2,lbnd),dpsic,'Wave',igk_k(:,iks))
                     !
                     !$acc parallel loop present(aux_r,phase,dpsic)
                     DO ir = 1,dffts_nnr
                        aux_r(ir) = aux_r(ir)+CONJG(psic(ir)*phase(ir))*dpsic(ir)
                     ENDDO
                     !$acc end parallel
                     !
                  ENDDO
               ENDIF
               !
            ENDIF
            !
            ! Sum up aux_r from band groups
            !
            IF (nbgrp > 1) THEN
               !$acc host_data use_device(aux_r)
               CALL mp_sum(aux_r,inter_bgrp_comm)
               !$acc end host_data
            ENDIF
            !
            ! The perturbation is in aux_r
            !
            IF(gamma_only) THEN
               CALL single_fwfft_gamma(dffts,npwq,npwqx,aux_r,aux_g,TRIM(fftdriver))
            ELSE
               CALL single_fwfft_k(dffts,npwq,npwqx,aux_r,aux_g,'Wave',igq_q(:,iq))
            ENDIF
            !
            !$acc update host(aux_g)
            !
            DO ig = 1,npwq ! pert acts only on body
               dng(ig,ipert) = dng(ig,ipert)+2._DP*k_grid%weight(iks_g)*aux_g(ig)/omega
            ENDDO
            !
            CALL update_bar_type( barra, 'dfpt', 1 )
            !
         ENDDO ! ipert
         !
         IF (l_frac_occ) THEN
            !$acc exit data delete(psi_dvpsi)
            DEALLOCATE(psi_dvpsi)
         ENDIF
         !
         IF( m == 0 ) CALL update_bar_type( barra, 'dfpt', 1 )
         !
      ENDDO ! K-POINT and SPIN
      !
      IF ( gamma_only ) THEN
         IF ( gstart == 2 ) dng(1,1:m) = zero
      ELSE
         IF ( gstart == 2 .AND. q_grid%l_pIsGamma(iq) ) dng(1,1:m) = zero
      ENDIF
      !
      CALL mp_sum(dng,inter_pool_comm)
      !
      !$acc exit data delete(eprec,eprec_loc,et_loc,dvpsi,dpsi,aux_r,aux_g)
      DEALLOCATE(eprec)
      DEALLOCATE(eprec_loc)
      DEALLOCATE(et_loc)
      DEALLOCATE(dvpsi)
      DEALLOCATE(dpsi)
      DEALLOCATE(aux_r)
      DEALLOCATE(aux_g)
      !
      IF (.NOT. gamma_only) THEN
         !$acc exit data delete(dpsic,phase,evckmq)
         DEALLOCATE(dpsic)
         DEALLOCATE(phase)
         DEALLOCATE(evckmq)
      ENDIF
      !
#if defined(__CUDA)
      CALL deallocate_gpu()
      CALL deallocate_linsolve_gpu()
#endif
      !
      CALL stop_bar_type( barra, 'dfpt' )
      !
    END SUBROUTINE
    !
END MODULE
