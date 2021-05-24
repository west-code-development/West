!
! Copyright (C) 2015-2017 M. Govoni
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
      USE wavefunctions_module,  ONLY : evc,psic
      USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
      USE mp_global,             ONLY : inter_image_comm,inter_pool_comm,my_image_id,inter_bgrp_comm
      USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
      USE fft_at_gamma,          ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
      USE buffers,               ONLY : get_buffer
      USE noncollin_module,      ONLY : noncolin,npol
      USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
      USE pwcom,                 ONLY : current_spin,nelup,neldw,isk,xk,npw,npwx,lsda,&
                                      & current_k,ngk,igk_k
      USE cell_base,             ONLY : omega
      USE control_flags,         ONLY : gamma_only
      USE uspp,                  ONLY : nkb, vkb
      USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,npwqx,npwq,igq_q,fftdriver
      USE io_push,               ONLY : io_push_title
      USE mp_world,              ONLY : world_comm
      USE types_bz_grid,         ONLY : k_grid, q_grid, compute_phase
      USE class_idistribute,     ONLY : idistribute
      USE distribution_center,   ONLY : occband
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
      INTEGER :: ipert, ig, ir, ibnd, ibnd2, lbnd, iks, ikqs, ikq, ik, is
      INTEGER :: nbndval, ierr
      INTEGER :: npwkq
      !
      REAL(DP) :: g0(3)
      REAL(DP) :: anorm
      REAL(DP), ALLOCATABLE :: eprec(:)
      REAL(DP), ALLOCATABLE :: eprec_loc(:)
      REAL(DP), ALLOCATABLE :: et_loc(:)
      !
      COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:),dpsi(:,:)
      COMPLEX(DP), ALLOCATABLE :: aux_r(:),aux_g(:)
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
      occband = idistribute()
      !
      dng=0.0_DP
      !
      CALL start_bar_type( barra, 'dfpt', MAX(m,1) * k_grid%nps )
      !
      IF (.NOT.gamma_only) THEN
         ALLOCATE( evckmq(npwx*npol,nbnd) )
         ALLOCATE( phase(dffts%nnr) )
      ENDIF
      !
      DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
         !
         ik = k_grid%ip(iks)
         is = k_grid%is(iks)
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
         !
         CALL occband%init( nbndval, 'b', 'occband', .FALSE. )
         !
         ! ... Number of G vectors for PW expansion of wfs at k
         !
         npw = ngk(iks)
         !
         ! ... Read wavefuctions at k in G space, for all bands, and store them in evc
         !
         IF(k_grid%nps>1) THEN
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
         ALLOCATE(eprec_loc(occband%nloc))
         ALLOCATE(et_loc(occband%nloc))
         CALL set_eprec(nbndval,evc(1,1),eprec)
         !
         DO lbnd = 1,occband%nloc
            ibnd = occband%l2g(lbnd)
            eprec_loc(lbnd) = eprec(ibnd)
            et_loc(lbnd) = et(ibnd,ikqs)
         ENDDO
         !
         ALLOCATE(dvpsi(npwx*npol,occband%nloc))
         ALLOCATE(dpsi(npwx*npol,occband%nloc))
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
            dvpsi=0._DP
            dpsi =0._DP
            !
            IF(gamma_only) THEN
               !
               ! double bands @ gamma
               DO lbnd = 1,occband%nloc-MOD(occband%nloc,2),2
                  !
                  ibnd = occband%l2g(lbnd)
                  ibnd2 = occband%l2g(lbnd+1)
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
               IF( MOD(occband%nloc,2) == 1 ) THEN
                  !
                  lbnd = occband%nloc
                  ibnd = occband%l2g(lbnd)
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
               DO lbnd = 1,occband%nloc
                  !
                  ibnd = occband%l2g(lbnd)
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
                  DO lbnd = 1,occband%nloc
                     !
                     ibnd = occband%l2g(lbnd)
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
            DEALLOCATE( aux_g )
            DEALLOCATE( aux_r )
            !
            ! - P_c | dvpsi >
            !
            CALL apply_alpha_pc_to_m_wfcs( nbndval, occband%nloc, dvpsi, (-1._DP,0._DP) )
            !
            CALL precondition_m_wfcts( occband%nloc, dvpsi, dpsi, eprec_loc )
            !
            IF( l_dost) THEN
               !
               ! The Sternheimer operator is (H_k - E_(k-q) + alpha * P_v)
               ! The Hamiltonian is evaluated at the k-point current_k in h_psi
               ! (see also PHonon/PH/cch_psi_all.f90, where H_(k+q) is evaluated)
               !
               CALL linsolve_sternheimer_m_wfcts (nbndval, occband%nloc, dvpsi, dpsi, et_loc, eprec_loc, tr2, ierr )
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
               DO lbnd = 1,occband%nloc
                  !
                  ibnd = occband%l2g(lbnd)
                  !
                  CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),dpsi(1,lbnd),psic,'Wave')
                  DO CONCURRENT (ir=1:dffts%nnr)
                     aux_r(ir) = aux_r(ir) + CMPLX( REAL( psic(ir),KIND=DP) * DIMAG( psic(ir)) , 0.0_DP, KIND=DP)
                  ENDDO
                  !
               ENDDO
               !
            ELSE
               !
               ALLOCATE( dpsic(dffts%nnr) )
               !
               DO lbnd = 1,occband%nloc
                  !
                  ibnd = occband%l2g(lbnd)
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
                  DO lbnd = 1,occband%nloc
                     !
                     ibnd = occband%l2g(lbnd)
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
               dng(ig,ipert) = dng(ig,ipert) + 2._DP * k_grid%weight(iks) * aux_g(ig) / omega
            ENDDO
            !
            DEALLOCATE( aux_g )
            DEALLOCATE( aux_r )
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
END MODULE
!!-----------------------------------------------------------------------
!SUBROUTINE dfpt (m,dvg,dng,tr2)
!  !-----------------------------------------------------------------------
!  !
!  USE kinds,                 ONLY : DP
!  USE io_global,             ONLY : stdout
!  USE wvfct,                 ONLY : nbnd,g2kin,et
!  USE fft_base,              ONLY : dfftp,dffts
!  USE gvect,                 ONLY : nl,gstart,g,ngm
!  USE wavefunctions_module,  ONLY : evc,psic
!  USE gvecw,                 ONLY : gcutw
!  USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
!  USE mp_global,             ONLY : inter_image_comm,inter_pool_comm,my_image_id
!  USE fft_at_gamma,          ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
!  USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
!  USE buffers,               ONLY : get_buffer
!  USE noncollin_module,      ONLY : noncolin,npol
!  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
!  USE pwcom,                 ONLY : current_spin,wk,nks,nelup,neldw,isk,igk_k,xk,npw,npwx,lsda,nkstot,&
!                                  & current_k,ngk
!  USE cell_base,             ONLY : tpiba2,omega
!  USE control_flags,         ONLY : gamma_only, io_level
!  USE io_files,              ONLY : tmp_dir, nwordwfc, iunwfc, diropn
!  USE uspp,                  ONLY : nkb, vkb, okvan
!  USE constants,             ONLY : e2,fpi
!  USE westcom,               ONLY : npwq,nbnd_occ,iuwfc,lrwfc,npwqx,fftdriver
!  USE io_push,               ONLY : io_push_title
!  USE mp_world,              ONLY : mpime,world_comm
!  USE types_coulomb,         ONLY : pot3D
!  !
!  IMPLICIT NONE
!  !
!  ! I/O
!  !
!  INTEGER,INTENT(IN) :: m
!  COMPLEX(DP),INTENT(IN) :: dvg(npwqx,m)
!  COMPLEX(DP),INTENT(OUT) :: dng(npwqx,m)
!  REAL(DP),INTENT(IN) :: tr2
!  !
!  ! Workspace
!  !
!  INTEGER :: ipert, ig, ir, ibnd, iks
!  INTEGER :: nbndval, ierr
!  !
!  REAL(DP) :: anorm, prod
!  REAL(DP),ALLOCATABLE :: eprec(:)
!  !
!  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:),dpsi(:,:)
!  COMPLEX(DP),ALLOCATABLE :: aux_r(:),aux_g(:)
!  COMPLEX(DP),ALLOCATABLE :: dpsic(:)
!  !
!  TYPE(bar_type) :: barra
!  !
!  LOGICAL :: conv_dfpt
!  LOGICAL :: exst,exst_mem
!  LOGICAL :: l_dost
!  !
!  CHARACTER(LEN=512) :: title
!  !
!  CALL mp_barrier( world_comm )
!  !
!  CALL report_dynamical_memory()
!  !
!  l_dost = ( tr2 >= 0._DP )
!  !
!  IF( l_dost ) THEN
!     WRITE(title,'(a,es14.6)') "Sternheimer eq. solver... with threshold = ", tr2
!  ELSE
!     WRITE(title,'(a,es14.6)') "Sternheimer eq. solver... with lite-solver"
!  ENDIF
!  CALL io_push_title(TRIM(ADJUSTL(title)))
!  !
!  dng=0.0_DP
!  !
!  CALL start_bar_type( barra, 'dfpt', MAX(m,1) * nks )
!  !IF(nks>1) CALL diropn(iuwfc,'wfc',lrwfc,exst)
!  !
!  DO iks = 1, nks  ! KPOINT-SPIN LOOP
!     !
!     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
!     !
!     current_k = iks
!     IF ( lsda ) current_spin = isk(iks)
!     call g2_kin( iks )
!     !
!     ! ... More stuff needed by the hamiltonian: nonlocal projectors
!     !
!     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
!     !npw = ngk(iks)
!     !
!     ! ... read in wavefunctions from the previous iteration
!     !
!     IF(nks>1) THEN
!        !iuwfc = 20
!        !lrwfc = nbnd * npwx * npol
!        !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
!        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
!        !CALL mp_bcast(evc,0,inter_image_comm)
!        !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
!        CALL mp_bcast(evc,0,inter_image_comm)
!     ENDIF
!!     !
!!     ! ... Needed for LDA+U
!!     !
!!     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
!!          CALL get_buffer ( wfcU, nwordwfcU, iunhub, iks )
!!     !
!!     current_k = iks
!!     current_spin = isk(iks)
!!     !
!!     CALL gk_sort(xk(1,iks),ngm,g,gcutw,npw,igk,g2kin)
!!     g2kin=g2kin*tpiba2
!!     !
!!     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
!!     !
!!     !
!!     CALL init_us_2 (npw, igk, xk (1, iks), vkb)
!     !
!     nbndval = nbnd_occ(iks)
!     IF(nbndval==0) THEN
!        CALL update_bar_type( barra,'dfpt', MAX(1,m) )
!        CYCLE
!     ENDIF
!     !
!     ALLOCATE(eprec(nbndval))
!     CALL set_eprec(nbndval,evc(1,1),eprec)
!     !
!     ALLOCATE(dvpsi(npwx*npol,nbndval))
!     ALLOCATE(dpsi(npwx*npol,nbndval))
!     !
!     DO ipert = 1, m
!        !
!        ALLOCATE(aux_g(npwqx))
!        ALLOCATE(aux_r(dffts%nnr))
!        !
!        aux_g = 0._DP
!        aux_r = 0._DP
!        !
!        DO ig = 1, npwq  ! perturbation acts only on body
!           aux_g(ig) = dvg(ig,ipert) * pot3D%sqvc(ig)
!        ENDDO
!        !
!        IF(gamma_only) THEN
!          CALL single_invfft_gamma(dffts,npwq,npwqx,aux_g,aux_r,TRIM(fftdriver))
!        ELSE
!          CALL single_invfft_k(dffts,npwq,npwqx,aux_g,aux_r,TRIM(fftdriver)) ! no igk
!        ENDIF
!        !
!        ! The perturbation is in aux_r
!        !
!        dvpsi=0._DP
!        dpsi =0._DP
!        !
!        IF(gamma_only) THEN
!           !
!           ! double bands @ gamma
!           DO ibnd=1,nbndval-MOD(nbndval,2),2
!              !
!              CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,ibnd+1),psic,'Wave')
!              DO ir=1,dffts%nnr
!                 psic(ir) = psic(ir) * REAL(aux_r(ir),KIND=DP)
!              ENDDO
!              CALL double_fwfft_gamma(dffts,npw,npwx,psic,dvpsi(1,ibnd),dvpsi(1,ibnd+1),'Wave')
!              !
!           ENDDO
!           !
!           ! single band @ gamma
!           IF( MOD(nbndval,2) == 1 ) THEN
!              ibnd=nbndval
!              !
!              CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),psic,'Wave')
!              DO ir=1,dffts%nnr
!                 psic(ir) = CMPLX( REAL(psic(ir),KIND=DP) * REAL(aux_r(ir),KIND=DP), 0._DP, KIND=DP)
!              ENDDO
!              CALL single_fwfft_gamma(dffts,npw,npwx,psic,dvpsi(1,ibnd),'Wave')
!              !
!           ENDIF
!           !
!        ELSE
!           !
!           ! only single bands
!           DO ibnd=1,nbndval
!              !
!              CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,current_k))
!              DO ir=1,dffts%nnr
!                 psic(ir) = psic(ir) * aux_r(ir)
!              ENDDO
!              CALL single_fwfft_k(dffts,npw,npwx,psic,dvpsi(1,ibnd),'Wave',igk_k(1,current_k))
!              !
!           ENDDO
!           !
!           IF(npol==2) THEN
!              DO ibnd=1,nbndval
!                 !
!                 CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ibnd),psic,'Wave',igk_k(1,current_k))
!                 DO ir=1,dffts%nnr
!                    psic(ir) = psic(ir) * aux_r(ir)
!                 ENDDO
!                 CALL single_fwfft_k(dffts,npw,npwx,psic,dvpsi(npwx+1,ibnd),'Wave',igk_k(1,current_k))
!                 !
!              ENDDO
!           ENDIF
!           !
!        ENDIF
!        !
!        DEALLOCATE(aux_g)
!        DEALLOCATE(aux_r)
!        !
!        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,dvpsi,(-1._DP,0._DP))
!        !
!        CALL precondition_m_wfcts( nbndval, dvpsi, dpsi, eprec )
!        !
!        IF( l_dost) THEN
!           !
!           CALL linsolve_sternheimer_m_wfcts (nbndval, nbndval, dvpsi, dpsi, et(1,iks), eprec, tr2, ierr )
!           !
!           IF(ierr/=0) THEN
!              WRITE(stdout, '(7X,"** WARNING : PERT ",i8," not converged, ierr = ",i8)') ipert,ierr
!           ENDIF
!           !
!        ENDIF
!        !
!        ALLOCATE(aux_r(dffts%nnr))
!        !
!        aux_r=0._DP
!        !
!        IF(gamma_only) THEN
!           !
!           ! double band @ gamma
!           DO ibnd=1,nbndval
!              !
!              CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),dpsi(1,ibnd),psic,'Wave')
!              DO ir=1,dffts%nnr
!                 prod =  REAL( psic(ir),KIND=DP) * DIMAG( psic(ir))
!                 aux_r(ir) = aux_r(ir) + CMPLX( prod, 0.0_DP, KIND=DP)
!              ENDDO
!              !
!           ENDDO
!           !
!        ELSE
!           !
!           ALLOCATE(dpsic(dffts%nnr))
!           !
!           ! only single bands
!           DO ibnd=1,nbndval
!              !
!              CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,current_k))
!              CALL single_invfft_k(dffts,npw,npwx,dpsi(1,ibnd),dpsic,'Wave',igk_k(1,current_k))
!              DO ir=1,dffts%nnr
!                 aux_r(ir) = aux_r(ir) + DCONJG(psic(ir))*dpsic(ir)
!              ENDDO
!              !
!           ENDDO
!           !
!           IF(npol==2) THEN
!              DO ibnd=1,nbndval
!                 !
!                 CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ibnd),psic,'Wave',igk_k(1,current_k))
!                 CALL single_invfft_k(dffts,npw,npwx,dpsi(npwx+1,ibnd),dpsic,'Wave',igk_k(1,current_k))
!                 DO ir=1,dffts%nnr
!                    aux_r(ir) = aux_r(ir) + DCONJG(psic(ir))*dpsic(ir)
!                 ENDDO
!                 !
!              ENDDO
!           ENDIF
!           !
!           DEALLOCATE(dpsic)
!           !
!        ENDIF
!        !
!        ! The perturbation is in aux_r
!        !
!        ALLOCATE(aux_g(npwqx))
!        IF(gamma_only) THEN
!           CALL single_fwfft_gamma(dffts,npwq,npwqx,aux_r,aux_g,TRIM(fftdriver))
!        ELSE
!           CALL single_fwfft_k(dffts,npwq,npwqx,aux_r,aux_g,TRIM(fftdriver)) ! no igk
!        ENDIF
!        !
!        DO ig=1,npwq ! pert acts only on body
!           dng(ig,ipert) = dng(ig,ipert) + 2._DP * wk(iks) * aux_g(ig) * pot3D%sqvc(ig) / omega
!        ENDDO
!        !
!        DEALLOCATE(aux_g)
!        DEALLOCATE(aux_r)
!        !
!        CALL update_bar_type( barra,'dfpt', 1 )
!        !
!     ENDDO ! ipert
!     !
!     IF( m==0 ) CALL update_bar_type( barra,'dfpt', 1 )
!     !
!     DEALLOCATE(eprec)
!     DEALLOCATE(dpsi)
!     DEALLOCATE(dvpsi)
!     !
!  ENDDO ! K-POINT and SPIN
!  !
!  IF( gstart==2 ) dng(1,1:m) = CMPLX( 0._DP, 0._DP, KIND=DP )
!  !
!  CALL mp_sum(dng,inter_pool_comm)
!  !
!  CALL mp_barrier( world_comm )
!  !
!  CALL stop_bar_type( barra, 'dfpt' )
!  !
!!  CALL close_buffer(iuwfc,'delete')
!  !IF( nks > 1 ) THEN
!  !   IF ( exst ) THEN
!  !      CLOSE(unit=iuwfc,status='keep')
!  !   ELSE
!  !      CLOSE(unit=iuwfc,status='delete')
!  !   ENDIF
!  !ENDIF
!  !
!END SUBROUTINE
