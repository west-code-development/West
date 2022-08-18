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
SUBROUTINE calc_vxc( sigma_vxcl, sigma_vxcnl )
  !-----------------------------------------------------------------------
  !
  ! store in sigma_vxc(n,iks) = < qp_bands(n),iks | V_vxc | qp_bands(n),iks >     n = 1,n_bands
  !
  ! IF (l_enable_off_diagonal .AND. l_full) store in
  ! sigma_vxc_full(ijpmap(m,n),iks) = < qp_bands(m),iks | V_vxc | qp_bands(n),iks >     n,m = 1,n_bands & m <= n
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout
  USE scf,                  ONLY : rho,rho_core,rhog_core
  USE gvect,                ONLY : gstart
  USE fft_base,             ONLY : dfftp,dffts
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,lsda,igk_k,nspin,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k
  USE wavefunctions,        ONLY : evc,psic
  USE westcom,              ONLY : qp_bands,n_bands,iuwfc,lrwfc,l_enable_off_diagonal,sigma_vxcl_full,&
                                 & sigma_vxcnl_full,ijpmap
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : vkb,nkb
  USE uspp_init,            ONLY : init_us_2
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE xc_lib,               ONLY : xclib_dft_is
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : kpt_pool
  USE exx,                  ONLY : vexx
  USE types_bz_grid,        ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP),INTENT(OUT) :: sigma_vxcl( n_bands, k_grid%nps )
  REAL(DP),INTENT(OUT) :: sigma_vxcnl( n_bands, k_grid%nps )
  !
  ! Workspace
  !
  REAL(DP) :: etxc_
  REAL(DP) :: vtxc_
  REAL(DP), ALLOCATABLE :: vxc(:,:)
  INTEGER :: ib,ir,iks,iks_g,jb_glob,ipair
  COMPLEX(DP) :: braket
  REAL(DP) :: nnr
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP), ALLOCATABLE :: xpsi(:,:),vxpsi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic1(:)
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  TYPE(idistribute) :: gwbnd
  !
  ALLOCATE( vxc(dfftp%nnr,nspin) )
  !
  IF (l_enable_off_diagonal) ALLOCATE( psic1 (dfftp%nnr) )
  !
  WRITE(stdout,'(5x,a)') ' '
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') 'Vxc'
  CALL io_push_bar()
  !
  gwbnd = idistribute()
  CALL gwbnd%init(n_bands,'i','n_bands',.FALSE.)
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc_, vtxc_, vxc )
  !
  sigma_vxcl = 0._DP
  sigma_vxcnl = 0._DP
  !
  nnr = REAL( dfftp%nr1*dfftp%nr2*dfftp%nr3, KIND=DP )
  !
  barra_load = kpt_pool%nloc
  CALL start_bar_type( barra, 'sigmavxc', barra_load )
  !
  ! LOOP
  !
  DO iks = 1, kpt_pool%nloc ! KPOINT-SPIN
     !
     iks_g = kpt_pool%l2g(iks)
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     ! NON-HYBRID CONTRIBUTION TO VXC
     !
     IF( gwbnd%nloc>0 ) THEN
        !
        IF(gamma_only) THEN
           !
           DO ib = 1, gwbnd%nloc
              CALL single_invfft_gamma(dffts,npw,npwx,evc(:,qp_bands(gwbnd%l2g(ib))),psic,'Wave')
              !
              DO jb_glob = 1, n_bands
                 !
                 braket = 0._DP
                 !
                 IF (l_enable_off_diagonal) ipair = ijpmap(jb_glob,gwbnd%l2g(ib))
                 !
                 IF (l_enable_off_diagonal .AND. jb_glob < gwbnd%l2g(ib)) THEN
                    CALL single_invfft_gamma(dffts,npw,npwx,evc(:,qp_bands(jb_glob)),psic1,'Wave')
                    DO ir = 1, dfftp%nnr
                       braket = braket + psic(ir)*CONJG(psic1(ir)) * vxc(ir,current_spin)
                    ENDDO
                    sigma_vxcl_full(ipair,iks_g) = REAL(braket,KIND=DP) / nnr
                 ELSEIF ( jb_glob == gwbnd%l2g(ib) ) THEN
                    DO ir = 1, dfftp%nnr
                       braket = braket + psic(ir)*CONJG(psic(ir)) * vxc(ir,current_spin)
                    ENDDO
                    sigma_vxcl(gwbnd%l2g(ib),iks_g)&
                    &= REAL(braket,KIND=DP) / nnr
                    IF (l_enable_off_diagonal) sigma_vxcl_full(ipair,iks_g) = REAL(braket,KIND=DP) / nnr
                 ENDIF
                 !
              ENDDO
              !
           ENDDO
           !
        ELSE
           !
           DO ib = 1, gwbnd%nloc
              CALL single_invfft_k(dffts,npw,npwx,evc(:,qp_bands(gwbnd%l2g(ib))),psic,'Wave',&
              & igk_k(:,current_k))
              braket = 0._DP
              DO ir = 1, dfftp%nnr
                 braket = braket + psic(ir)*CONJG(psic(ir)) * vxc(ir,current_spin)
              ENDDO
              sigma_vxcl(gwbnd%l2g(ib),iks_g) = REAL(braket,KIND=DP) / nnr
           ENDDO
           !
           IF(noncolin) THEN
              !
              DO ib = 1,gwbnd%nloc
                 CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx:npwx*2,qp_bands(gwbnd%l2g(ib))),&
                 & psic,'Wave',igk_k(:,current_k))
                 braket = 0._DP
                 DO ir = 1, dfftp%nnr
                    braket = braket + psic(ir)*CONJG(psic(ir)) * vxc(ir,current_spin)
                 ENDDO
                 sigma_vxcl(gwbnd%l2g(ib),iks_g) = &
                 & sigma_vxcl(gwbnd%l2g(ib),iks_g) + REAL(braket,KIND=DP) / nnr
              ENDDO
              !
           ENDIF
           !
        ENDIF
        !
     ENDIF
     !
     ! HYBRID CONTRIBUTION TO VXC
     !
     IF(  xclib_dft_is('hybrid') ) THEN
        !
        IF( gwbnd%nloc>0 ) THEN
           !
           ALLOCATE( vxpsi(npwx*npol,gwbnd%nloc ) )
           ALLOCATE(  xpsi(npwx*npol,gwbnd%nloc ) )
           !
           xpsi = 0._DP
           DO ib=1,gwbnd%nloc
              xpsi(:,ib) = evc(:,qp_bands(gwbnd%l2g(ib)))
           ENDDO
           vxpsi = 0._DP
           CALL vexx( npwx, npw, gwbnd%nloc, xpsi, vxpsi )
           !
           IF( gamma_only ) THEN
              !
              DO ib = 1,gwbnd%nloc
                 !
                 DO jb_glob = 1, n_bands
                    !
                    IF (l_enable_off_diagonal) ipair = ijpmap(jb_glob,gwbnd%l2g(ib))
                    !
                    IF (l_enable_off_diagonal .AND. jb_glob < gwbnd%l2g(ib)) THEN
                       braket = 2._DP * REAL( ZDOTC( npw, evc(1,qp_bands(jb_glob)),1,vxpsi(1,ib),1), KIND=DP )
                       IF(gstart==2) braket = braket - REAL( evc(1,qp_bands(jb_glob)), KIND=DP) * REAL( vxpsi(1,ib), KIND=DP)
                       sigma_vxcnl_full(ipair,iks_g) = REAL( braket, KIND=DP )
                    ELSEIF ( jb_glob == gwbnd%l2g(ib) ) THEN
                       braket = 2._DP * DDOT( 2*npw, xpsi(1,ib), 1, vxpsi(1,ib), 1)
                       IF(gstart==2) braket = braket - REAL( xpsi(1,ib), KIND=DP) * REAL( vxpsi(1,ib), KIND=DP)
                       sigma_vxcnl(gwbnd%l2g(ib),iks_g) = REAL( braket, KIND=DP )
                       IF (l_enable_off_diagonal) sigma_vxcnl_full(ipair,iks_g) = REAL( braket, KIND=DP )
                    ENDIF
                    !
                 ENDDO
                 !
              ENDDO
              !
           ELSE
              !
              DO ib = 1,gwbnd%nloc
                 braket = ZDOTC( npw, xpsi(1,ib),1,vxpsi(1,ib),1)
                 sigma_vxcnl(gwbnd%l2g(ib),iks_g) = REAL( braket, KIND=DP )
              ENDDO
              !
              IF(noncolin) THEN
                 !
                 DO ib = 1, gwbnd%nloc
                    braket = ZDOTC( npw, xpsi(1+npwx,ib),1,vxpsi(1+npwx,ib),1)
                    sigma_vxcnl(gwbnd%l2g(ib),iks_g) = &
                    & sigma_vxcnl(gwbnd%l2g(ib),iks_g) + REAL( braket, KIND=DP )
                 ENDDO
              ENDIF
              !
           ENDIF
           !
           DEALLOCATE( vxpsi )
           DEALLOCATE(  xpsi )
           !
        ENDIF
        !
     ENDIF
     !
     CALL update_bar_type( barra, 'sigmavxc', 1 )
     !
  ENDDO
  !
  CALL stop_bar_type( barra, 'sigmavxc' )
  !
  CALL mp_sum( sigma_vxcl, intra_bgrp_comm )
  CALL mp_sum( sigma_vxcnl, intra_bgrp_comm )
  CALL mp_sum( sigma_vxcl, inter_pool_comm )
  CALL mp_sum( sigma_vxcnl, inter_pool_comm )
  CALL mp_sum( sigma_vxcl, inter_image_comm )
  CALL mp_sum( sigma_vxcnl, inter_image_comm )
  !
  IF (l_enable_off_diagonal) THEN
     CALL mp_sum( sigma_vxcl_full, intra_bgrp_comm )
     CALL mp_sum( sigma_vxcnl_full, intra_bgrp_comm )
     CALL mp_sum( sigma_vxcl_full, inter_image_comm )
     CALL mp_sum( sigma_vxcnl_full, inter_image_comm )
  ENDIF
  !
  DEALLOCATE( vxc )
  IF (l_enable_off_diagonal) DEALLOCATE( psic1 )
  !
END SUBROUTINE
