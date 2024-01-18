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
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,lsda,igk_k,nspin,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k
  USE westcom,              ONLY : qp_bands,n_bands,iuwfc,lrwfc,l_enable_off_diagonal,&
                                 & sigma_vxcl_full,sigma_vxcnl_full,ijpmap
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE xc_lib,               ONLY : xclib_dft_is
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : kpt_pool
  USE exx,                  ONLY : use_ace,vexx,vexxace_gamma,vexxace_k
  USE types_bz_grid,        ONLY : k_grid
  USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu
#endif
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
  REAL(DP) :: ee_
  REAL(DP), ALLOCATABLE :: vxc(:,:)
  INTEGER :: ib,ir,iks,iks_g,is,jb_glob,ipair
  COMPLEX(DP) :: braket
  REAL(DP) :: nnr
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  INTEGER :: dfftp_nnr
  COMPLEX(DP), ALLOCATABLE :: xpsi(:,:),vxpsi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic1(:)
  !$acc declare device_resident(psic1)
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  TYPE(idistribute) :: gwbnd
  !
  ALLOCATE( vxc(dfftp%nnr,nspin) )
  !
  IF (l_enable_off_diagonal) ALLOCATE( psic1 (dfftp%nnr) )
  !
  WRITE(stdout,*)
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') 'Vxc'
  CALL io_push_bar()
  !
  gwbnd = idistribute()
  CALL gwbnd%init(n_bands,'i','n_bands',.FALSE.)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  dfftp_nnr = dfftp%nnr
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc_, vtxc_, vxc )
  !
  !$acc enter data copyin(vxc)
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
     is = k_grid%is(iks_g)
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     npw = ngk(iks)
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     ! NON-HYBRID CONTRIBUTION TO VXC
     !
     IF( gwbnd%nloc>0 ) THEN
        !
        IF(gamma_only) THEN
           !
           DO ib = 1, gwbnd%nloc
              CALL single_invfft_gamma(dffts,npw,npwx,evc(:,qp_bands(gwbnd%l2g(ib),is)),psic,'Wave')
              !
              DO jb_glob = 1, n_bands
                 !
                 braket = 0._DP
                 !
                 IF (l_enable_off_diagonal) ipair = ijpmap(jb_glob,gwbnd%l2g(ib))
                 !
                 IF (l_enable_off_diagonal .AND. jb_glob < gwbnd%l2g(ib)) THEN
                    CALL single_invfft_gamma(dffts,npw,npwx,evc(:,qp_bands(jb_glob,is)),psic1,'Wave')
                    !$acc parallel loop reduction(+:braket) present(psic1,vxc) copy(braket)
                    DO ir = 1, dfftp_nnr
                       braket = braket + psic(ir) * CONJG(psic1(ir)) * vxc(ir,current_spin)
                    ENDDO
                    !$acc end parallel
                    sigma_vxcl_full(ipair,iks_g) = REAL(braket,KIND=DP) / nnr
                 ELSEIF ( jb_glob == gwbnd%l2g(ib) ) THEN
                    !$acc parallel loop reduction(+:braket) present(vxc) copy(braket)
                    DO ir = 1, dfftp_nnr
                       braket = braket + psic(ir) * CONJG(psic(ir)) * vxc(ir,current_spin)
                    ENDDO
                    !$acc end parallel
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
              CALL single_invfft_k(dffts,npw,npwx,evc(:,qp_bands(gwbnd%l2g(ib),is)),psic,'Wave',&
              & igk_k(:,current_k))
              braket = 0._DP
              !$acc parallel loop reduction(+:braket) present(vxc) copy(braket)
              DO ir = 1, dfftp_nnr
                 braket = braket + psic(ir) * CONJG(psic(ir)) * vxc(ir,current_spin)
              ENDDO
              !$acc end parallel
              sigma_vxcl(gwbnd%l2g(ib),iks_g) = REAL(braket,KIND=DP) / nnr
           ENDDO
           !
           IF(noncolin) THEN
              !
              DO ib = 1,gwbnd%nloc
                 CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx:npwx*2,qp_bands(gwbnd%l2g(ib),is)),&
                 & psic,'Wave',igk_k(:,current_k))
                 braket = 0._DP
                 !$acc parallel loop reduction(+:braket) present(vxc) copy(braket)
                 DO ir = 1, dfftp_nnr
                    braket = braket + psic(ir) * CONJG(psic(ir)) * vxc(ir,current_spin)
                 ENDDO
                 !$acc end parallel
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
     IF( xclib_dft_is('hybrid') ) THEN
        !
        IF( gwbnd%nloc>0 ) THEN
           !
           ALLOCATE( vxpsi(npwx*npol,gwbnd%nloc ) )
           ALLOCATE(  xpsi(npwx*npol,gwbnd%nloc ) )
           !
           xpsi = 0._DP
           DO ib=1,gwbnd%nloc
              xpsi(:,ib) = evc(:,qp_bands(gwbnd%l2g(ib),is))
           ENDDO
           !
           vxpsi = 0._DP
           IF( use_ace ) THEN
              IF( gamma_only ) THEN
                 CALL vexxace_gamma( npwx, gwbnd%nloc, xpsi, ee_, vxpsi )
              ELSE
                 CALL vexxace_k( npwx, gwbnd%nloc, xpsi, ee_, vxpsi )
              ENDIF
           ELSE
              CALL vexx( npwx, npw, gwbnd%nloc, xpsi, vxpsi )
           ENDIF
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
                       braket = 2._DP * REAL( ZDOTC( npw, evc(1,qp_bands(jb_glob,is)),1,vxpsi(1,ib),1), KIND=DP )
                       IF(gstart==2) braket = braket - REAL( evc(1,qp_bands(jb_glob,is)), KIND=DP) * REAL( vxpsi(1,ib), KIND=DP)
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
     CALL mp_sum( sigma_vxcl_full, inter_pool_comm )
     CALL mp_sum( sigma_vxcnl_full, inter_pool_comm )
     CALL mp_sum( sigma_vxcl_full, inter_image_comm )
     CALL mp_sum( sigma_vxcnl_full, inter_image_comm )
  ENDIF
  !
  !$acc exit data delete(vxc)
  DEALLOCATE( vxc )
  IF (l_enable_off_diagonal) DEALLOCATE( psic1 )
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
