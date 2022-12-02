!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE wbse_bse_kernel(current_spin,nbndval_k,evc1,bse_kd1)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE pwcom,                ONLY : npwx,nks
  USE westcom,              ONLY : nbndval0x
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin,nbndval_k
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,nbndval0x,nks)
  COMPLEX(DP), INTENT(INOUT) :: bse_kd1(npwx,nbndval0x)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('bse_kernel')
#else
  CALL start_clock('bse_kernel')
#endif
  !
  IF(gamma_only) THEN
     CALL bse_kernel_finite_field_gamma(current_spin,nbndval_k,evc1,bse_kd1)
  ELSE
     CALL errore('wbse_bse_kernel','Only Gamma is supported',1)
  ENDIF
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('bse_kernel')
#else
  CALL stop_clock('bse_kernel')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE bse_kernel_finite_field_gamma(current_spin,nbndval_k,evc1,bse_kd1)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE fft_base,              ONLY : dffts
  USE mp,                    ONLY : mp_sum
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,double_invfft_gamma,&
                                  & double_fwfft_gamma
  USE mp_global,             ONLY : inter_image_comm,inter_bgrp_comm
  USE pwcom,                 ONLY : npw,npwx,nks,isk,ngk
  USE westcom,               ONLY : l_lanczos,nbnd_occ,nbndval0x,l_local_repr,u_matrix,idx_matrix,&
                                  & n_bse_idx
  USE distribution_center,   ONLY : aband,bandpair
  USE wbse_io,               ONLY : read_bse_pots_g2r,read_bse_pots_g2g
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : psic=>psic_d
  USE west_gpu,              ONLY : raux1,raux2,caux1,caux2,kd1_ij,gaux
  USE cublas
#else
  USE wavefunctions,         ONLY : psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin,nbndval_k
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,nbndval0x,nks)
  COMPLEX(DP), INTENT(INOUT) :: bse_kd1(npwx,nbndval0x)
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,do_idx,n_done,ig,ir
  INTEGER :: ibnd_index,jbnd_index,il1,ig1
  INTEGER :: nbndval_q,current_spin_ikq,ikq
  INTEGER :: dffts_nnr
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: raux1(:),raux2(:)
  COMPLEX(DP), ALLOCATABLE :: caux1(:,:),caux2(:,:),kd1_ij(:,:),gaux(:)
#endif
  COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
  COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
  COMPLEX(DP), PARAMETER :: mone = (-1._DP,0._DP)
  !
  dffts_nnr = dffts%nnr
  !
#if !defined(__CUDA)
  IF(l_local_repr) THEN
     ALLOCATE(caux1(npwx,nbndval0x))
  ENDIF
  ALLOCATE(caux2(npwx,nbndval0x))
  do_idx = MAXVAL(n_bse_idx)
  ALLOCATE(raux1(dffts%nnr))
  IF(.NOT. l_lanczos) THEN
     ALLOCATE(kd1_ij(npwx,do_idx))
     ALLOCATE(raux2(dffts%nnr))
  ENDIF
  ALLOCATE(gaux(npwx))
#endif
  !
  DO ikq = 1, nks
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     nbndval_q = nbnd_occ(ikq)
     npw = ngk(ikq)
     !
     IF(l_local_repr) THEN
        !$acc host_data use_device(evc1,u_matrix,caux1)
        CALL ZGEMM('N','N',npw,nbndval0x,nbndval0x,one,evc1(:,:,ikq),npwx,&
        & u_matrix(:,:,current_spin),nbndval0x,zero,caux1,npwx)
        !$acc end host_data
     ENDIF
     !
     do_idx = n_bse_idx(current_spin)
     !
     IF(.NOT. l_lanczos) THEN
        !$acc kernels present(kd1_ij)
        kd1_ij(:,:) = zero
        !$acc end kernels
     ENDIF
     !
     !$acc kernels present(caux2)
     caux2 = zero
     !$acc end kernels
     !
     DO il1 = 1, bandpair%nlocx
        !
        ig1 = bandpair%l2g(il1) ! global index of n_total
        !
        IF(ig1 < 1 .OR. ig1 > do_idx) CYCLE
        !
        ibnd_index = idx_matrix(ig1,1,current_spin)
        jbnd_index = idx_matrix(ig1,2,current_spin)
        !
        IF(l_lanczos) THEN
           !
           ! READ response at iq,ik,ispin
           !
           CALL read_bse_pots_g2r(raux1,ibnd_index,jbnd_index,current_spin)
           !
           IF(l_local_repr) THEN
              !$acc host_data use_device(caux1)
              CALL single_invfft_gamma(dffts,npw,npwx,caux1(:,jbnd_index),psic,'Wave')
              !$acc end host_data
           ELSE
              !$acc host_data use_device(evc1)
              CALL single_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),psic,'Wave')
              !$acc end host_data
           ENDIF
           !
           !$acc parallel loop present(raux1)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*raux1(ir),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(gaux)
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,gaux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux2,gaux)
           DO ig = 1, npw
              caux2(ig,ibnd_index) = caux2(ig,ibnd_index)+gaux(ig)
           ENDDO
           !$acc end parallel
           !
        ELSE
           !
           CALL read_bse_pots_g2g(gaux,ibnd_index,jbnd_index,current_spin)
           !
           !$acc update device(gaux)
           !
           !$acc parallel loop present(kd1_ij,gaux)
           DO ig = 1, npw
              kd1_ij(ig,ig1) = gaux(ig)
           ENDDO
           !$acc end parallel
           !
        ENDIF
        !
     ENDDO
     !
     IF(l_lanczos) THEN
        !$acc update host(caux2)
        CALL mp_sum(caux2,inter_image_comm)
        !$acc update device(caux2)
     ELSE
        !$acc update host(kd1_ij)
        CALL mp_sum(kd1_ij,inter_image_comm)
        !$acc update device(kd1_ij)
     ENDIF
     !
     ! LOOP OVER BANDS AT KPOINT
     !
     IF(.NOT. l_lanczos) THEN
        !
        ! Davidson method
        !
        ibnd = 0
        jbnd = 0
        !
        DO il1 = 1, aband%nloc, 2
           !
           ibnd = aband%l2g(il1)
           jbnd = aband%l2g(il1+1)
           !
           !$acc kernels present(raux1,raux2)
           raux1(:) = 0._DP
           raux2(:) = 0._DP
           !$acc end kernels
           !
           n_done = 0
           !
           ! LOOP OVER BANDS AT QPOINT
           !
           DO ig1 = 1, do_idx
              !
              IF(n_done > MIN(do_idx,(nbndval_q+nbndval_k))) CYCLE
              !
              ibnd_index = idx_matrix(ig1,1,current_spin)
              jbnd_index = idx_matrix(ig1,2,current_spin)
              !
              IF(ibnd_index == ibnd) THEN
                 n_done = n_done+1
                 !
                 IF(l_local_repr) THEN
                    !$acc host_data use_device(caux1,kd1_ij)
                    CALL double_invfft_gamma(dffts,npw,npwx,caux1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
                    !$acc end host_data
                 ELSE
                    !$acc host_data use_device(evc1,kd1_ij)
                    CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),kd1_ij(:,ig1),psic,'Wave')
                    !$acc end host_data
                 ENDIF
                 !
                 !$acc parallel loop present(raux1)
                 DO ir = 1, dffts_nnr
                    raux1(ir) = raux1(ir)+REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
              IF(ibnd_index == jbnd) THEN
                 n_done = n_done+1
                 !
                 IF(l_local_repr) THEN
                    !$acc host_data use_device(caux1,kd1_ij)
                    CALL double_invfft_gamma(dffts,npw,npwx,caux1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
                    !$acc end host_data
                 ELSE
                    !$acc host_data use_device(evc1,kd1_ij)
                    CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),kd1_ij(:,ig1),psic,'Wave')
                    !$acc end host_data
                 ENDIF
                 !
                 !$acc parallel loop present(raux2)
                 DO ir = 1, dffts_nnr
                    raux2(ir) = raux2(ir)+REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ENDDO
           !
           ! Back to reciprocal space
           !
           IF(il1 < aband%nloc) THEN
              !$acc parallel loop present(raux2)
              DO ir = 1, dffts_nnr
                 psic(ir) = CMPLX(raux1(ir),raux2(ir),KIND=DP)
              ENDDO
              !$acc end parallel
           ELSE
              !$acc parallel loop present(raux2)
              DO ir = 1, dffts_nnr
                 psic(ir) = CMPLX(raux1(ir),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
           !$acc host_data use_device(caux2)
           IF(il1 < aband%nloc) THEN
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,caux2(:,ibnd),caux2(:,jbnd),'Wave')
           ELSE
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,caux2(:,ibnd),'Wave')
           ENDIF
           !$acc end host_data
           !
        ENDDO
        !
        !$acc update host(caux2)
        CALL mp_sum(caux2,inter_bgrp_comm)
        !$acc update device(caux2)
        !
     ENDIF
     !
     IF(l_local_repr) THEN
        !
        ! U^{+}(\xi)
        !
        !$acc host_data use_device(caux2,u_matrix,bse_kd1)
        CALL ZGEMM('N','C',npw,nbndval0x,nbndval0x,mone,caux2,npwx,u_matrix(:,:,current_spin),&
        & nbndval0x,one,bse_kd1,npwx)
        !$acc end host_data
        !
     ELSE
        !
        !$acc parallel loop collapse(2) present(bse_kd1,caux2)
        DO ibnd = 1, nbndval0x
           DO ig = 1, npw
              bse_kd1(ig,ibnd) = bse_kd1(ig,ibnd)-caux2(ig,ibnd)
           ENDDO
        ENDDO
        !$acc end parallel
        !
     ENDIF
     !
  ENDDO
  !
#if !defined(__CUDA)
  IF(l_local_repr) THEN
     DEALLOCATE(caux1)
  ENDIF
  DEALLOCATE(caux2)
  DEALLOCATE(raux1)
  IF(.NOT. l_lanczos) THEN
     DEALLOCATE(kd1_ij)
     DEALLOCATE(raux2)
  ENDIF
  DEALLOCATE(gaux)
#endif
  !
END SUBROUTINE
