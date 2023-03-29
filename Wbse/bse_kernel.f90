!
! Copyright (C) 2015-2023 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE bse_kernel_gamma(current_spin,evc1,bse_kd1)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE fft_base,              ONLY : dffts
  USE mp,                    ONLY : mp_sum
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,double_invfft_gamma,double_fwfft_gamma
  USE mp_global,             ONLY : inter_bgrp_comm,nbgrp,my_bgrp_id
  USE pwcom,                 ONLY : npw,npwx,nks,isk,ngk
  USE westcom,               ONLY : nbndval0x,n_trunc_bands,l_local_repr,u_matrix,idx_matrix,&
                                  & n_bse_idx
  USE distribution_center,   ONLY : aband
  USE wbse_io,               ONLY : read_bse_pots_g
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : psic=>psic_d
  USE west_gpu,              ONLY : raux1,raux2,caux1,caux2,gaux
  USE cublas
#else
  USE wavefunctions,         ONLY : psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,aband%nlocx,nks)
  COMPLEX(DP), INTENT(INOUT) :: bse_kd1(npwx,aband%nlocx)
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,my_ibnd,my_jbnd,lbnd,ipair
  INTEGER :: ig,ir
  INTEGER :: nbnd_do,do_idx,current_spin_ikq,ikq
  INTEGER :: dffts_nnr,aband_nloc
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: raux1(:),raux2(:)
  COMPLEX(DP), ALLOCATABLE :: caux1(:,:),caux2(:,:),gaux(:)
#endif
  COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
  COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
  COMPLEX(DP), PARAMETER :: mone = (-1._DP,0._DP)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('bse_kernel')
#else
  CALL start_clock('bse_kernel')
#endif
  !
  nbnd_do = nbndval0x-n_trunc_bands
  dffts_nnr = dffts%nnr
  aband_nloc = aband%nloc
  !
#if !defined(__CUDA)
  ALLOCATE(raux1(dffts%nnr))
  ALLOCATE(raux2(dffts%nnr))
  ALLOCATE(caux1(npwx,nbnd_do))
  ALLOCATE(caux2(npwx,nbnd_do))
  ALLOCATE(gaux(npwx))
#endif
  !
  DO ikq = 1, nks
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     npw = ngk(ikq)
     !
     IF(l_local_repr) THEN
        !
        !$acc host_data use_device(evc1,u_matrix,caux1)
        CALL ZGEMM('N','N',npw,nbnd_do,aband_nloc,one,evc1(:,:,ikq),npwx,&
        & u_matrix(:,:,current_spin),aband_nloc,zero,caux1,npwx)
        !$acc end host_data
        !
     ELSE
        !
        !$acc kernels present(caux1)
        caux1(:,:) = zero
        !$acc end kernels
        !
        !$acc parallel loop collapse(2) present(caux1,evc1)
        DO lbnd = 1, aband_nloc
           !
           ! ibnd = aband%l2g(lbnd)
           !
           DO ig = 1, npw
              ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
              caux1(ig,ibnd) = evc1(ig,lbnd,ikq)
           ENDDO
           !
        ENDDO
        !$acc end parallel
        !
     ENDIF
     !
     !$acc update host(caux1)
     CALL mp_sum(caux1,inter_bgrp_comm)
     !$acc update device(caux1)
     !
     ! LOOP OVER BANDS AT KPOINT
     !
     !$acc kernels present(caux2)
     caux2(:,:) = zero
     !$acc end kernels
     !
     do_idx = n_bse_idx(current_spin)
     !
     DO lbnd = 1, aband%nloc, 2
        !
        my_ibnd = aband%l2g(lbnd)
        my_jbnd = aband%l2g(lbnd+1)
        !
        !$acc kernels present(raux1,raux2)
        raux1(:) = 0._DP
        raux2(:) = 0._DP
        !$acc end kernels
        !
        ! LOOP OVER BANDS AT QPOINT
        !
        DO ipair = 1, do_idx
           !
           ibnd = idx_matrix(ipair,1,current_spin)-n_trunc_bands
           jbnd = idx_matrix(ipair,2,current_spin)-n_trunc_bands
           !
           IF(ibnd == my_ibnd .OR. ibnd == my_jbnd) THEN
              !
              CALL read_bse_pots_g(gaux,ibnd,jbnd,current_spin)
              !
              !$acc update device(gaux)
              !
              !$acc host_data use_device(caux1,gaux)
              CALL double_invfft_gamma(dffts,npw,npwx,caux1(:,jbnd),gaux,psic,'Wave')
              !$acc end host_data
              !
              IF(ibnd == my_ibnd) THEN
                 !$acc parallel loop present(raux1)
                 DO ir = 1, dffts_nnr
                    raux1(ir) = raux1(ir)+REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
              IF(ibnd == my_jbnd) THEN
                 !$acc parallel loop present(raux2)
                 DO ir = 1, dffts_nnr
                    raux2(ir) = raux2(ir)+REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ENDIF
           !
        ENDDO
        !
        ! Back to reciprocal space
        !
        IF(lbnd < aband%nloc) THEN
           !$acc parallel loop present(raux1,raux2)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(raux1(ir),raux2(ir),KIND=DP)
           ENDDO
           !$acc end parallel
        ELSE
           !$acc parallel loop present(raux1)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(raux1(ir),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        IF(lbnd < aband%nloc) THEN
           !$acc host_data use_device(caux2)
           CALL double_fwfft_gamma(dffts,npw,npwx,psic,caux2(:,my_ibnd),caux2(:,my_jbnd),'Wave')
           !$acc end host_data
        ELSE
           !$acc host_data use_device(caux2)
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,caux2(:,my_ibnd),'Wave')
           !$acc end host_data
        ENDIF
        !
     ENDDO
     !
     !$acc update host(caux2)
     CALL mp_sum(caux2,inter_bgrp_comm)
     !$acc update device(caux2)
     !
     IF(l_local_repr) THEN
        !
        ! U^{+}(\xi)
        !
        !$acc host_data use_device(caux2,u_matrix,bse_kd1)
        CALL ZGEMM('N','C',npw,aband_nloc,nbnd_do,mone,caux2,npwx,u_matrix(:,:,current_spin),&
        & aband_nloc,one,bse_kd1,npwx)
        !$acc end host_data
        !
     ELSE
        !
        !$acc parallel loop collapse(2) present(bse_kd1,caux2)
        DO lbnd = 1, aband_nloc
           !
           ! ibnd = aband%l2g(lbnd)
           !
           DO ig = 1, npw
              ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
              bse_kd1(ig,lbnd) = bse_kd1(ig,lbnd)-caux2(ig,ibnd)
           ENDDO
        ENDDO
        !$acc end parallel
        !
     ENDIF
     !
  ENDDO
  !
#if !defined(__CUDA)
  DEALLOCATE(raux1)
  DEALLOCATE(raux2)
  DEALLOCATE(caux1)
  DEALLOCATE(caux2)
  DEALLOCATE(gaux)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('bse_kernel')
#else
  CALL stop_clock('bse_kernel')
#endif
  !
END SUBROUTINE
