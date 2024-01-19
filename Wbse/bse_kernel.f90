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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE bse_kernel_gamma(current_spin,evc1,bse_k1d,sf)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE fft_base,              ONLY : dffts
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,double_invfft_gamma,double_fwfft_gamma
  USE pwcom,                 ONLY : npw,npwx,isk,ngk
  USE exx,                   ONLY : exxalfa
  USE westcom,               ONLY : nbndval0x,n_trunc_bands,l_local_repr,u_matrix,idx_matrix,&
                                  & n_bse_idx,l_hybrid_tddft
  USE distribution_center,   ONLY : kpt_pool,band_group
  USE wbse_io,               ONLY : read_bse_pots_g
  USE wbse_bgrp,             ONLY : gather_bands
  USE west_mp,               ONLY : west_mp_wait
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : psic=>psic_d
  USE west_gpu,              ONLY : raux1,raux2,caux1,caux2,caux3,gaux
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
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,nbndval0x-n_trunc_bands)
  COMPLEX(DP), INTENT(INOUT) :: bse_k1d(npwx,band_group%nlocx)
  LOGICAL, INTENT(IN) :: sf
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,my_ibnd,my_jbnd,lbnd,ipair
  INTEGER :: ig,ir
  INTEGER :: nbnd_do,do_idx,current_spin_ikq,ikq,ikq_g,ikq_do
  INTEGER :: dffts_nnr,band_group_nloc,band_group_myoffset
  INTEGER :: req
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: raux1(:),raux2(:)
  COMPLEX(DP), ALLOCATABLE :: caux1(:,:),caux2(:,:),caux3(:,:),gaux(:)
#endif
  COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
  COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
#if defined(__CUDA)
  CALL start_clock_gpu('bse_kernel')
#else
  CALL start_clock('bse_kernel')
#endif
  !
  nbnd_do = nbndval0x-n_trunc_bands
  dffts_nnr = dffts%nnr
  band_group_nloc = band_group%nloc
  band_group_myoffset = band_group%myoffset
  !
#if !defined(__CUDA)
  ALLOCATE(raux1(dffts%nnr))
  ALLOCATE(raux2(dffts%nnr))
  ALLOCATE(caux1(npwx,nbnd_do))
  ALLOCATE(caux2(npwx,band_group%nlocx))
  IF(l_local_repr) ALLOCATE(caux3(npwx,nbnd_do))
  ALLOCATE(gaux(npwx))
#endif
  !
  DO ikq = 1, kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     npw = ngk(ikq)
     !
     IF(sf) THEN
        !
        ! spin-flip cannot be parallelized over spin
        ! ikq_do == ikq_g == local (flipped) spin index == global (flipped) spin index
        !
        ikq_do = flks(ikq)
        ikq_g = ikq_do
        !
     ELSE
        !
        ! spin-conserving can be parallelized over spin
        ! ikq_do == local spin index
        ! ikq_g == global spin index
        !
        ikq_do = ikq
        ikq_g = kpt_pool%l2g(ikq)
        !
     ENDIF
     !
     IF(l_local_repr) THEN
        !$acc host_data use_device(evc1,u_matrix,caux1)
        CALL ZGEMM('N','N',npw,nbnd_do,nbnd_do,one,evc1,npwx,u_matrix(1,1,ikq_do),nbnd_do,zero,&
        & caux1,npwx)
        !$acc end host_data
     ELSE
        !$acc kernels present(caux1,evc1)
        caux1(:,:) = evc1
        !$acc end kernels
     ENDIF
     !
     ! LOOP OVER BANDS AT KPOINT
     !
     !$acc kernels present(caux2)
     caux2(:,:) = zero
     !$acc end kernels
     !
     do_idx = n_bse_idx(ikq_do)
     !
     DO lbnd = 1, band_group%nloc, 2
        !
        my_ibnd = band_group%l2g(lbnd)
        my_jbnd = band_group%l2g(lbnd+1)
        IF(lbnd == band_group%nloc) my_jbnd = -1
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
           ibnd = idx_matrix(ipair,1,ikq_do)-n_trunc_bands
           jbnd = idx_matrix(ipair,2,ikq_do)-n_trunc_bands
           !
           IF(ibnd == my_ibnd .OR. ibnd == my_jbnd) THEN
              !
              ! read_bse_pots_g uses global spin index
              !
              CALL read_bse_pots_g(gaux,ibnd,jbnd,ikq_g)
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
        IF(lbnd < band_group%nloc) THEN
           !
           !$acc parallel loop present(raux1,raux2)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(raux1(ir),raux2(ir),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(caux2)
           CALL double_fwfft_gamma(dffts,npw,npwx,psic,caux2(:,lbnd),caux2(:,lbnd+1),'Wave')
           !$acc end host_data
           !
        ELSE
           !
           !$acc parallel loop present(raux1)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(raux1(ir),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(caux2)
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,caux2(:,lbnd),'Wave')
           !$acc end host_data
           !
        ENDIF
        !
     ENDDO
     !
     IF(l_hybrid_tddft) THEN
        !
        ! multiply the fraction of exx for hybrid tddft calculations
        !
        !$acc kernels present(caux2)
        caux2(:,:) = caux2*exxalfa
        !$acc end kernels
        !
     ENDIF
     !
     CALL gather_bands(caux2,caux1,req)
     CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
     !$acc update device(caux1)
#endif
     !
     IF(l_local_repr) THEN
        !
        !$acc kernels present(caux3,caux1)
        caux3(:,:) = caux1
        !$acc end kernels
        !
        !$acc host_data use_device(caux3,u_matrix,caux1)
        CALL ZGEMM('N','C',npw,nbnd_do,nbnd_do,one,caux3,npwx,u_matrix(1,1,ikq_do),nbnd_do,zero,&
        & caux1,npwx)
        !$acc end host_data
        !
     ENDIF
     !
     !$acc parallel loop collapse(2) present(bse_k1d,caux1)
     DO lbnd = 1, band_group_nloc
        !
        ! ibnd = band_group%l2g(lbnd)
        !
        DO ig = 1, npw
           bse_k1d(ig,lbnd) = bse_k1d(ig,lbnd)-caux1(ig,band_group_myoffset+lbnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
#if !defined(__CUDA)
  DEALLOCATE(raux1)
  DEALLOCATE(raux2)
  DEALLOCATE(caux1)
  DEALLOCATE(caux2)
  IF(l_local_repr) DEALLOCATE(caux3)
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
