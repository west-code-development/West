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
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE get_brak_hyper_parallel(dvpsi,NRHS,NLSTEPS,x,brak,idistr)
  !-----------------------------------------------------------------------
  !
  ! ... brak = < dvpsi | x >
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol
  USE class_idistribute,    ONLY : idistribute
#if defined(__CUDA)
  USE mp,                   ONLY : mp_sum,mp_waitall
  USE west_mp,              ONLY : mp_circular_shift_left_begin
  USE west_gpu,             ONLY : tmp=>tmp_r3,dvpsi_h,memcpy_H2D
#else
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute), INTENT(IN) :: idistr
  COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol,idistr%nlocx)
  INTEGER, INTENT(IN) :: NRHS,NLSTEPS
  COMPLEX(DP), INTENT(IN) :: x(npwx*npol,NRHS,NLSTEPS)
  REAL(DP), INTENT(OUT) :: brak(idistr%nglob,NLSTEPS,NRHS)
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i1_glob,il
  INTEGER :: icycl,idx,nblock_i
#if defined(__CUDA)
  INTEGER :: reqs(2)
#else
  REAL(DP), ALLOCATABLE :: tmp(:,:,:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('brak')
#else
  CALL start_clock('brak')
  !
  ALLOCATE(tmp(idistr%nlocx,NRHS,NLSTEPS))
#endif
  !
  ! -------------------
  ! BEGIN
  ! brak = < dvpsi, x >
  ! -------------------
  !
  ! Initialize to zero
  !
  !$acc kernels present(brak)
  brak(:,:,:) = 0.0_DP
  !$acc end kernels
  !
  DO icycl = 0,nimage-1
#if defined(__CUDA)
     !
     ! Cycle dvpsi (start)
     !
     IF(MOD(icycl,2) == 0) THEN
        CALL mp_circular_shift_left_begin(dvpsi,dvpsi_h,icycl,inter_image_comm,reqs)
        !
        !$acc update device(dvpsi)
     ELSE
        CALL mp_circular_shift_left_begin(dvpsi_h,dvpsi,icycl,inter_image_comm,reqs)
        !
        CALL memcpy_H2D(dvpsi,dvpsi_h,npwx*npol*idistr%nlocx)
     ENDIF
#endif
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
#if defined(__CUDA)
     !$acc host_data use_device(dvpsi,x,tmp)
     CALL glbrak_gamma_gpu(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
     !$acc end host_data
#else
     CALL glbrak_gamma(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
#endif
     !
     !$acc parallel loop collapse(3) present(brak,tmp)
     DO il = 1,NLSTEPS
        DO i2 = 1,NRHS
           DO i1 = 1,nblock_i
              !
              ! i1_glob = idistr%l2g(i1,idx)
              !
              i1_glob = nimage*(i1-1)+idx+1
              brak(i1_glob,il,i2) = tmp(i1,i2,il)
              !
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! Cycle dvpsi (end)
     !
#if defined(__CUDA)
     CALL mp_waitall(reqs)
#else
     CALL mp_circular_shift_left(dvpsi,icycl,inter_image_comm)
#endif
     !
  ENDDO
  !
  !$acc update self(brak)
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('brak')
#else
  DEALLOCATE(tmp)
  !
  CALL stop_clock('brak')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE get_brak_hyper_parallel_complex(dvpsi,NRHS,NLSTEPS,x,brak,idistr)
  !-----------------------------------------------------------------------
  !
  ! ... brak = < dvpsi | x >
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol
  USE class_idistribute,    ONLY : idistribute
#if defined(__CUDA)
  USE mp,                   ONLY : mp_sum,mp_waitall
  USE west_mp,              ONLY : mp_circular_shift_left_begin
  USE west_gpu,             ONLY : tmp=>tmp_c3,dvpsi_h,memcpy_H2D
#else
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute), INTENT(IN) :: idistr
  COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol,idistr%nlocx)
  INTEGER, INTENT(IN) :: NRHS,NLSTEPS
  COMPLEX(DP), INTENT(IN) :: x(npwx*npol,NRHS,NLSTEPS)
  COMPLEX(DP), INTENT(OUT) :: brak(idistr%nglob,NLSTEPS,NRHS)
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i1_glob,il
  INTEGER :: icycl,idx,nblock_i
#if defined(__CUDA)
  INTEGER :: reqs(2)
#else
  COMPLEX(DP), ALLOCATABLE :: tmp(:,:,:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('brak')
#else
  CALL start_clock('brak')
  !
  ALLOCATE(tmp(idistr%nlocx,NRHS,NLSTEPS))
#endif
  !
  ! -------------------
  ! BEGIN
  ! brak = < dvpsi, x >
  ! -------------------
  !
  ! Initialize to zero
  !
  !$acc kernels present(brak)
  brak(:,:,:) = 0.0_DP
  !$acc end kernels
  !
  DO icycl = 0,nimage-1
#if defined(__CUDA)
     !
     ! Cycle dvpsi (start)
     !
     IF(MOD(icycl,2) == 0) THEN
        CALL mp_circular_shift_left_begin(dvpsi,dvpsi_h,icycl,inter_image_comm,reqs)
        !
        !$acc update device(dvpsi)
     ELSE
        CALL mp_circular_shift_left_begin(dvpsi_h,dvpsi,icycl,inter_image_comm,reqs)
        !
        CALL memcpy_H2D(dvpsi,dvpsi_h,npwx*npol*idistr%nlocx)
     ENDIF
#endif
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
#if defined(__CUDA)
     !$acc host_data use_device(dvpsi,x,tmp)
     CALL glbrak_k_gpu(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
     !$acc end host_data
#else
     CALL glbrak_k(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
#endif
     !
     !$acc parallel loop collapse(3) present(brak,tmp)
     DO il = 1,NLSTEPS
        DO i2 = 1,NRHS
           DO i1 = 1,nblock_i
              !
              ! i1_glob = idistr%l2g(i1,idx)
              !
              i1_glob = nimage*(i1-1)+idx+1
              brak(i1_glob,il,i2) = tmp(i1,i2,il)
              !
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! Cycle dvpsi (end)
     !
#if defined(__CUDA)
     CALL mp_waitall(reqs)
#else
     CALL mp_circular_shift_left(dvpsi,icycl,inter_image_comm)
#endif
     !
  ENDDO
  !
  !$acc update self(brak)
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('brak')
#else
  DEALLOCATE(tmp)
  !
  CALL stop_clock('brak')
#endif
  !
END SUBROUTINE
