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
SUBROUTINE get_brak_hyper_parallel(dvpsi,NRHS,NLSTEPS,x,brak,idistr)
  !-----------------------------------------------------------------------
  !
  ! ... brak = < dvpsi | x >
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol
  USE class_idistribute,    ONLY : idistribute
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
  REAL(DP), ALLOCATABLE :: tmp(:,:,:)
  !
  CALL start_clock("brak")
  !
  ALLOCATE(tmp(idistr%nlocx,NRHS,NLSTEPS))
  !
  ! -------------------
  ! BEGIN
  ! brak = < dvpsi, x >
  ! -------------------
  !
  ! Initialize to zero
  !
  brak = 0.0_DP
  tmp = 0.0_DP
  !
  DO icycl = 0,nimage-1
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     CALL glbrak_gamma(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
     !
     DO i1 = 1,nblock_i
        !
        i1_glob = idistr%l2g(i1,idx)
        !
        DO il = 1,NLSTEPS
           DO i2 = 1,NRHS
              brak(i1_glob,il,i2) = tmp(i1,i2,il)
           ENDDO
        ENDDO
     ENDDO
     !
     ! Cycle dvpsi
     !
     CALL mp_circular_shift_left(dvpsi,icycl,inter_image_comm)
     !
  ENDDO
  !
  ! Syncronize brak
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
  DEALLOCATE(tmp)
  !
  CALL stop_clock("brak")
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
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol
  USE class_idistribute,    ONLY : idistribute
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
  COMPLEX(DP), ALLOCATABLE :: tmp(:,:,:)
  !
  CALL start_clock("brak")
  !
  ALLOCATE(tmp(idistr%nlocx,NRHS,NLSTEPS))
  !
  ! -------------------
  ! BEGIN
  ! brak = < dvpsi, x >
  ! -------------------
  !
  ! Initialize to zero
  !
  brak = 0.0_DP
  tmp = 0.0_DP
  !
  DO icycl = 0,nimage-1
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     CALL glbrak_k(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
     !
     DO i1 = 1,nblock_i
        !
        i1_glob = idistr%l2g(i1,idx)
        !
        DO il = 1,NLSTEPS
           DO i2 = 1,NRHS
              brak(i1_glob,il,i2) = tmp(i1,i2,il)
           ENDDO
        ENDDO
     ENDDO
     !
     ! Cycle dvpsi
     !
     CALL mp_circular_shift_left(dvpsi,icycl,inter_image_comm)
     !
  ENDDO
  !
  ! Syncronize brak
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
  DEALLOCATE(tmp)
  !
  CALL stop_clock("brak")
  !
END SUBROUTINE
!
#if defined(__CUDA)
!-----------------------------------------------------------------------
SUBROUTINE get_brak_hyper_parallel_gpu(dvpsi,NRHS,NLSTEPS,x,brak,idistr)
  !-----------------------------------------------------------------------
  !
  ! ... brak = < dvpsi | x >
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_waitall
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol
  USE class_idistribute,    ONLY : idistribute
  USE west_mp,              ONLY : mp_circular_shift_left_begin
  USE west_gpu,             ONLY : tmp_r3,dvpsi_h,memcpy_H2D
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
  INTEGER :: nblock_i
  INTEGER :: icycl,idx
  INTEGER :: reqs(2)
  !
  CALL start_clock_gpu("brak")
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
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     !$acc host_data use_device(dvpsi,x,tmp_r3)
     CALL glbrak_gamma_gpu(dvpsi,x,tmp_r3,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
     !$acc end host_data
     !
     !$acc parallel loop collapse(3) present(brak,tmp_r3)
     DO i1 = 1,nblock_i
        DO il = 1,NLSTEPS
           DO i2 = 1,NRHS
              !
              ! i1_glob = idistr%l2g(i1,idx)
              !
              i1_glob = nimage*(i1-1)+idx+1
              brak(i1_glob,il,i2) = tmp_r3(i1,i2,il)
              !
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! Cycle dvpsi (end)
     !
     CALL mp_waitall(reqs)
     !
  ENDDO
  !
  !$acc update self(brak)
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
  CALL stop_clock_gpu("brak")
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE get_brak_hyper_parallel_complex_gpu(dvpsi,NRHS,NLSTEPS,x,brak,idistr)
  !-----------------------------------------------------------------------
  !
  ! ... brak = < dvpsi | x >
  !
  USE kinds,                ONLY : DP,sgl
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_waitall
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol
  USE class_idistribute,    ONLY : idistribute
  USE west_mp,              ONLY : mp_circular_shift_left_begin
  USE west_gpu,             ONLY : tmp_c3,dvpsi_h,memcpy_H2D
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
  INTEGER :: nblock_i
  INTEGER :: icycl,idx
  INTEGER :: reqs(2)
  !
  CALL start_clock_gpu("brak")
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
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     !$acc host_data use_device(dvpsi,x,tmp_c3)
     CALL glbrak_k_gpu(dvpsi,x,tmp_c3,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
     !$acc end host_data
     !
     !$acc parallel loop collapse(3) present(brak,tmp_c3)
     DO i1 = 1,nblock_i
        DO il = 1,NLSTEPS
           DO i2 = 1,NRHS
              !
              ! i1_glob = idistr%l2g(i1,idx)
              !
              i1_glob = nimage*(i1-1)+idx+1
              brak(i1_glob,il,i2) = tmp_c3(i1,i2,il)
              !
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! Cycle dvpsi (end)
     !
     CALL mp_waitall(reqs)
     !
  ENDDO
  !
  !$acc update self(brak)
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
  CALL stop_clock_gpu("brak")
  !
END SUBROUTINE
#endif
