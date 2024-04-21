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
  USE mp,                   ONLY : mp_sum,mp_waitall
  USE west_mp,              ONLY : west_mp_circ_shift_start
#if defined(__CUDA)
  USE west_gpu,             ONLY : tmp=>tmp_r3,dvpsi2
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
  INTEGER :: reqs(2)
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: tmp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dvpsi2(:,:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('brak')
#else
  CALL start_clock('brak')
#endif
  !
#if !defined(__CUDA)
  ALLOCATE(tmp(idistr%nlocx,NRHS,NLSTEPS))
  ALLOCATE(dvpsi2(npwx*npol,idistr%nlocx))
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
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     ! Cycle dvpsi (start)
     !
     IF(MOD(icycl,2) == 0) THEN
        !
        ! Start sending dvpsi to dvpsi2 for next round
        !
        CALL west_mp_circ_shift_start(dvpsi,dvpsi2,icycl,inter_image_comm,reqs)
        !
        ! Use dvpsi this round
        !
        !$acc update device(dvpsi)
        !
        CALL glbrak_gamma(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
        !
     ELSE
        !
        ! Start sending dvpsi2 to dvpsi for next round
        !
        CALL west_mp_circ_shift_start(dvpsi2,dvpsi,icycl,inter_image_comm,reqs)
        !
        ! Use dvpsi2 this round
        !
        !$acc update device(dvpsi2)
        !
        CALL glbrak_gamma(dvpsi2,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
        !
     ENDIF
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
     CALL mp_waitall(reqs)
     !
  ENDDO
  !
  !$acc update host(brak)
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
#if !defined(__CUDA)
  DEALLOCATE(tmp)
  DEALLOCATE(dvpsi2)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('brak')
#else
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
  USE mp,                   ONLY : mp_sum,mp_waitall
  USE west_mp,              ONLY : west_mp_circ_shift_start
#if defined(__CUDA)
  USE west_gpu,             ONLY : tmp=>tmp_c3,dvpsi2
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
  INTEGER :: reqs(2)
#if !defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE :: tmp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dvpsi2(:,:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('brak')
#else
  CALL start_clock('brak')
#endif
  !
#if !defined(__CUDA)
  ALLOCATE(tmp(idistr%nlocx,NRHS,NLSTEPS))
  ALLOCATE(dvpsi2(npwx*npol,idistr%nlocx))
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
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     ! Cycle dvpsi (start)
     !
     IF(MOD(icycl,2) == 0) THEN
        !
        ! Start sending dvpsi to dvpsi2 for next round
        !
        CALL west_mp_circ_shift_start(dvpsi,dvpsi2,icycl,inter_image_comm,reqs)
        !
        ! Use dvpsi this round
        !
        !$acc update device(dvpsi)
        !
        CALL glbrak_k(dvpsi,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
        !
     ELSE
        !
        ! Start sending dvpsi2 to dvpsi for next round
        !
        CALL west_mp_circ_shift_start(dvpsi2,dvpsi,icycl,inter_image_comm,reqs)
        !
        ! Use dvpsi2 this round
        !
        !$acc update device(dvpsi2)
        !
        CALL glbrak_k(dvpsi2,x,tmp,npw,npwx,nblock_i,NLSTEPS*NRHS,idistr%nlocx,npol)
        !
     ENDIF
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
     CALL mp_waitall(reqs)
     !
  ENDDO
  !
  !$acc update host(brak)
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
#if !defined(__CUDA)
  DEALLOCATE(tmp)
  DEALLOCATE(dvpsi2)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('brak')
#else
  CALL stop_clock('brak')
#endif
  !
END SUBROUTINE
