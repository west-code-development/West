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
  USE gvect,                ONLY : gstart
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
     CALL DGEMM('T','N',nblock_i,NLSTEPS*NRHS,2*npw,2.0_DP,dvpsi,2*npwx*npol,x,2*npwx*npol,&
     & 0.0_DP,tmp,idistr%nlocx)
     IF(gstart == 2) THEN
        CALL DGER(nblock_i,NLSTEPS*NRHS,-1.0_DP,dvpsi,2*npwx*npol,x,2*npwx*npol,tmp,idistr%nlocx)
     ENDIF
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
  USE noncollin_module,     ONLY : noncolin,npol
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
  COMPLEX(DP), PARAMETER :: one = (1.0_DP,0.0_DP)
  COMPLEX(DP), PARAMETER :: zero = (0.0_DP,0.0_DP)
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
     CALL ZGEMM('C','N',nblock_i,NLSTEPS*NRHS,npw,one,dvpsi,npwx*npol,x,npwx*npol,zero,&
     & tmp,idistr%nlocx)
     IF(noncolin) THEN
        CALL ZGEMM('C','N',nblock_i,NLSTEPS*NRHS,npw,one,dvpsi(1+npwx,1),npwx*npol,&
        & x(1+npwx,1,1),npwx*npol,one,tmp,idistr%nlocx)
     ENDIF
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
