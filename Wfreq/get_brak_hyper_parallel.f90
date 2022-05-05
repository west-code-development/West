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
  COMPLEX(DP) :: dvpsi(npwx*npol,idistr%nlocx)
  INTEGER, INTENT(IN) :: NRHS,NLSTEPS
  COMPLEX(DP), INTENT(IN) :: x(npwx*npol,NRHS,NLSTEPS)
  REAL(DP), INTENT(OUT) :: brak(idistr%nglob,NLSTEPS,NRHS)
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i1_glob,il
  INTEGER :: icycl,idx,nblock_i
  !
  REAL(DP), EXTERNAL :: DDOT
  !
  CALL start_clock("brak")
  !
  ! -------------------
  ! BEGIN
  ! brak = < dvpsi, x >
  ! -------------------
  !
  ! Initialize to zero
  !
  brak = 0.0_DP
  !
  DO icycl = 0,nimage-1
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     DO i1 = 1,nblock_i
        !
        i1_glob = idistr%l2g(i1,idx)
        !
        DO il = 1,NLSTEPS
           DO i2 = 1,NRHS
              brak(i1_glob,il,i2) = 2.0_DP*DDOT(2*npw,dvpsi(1,i1),1,x(1,i2,il),1)
              IF(gstart == 2) brak(i1_glob,il,i2) = brak(i1_glob,il,i2)-REAL(dvpsi(1,i1),KIND=DP)*REAL(x(1,i2,il),KIND=DP)
           ENDDO
        ENDDO
     ENDDO
     !
     ! Cycle the dvpsi array
     !
     CALL mp_circular_shift_left(dvpsi,icycl,inter_image_comm)
     !
  ENDDO
  !
  ! Syncronize brak
  !
  CALL mp_sum(brak,intra_bgrp_comm)
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
  USE gvect,                ONLY : gstart
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : noncolin,npol
  USE class_idistribute,    ONLY : idistribute
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute), INTENT(IN) :: idistr
  COMPLEX(DP) :: dvpsi(npwx*npol,idistr%nlocx)
  INTEGER, INTENT(IN) :: NRHS,NLSTEPS
  COMPLEX(DP), INTENT(IN) :: x(npwx*npol,NRHS,NLSTEPS)
  COMPLEX(DP), INTENT(OUT) :: brak(idistr%nglob,NLSTEPS,NRHS)
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i1_glob,il
  INTEGER :: icycl,idx,nblock_i
  !
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock("brak")
  !
  ! -------------------
  ! BEGIN
  ! brak = < dvpsi, x >
  ! -------------------
  !
  ! Initialize to zero
  !
  brak = 0.0_DP
  !
  DO icycl = 0,nimage-1
     !
     idx = MOD(my_image_id+icycl,nimage)
     nblock_i = idistr%nglob/nimage
     IF(idx < MOD(idistr%nglob,nimage)) nblock_i = nblock_i+1
     !
     DO i1 = 1,nblock_i
        !
        i1_glob = idistr%l2g(i1,idx)
        !
        DO il = 1,NLSTEPS
           DO i2 = 1,NRHS
              brak(i1_glob,il,i2) = ZDOTC(npw,dvpsi(1,i1),1,x(1,i2,il),1)
           ENDDO
        ENDDO
        IF(noncolin) THEN
           DO il = 1,NLSTEPS
              DO i2 = 1,NRHS
                 brak(i1_glob,il,i2) = brak(i1_glob,il,i2)+ZDOTC(npw,dvpsi(1+npwx,i1),1,x(1+npwx,i2,il),1)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     ! Cycle the dvpsi array
     !
     CALL mp_circular_shift_left(dvpsi,icycl,inter_image_comm)
     !
  ENDDO
  !
  ! Syncronize brak
  !
  CALL mp_sum(brak,intra_bgrp_comm)
  !
  CALL stop_clock("brak")
  !
END SUBROUTINE
