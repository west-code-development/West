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
MODULE west_mp
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP,i8b
  USE parallel_include
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: mp_alltoallv
  !
  INTERFACE mp_alltoallv
    MODULE PROCEDURE mp_alltoallv_i4_1d, mp_alltoallv_i8_1d, mp_alltoallv_r8_1d, &
                     mp_alltoallv_r8_2d, mp_alltoallv_c16_1d, mp_alltoallv_c16_2d
  END INTERFACE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_alltoallv_i4_1d(send_buf,send_count,send_displ,recv_buf,recv_count,recv_displ,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: send_buf(:)
      INTEGER, INTENT(IN) :: send_count(:)
      INTEGER, INTENT(IN) :: send_displ(:)
      INTEGER, INTENT(OUT) :: recv_buf(:)
      INTEGER, INTENT(IN) :: recv_count(:)
      INTEGER, INTENT(IN) :: recv_displ(:)
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      !
      CALL MPI_ALLTOALLV(send_buf,send_count,send_displ,MPI_INTEGER,recv_buf,recv_count,&
      & recv_displ,MPI_INTEGER,comm,ierr)
      !
    END SUBROUTINE mp_alltoallv_i4_1d
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_alltoallv_i8_1d(send_buf,send_count,send_displ,recv_buf,recv_count,recv_displ,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER(i8b), INTENT(IN) :: send_buf(:)
      INTEGER, INTENT(IN) :: send_count(:)
      INTEGER, INTENT(IN) :: send_displ(:)
      INTEGER(i8b), INTENT(OUT) :: recv_buf(:)
      INTEGER, INTENT(IN) :: recv_count(:)
      INTEGER, INTENT(IN) :: recv_displ(:)
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      !
      CALL MPI_ALLTOALLV(send_buf,send_count,send_displ,MPI_INTEGER8,recv_buf,recv_count,&
      & recv_displ,MPI_INTEGER8,comm,ierr)
      !
    END SUBROUTINE mp_alltoallv_i8_1d
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_alltoallv_r8_1d(send_buf,send_count,send_displ,recv_buf,recv_count,recv_displ,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      REAL(DP), INTENT(IN) :: send_buf(:)
      INTEGER, INTENT(IN) :: send_count(:)
      INTEGER, INTENT(IN) :: send_displ(:)
      REAL(DP), INTENT(OUT) :: recv_buf(:)
      INTEGER, INTENT(IN) :: recv_count(:)
      INTEGER, INTENT(IN) :: recv_displ(:)
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      !
      CALL MPI_ALLTOALLV(send_buf,send_count,send_displ,MPI_DOUBLE_PRECISION,recv_buf,recv_count,&
      & recv_displ,MPI_DOUBLE_PRECISION,comm,ierr)
      !
    END SUBROUTINE mp_alltoallv_r8_1d
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_alltoallv_r8_2d(send_buf,send_count,send_displ,recv_buf,recv_count,recv_displ,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      REAL(DP), INTENT(IN) :: send_buf(:,:)
      INTEGER, INTENT(IN) :: send_count(:)
      INTEGER, INTENT(IN) :: send_displ(:)
      REAL(DP), INTENT(OUT) :: recv_buf(:,:)
      INTEGER, INTENT(IN) :: recv_count(:)
      INTEGER, INTENT(IN) :: recv_displ(:)
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      !
      CALL MPI_ALLTOALLV(send_buf,send_count,send_displ,MPI_DOUBLE_PRECISION,recv_buf,recv_count,&
      & recv_displ,MPI_DOUBLE_PRECISION,comm,ierr)
      !
    END SUBROUTINE mp_alltoallv_r8_2d
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_alltoallv_c16_1d(send_buf,send_count,send_displ,recv_buf,recv_count,recv_displ,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: send_buf(:)
      INTEGER, INTENT(IN) :: send_count(:)
      INTEGER, INTENT(IN) :: send_displ(:)
      COMPLEX(DP), INTENT(OUT) :: recv_buf(:)
      INTEGER, INTENT(IN) :: recv_count(:)
      INTEGER, INTENT(IN) :: recv_displ(:)
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      !
      CALL MPI_ALLTOALLV(send_buf,send_count,send_displ,MPI_DOUBLE_COMPLEX,recv_buf,recv_count,&
      & recv_displ,MPI_DOUBLE_COMPLEX,comm,ierr)
      !
    END SUBROUTINE mp_alltoallv_c16_1d
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_alltoallv_c16_2d(send_buf,send_count,send_displ,recv_buf,recv_count,recv_displ,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: send_buf(:,:)
      INTEGER, INTENT(IN) :: send_count(:)
      INTEGER, INTENT(IN) :: send_displ(:)
      COMPLEX(DP), INTENT(OUT) :: recv_buf(:,:)
      INTEGER, INTENT(IN) :: recv_count(:)
      INTEGER, INTENT(IN) :: recv_displ(:)
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      !
      CALL MPI_ALLTOALLV(send_buf,send_count,send_displ,MPI_DOUBLE_COMPLEX,recv_buf,recv_count,&
      & recv_displ,MPI_DOUBLE_COMPLEX,comm,ierr)
      !
    END SUBROUTINE mp_alltoallv_c16_2d
    !
END MODULE
