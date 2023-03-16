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
MODULE west_mp
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP,i8b,sgl
  USE parallel_include
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: mp_alltoallv
  PUBLIC :: mp_circular_shift_left_begin
  PUBLIC :: mp_circular_shift_left_c16_4d
  PUBLIC :: mp_root_sum_c16_3d
  !
  INTERFACE mp_alltoallv
    MODULE PROCEDURE mp_alltoallv_i4_1d, mp_alltoallv_i8_1d, mp_alltoallv_r8_1d, &
                     mp_alltoallv_r8_2d, mp_alltoallv_c16_1d, mp_alltoallv_c16_2d
  END INTERFACE
  !
  INTERFACE mp_circular_shift_left_begin
    MODULE PROCEDURE mp_circular_shift_left_begin_c8_2d, mp_circular_shift_left_begin_c16_2d
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
    END SUBROUTINE
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
    END SUBROUTINE
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
    END SUBROUTINE
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
    END SUBROUTINE
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
    END SUBROUTINE
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
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_circular_shift_left_begin_c8_2d(send_buf,recv_buf,itag,comm,requests)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(sgl), INTENT(IN) :: send_buf(:,:)
      COMPLEX(sgl), INTENT(OUT) :: recv_buf(:,:)
      INTEGER, INTENT(IN) :: itag
      INTEGER, INTENT(IN) :: comm
      INTEGER, INTENT(OUT) :: requests(2)
      !
      ! Workspace
      !
      INTEGER :: ierr
      INTEGER :: nproc
      INTEGER :: mpime
      INTEGER :: sour
      INTEGER :: dest
      !
      CALL MPI_COMM_SIZE(comm,nproc,ierr)
      CALL MPI_COMM_RANK(comm,mpime,ierr)
      !
      sour = MOD(mpime+1,nproc)
      dest = MOD(mpime-1+nproc,nproc)
      !
      CALL MPI_IRECV(recv_buf,SIZE(recv_buf),MPI_COMPLEX,sour,itag,comm,requests(1),ierr)
      CALL MPI_ISEND(send_buf,SIZE(send_buf),MPI_COMPLEX,dest,itag,comm,requests(2),ierr)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_circular_shift_left_begin_c16_2d(send_buf,recv_buf,itag,comm,requests)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: send_buf(:,:)
      COMPLEX(DP), INTENT(OUT) :: recv_buf(:,:)
      INTEGER, INTENT(IN) :: itag
      INTEGER, INTENT(IN) :: comm
      INTEGER, INTENT(OUT) :: requests(2)
      !
      ! Workspace
      !
      INTEGER :: ierr
      INTEGER :: nproc
      INTEGER :: mpime
      INTEGER :: sour
      INTEGER :: dest
      !
      CALL MPI_COMM_SIZE(comm,nproc,ierr)
      CALL MPI_COMM_RANK(comm,mpime,ierr)
      !
      sour = MOD(mpime+1,nproc)
      dest = MOD(mpime-1+nproc,nproc)
      !
      CALL MPI_IRECV(recv_buf,SIZE(recv_buf),MPI_DOUBLE_COMPLEX,sour,itag,comm,requests(1),ierr)
      CALL MPI_ISEND(send_buf,SIZE(send_buf),MPI_DOUBLE_COMPLEX,dest,itag,comm,requests(2),ierr)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_circular_shift_left_c16_4d(buf,itag,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(INOUT) :: buf(:,:,:,:)
      INTEGER, INTENT(IN) :: itag
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      INTEGER :: nproc
      INTEGER :: mpime
      INTEGER :: sour
      INTEGER :: dest
      INTEGER :: istat(MPI_STATUS_SIZE)
      !
      CALL MPI_COMM_SIZE(comm,nproc,ierr)
      CALL MPI_COMM_RANK(comm,mpime,ierr)
      !
      sour = MOD(mpime+1,nproc)
      dest = MOD(mpime-1+nproc,nproc)
      !
      CALL MPI_SENDRECV_REPLACE(buf,SIZE(buf),MPI_DOUBLE_COMPLEX,dest,itag,sour,itag,comm,istat,ierr)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_root_sum_c16_3d(buf,root,comm)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(INOUT) :: buf(:,:,:)
      INTEGER, INTENT(IN) :: root
      INTEGER, INTENT(IN) :: comm
      !
      ! Workspace
      !
      INTEGER :: ierr
      INTEGER :: mpime
      !
      CALL MPI_COMM_RANK(comm,mpime,ierr)
      !
      IF(mpime == root) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,buf,SIZE(buf),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
      ELSE
         CALL MPI_REDUCE(buf,MPI_IN_PLACE,SIZE(buf),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
      ENDIF
      !
    END SUBROUTINE
    !
END MODULE
