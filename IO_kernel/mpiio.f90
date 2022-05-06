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
MODULE mod_mpiio
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_write_dmsg_at(file_name,dmsg,nlen,myoffset)
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP,i8b
      USE mp_global,   ONLY : me_bgrp,inter_image_comm
      USE parallel_include
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen
      INTEGER(i8b),INTENT(IN) :: myoffset
      REAL(DP),INTENT(IN) :: dmsg(*)
      !
      ! Workspace
      !
      INTEGER :: fh
      INTEGER :: ierr
      INTEGER :: stat(MPI_STATUS_SIZE)
      REAL(DP),PARAMETER :: a = 1._DP
      INTEGER(MPI_OFFSET_KIND) :: mp_offset
      !
      IF(me_bgrp == 0) THEN
         !
         CALL MPI_FILE_OPEN(inter_image_comm,file_name,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
         mp_offset = myoffset*SIZEOF(a)
         CALL MPI_FILE_WRITE_AT_ALL(fh,mp_offset,dmsg,nlen,MPI_DOUBLE_PRECISION,STAT,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_write_zmsg_at(file_name,zmsg,nlen,myoffset)
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP,i8b
      USE parallel_include
      USE ISO_C_BINDING
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen
      INTEGER(i8b),INTENT(IN) :: myoffset
      COMPLEX(DP),INTENT(IN),TARGET :: zmsg(*)
      !
      ! Workspace
      !
      INTEGER(MPI_OFFSET_KIND) :: mp_offset
      REAL(DP),POINTER :: p(:)
      !
      mp_offset = myoffset*2
      CALL C_F_POINTER(C_LOC(zmsg),p,[2*nlen])
      CALL mp_write_dmsg_at(file_name,p,2*nlen,mp_offset)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_read_dmsg_at(file_name,dmsg,nlen,myoffset)
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP,i8b
      USE mp_global,   ONLY : me_bgrp,intra_bgrp_comm,inter_image_comm
      USE mp,          ONLY : mp_bcast
      USE parallel_include
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen
      INTEGER(i8b),INTENT(IN) :: myoffset
      REAL(DP),INTENT(OUT) :: dmsg(*)
      !
      ! Workspace
      !
      INTEGER :: fh
      INTEGER :: ierr
      INTEGER :: stat(MPI_STATUS_SIZE)
      REAL(DP),PARAMETER :: a = 1._DP
      INTEGER(MPI_OFFSET_KIND) :: mp_offset
      !
      IF(me_bgrp == 0) THEN
         !
         CALL MPI_FILE_OPEN(inter_image_comm,file_name,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
         mp_offset = myoffset*SIZEOF(a)
         CALL MPI_FILE_READ_AT_ALL(fh,mp_offset,dmsg,nlen,MPI_DOUBLE_PRECISION,STAT,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
         !
      ENDIF
      !
      CALL mp_bcast(dmsg(1:nlen),0,intra_bgrp_comm)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_read_zmsg_at(file_name,zmsg,nlen,myoffset)
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP,i8b
      USE parallel_include
      USE ISO_C_BINDING
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen
      INTEGER(i8b),INTENT(IN) :: myoffset
      COMPLEX(DP),INTENT(OUT),TARGET :: zmsg(*)
      !
      ! Workspace
      !
      INTEGER(MPI_OFFSET_KIND) :: mp_offset
      REAL(DP),POINTER :: p(:)
      !
      mp_offset = myoffset*2
      CALL C_F_POINTER(C_LOC(zmsg),p,[2*nlen])
      CALL mp_read_dmsg_at(file_name,p,2*nlen,mp_offset)
      !
    END SUBROUTINE
END MODULE
