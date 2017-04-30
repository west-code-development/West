!
! Copyright (C) 2015-2016 M. Govoni 
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
  SAVE
  !
  INTEGER :: io_comm
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE set_io_comm( )
      !-----------------------------------------------------------------------
      !
      USE mp,        ONLY : mp_sum,mp_barrier
      USE mp_world,  ONLY : mpime,world_comm
      USE mp_images, ONLY : my_image_id,nimage
      USE mp_bands,  ONLY : me_bgrp
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      INTEGER :: old_group, new_group 
      INTEGER :: ierr
      INTEGER,ALLOCATABLE :: will_use(:)
      !
      ALLOCATE( will_use(0:nimage-1) )
      !
      will_use = 0 
      !
      IF(me_bgrp==0) will_use(my_image_id) = mpime 
      !
      CALL mp_sum(will_use,world_comm)
      !
      CALL MPI_COMM_GROUP(world_comm,old_group,ierr)
      CALL MPI_GROUP_INCL(old_group,nimage,will_use,new_group,ierr)
      CALL MPI_COMM_CREATE(world_comm,new_group,io_comm,ierr)
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_master_creates_and_preallocates( file_name, ntot )  
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP
      USE mp_world,    ONLY : mpime,world_comm
      USE mp_images,   ONLY : nimage
      USE mp,          ONLY : mp_barrier
      USE parallel_include
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: ntot
      !
      ! Workspace
      !
      INTEGER(KIND=MPI_OFFSET_KIND) :: file_size
      INTEGER :: fh 
      INTEGER :: ierr
      REAL(DP),PARAMETER :: a=1._DP
      !
      IF( mpime==0 ) THEN
         !
         file_size = SIZEOF( a ) * ntot
         CALL MPI_FILE_OPEN(MPI_COMM_SELF, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr) 
         CALL MPI_FILE_SET_SIZE (fh, file_size, ierr)
         CALL MPI_FILE_CLOSE( fh, ierr)
         !
      ENDIF
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_write_dmsg_at( file_name, dmsg, nlen, myoffset )  
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP
      USE mp_global,   ONLY : me_bgrp,world_comm
      USE mp_images,   ONLY : my_image_id
      USE mp,          ONLY : mp_barrier
      USE parallel_include
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen, myoffset
      REAL(DP),INTENT(IN) :: dmsg(*)
      !
      ! Workspace
      !
      INTEGER :: fh 
      INTEGER :: ierr
      INTEGER :: stat(MPI_STATUS_SIZE)
      REAL(DP),PARAMETER :: a=1._DP
      INTEGER(KIND=MPI_OFFSET_KIND) :: mp_offset 
      !
      IF(me_bgrp==0) THEN
         !
         CALL MPI_FILE_OPEN(io_comm, file_name, MPI_MODE_WRONLY , MPI_INFO_NULL, fh, ierr)
         mp_offset = SIZEOF(a)*myoffset
         CALL MPI_FILE_WRITE_AT(fh, mp_offset, dmsg, nlen, MPI_DOUBLE_PRECISION, STAT, ierr)
         CALL MPI_FILE_CLOSE( fh, ierr ) 
         !
      ENDIF
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_write_zmsg_at( file_name, zmsg, nlen, myoffset )  
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP
      USE mp_global,   ONLY : me_bgrp,world_comm
      USE mp_images,   ONLY : my_image_id
      USE mp,          ONLY : mp_barrier
      USE parallel_include
      USE ISO_C_BINDING 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen, myoffset
      COMPLEX(DP),INTENT(IN),TARGET :: zmsg(*)
      !
      ! Workspace
      !
      REAL(DP),POINTER :: p(:)
      !
      CALL C_F_POINTER(C_LOC(zmsg),p,[2*nlen])
      CALL mp_write_dmsg_at( file_name, p, 2*nlen, 2*myoffset)
      !
    END SUBROUTINE
    !
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_read_dmsg_at( file_name, dmsg, nlen, myoffset )  
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP
      USE mp_global,   ONLY : me_bgrp,world_comm
      USE mp_world,    ONLY : mpime
      USE mp_images,   ONLY : my_image_id
      USE mp_bands,    ONLY : intra_bgrp_comm
      USE mp,          ONLY : mp_barrier,mp_bcast
      USE parallel_include
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen, myoffset
      REAL(DP),INTENT(OUT) :: dmsg(*)
      !
      ! Workspace
      !
      INTEGER :: fh 
      INTEGER :: ierr
      INTEGER :: stat(MPI_STATUS_SIZE)
      REAL(DP),PARAMETER :: a=1._DP
      INTEGER(KIND=MPI_OFFSET_KIND) :: mp_offset 
      INTEGER :: mys,myr
      !
      IF(me_bgrp==0) THEN
         CALL MPI_COMM_SIZE(io_comm, mys, ierr)
         CALL MPI_COMM_RANK(io_comm, myr, ierr)
         !
         CALL MPI_FILE_OPEN(io_comm, file_name, MPI_MODE_RDONLY , MPI_INFO_NULL, fh, ierr)
         mp_offset = SIZEOF(a)*myoffset
         CALL MPI_FILE_READ_AT(fh, mp_offset, dmsg, nlen, MPI_DOUBLE_PRECISION, STAT, ierr)
         CALL MPI_FILE_CLOSE( fh, ierr ) 
         !
      ENDIF
      !
      CALL mp_bcast( dmsg(1:nlen), 0, intra_bgrp_comm )
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mp_read_zmsg_at( file_name, zmsg, nlen, myoffset )  
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP
      USE mp_global,   ONLY : me_bgrp,world_comm
      USE mp_world,    ONLY : mpime
      USE mp_images,   ONLY : my_image_id
      USE mp_bands,    ONLY : intra_bgrp_comm
      USE mp,          ONLY : mp_barrier,mp_bcast
      USE parallel_include
      USE ISO_C_BINDING
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      CHARACTER(*),INTENT(IN) :: file_name
      INTEGER,INTENT(IN) :: nlen, myoffset
      COMPLEX(DP),INTENT(OUT),TARGET :: zmsg(*)
      !
      ! Workspace
      !
      REAL(DP),POINTER :: p(:)
      !
      CALL C_F_POINTER(C_LOC(zmsg),p,[2*nlen])
      CALL mp_read_dmsg_at( file_name, p, 2*nlen, 2*myoffset)
      !
    END SUBROUTINE
!
!
!
!SUBROUTINE mp_write_d2_at( file_name, io_comm, msg, n1, nloc, myoffset )  
!   !
!   USE kinds,       ONLY : DP
!   USE mp_global,   ONLY : me_bgrp,world_comm
!   USE mp_images,   ONLY : my_image_id
!   USE mp,          ONLY : mp_barrier
!   USE parallel_include
!   !
!   IMPLICIT NONE
!   !
!   CHARACTER(*) :: file_name
!   INTEGER :: io_comm, n1, nloc, myoffset
!   REAL(DP) :: msg(n1,nloc)
!   !
!   INTEGER :: fh 
!   INTEGER :: ierr
!   INTEGER :: stat(MPI_STATUS_SIZE)
!   REAL(DP),PARAMETER :: a=1._DP
!   INTEGER(KIND=MPI_OFFSET_KIND) :: mp_offset 
!   !
!   IF(me_bgrp==0) THEN
!      !
!      CALL MPI_FILE_OPEN(io_comm, file_name, MPI_MODE_WRONLY , MPI_INFO_NULL, fh, ierr)
!      mp_offset = SIZEOF(a)*n1*myoffset
!      CALL MPI_FILE_WRITE_AT(fh, mp_offset, msg, n1*nloc, MPI_DOUBLE_PRECISION, STAT, ierr)
!      CALL MPI_FILE_CLOSE( fh, ierr ) 
!      !
!   ENDIF
!   !
!   CALL mp_barrier(world_comm)
!   !
!END SUBROUTINE
!
!
!
!SUBROUTINE mp_write_d3_at( file_name, io_comm, msg, n1, n2, nloc, myoffset )  
!   !
!   USE kinds,       ONLY : DP
!   USE mp_global,   ONLY : me_bgrp,world_comm
!   USE mp_images,   ONLY : my_image_id
!   USE mp,          ONLY : mp_barrier
!   USE parallel_include
!   !
!   IMPLICIT NONE
!   !
!   CHARACTER(*) :: file_name
!   INTEGER :: io_comm, n1, n2, nloc, myoffset
!   REAL(DP) :: msg(n1,nloc)
!   !
!   INTEGER :: fh 
!   INTEGER :: ierr
!   INTEGER :: stat(MPI_STATUS_SIZE)
!   REAL(DP),PARAMETER :: a=1._DP
!   INTEGER(KIND=MPI_OFFSET_KIND) :: mp_offset 
!   !
!   IF(me_bgrp==0) THEN
!      !
!      CALL MPI_FILE_OPEN(io_comm, file_name, MPI_MODE_WRONLY , MPI_INFO_NULL, fh, ierr)
!      mp_offset = SIZEOF(a)*n1*n2*myoffset
!      CALL MPI_FILE_WRITE_AT(fh, mp_offset, msg, n1*n2*nloc, MPI_DOUBLE_PRECISION, STAT, ierr)
!      CALL MPI_FILE_CLOSE( fh, ierr ) 
!      !
!   ENDIF
!   !
!   CALL mp_barrier(world_comm)
!   !
!END SUBROUTINE
!
!
!
!
!SUBROUTINE mp_read_d2_at( file_name, io_comm, msg, n1, nloc, myoffset )  
!   !
!   USE kinds,       ONLY : DP
!   USE mp_global,   ONLY : me_bgrp,world_comm
!   USE mp_world,    ONLY : mpime
!   USE mp_images,   ONLY : my_image_id
!   USE mp_bands,    ONLY : intra_bgrp_comm
!   USE mp,          ONLY : mp_barrier,mp_bcast
!   USE parallel_include
!   !
!   IMPLICIT NONE
!   !
!   CHARACTER(*) :: file_name
!   INTEGER :: io_comm, n1, nloc, myoffset
!   REAL(DP) :: msg(n1,nloc)
!   !
!   INTEGER :: fh 
!   INTEGER :: ierr
!   INTEGER :: stat(MPI_STATUS_SIZE)
!   REAL(DP),PARAMETER :: a=1._DP
!   INTEGER(KIND=MPI_OFFSET_KIND) :: mp_offset 
!   INTEGER :: mys,myr
!   !
!   IF(me_bgrp==0) THEN
!      CALL MPI_COMM_SIZE(io_comm, mys, ierr)
!      CALL MPI_COMM_RANK(io_comm, myr, ierr)
!      !
!      CALL MPI_FILE_OPEN(io_comm, file_name, MPI_MODE_RDONLY , MPI_INFO_NULL, fh, ierr)
!      mp_offset = SIZEOF(a)*n1*myoffset
!      CALL MPI_FILE_READ_AT(fh, mp_offset, msg, n1*nloc, MPI_DOUBLE_PRECISION, STAT, ierr)
!      CALL MPI_FILE_CLOSE( fh, ierr ) 
!      !
!   ENDIF
!   !
!   CALL mp_bcast( msg, 0, intra_bgrp_comm )
!   CALL mp_barrier(world_comm)
!   !
!END SUBROUTINE
!
!
!
!
!SUBROUTINE mp_read_d3_at( file_name, io_comm, msg, n1, n2, nloc, myoffset )  
!   !
!   USE kinds,       ONLY : DP
!   USE mp_global,   ONLY : me_bgrp,world_comm
!   USE mp_world,    ONLY : mpime
!   USE mp_images,   ONLY : my_image_id
!   USE mp_bands,    ONLY : intra_bgrp_comm
!   USE mp,          ONLY : mp_barrier,mp_bcast
!   USE parallel_include
!   !
!   IMPLICIT NONE
!   !
!   CHARACTER(*) :: file_name
!   INTEGER :: io_comm, n1, n2, nloc, myoffset
!   REAL(DP) :: msg(n1,n2,nloc)
!   !
!   INTEGER :: fh 
!   INTEGER :: ierr
!   INTEGER :: stat(MPI_STATUS_SIZE)
!   REAL(DP),PARAMETER :: a=1._DP
!   INTEGER(KIND=MPI_OFFSET_KIND) :: mp_offset 
!   !
!   IF(me_bgrp==0) THEN
!      !
!      CALL MPI_FILE_OPEN(io_comm, file_name, MPI_MODE_RDONLY , MPI_INFO_NULL, fh, ierr)
!      mp_offset = SIZEOF(a)*n1*n2*myoffset
!      CALL MPI_FILE_READ_AT(fh, mp_offset, msg, n1*n2*nloc, MPI_DOUBLE_PRECISION, STAT, ierr)
!      CALL MPI_FILE_CLOSE( fh, ierr ) 
!      !
!   ENDIF
!   !
!   CALL mp_bcast( msg, 0, intra_bgrp_comm )
!   CALL mp_barrier(world_comm)
!   !
!END SUBROUTINE
   !
END MODULE
