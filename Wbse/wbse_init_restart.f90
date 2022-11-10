!
! Copyright (C) 2015-2022 M. Govoni
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
MODULE wbse_init_restart
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    SUBROUTINE wbse_index_matrix_write(fname,size_list,size_column,index_matrix)
      !
      USE kinds,               ONLY : i8b
      USE mp_world,            ONLY : world_comm
      USE io_global,           ONLY : ionode
      USE mp,                  ONLY : mp_barrier
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,&
                                    & HD_ID_DIMENSION
      USE base64_module,       ONLY : islittleendian
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: size_list,size_column
      INTEGER,INTENT(IN) :: index_matrix(size_list,size_column)
      CHARACTER(LEN=*),INTENT(IN) :: fname
      !
      ! Workspace
      !
      INTEGER :: iun
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! CREATE THE STATUS FILE
      !
      IF(ionode) THEN
         !
         header = 0
         header(HD_ID_VERSION) = HD_VERSION
         IF(islittleendian()) THEN
            header(HD_ID_LITTLE_ENDIAN) = 1
         ENDIF
         header(HD_ID_DIMENSION) = size_list
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
         offset = 1
         WRITE(iun,POS=offset) header
         offset = offset+HD_LENGTH*SIZEOF(header(1))
         WRITE(iun,POS=offset) index_matrix(1:size_list,1:size_column)
         CLOSE(iun)
         !
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_index_matrix_read(fname,size_list0,size_list1,size_column,index_matrix)
      !
      USE kinds,               ONLY : i8b
      USE io_global,           ONLY : ionode
      USE mp,                  ONLY : mp_bcast
      USE mp_global,           ONLY : intra_image_comm
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,&
                                    & HD_ID_DIMENSION
      USE base64_module,       ONLY : islittleendian
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: size_list0,size_column
      INTEGER,INTENT(OUT) :: size_list1
      INTEGER,INTENT(OUT) :: index_matrix(size_list0,size_column)
      CHARACTER(LEN=*),INTENT(IN) :: fname
      !
      ! Workspace
      !
      INTEGER :: ierr,iun
      INTEGER,ALLOCATABLE :: index_matrix_tmp(:,:)
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      IF(ionode) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
         IF(ierr /= 0) THEN
            CALL errore('wbse_index_matrix_read','Cannot read file:'//TRIM(fname),1)
         ENDIF
         !
         offset = 1
         READ(iun,POS=offset) header
         IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
            CALL errore('wbse_index_matrix_read','Unknown file format:'//TRIM(fname),1)
         ENDIF
         IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
            .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
            CALL errore('wbse_index_matrix_read','Endianness mismatch:'//TRIM(fname),1)
         ENDIF
         !
         size_list1 = header(HD_ID_DIMENSION)
         !
      ENDIF
      !
      CALL mp_bcast(size_list1,0,intra_image_comm)
      !
      ALLOCATE(index_matrix_tmp(size_list1,size_column))
      !
      IF(ionode) THEN
         !
         READ(iun,POS=offset) index_matrix_tmp(1:size_list1,1:size_column)
         CLOSE(iun)
         !
      ENDIF
      !
      CALL mp_bcast(index_matrix_tmp,0,intra_image_comm)
      !
      index_matrix(:,:) = 0
      index_matrix(1:size_list1,1:size_column) = index_matrix_tmp(1:size_list1,1:size_column)
      !
      DEALLOCATE(index_matrix_tmp)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_status_restart_write(fname,size_list,restart_matrix)
      !
      USE kinds,               ONLY : i8b
      USE mp_world,            ONLY : world_comm
      USE io_global,           ONLY : ionode
      USE mp,                  ONLY : mp_barrier
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,&
                                    & HD_ID_DIMENSION
      USE base64_module,       ONLY : islittleendian
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: size_list
      INTEGER,INTENT(IN) :: restart_matrix(size_list)
      CHARACTER(LEN=*),INTENT(IN) :: fname
      !
      ! Workspace
      !
      INTEGER :: iun
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! CREATE THE STATUS FILE
      !
      IF(ionode) THEN
         !
         header = 0
         header(HD_ID_VERSION) = HD_VERSION
         IF(islittleendian()) THEN
            header(HD_ID_LITTLE_ENDIAN) = 1
         ENDIF
         header(HD_ID_DIMENSION) = size_list
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
         offset = 1
         WRITE(iun,POS=offset) header
         offset = offset+HD_LENGTH*SIZEOF(header(1))
         WRITE(iun,POS=offset) restart_matrix(1:size_list)
         CLOSE(iun)
         !
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_status_restart_read(fname,size_list,restart_matrix,done_calc)
      !
      USE kinds,               ONLY : i8b
      USE io_global,           ONLY : ionode
      USE mp,                  ONLY : mp_bcast
      USE mp_global,           ONLY : intra_image_comm
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,&
                                    & HD_ID_DIMENSION
      USE base64_module,       ONLY : islittleendian
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: size_list
      INTEGER,INTENT(OUT) :: restart_matrix(size_list)
      CHARACTER(LEN=*),INTENT(IN) :: fname
      LOGICAL,INTENT(OUT) :: done_calc
      !
      ! Workspace
      !
      INTEGER :: ierr,iun
      INTEGER :: size_list_tmp
      INTEGER,ALLOCATABLE :: restart_matrix_tmp(:)
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      IF(ionode) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
         IF(ierr /= 0) THEN
            CALL errore('wbse_status_restart_read','Cannot read file:'//TRIM(fname),1)
         ENDIF
         !
         offset = 1
         READ(iun,POS=offset) header
         IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
            CALL errore('wbse_status_restart_read','Unknown file format:'//TRIM(fname),1)
         ENDIF
         IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
            .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
            CALL errore('wbse_status_restart_read','Endianness mismatch:'//TRIM(fname),1)
         ENDIF
         !
         size_list_tmp = header(HD_ID_DIMENSION)
         !
      ENDIF
      !
      CALL mp_bcast(size_list_tmp,0,intra_image_comm)
      !
      ALLOCATE(restart_matrix_tmp(size_list_tmp))
      !
      IF(ionode) THEN
         !
         READ(iun,POS=offset) restart_matrix_tmp(1:size_list_tmp)
         CLOSE(iun)
         !
      ENDIF
      !
      CALL mp_bcast(restart_matrix_tmp,0,intra_image_comm)
      !
      restart_matrix(:) = 0
      restart_matrix(1:size_list_tmp) = restart_matrix_tmp(1:size_list_tmp)
      !
      DEALLOCATE(restart_matrix_tmp)
      !
      IF(ANY(restart_matrix(1:size_list) == 0)) THEN
         done_calc = .FALSE.
      ELSE
         done_calc = .TRUE.
      ENDIF
      !
    END SUBROUTINE
    !
END MODULE
