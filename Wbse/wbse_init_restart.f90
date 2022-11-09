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
  USE iotk_module
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wbse_index_matrix_write
  PUBLIC :: wbse_index_matrix_read
  PUBLIC :: wbse_status_restart_write
  PUBLIC :: wbse_status_restart_read
  !
  INTERFACE wbse_index_matrix_write
     MODULE PROCEDURE wbse_index_matrix_write3d
  END INTERFACE
  INTERFACE wbse_index_matrix_read
     MODULE PROCEDURE wbse_index_matrix_read3d
  END INTERFACE
  INTERFACE wbse_status_restart_write
     MODULE PROCEDURE wbse_status_restart_write3d
  END INTERFACE
  INTERFACE wbse_status_restart_read
     MODULE PROCEDURE wbse_status_restart_read3d
  END INTERFACE
  !
  CONTAINS
    !
    SUBROUTINE wbse_index_matrix_write3d(filename,size_list,size_column,index_matrix)
      !
      USE mp_world,             ONLY : root,world_comm
      USE io_global,            ONLY : ionode
      USE mp,                   ONLY : mp_barrier,mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: size_list,size_column
      REAL(DP), INTENT(IN) :: index_matrix(size_list,size_column)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !
      ! Workspace
      !
      INTEGER :: ierr, iun
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! CREATE THE STATUS FILE
      !
      IF(ionode) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit(iun, ierr)
         CALL iotk_open_write(iun, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
         !
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      CALL errore('wbse_index_matrix_write', 'cannot open restart file for writing', ierr)
      !
      IF(ionode) THEN
         CALL iotk_write_begin(iun, 'WBSE_INDEX_SIZE')
         CALL iotk_write_dat(iun, 'size_list', size_list)
         CALL iotk_write_dat(iun, 'size_column', size_column)
         CALL iotk_write_end(iun, 'WBSE_INDEX_SIZE')
         !
         CALL iotk_write_begin(iun, 'WBSE_INIT_INDEX')
         CALL iotk_write_dat(iun, 'index_matrix', index_matrix(:,1:size_column))
         CALL iotk_write_end(iun, 'WBSE_INIT_INDEX')
      ENDIF
      !
      IF(ionode) CALL iotk_close_write(iun)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_index_matrix_read3d(filename,size_list0,size_list1,size_column,index_matrix)
      !
      USE io_global,            ONLY : ionode
      USE mp_world,             ONLY : world_comm,root
      USE mp,                   ONLY : mp_bcast
      USE mp_global,            ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: size_list0,size_column
      INTEGER, INTENT(OUT) :: size_list1
      REAL(DP), INTENT(OUT) :: index_matrix(size_list0,size_column)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !
      ! Workspace
      !
      INTEGER :: ierr,iun
      INTEGER :: size_list_tmp,size_column_tmp
      REAL(DP), ALLOCATABLE :: index_matrix_tmp(:,:)
      !
      ierr = 0
      IF(ionode) THEN
         CALL iotk_free_unit(iun, ierr)
         CALL iotk_open_read(iun, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      IF(ierr /=0) CALL errore('wbse_index_matrix_read', 'cannot open restart file for reading', ierr)
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iun, 'WBSE_INDEX_SIZE')
         CALL iotk_scan_dat(iun, 'size_list', size_list_tmp)
         CALL iotk_scan_dat(iun, 'size_column', size_column_tmp)
         CALL iotk_scan_end(iun, 'WBSE_INDEX_SIZE')
      ENDIF
      !
      CALL mp_bcast(size_list_tmp, 0, intra_image_comm)
      CALL mp_bcast(size_column_tmp, 0, intra_image_comm)
      !
      ALLOCATE(index_matrix_tmp(size_list_tmp, size_column_tmp))
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iun, 'WBSE_INIT_INDEX')
         CALL iotk_scan_dat(iun, 'index_matrix', index_matrix_tmp(1:size_list_tmp,1:size_column_tmp))
         CALL iotk_scan_end(iun, 'WBSE_INIT_INDEX')
         CALL iotk_close_read(iun)
      ENDIF
      !
      CALL mp_bcast(index_matrix_tmp, 0, intra_image_comm)
      !
      size_list1 = size_list_tmp
      index_matrix(:,:) = 0._DP
      index_matrix(1:size_list_tmp,1:size_column_tmp) = index_matrix_tmp(1:size_list_tmp,1:size_column_tmp)
      !
      DEALLOCATE(index_matrix_tmp)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_status_restart_write3d(filename,size_list,restart_matrix,done_calc)
      !
      USE mp_world,             ONLY : root,world_comm
      USE io_global,            ONLY : ionode
      USE mp,                   ONLY : mp_barrier,mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: size_list
      REAL(DP), INTENT(IN) :: restart_matrix(size_list)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      LOGICAL, INTENT(IN) :: done_calc
      !
      ! Workspace
      !
      INTEGER :: ierr, iun
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! CREATE THE STATUS FILE
      !
      IF(ionode) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit(iun, ierr)
         CALL iotk_open_write(iun, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
         !
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      CALL errore('wbse_status_restart_write', 'cannot open restart file for writing', ierr)
      !
      IF(ionode) THEN
         CALL iotk_write_begin(iun, 'STATUS_OF_CALCULATION')
         CALL iotk_write_dat(iun, 'status', done_calc)
         CALL iotk_write_end(iun, 'STATUS_OF_CALCULATION')
         !
         CALL iotk_write_begin(iun, 'WBSE_INIT_SIZE')
         CALL iotk_write_dat(iun, 'size_list', size_list)
         CALL iotk_write_end(iun, 'WBSE_INIT_SIZE')
         !
         CALL iotk_write_begin(iun, 'WBSE_INIT_RESTART')
         CALL iotk_write_dat(iun, 'restart_matrix', restart_matrix(:))
         CALL iotk_write_end(iun, 'WBSE_INIT_RESTART')
         !
         CALL iotk_close_write(iun)
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_status_restart_read3d(filename,size_list,restart_matrix,done_calc)
      !
      USE io_global,            ONLY : ionode
      USE mp_world,             ONLY : world_comm,root
      USE mp,                   ONLY : mp_bcast
      USE mp_global,            ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: size_list
      REAL(DP), INTENT(OUT) :: restart_matrix(size_list)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      LOGICAL, INTENT(OUT) :: done_calc
      !
      ! Workspace
      !
      INTEGER :: ierr,iun
      INTEGER :: size_list_tmp
      REAL(DP), ALLOCATABLE :: restart_matrix_tmp(:)
      !
      ierr = 0
      IF(ionode) THEN
         CALL iotk_free_unit(iun, ierr)
         CALL iotk_open_read(iun, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      IF(ierr /=0) CALL errore('wbse_status_restart_read', 'cannot open restart file for reading', ierr)
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iun, 'STATUS_OF_CALCULATION')
         CALL iotk_scan_dat(iun, 'status', done_calc)
         CALL iotk_scan_end(iun, 'STATUS_OF_CALCULATION')
         !
         CALL iotk_scan_begin(iun, 'WBSE_INIT_SIZE')
         CALL iotk_scan_dat(iun, 'size_list', size_list_tmp)
         CALL iotk_scan_end(iun, 'WBSE_INIT_SIZE')
      ENDIF
      !
      CALL mp_bcast(done_calc, 0, intra_image_comm)
      CALL mp_bcast(size_list_tmp, 0, intra_image_comm)
      !
      ALLOCATE(restart_matrix_tmp(size_list_tmp))
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iun, 'WBSE_INIT_RESTART')
         CALL iotk_scan_dat(iun, 'restart_matrix', restart_matrix_tmp(:))
         CALL iotk_scan_end(iun, 'WBSE_INIT_RESTART')
         !
         CALL iotk_close_read(iun)
      ENDIF
      !
      CALL mp_bcast(restart_matrix_tmp, 0, intra_image_comm)
      !
      restart_matrix(:) = 0._DP
      restart_matrix(1:size_list_tmp) = restart_matrix_tmp(1:size_list_tmp)
      !
      DEALLOCATE(restart_matrix_tmp)
      !
    END SUBROUTINE
    !
END MODULE
