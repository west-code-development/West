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
  USE io_files,  ONLY : tmp_dir
  !
  IMPLICIT NONE
  !
  PUBLIC :: wbse_init_restart_para_read
  PUBLIC :: wbse_init_restart_para_write
  !
  INTERFACE wbse_index_matrix_write
     MODULE PROCEDURE wbse_index_matrix_write3d
  END INTERFACE
  INTERFACE wbse_index_matrix_read
     MODULE PROCEDURE wbse_index_matrix_read3d
  END INTERFACE
  INTERFACE wbse_stat_restart_read
     MODULE PROCEDURE wbse_stat_restart_read3d
  END INTERFACE
  INTERFACE wbse_stat_restart_write
     MODULE PROCEDURE wbse_stat_restart_write3d
  END INTERFACE
  INTERFACE wbse_pdep_coeffie_read
     MODULE PROCEDURE wbse_pdep_coeffie_read3d
  END INTERFACE
  INTERFACE wbse_pdep_coeffie_write
     MODULE PROCEDURE wbse_pdep_coeffie_write3d
  END INTERFACE
  !
  CONTAINS
    !
    SUBROUTINE wbse_index_matrix_write3d(filename,size_list,size_column,index_matrix)
      !
      USE mp_world,             ONLY : root,world_comm
      USE io_global,            ONLY : ionode
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
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
      INTEGER :: ierr, iunout
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
         CALL iotk_free_unit(iunout, ierr)
         CALL iotk_open_write(iunout, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
         !
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      CALL errore('wbse_index_matrix_write', 'cannot open restart file for writing', ierr)
      !
      IF(ionode) THEN
         CALL iotk_write_begin(iunout, 'WBSE_INDEX_SIZE')
         CALL iotk_write_dat(iunout, 'size_list', size_list)
         CALL iotk_write_dat(iunout, 'size_column', size_column)
         CALL iotk_write_end(iunout, 'WBSE_INDEX_SIZE')
         !
         CALL iotk_write_begin(iunout, 'WBSE_INIT_INDEX')
         CALL iotk_write_dat(iunout, 'index_matrix', index_matrix(:,1:size_column))
         CALL iotk_write_end(iunout, 'WBSE_INIT_INDEX')
      ENDIF
      !
      IF(ionode) CALL iotk_close_write(iunout)
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
      INTEGER :: ierr,iunout
      INTEGER :: size_list_tmp,size_column_tmp
      REAL(DP), ALLOCATABLE :: index_matrix_tmp(:,:)
      !
      ierr = 0
      IF(ionode) THEN
         CALL iotk_free_unit(iunout, ierr)
         CALL iotk_open_read(iunout, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      IF(ierr /=0) CALL errore('wbse_index_matrix_read', 'cannot open restart file for reading', ierr)
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iunout, 'WBSE_INDEX_SIZE')
         CALL iotk_scan_dat(iunout, 'size_list', size_list_tmp)
         CALL iotk_scan_dat(iunout, 'size_column', size_column_tmp)
         CALL iotk_scan_end(iunout, 'WBSE_INDEX_SIZE')
      ENDIF
      !
      CALL mp_bcast(size_list_tmp, 0, intra_image_comm)
      CALL mp_bcast(size_column_tmp, 0, intra_image_comm)
      !
      ALLOCATE(index_matrix_tmp(size_list_tmp, size_column_tmp))
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iunout, 'WBSE_INIT_INDEX')
         CALL iotk_scan_dat(iunout, 'index_matrix', index_matrix_tmp(1:size_list_tmp,1:size_column_tmp))
         CALL iotk_scan_end(iunout, 'WBSE_INIT_INDEX')
         CALL iotk_close_read(iunout)
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
    SUBROUTINE wbse_pdep_coeffie_write3d(filename,size_list,alpha_ija_vx,alpha_ija_vc)
      !
      USE mp_world,             ONLY : root,world_comm
      USE io_global,            ONLY : ionode
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: size_list
      REAL(DP), INTENT(IN) :: alpha_ija_vx(size_list)
      REAL(DP), INTENT(IN) :: alpha_ija_vc(size_list)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !
      ! Workspace
      !
      INTEGER :: ierr, iunout
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
         CALL iotk_free_unit(iunout, ierr)
         CALL iotk_open_write(iunout, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
         !
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      CALL errore('wbse_pdep_coeffie_write', 'cannot open restart file for writing', ierr)
      !
      IF(ionode) THEN
         CALL iotk_write_begin(iunout, 'WBSE_PDEP_COEFFS_SIZE')
         CALL iotk_write_dat(iunout, 'size_list', size_list)
         CALL iotk_write_end(iunout, 'WBSE_PDEP_COEFFS_SIZE')
         !
         CALL iotk_write_begin(iunout, 'WBSE_PDEP_COEFFS_VX')
         CALL iotk_write_dat(iunout, 'coeff_matrix', alpha_ija_vx(:))
         CALL iotk_write_end(iunout, 'WBSE_PDEP_COEFFS_VX')
         !
         CALL iotk_write_begin(iunout, 'WBSE_PDEP_COEFFS_VC')
         CALL iotk_write_dat(iunout, 'coeff_matrix', alpha_ija_vc(:))
         CALL iotk_write_end(iunout, 'WBSE_PDEP_COEFFS_VC')
      ENDIF
      !
      IF(ionode) CALL iotk_close_write(iunout)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_pdep_coeffie_read3d(filename,size_list,alpha_ija_vx,alpha_ija_vc)
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
      REAL(DP), INTENT(OUT) :: alpha_ija_vx(size_list)
      REAL(DP), INTENT(OUT) :: alpha_ija_vc(size_list)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !
      ! Workspace
      !
      INTEGER :: ierr,iunout
      INTEGER :: size_list_tmp
      REAL(DP), ALLOCATABLE :: tmp_distr_vx(:)
      REAL(DP), ALLOCATABLE :: tmp_distr_vc(:)
      !
      ierr = 0
      IF(ionode) THEN
         CALL iotk_free_unit(iunout, ierr)
         CALL iotk_open_read(iunout, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      IF(ierr /=0) CALL errore('wbse_pdep_coeffie_read', 'cannot open restart file for reading', ierr)
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iunout, 'WBSE_PDEP_COEFFS_SIZE')
         CALL iotk_scan_dat(iunout, 'size_list', size_list_tmp)
         CALL iotk_scan_end(iunout, 'WBSE_PDEP_COEFFS_SIZE')
      ENDIF
      !
      CALL mp_bcast(size_list_tmp, 0, intra_image_comm)
      !
      ALLOCATE(tmp_distr_vx(size_list_tmp))
      ALLOCATE(tmp_distr_vc(size_list_tmp))
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iunout, 'WBSE_PDEP_COEFFS_VX')
         CALL iotk_scan_dat(iunout, 'coeff_matrix', tmp_distr_vx(:))
         CALL iotk_scan_end(iunout, 'WBSE_PDEP_COEFFS_VX')
      ENDIF
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iunout, 'WBSE_PDEP_COEFFS_VC')
         CALL iotk_scan_dat(iunout, 'coeff_matrix', tmp_distr_vc(:))
         CALL iotk_scan_end(iunout, 'WBSE_PDEP_COEFFS_VC')
      ENDIF
      !
      IF(ionode) CALL iotk_close_read(iunout)
      !
      CALL mp_bcast(tmp_distr_vx, 0, intra_image_comm)
      CALL mp_bcast(tmp_distr_vc, 0, intra_image_comm)
      !
      alpha_ija_vx(:) = 0._DP
      alpha_ija_vc(:) = 0._DP
      alpha_ija_vx(1:size_list_tmp) = tmp_distr_vx(1:size_list_tmp)
      alpha_ija_vc(1:size_list_tmp) = tmp_distr_vc(1:size_list_tmp)
      !
      DEALLOCATE(tmp_distr_vx)
      DEALLOCATE(tmp_distr_vc)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_stat_restart_write3d(filename,size_list,restart_matrix,done_calc)
      !
      USE mp_world,             ONLY : root,world_comm
      USE io_global,            ONLY : ionode
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
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
      INTEGER :: ierr, iunout
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
         CALL iotk_free_unit(iunout, ierr)
         CALL iotk_open_write(iunout, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
         !
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      CALL errore('wbse_stat_restart_write', 'cannot open restart file for writing', ierr)
      !
      IF(ionode) THEN
         CALL iotk_write_begin(iunout, 'STATUS_OF_CALCULATION')
         CALL iotk_write_dat(iunout, 'status', done_calc)
         CALL iotk_write_end(iunout, 'STATUS_OF_CALCULATION')
         !
         CALL iotk_write_begin(iunout, 'WBSE_INIT_SIZE')
         CALL iotk_write_dat(iunout, 'size_list', size_list)
         CALL iotk_write_end(iunout, 'WBSE_INIT_SIZE')
         !
         CALL iotk_write_begin(iunout, 'WBSE_INIT_RESTART')
         CALL iotk_write_dat(iunout, 'restart_matrix', restart_matrix(:))
         CALL iotk_write_end(iunout, 'WBSE_INIT_RESTART')
         !
         CALL iotk_close_write(iunout)
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_stat_restart_read3d(filename,size_list,restart_matrix,done_calc)
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
      INTEGER :: ierr,iunout
      INTEGER :: size_list_tmp
      REAL(DP), ALLOCATABLE :: restart_matrix_tmp(:)
      !
      ierr = 0
      IF(ionode) THEN
         CALL iotk_free_unit(iunout, ierr)
         CALL iotk_open_read(iunout, FILE=TRIM(filename), BINARY=.FALSE., IERR=ierr)
      ENDIF
      !
      CALL mp_bcast(ierr, root, world_comm)
      !
      IF(ierr /=0) CALL errore('wbse_stat_restart_read', 'cannot open restart file for reading', ierr)
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iunout, 'STATUS_OF_CALCULATION')
         CALL iotk_scan_dat(iunout, 'status', done_calc)
         CALL iotk_scan_end(iunout, 'STATUS_OF_CALCULATION')
         !
         CALL iotk_scan_begin(iunout, 'WBSE_INIT_SIZE')
         CALL iotk_scan_dat(iunout, 'size_list', size_list_tmp)
         CALL iotk_scan_end(iunout, 'WBSE_INIT_SIZE')
      ENDIF
      !
      CALL mp_bcast(done_calc, 0, intra_image_comm)
      CALL mp_bcast(size_list_tmp, 0, intra_image_comm)
      !
      ALLOCATE(restart_matrix_tmp(size_list_tmp))
      !
      IF(ionode) THEN
         CALL iotk_scan_begin(iunout, 'WBSE_INIT_RESTART')
         CALL iotk_scan_dat(iunout, 'restart_matrix', restart_matrix_tmp(:))
         CALL iotk_scan_end(iunout, 'WBSE_INIT_RESTART')
         !
         CALL iotk_close_read(iunout)
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
    SUBROUTINE wbse_init_restart_para_write(size_list,alpha_ija_vx, alpha_ija_vc)
      !
      USE mp_world,             ONLY : world_comm
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
      USE mp_global,            ONLY : my_pool_id,my_bgrp_id,me_bgrp,root_bgrp
      USE distribution_center,  ONLY : bseparal
      USE westcom,              ONLY : wbse_init_save_dir
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: size_list
      REAL(DP), INTENT(IN) :: alpha_ija_vx(size_list)
      REAL(DP), INTENT(IN) :: alpha_ija_vc(size_list)
      !
      ! Workspace
      !
      INTEGER :: ierr, iunout, image_id
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=6) :: my_label
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      image_id = bseparal%mylevelid
      WRITE(my_label,'(i6.6)') image_id
      filename = TRIM(wbse_init_save_dir)//'aux_imageid_'//my_label//'.dat'
      !
      IF(my_pool_id /= 0) RETURN
      IF(my_bgrp_id /= 0) RETURN
      !
      ! Resume all components
      !
      IF(me_bgrp == root_bgrp) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit(iunout, ierr)
         CALL iotk_open_write(iunout, FILE=TRIM(filename), BINARY=.FALSE.)
         !
         CALL iotk_write_begin(iunout, 'WBSE_PDEP_COEFFS_SIZE')
         CALL iotk_write_dat(iunout, 'size_list', size_list)
         CALL iotk_write_end(iunout, 'WBSE_PDEP_COEFFS_SIZE')
         !
         CALL iotk_write_begin(iunout, 'WBSE_PDEP_COEFFS_VX')
         CALL iotk_write_dat(iunout, 'coeff_matrix', alpha_ija_vx(:))
         CALL iotk_write_end(iunout, 'WBSE_PDEP_COEFFS_VX')
         !
         CALL iotk_write_begin(iunout, 'WBSE_PDEP_COEFFS_VC')
         CALL iotk_write_dat(iunout, 'coeff_matrix', alpha_ija_vc(:))
         CALL iotk_write_end(iunout, 'WBSE_PDEP_COEFFS_VC')
         !
         CALL iotk_close_write(iunout)
         !
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE wbse_init_restart_para_read(size_list,alpha_ija_vx,alpha_ija_vc)
      !
      USE mp_world,             ONLY : world_comm
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_global,            ONLY : my_pool_id,my_bgrp_id,me_bgrp,root_bgrp,inter_bgrp_comm,&
                                     & inter_pool_comm
      USE distribution_center,  ONLY : bseparal
      USE westcom,              ONLY : wbse_init_save_dir
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: size_list
      REAL(DP), INTENT(OUT) :: alpha_ija_vx(size_list)
      REAL(DP), INTENT(OUT) :: alpha_ija_vc(size_list)
      !
      ! Workspace
      !
      INTEGER :: ierr,iunout,image_id
      INTEGER :: size_list_tmp
      REAL(DP), ALLOCATABLE :: alpha_ija_vx_tmp(:)
      REAL(DP), ALLOCATABLE :: alpha_ija_vc_tmp(:)
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=6) :: my_label
      !
      CALL mp_barrier(world_comm)
      !
      image_id = bseparal%mylevelid
      WRITE(my_label,'(i6.6)') image_id
      filename = TRIM(wbse_init_save_dir)//'aux_imageid_'//my_label//'.dat'
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp == root_bgrp) THEN
            !
            ! ... open XML descriptor
            !
            CALL iotk_free_unit(iunout, ierr)
            CALL iotk_open_read(iunout, FILE=TRIM(filename), BINARY=.TRUE., IERR=ierr)
            !
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(ierr,0,inter_bgrp_comm)
      !
      CALL errore('wbse_init_restart_para_read ', 'cannot open restart file for reading', ierr)
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp == root_bgrp) THEN
            CALL iotk_scan_begin(iunout, 'WBSE_PDEP_COEFFS_SIZE')
            CALL iotk_scan_dat(iunout, 'size_list', size_list_tmp)
            CALL iotk_scan_end(iunout, 'WBSE_PDEP_COEFFS_SIZE')
         ENDIF
         !
      ENDIF
      CALL mp_bcast(size_list_tmp,0,inter_bgrp_comm)
      CALL mp_bcast(size_list_tmp,0,inter_pool_comm)
      !
      ALLOCATE(alpha_ija_vx_tmp(size_list_tmp))
      ALLOCATE(alpha_ija_vc_tmp(size_list_tmp))
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp == root_bgrp) THEN
            CALL iotk_scan_begin(iunout, 'WBSE_PDEP_COEFFS_VX')
            CALL iotk_scan_dat(iunout, 'coeff_matrix', alpha_ija_vx_tmp(:))
            CALL iotk_scan_end(iunout, 'WBSE_PDEP_COEFFS_VX')
            !
            CALL iotk_scan_begin(iunout, 'WBSE_PDEP_COEFFS_VC')
            CALL iotk_scan_dat(iunout, 'coeff_matrix', alpha_ija_vc_tmp(:))
            CALL iotk_scan_end(iunout, 'WBSE_PDEP_COEFFS_VC')
            !
            CALL iotk_close_read(iunout)
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(alpha_ija_vx_tmp,0,inter_bgrp_comm)
      CALL mp_bcast(alpha_ija_vx_tmp,0,inter_pool_comm)
      CALL mp_bcast(alpha_ija_vc_tmp,0,inter_bgrp_comm)
      CALL mp_bcast(alpha_ija_vc_tmp,0,inter_pool_comm)
      !
      alpha_ija_vx = 0._DP
      alpha_ija_vc = 0._DP
      alpha_ija_vx(1:size_list_tmp) = alpha_ija_vx_tmp(1:size_list_tmp)
      alpha_ija_vc(1:size_list_tmp) = alpha_ija_vc_tmp(1:size_list_tmp)
      !
      DEALLOCATE(alpha_ija_vx_tmp, alpha_ija_vc_tmp)
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
END MODULE
