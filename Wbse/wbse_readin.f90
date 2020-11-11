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
SUBROUTINE wbse_readin()
  !-----------------------------------------------------------------------
  !
  USE pwcom
  USE westcom
  USE ions_base,        ONLY : nat
  USE uspp,             ONLY : okvan
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : nproc,mpime,root
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: iunit =5, ios
  !
  CALL start_clock('wbse_readin')
  !
  !TODO: change to new vest fetch namelist
  CALL fetch_input_yml(3,(/1,6,7/),.TRUE.,.FALSE.)
  !CALL wbse_fetch_namelist(3,(/1,2,3/),.TRUE.)
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  CALL read_pwout( )
  !
  ! PW checks
  !
  IF (domag) CALL errore('wbse_readin','domag version not available',1)
  IF (okvan) CALL errore('wbse_readin','ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) CALL errore('wbse_readin', 'double grid not implemented',1)
  !
  CALL stop_clock('wbse_readin')
  !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_init_readin()
  !-----------------------------------------------------------------------
  !
  USE pwcom
  USE westcom
  USE ions_base,        ONLY : nat
  USE uspp,             ONLY : okvan
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : nproc,mpime,root
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: iunit =5, ios
  !
  CALL start_clock('wbse_init_readin')
  !
  !TODO: change to new west version of fech namelist
  CALL fetch_input_yml(3,(/1,6,8/),.TRUE.,.FALSE.)
  !CALL wbse_fetch_namelist(3,(/1,2,4/))
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  CALL read_pwout( )
  !
  ! PW checks
  !
  IF (domag) CALL errore('wbse_readin','domag version not available',1)
  IF (okvan) CALL errore('wbse_readin','ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) CALL errore('wbse_readin', 'double grid not implemented',1)
  !
  CALL stop_clock('wbse_init_readin')
  !
END SUBROUTINE
