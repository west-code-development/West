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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_init_readin()
  !-----------------------------------------------------------------------
  !
  USE gvecs,            ONLY : doublegrid
  USE uspp,             ONLY : okvan
  USE noncollin_module, ONLY : domag
  USE mp_global,        ONLY : npool
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: needwf
  !
  CALL start_clock('wbse_init_readin')
  !
  ! READ INPUT_WEST
  !
  CALL fetch_input_yml(1,(/1/),.TRUE.)
  !
  ! read the input file produced by the pwscf program
  ! allocate memory and recalculate what is needed
  !
  needwf = .TRUE.
  CALL read_file_new(needwf)
  !
  ! PW checks
  !
  IF(domag) CALL errore('wbse_init_readin','domag version not available',1)
  IF(okvan) CALL errore('wbse_init_readin','ultrasoft pseudopotential not implemented',1)
  IF(doublegrid) CALL errore('wbse_init_readin','double grid not implemented',1)
  IF(npool > 1) CALL errore('wbse_init_readin','pools not implemented',1)
  !
  ! READ other sections of the input file
  !
  CALL fetch_input_yml(2,(/5,6/),.TRUE.)
  !
  CALL stop_clock('wbse_init_readin')
  !
END SUBROUTINE
