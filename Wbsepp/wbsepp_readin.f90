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
SUBROUTINE wbsepp_readin()
  !-----------------------------------------------------------------------
  !
  USE gvecs,            ONLY : doublegrid
  USE uspp,             ONLY : okvan
  USE noncollin_module, ONLY : domag
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: needwf
  !
  CALL start_clock('wbsepp_readin')
  !
  CALL fetch_input_yml(2,(/1,8/),.TRUE.,.FALSE.)
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  needwf = .TRUE.
  CALL read_file_new(needwf)
  !
  ! PW checks
  !
  IF(domag) CALL errore('wbse_readin','domag version not available',1)
  IF(okvan) CALL errore('wbse_readin','ultrasoft pseudopotential not implemented',1)
  IF(doublegrid) CALL errore('wbse_readin','double grid not implemented',1)
  !
  CALL stop_clock('wbsepp_readin')
  !
END SUBROUTINE
