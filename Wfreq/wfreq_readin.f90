!
! Copyright (C) 2015-2024 M. Govoni
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
SUBROUTINE wfreq_readin()
  !-----------------------------------------------------------------------
  !
  USE gvecs,            ONLY : doublegrid
  USE uspp,             ONLY : okvan
  USE mp_global,        ONLY : npool
  USE pwcom,            ONLY : nkstot,lsda
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: needwf
  INTEGER :: nkpt
  !
  CALL start_clock('wfreq_readin')
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
  ! READ other sections of the input file
  !
  CALL fetch_input_yml(2,(/2,3/),.TRUE.)
  !
  ! checks
  !
  IF(lsda) THEN
     nkpt = nkstot/2
  ELSE
     nkpt = nkstot
  ENDIF
  !
  IF(okvan) CALL errore('wfreq_readin','ultrasoft pseudopotential not implemented',1)
  IF(doublegrid) CALL errore('wfreq_readin','double grid not implemented',1)
  IF(npool > 1 .AND. nkpt > 1) &
  & CALL errore('wfreq_readin','pools only implemented for spin, not k-points',1)
  !
  CALL stop_clock('wfreq_readin')
  !
END SUBROUTINE
