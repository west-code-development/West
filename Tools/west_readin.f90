!
! Copyright (C) 2015-2025 M. Govoni
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
SUBROUTINE west_readin(code)
  !-----------------------------------------------------------------------
  !
  USE gvecs,            ONLY : doublegrid
  USE uspp,             ONLY : okvan
  USE mp_global,        ONLY : npool,nbgrp
  USE pwcom,            ONLY : nkstot,lsda
  USE symm_base,        ONLY : nosym
  USE control_flags,    ONLY : noinv
  USE westcom,          ONLY : l_spin_flip
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  CHARACTER(*),INTENT(IN) :: code
  !
  ! Workspace
  !
  LOGICAL :: needwf
  INTEGER :: nkpt
  !
  CALL start_clock('west_readin')
  !
  ! Read west_control
  !
  CALL fetch_input_yml(1,(/1/),.TRUE.)
  !
  ! Read pwscf
  !
  needwf = .TRUE.
  CALL read_file_new(needwf)
  !
  ! Read other input sections
  !
  SELECT CASE(TRIM(code))
  CASE('WSTAT')
     CALL fetch_input_yml(2,(/2,5/),.TRUE.)
  CASE('WFREQ')
     CALL fetch_input_yml(2,(/2,3/),.TRUE.)
  CASE('WESTPP')
     CALL fetch_input_yml(2,(/2,4/),.TRUE.)
  CASE('WBSE_INIT')
     CALL fetch_input_yml(2,(/5,6/),.TRUE.)
  CASE('WBSE')
     CALL fetch_input_yml(2,(/6,7/),.TRUE.)
  CASE DEFAULT
     CALL errore('west_readin','unknown code',1)
  END SELECT
  !
  ! General checks
  !
  IF(lsda) THEN
     nkpt = nkstot/2
  ELSE
     nkpt = nkstot
  ENDIF
  !
  IF(okvan) CALL errore('west_readin','ultrasoft pseudopotential not implemented',1)
  IF(doublegrid) CALL errore('west_readin','double grid not implemented',1)
  IF(nkpt > 1) THEN
     IF(npool > 1) CALL errore('west_readin','pools only implemented for spin, not k-points',1)
     IF(.NOT. nosym) CALL errore('west_readin','pwscf flags nosym, noinv required for k-points',1)
     IF(.NOT. noinv) CALL errore('west_readin','pwscf flags nosym, noinv required for k-points',1)
  ENDIF
  !
  ! Code specific checks
  !
  SELECT CASE(TRIM(code))
  CASE('WESTPP')
     IF(nbgrp > 1) CALL errore('west_readin','band groups not implemented for westpp',1)
     IF(npool > 1) CALL errore('west_readin','pools not implemented for westpp',1)
  CASE('WBSE_INIT','WBSE')
     IF(npool > 1 .AND. l_spin_flip) CALL errore('west_readin','pools not implemented for spin flip',1)
  END SELECT
  !
  CALL stop_clock('west_readin')
  !
END SUBROUTINE
