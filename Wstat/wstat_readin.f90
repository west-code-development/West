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
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE wstat_readin()
  !-----------------------------------------------------------------------
  !
  USE gvecs,            ONLY : doublegrid
  USE uspp,             ONLY : okvan
  USE mp_global,        ONLY : nbgrp,npool
  USE xc_lib,           ONLY : xclib_dft_is
  USE pwcom,            ONLY : nkstot,lsda
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: needwf
  INTEGER :: nkpt
  !
  CALL start_clock('wstat_readin')
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
  IF(lsda) THEN
     nkpt = nkstot/2
  ELSE
     nkpt = nkstot
  ENDIF
  !
  IF(okvan) CALL errore('wstat_readin','ultrasoft pseudopotential not implemented',1)
  IF(doublegrid) CALL errore('wstat_readin','double grid not implemented',1)
  IF(nbgrp > 1 .AND. xclib_dft_is('hybrid')) CALL errore('wstat_readin','band groups not implemented for EXX',1)
  IF(npool > 1 .AND. nkpt > 1) CALL errore('wstat_readin','pools only implemented for spin, not k-points',1)
  !
  ! READ other sections of the input file
  !
  CALL fetch_input_yml(2,(/2,5/),.TRUE.)
  !
  CALL stop_clock('wstat_readin')
  !
END SUBROUTINE
