!
! Copyright (C) 2015-2021 M. Govoni
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
  USE uspp,             ONLY : okvan
  USE gvecs,            ONLY : doublegrid
  USE spin_orb,         ONLY : domag
  USE mp_bands,         ONLY : nbgrp
  USE funct,            ONLY : dft_is_hybrid
  USE westcom,          ONLY : l_kinetic_only
  !
  IMPLICIT NONE
  !
  CALL start_clock('wstat_readin')
  !
  ! READ INPUT_WEST
  !
  CALL fetch_input_yml(1,(/1/),.TRUE.,.FALSE.)
  !
  ! read the input file produced by the pwscf program
  ! allocate memory and recalculate what is needed
  !
  CALL read_pwout( )
  !
  ! PW checks
  !
  IF (domag) CALL errore('wstat_readin','domag version not available',1)
  IF (okvan) CALL errore('wstat_readin','ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) CALL errore('wstat_readin','double grid not implemented',1)
  IF (nbgrp > 1) THEN
     IF (dft_is_hybrid()) CALL errore('wstat_readin','band groups not implemented for hybrids',1)
     IF (l_kinetic_only) CALL errore('wstat_readin','band groups not implemented for kinetic only',1)
  ENDIF
  !
  ! READ other sections of the input file
  !
  CALL fetch_input_yml(2,(/2,5/),.TRUE.,.FALSE.)
  !
  CALL stop_clock('wstat_readin')
  !
END SUBROUTINE
