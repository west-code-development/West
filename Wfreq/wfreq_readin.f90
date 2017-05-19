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
SUBROUTINE wfreq_readin()
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
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: iunit =5, ios
  !
  CALL start_clock('wfreq_readin')
  !
  !CALL fetch_namelist(3,(/1,2,3/))
  CALL fetch_input(3,(/1,2,3/))
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  !CALL read_file( )
  CALL read_pwout( )
  !
  ! PW checks
  !
  IF (domag) CALL errore('wfreq_readin','domag version not available',1)
  IF (okvan) CALL errore('wfreq_readin','ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) CALL errore('wfreq_readin', 'double grid not implemented',1)
  !
  CALL stop_clock('wfreq_readin')
  !
END SUBROUTINE
