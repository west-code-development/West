!
! Copyright (C) 2015-2017 M. Govoni 
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
SUBROUTINE westpp_readin()
  !-----------------------------------------------------------------------
  !
  USE ions_base,        ONLY : nat
  USE uspp,             ONLY : okvan
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast
  USE gvecs,            ONLY : doublegrid
  USE spin_orb,         ONLY : domag
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iunit =5, ios
  !
  CALL start_clock('westpp_readin')
  !
  !CALL fetch_namelist(3,(/1,2,4/))
  !
  ! READ INPUT_WEST
  !
  CALL fetch_input(1,(/1/),.TRUE.)
  !CALL fetch_input(3,(/1,2,4/),.TRUE.)
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  !CALL read_file( )
  CALL read_pwout( )
  !
  ! PW checks
  !
  IF (domag) CALL errore('westpp_readin','domag version not available',1)
  IF (okvan) CALL errore('westpp_readin','ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) CALL errore('westpp_readin', 'double grid not implemented',1)
  !
  ! READ other sections of the input file
  !
  CALL fetch_input(2,(/2,4/),.TRUE.)
  !
  CALL stop_clock('westpp_readin')
  !
END SUBROUTINE
