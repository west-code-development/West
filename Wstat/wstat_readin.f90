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
SUBROUTINE wstat_readin()
  !-----------------------------------------------------------------------
  !
  USE ions_base,        ONLY : nat
  USE uspp,             ONLY : okvan
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : nproc,mpime,root
  USE gvecs,            ONLY : doublegrid
  USE spin_orb,         ONLY : domag
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iunit =5, ios
  !
  CALL start_clock('wstat_readin')
  !
  !CALL fetch_namelist(2,(/1,2/))
  !
  ! READ INPUT_WEST
  !
  CALL fetch_input_yml(1,(/1/),.TRUE.,.FALSE.)
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  !CALL read_file( )
  CALL read_pwout( )
  !
  ! PW checks
  !
  IF (domag) CALL errore('wstat_readin','domag version not available',1)
  IF (okvan) CALL errore('wstat_readin','ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) CALL errore('wstat_readin', 'double grid not implemented',1)
  !
  ! READ other sections of the input file
  !
  CALL fetch_input_yml(1,(/2/),.TRUE.,.FALSE.)
  !
  CALL stop_clock('wstat_readin')
  !
END SUBROUTINE
