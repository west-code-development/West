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
SUBROUTINE set_dirs( )
  !-----------------------------------------------------------------------
  !
  USE io_files,               ONLY : tmp_dir
  USE westcom,                ONLY : west_prefix, wstat_save_dir, wstat_restart_dir
  !
  IMPLICIT NONE
  !
  wstat_save_dir    = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
  wstat_restart_dir = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.restart'
  !
END SUBROUTINE
