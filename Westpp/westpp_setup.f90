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
SUBROUTINE westpp_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,               ONLY : westpp_save_dir
  !
  IMPLICIT NONE
  !
  CALL do_setup()
  !
  CALL set_npwq()
  !
  CALL set_nbndocc()
  !
  CALL my_mkdir(westpp_save_dir)
  !
END SUBROUTINE
