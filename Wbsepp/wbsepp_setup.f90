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
SUBROUTINE wbsepp_setup
  !-----------------------------------------------------------------------
  !
  USE types_coulomb,          ONLY : pot3D
  USE westcom,                ONLY : west_prefix,nbnd_occ,l_use_ecutrho,nbndval0x
  USE kinds,                  ONLY : DP
  USE io_files,               ONLY : tmp_dir
  USE westcom,                ONLY : wbse_save_dir
  !
  IMPLICIT NONE
  !
  CALL do_setup()
  !
  l_use_ecutrho = .FALSE.
  !
  CALL set_npwq()
  !
  CALL pot3D%init('Rho',.FALSE.,'gb')
  !
  CALL set_nbndocc()
  !
  nbndval0x = nbnd_occ(1)
  !
  wbse_save_dir = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.save'
  !
END SUBROUTINE
