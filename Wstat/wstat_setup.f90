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
SUBROUTINE wstat_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,n_pdep_basis,n_pdep_eigen,&
                                   & n_pdep_times,wstat_save_dir
  USE kinds,                  ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  CALL set_npwq()
  !
  CALL set_nbndocc()
  !
  CALL my_mkdir(wstat_save_dir)
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
END SUBROUTINE
