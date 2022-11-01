!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE wbse_setup(code)
  !-----------------------------------------------------------------------
  !
  USE bse_module,       ONLY : bse_calc
  USE westcom,          ONLY : alphapv_dfpt,n_pdep_basis,n_pdep_eigen,n_pdep_times,l_use_ecutrho,&
                             & wbse_init_save_dir,wbse_save_dir,nbndval0x,l_qp_correction,nbnd_occ
  USE kinds,            ONLY : DP
  USE types_coulomb,    ONLY : pot3D
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9), INTENT(IN):: code
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  CALL wbse_input()
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  l_use_ecutrho = .FALSE.
  !
  CALL set_npwq()
  CALL pot3D%init('Rho',.FALSE.,'gb')
  CALL pot3D%print_divergence()
  !
  CALL set_nbndocc()
  !
  nbndval0x = MAXVAL(nbnd_occ(:))
  !
  CALL west_dv_setup(bse_calc)
  !
  IF (TRIM(code) == 'WBSE') THEN
     CALL my_mkdir(wbse_save_dir)
  ELSE
     CALL my_mkdir(wbse_init_save_dir)
  ENDIF
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
  IF (l_qp_correction) THEN
     CALL read_qp_eigs()
  ENDIF
  !
  ! read ovl_matrix and u_matrix, and compute macroscopic term, if any
  !
  IF (bse_calc) THEN
     CALL bse_init()
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_input
  !-----------------------------------------------------------------------
  !
  USE bse_module,       ONLY : bse_calc,l_wannier_repr,ovl_thr
  USE westcom,          ONLY : wstat_calculation,n_pdep_times,n_pdep_eigen,n_pdep_maxiter,&
                             & n_pdep_read_from_file,trev_pdep,trev_pdep_rel,wbse_calculation,&
                             & n_liouville_times,n_liouville_eigen,n_liouville_maxiter,&
                             & n_liouville_read_from_file,trev_liouville_rel,trev_liouville,&
                             & l_bse_calculation,l_use_localise_repr,overlap_thr
  !
  IMPLICIT NONE
  !
  SELECT CASE(wbse_calculation)
  CASE('l','d','r','R')
     wstat_calculation = 'R'
  CASE('L','D','s','S')
     wstat_calculation = 'S'
  CASE DEFAULT
     wstat_calculation = wbse_calculation
  END SELECT
  !
  n_pdep_times = n_liouville_times
  n_pdep_eigen = n_liouville_eigen
  n_pdep_maxiter = n_liouville_maxiter
  n_pdep_read_from_file = n_liouville_read_from_file
  trev_pdep_rel = trev_liouville_rel
  trev_pdep = trev_liouville
  bse_calc = l_bse_calculation
  l_wannier_repr = l_use_localise_repr
  ovl_thr = overlap_thr
  !
END SUBROUTINE
