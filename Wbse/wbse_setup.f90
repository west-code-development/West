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
SUBROUTINE wbse_setup()
  !-----------------------------------------------------------------------
  !
  USE westcom,          ONLY : localization,l_use_localise_repr,l_use_bisection_thr,&
                             & macropol_calculation,l_macropol,solver,l_bse_calculation,&
                             & wbse_calculation,l_davidson,l_lanczos,qp_correction,l_qp_correction,&
                             & spin_excitation,l_bse_triplet,wstat_calculation,n_pdep_times,&
                             & n_pdep_eigen,n_pdep_basis,n_pdep_maxiter,n_pdep_read_from_file,&
                             & trev_pdep_rel,trev_pdep,n_liouville_times,n_liouville_eigen,&
                             & n_liouville_maxiter,n_liouville_read_from_file,trev_liouville_rel,&
                             & trev_liouville,alphapv_dfpt,l_use_ecutrho,nbndval0x,nbnd_occ,&
                             & wbse_save_dir
  USE kinds,            ONLY : DP
  USE types_coulomb,    ONLY : pot3D
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  SELECT CASE(TRIM(localization))
  CASE('N','n')
     l_use_localise_repr = .FALSE.
     l_use_bisection_thr = .FALSE.
  CASE('B','b')
     l_use_localise_repr = .TRUE.
     l_use_bisection_thr = .TRUE.
  END SELECT
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
  END SELECT
  !
  SELECT CASE(TRIM(solver))
  CASE('BSE','bse')
     l_bse_calculation = .TRUE.
  CASE('TDDFT','tddft')
     l_bse_calculation = .FALSE.
  END SELECT
  !
  SELECT CASE(wbse_calculation)
  CASE('D','d')
     l_davidson = .TRUE.
  CASE('L','l')
     l_lanczos  = .TRUE.
  END SELECT
  !
  IF (TRIM(qp_correction) == '') THEN
     l_qp_correction = .FALSE.
  ELSE
     l_qp_correction = .TRUE.
  ENDIF
  !
  SELECT CASE(spin_excitation)
  CASE('S','s','singlet')
     l_bse_triplet = .FALSE.
  CASE('T','t','triplet')
     l_bse_triplet = .TRUE.
  END SELECT
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
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  n_pdep_maxiter = n_liouville_maxiter
  n_pdep_read_from_file = n_liouville_read_from_file
  trev_pdep_rel = trev_liouville_rel
  trev_pdep = trev_liouville
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  l_use_ecutrho = .FALSE.
  !
  CALL set_npwq()
  !
  CALL pot3D%init('Rho',.FALSE.,'gb')
  !
  CALL set_nbndocc()
  !
  nbndval0x = MAXVAL(nbnd_occ(:))
  !
  CALL west_dv_setup(l_bse_calculation)
  !
  CALL my_mkdir(wbse_save_dir)
  !
  IF (l_qp_correction) THEN
     CALL read_qp_eigs()
  ENDIF
  !
  ! read ovl_matrix and u_matrix, and compute macroscopic term, if any
  !
  IF (l_bse_calculation) THEN
     CALL bse_start()
  ENDIF
  !
END SUBROUTINE
