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
SUBROUTINE wbse_setup(code)
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq,west_prefix,&
                                   & n_pdep_basis,n_pdep_eigen,n_pdep_times,l_use_ecutrho,&
                                   & wbse_init_save_dir, wbse_save_dir,&
                                   nbndval0x,l_qp_correction, nbnd_occ
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ngm,ngmx
  USE constants,              ONLY : e2,fpi
  USE cell_base,              ONLY : tpiba2
  USE control_flags,          ONLY : gamma_only
  !wbsecom combined into westcom
  !USE wbsecom,                ONLY : l_davidson, l_lanzcos, nbndval0x
  USE bse_module,             ONLY : bse_calc
  USE types_coulomb,      ONLY : pot3D
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: ig
  CHARACTER(LEN=9), INTENT(IN):: code
  !
  CALL do_setup ( )
  !
  CALL wbse_input( )
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  l_use_ecutrho = .false.
  !
  CALL set_npwq()
  !TODO: Dense?   store_sqvc(sqvc,ngm,2,isz)
  CALL pot3D%init('Dense',.FALSE.,'gb')
  CALL pot3D%print_divergence()
  !
  CALL set_nbndocc()
  !
  nbndval0x = maxval(nbnd_occ(:))
  !
  CALL west_dv_setup(bse_calc)
  !
  IF (TRIM(code) .eq. 'WBSE') THEN
      CALL my_mkdir(wbse_save_dir)
  ELSE
      CALL my_mkdir(wbse_init_save_dir)
  ENDIF
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
    IF (l_qp_correction) THEN
     !
     CALL read_qp_eigs ()
     !
  ENDIF
  !
  ! read ovl_matrix and u_matrix, and compute macroscopic term, if any
  !
  IF (bse_calc) THEN
     !
     CALL bse_init()
     !
  ENDIF
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_input
  !-----------------------------------------------------------------------
  !
  USE bse_module,    ONLY : bse_calc, &
                            l_wannier_repr, &
                            ovl_thr
  USE westcom,       ONLY : wstat_calculation, &
                            n_pdep_times, &
                            n_pdep_eigen, &
                            n_pdep_maxiter, &
                            n_pdep_read_from_file, &
                            trev_pdep, &
                            trev_pdep_rel, &
                            wbse_calculation, n_plep_times,n_plep_eigen,n_plep_maxiter,&
                            n_plep_read_from_file,trev_plep_rel,trev_plep, l_bse_calculation,&
                             l_use_localise_repr, overlap_thr
  !wbsecom combined into westcom
  !USE wbsecom
  !
  !
  IMPLICIT NONE
  !
  wstat_calculation = wbse_calculation
  n_pdep_times      = n_plep_times
  n_pdep_eigen      = n_plep_eigen
  n_pdep_maxiter    = n_plep_maxiter
  n_pdep_read_from_file = n_plep_read_from_file
  trev_pdep_rel         = trev_plep_rel
  trev_pdep             = trev_plep
  !
  bse_calc          = l_bse_calculation
  l_wannier_repr    = l_use_localise_repr
  ovl_thr           = overlap_thr
  !
END SUBROUTINE
