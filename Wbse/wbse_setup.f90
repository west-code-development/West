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
SUBROUTINE wbse_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq,west_prefix,&
                                   & n_pdep_basis,n_pdep_eigen,n_pdep_times,l_use_ecutrho,&
                                   & wstat_save_dir, wstat_restart_dir, nbnd_occ
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ngm,ngmx
  USE constants,              ONLY : e2,fpi
  USE cell_base,              ONLY : tpiba2
  USE io_files,               ONLY : tmp_dir
  USE control_flags,          ONLY : gamma_only
  USE wbsecom,                ONLY : l_davidson, l_lanzcos, nbndval0x 
  USE bse_module,             ONLY : bse_calc
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: ig
  !
  CALL do_setup ( ) 
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  CALL set_npwq()
  !
  CALL set_nbndocc()
  !
  nbndval0x = maxval(nbnd_occ(:))
  !
  CALL wbse_input( )
  !
  CALL west_dv_setup(bse_calc)
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
  IF (l_lanzcos) THEN
     !
     wstat_save_dir = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.lanzcos.save'
     !
  ELSEIF (l_davidson) THEN
     ! 
     wstat_save_dir = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.david.save'
     !
  ELSE
     !
     wstat_save_dir = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.init.save'
     !
  ENDIF
  !
  CALL my_mkdir( wstat_save_dir )
  !
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
                            trev_pdep_rel
  USE wbsecom
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
