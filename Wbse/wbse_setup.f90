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
SUBROUTINE wbse_setup
  !-----------------------------------------------------------------------
  !
  USE io_global,              ONLY : stdout
  USE westcom,                ONLY : alphapv_dfpt,npwq0,sqvc,wstat_dirname,west_prefix,nbnd_occ,&
                                   & n_pdep_basis,n_pdep_eigen,n_pdep_times,isz,l_use_ecutrho
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE pwcom,                  ONLY : npw,npwx
  USE kinds,                  ONLY : DP
  USE wavefunctions_module,   ONLY : evc
  USE gvect,                  ONLY : gstart,g,ig_l2g,ngm,ngmx
  USE constants,              ONLY : e2,fpi
  USE cell_base,              ONLY : tpiba2, alat
  USE io_files,               ONLY : tmp_dir
  USE wbsecom,                ONLY : l_bse_calculation, l_lanzcos, l_davidson
  USE wbsecom,                ONLY : nbndval0x,l_qp_correction 
  USE bse_module,             ONLY : bse_calc
  USE qbox_interface,         ONLY : l_load_qbox_wfc, load_qbox_wfc, qbox_ks_wfc_filename
  !
  IMPLICIT NONE
  !
  REAL(DP) :: q(3)
  REAL(DP) :: qq
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: ig
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
  CALL set_npwq0()
  !
  ALLOCATE(sqvc(ngm))
  !
  CALL store_sqvc(sqvc,ngm,2,isz)
  !
  CALL set_nbndocc()
  !
  nbndval0x = maxval(nbnd_occ(:))
  !
  CALL west_dv_setup(bse_calc)
  !
  IF (l_lanzcos) THEN
     !
     wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.lanzcos.save'
     !
  ELSEIF (l_davidson) THEN
     ! 
     wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.david.save'
     !
  ELSE
     !
     wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.init.save'
     !
  ENDIF
  !
  CALL my_mkdir( wstat_dirname )
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
  ! if l_load_qbox_wfc == TRUE, overwrite evc by qbox wfc
  !
  !IF (l_load_qbox_wfc) THEN
     !
  !   CALL load_qbox_wfc(evc(:,1:nbndval0x), qbox_ks_wfc_filename)
     !
  !ENDIF
  !
  ! read qp_energies from file, if any
  !  
  IF (l_qp_correction) THEN
     !
     CALL read_qp_eigs ()
     !CALL read_ks_wfc  ()
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
  !
END SUBROUTINE
!
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
