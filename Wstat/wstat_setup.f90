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
SUBROUTINE wstat_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq,west_prefix,&
                                   & n_pdep_basis,n_pdep_eigen,n_pdep_times,isz,l_use_ecutrho,&
                                   & wstat_save_dir, wstat_restart_dir
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ig_l2g,ngm,ngmx
  USE constants,              ONLY : e2,fpi
  USE cell_base,              ONLY : tpiba2
  USE io_files,               ONLY : tmp_dir
  USE control_flags,          ONLY : gamma_only
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
  CALL my_mkdir( wstat_save_dir )
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
END SUBROUTINE
