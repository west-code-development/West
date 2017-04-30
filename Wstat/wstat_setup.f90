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
SUBROUTINE wstat_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq0,sqvc,wstat_dirname,west_prefix,&
                                   & n_pdep_basis,n_pdep_eigen,n_pdep_times,isz,l_use_ecutrho
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE pwcom,                  ONLY : npw,npwx
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ig_l2g,ngm,ngmx
  USE constants,              ONLY : e2,fpi
  USE cell_base,              ONLY : tpiba2
  USE io_files,               ONLY : tmp_dir
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
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  CALL set_npwq0()
  !
  ALLOCATE(sqvc(npwq0))
  !
  CALL store_sqvc(sqvc,npwq0,1,isz)
  !
  CALL set_nbndocc()
  !
  wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
  CALL my_mkdir( wstat_dirname )
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
END SUBROUTINE 
