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
SUBROUTINE westpp_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq0,sqvc,west_prefix,westpp_dirname,&
                                   & n_imfreq,nbnd_occ,l_macropol,macropol_calculation,&
                                   & n_refreq,isz,qp_bandrange,westpp_calculation,westpp_n_pdep_eigen_to_use
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE pwcom,                  ONLY : npw,nbnd
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ig_l2g
  USE io_files,               ONLY : tmp_dir
  USE distribution_center,    ONLY : pert,macropert,ifr,rfr,aband
  USE class_idistribute,      ONLY : idistribute
  USE wavefunctions_module,   ONLY : evc
  USE mod_mpiio,              ONLY : set_io_comm
  USE pdep_db,                ONLY : pdep_db_read
  !
  IMPLICIT NONE
  !
  REAL(DP) :: q(3)
  REAL(DP) :: qq
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: ig, i
  !
  CALL do_setup ( ) 
  !
  CALL set_npwq0()
  !
  ALLOCATE(sqvc(npwq0))
  !
  CALL store_sqvc(sqvc,npwq0,1,isz)
  !
  CALL set_nbndocc()
  !
  westpp_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.westpp.save'
  CALL my_mkdir( westpp_dirname )
  !
  DO i = 1, 8
     IF( westpp_calculation(i:i) == 'r' .OR. westpp_calculation(i:i) == 'R' ) THEN 
        aband = idistribute()
        CALL aband%init(nbnd,'i','nbnd',.TRUE.)
     ENDIF  
     IF( westpp_calculation(i:i) == 'w' .OR. westpp_calculation(i:i) == 'W' ) THEN 
        aband = idistribute()
        CALL aband%init(nbnd,'i','nbnd',.TRUE.)
     ENDIF  
     IF( westpp_calculation(i:i) == 'e' .OR. westpp_calculation(i:i) == 'E' ) THEN
        pert = idistribute()
        CALL pert%init(westpp_n_pdep_eigen_to_use,'i','npdep',.TRUE.)
        CALL pdep_db_read(westpp_n_pdep_eigen_to_use)
     ENDIF 
     IF( westpp_calculation(i:i) == 's' .OR. westpp_calculation(i:i) == 'S' ) THEN
        pert = idistribute()
        CALL pert%init(westpp_n_pdep_eigen_to_use,'i','npdep',.TRUE.)
        CALL pdep_db_read(westpp_n_pdep_eigen_to_use)
     ENDIF 
  ENDDO
  !
  CALL set_io_comm( ) ! this defines communicator between heads of each image (me_bgrp==0) 
  !
END SUBROUTINE 
