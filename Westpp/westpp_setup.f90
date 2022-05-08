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
SUBROUTINE westpp_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : westpp_save_dir,westpp_calculation,westpp_n_pdep_eigen_to_use
  USE pwcom,                  ONLY : nbnd
  USE distribution_center,    ONLY : pert,aband
  USE class_idistribute,      ONLY : idistribute
  !
  IMPLICIT NONE
  !
  INTEGER :: i
  !
  CALL do_setup()
  !
  CALL set_npwq()
  !
  CALL set_nbndocc()
  !
  CALL my_mkdir(westpp_save_dir)
  !
  DO i = 1,8
     IF(westpp_calculation(i:i) == 'r' .OR. westpp_calculation(i:i) == 'R') THEN
        aband = idistribute()
        CALL aband%init(nbnd,'i','nbnd',.TRUE.)
     ENDIF
     IF(westpp_calculation(i:i) == 'w' .OR. westpp_calculation(i:i) == 'W') THEN
        aband = idistribute()
        CALL aband%init(nbnd,'i','nbnd',.TRUE.)
     ENDIF
     IF(westpp_calculation(i:i) == 'e' .OR. westpp_calculation(i:i) == 'E') THEN
        pert = idistribute()
        CALL pert%init(westpp_n_pdep_eigen_to_use,'i','npdep',.TRUE.)
     ENDIF
     IF(westpp_calculation(i:i) == 's' .OR. westpp_calculation(i:i) == 'S') THEN
        pert = idistribute()
        CALL pert%init(westpp_n_pdep_eigen_to_use,'i','npdep',.TRUE.)
     ENDIF
  ENDDO
  !
END SUBROUTINE
