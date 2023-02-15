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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_init_setup()
  !-----------------------------------------------------------------------
  !
  USE westcom,          ONLY : bse_method,l_pdep,localization,l_local_repr,l_use_ecutrho,&
                             & wbse_init_save_dir
  USE kinds,            ONLY : DP
  USE types_coulomb,    ONLY : pot3D
  USE mp_global,        ONLY : nbgrp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  SELECT CASE(TRIM(bse_method))
  CASE('PDEP','pdep')
     l_pdep = .TRUE.
  CASE('FF_QBOX','FF_Qbox','ff_qbox')
     l_pdep = .FALSE.
  END SELECT
  !
  IF(.NOT. l_pdep .AND. nbgrp > 1) CALL errore('wbse_init_setup','band groups not implemented for FF_Qbox',1)
  !
  SELECT CASE(TRIM(localization))
  CASE('N','n')
     l_local_repr = .FALSE.
  CASE('B','b','W','w')
     l_local_repr = .TRUE.
  END SELECT
  !
  l_use_ecutrho = .FALSE.
  !
  CALL set_npwq()
  !
  CALL pot3D%init('Wave',.FALSE.,'default')
  !
  CALL set_nbndocc()
  !
  CALL my_mkdir(wbse_init_save_dir)
  !
END SUBROUTINE
