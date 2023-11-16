!
! Copyright (C) 2015-2023 M. Govoni
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
  USE westcom,              ONLY : solver,l_bse,bse_method,l_pdep,localization,l_local_repr,&
                                 & l_use_ecutrho,wbse_init_save_dir,l_hybrid_tddft
  USE kinds,                ONLY : DP
  USE types_coulomb,        ONLY : pot3D
  USE mp_global,            ONLY : npool,nbgrp
  USE xc_lib,               ONLY : xclib_dft_is
  USE exx_base,             ONLY : exxdiv_treatment,erfc_scrlen
  USE pwcom,                ONLY : nkstot,nks
  USE distribution_center,  ONLY : kpt_pool
  USE class_idistribute,    ONLY : idistribute,IDIST_BLK
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  SELECT CASE(TRIM(solver))
  CASE('BSE','bse')
     l_bse = .TRUE.
  CASE('TDDFT','tddft')
     l_bse = .FALSE.
  END SELECT
  !
  ! ground state hybrid DFT + TDDFT -> TD-hybrid-DFT
  !
  IF((.NOT. l_bse) .AND. xclib_dft_is('hybrid')) THEN
     l_hybrid_tddft = .TRUE.
  ELSE
     l_hybrid_tddft = .FALSE.
  ENDIF
  !
  SELECT CASE(TRIM(bse_method))
  CASE('PDEP','pdep')
     l_pdep = .TRUE.
  CASE('FF_QBOX','FF_Qbox','ff_qbox')
     l_pdep = .FALSE.
  END SELECT
  !
  IF(.NOT. l_pdep) THEN
     IF(npool > 1) CALL errore('wbse_init_setup','pools not implemented for FF_Qbox',1)
     IF(nbgrp > 1) CALL errore('wbse_init_setup','band groups not implemented for FF_Qbox',1)
#if defined(__CUDA)
     CALL errore('wbse_init_setup','GPU not implemented for FF_Qbox',1)
#endif
  ENDIF
  !
  SELECT CASE(localization)
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
  IF(l_hybrid_tddft) THEN
     !
     IF(erfc_scrlen > 0._DP) THEN
        !
        ! HSE functional, mya = 1._DP, myb = -1._DP, mymu = erfc_scrlen
        !
        CALL pot3D%init('Rho',.FALSE.,exxdiv_treatment,mya=1._DP,myb=-1._DP,mymu=erfc_scrlen)
        !
     ELSE
        !
        ! PBE0 functional, mya = 1._DP, myb = 0._DP, mymu = 1._DP to avoid divergence
        !
        CALL pot3D%init('Rho',.FALSE.,exxdiv_treatment,mya=1._DP,myb=0._DP,mymu=1._DP)
        !
     ENDIF
     !
  ELSE
     !
     CALL pot3D%init('Wave',.FALSE.,'default')
     !
  ENDIF
  !
  !$acc enter data copyin(pot3D)
  !$acc enter data copyin(pot3D%sqvc)
  !
  CALL pot3D%print_divergence()
  !
  CALL set_nbndocc()
  !
  CALL my_mkdir(wbse_init_save_dir)
  !
  kpt_pool = idistribute()
  CALL kpt_pool%init(nkstot,'p','nkstot',.FALSE.,IDIST_BLK)
  !
  IF(kpt_pool%nloc /= nks) CALL errore('wbse_init_setup','unexpected kpt_pool init error',1)
  !
END SUBROUTINE
