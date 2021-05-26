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
FUNCTION get_alpha_pv()
  !-----------------------------------------------------------------------
  !
  USE pwcom,                  ONLY : nbnd,et
  USE kinds,                  ONLY : DP
  USE mp,                     ONLY : mp_min,mp_max
  USE mp_global,              ONLY : inter_pool_comm
  USE klist,                  ONLY : nks
  USE westcom,                ONLY : nbnd_occ
  USE types_bz_grid,          ONLY : k_grid
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: get_alpha_pv
  !
  INTEGER :: ibnd, iks
  REAL(DP) :: emin, emax
  REAL(DP) :: alpha_pv
  !
  CALL set_nbndocc()
  !
  ! Calculate ALPHA_PV
  !
  emin = et (1, 1)
  DO iks = 1, k_grid%nps
     DO ibnd = 1, nbnd
        emin = MIN (emin, et (ibnd, iks) )
     ENDDO
  ENDDO
  !
  CALL mp_min( emin, inter_pool_comm) 
  !
  emax = et (1, 1)
  DO iks = 1, k_grid%nps 
     DO ibnd = 1, nbnd_occ(iks)
        emax = MAX (emax, et (ibnd, iks) )
     ENDDO
  ENDDO
  !
  CALL mp_max( emax, inter_pool_comm)
  !
  alpha_pv = 2.0_DP * (emax - emin) 
  ! avoid zero value for alpha_pv
  alpha_pv = MAX (alpha_pv, 1.0d-2)
  !
  get_alpha_pv = CMPLX( alpha_pv, 0._DP, KIND=DP )
  !
END FUNCTION 
