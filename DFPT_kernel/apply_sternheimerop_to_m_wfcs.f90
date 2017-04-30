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
SUBROUTINE apply_sternheimerop_to_m_wfcs(nbndval, psi, hpsi, e, alpha, m)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : npw,npwx
  USE wvfct,                ONLY : nbnd
  USE becmod,               ONLY : bec_type, becp, calbec
  USE uspp,                 ONLY : nkb, vkb
  USE mp_global,            ONLY : intra_pool_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE westcom,              ONLY : l_kinetic_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nbndval, m
  ! input: the number of val bands
  ! input: the number of bands
  REAL(DP), INTENT(IN) :: e (m)
  COMPLEX(DP), INTENT(IN) :: alpha
  ! input: the eigenvalue
  COMPLEX(DP), INTENT(IN)  :: psi (npwx*npol,m)
  COMPLEX(DP), INTENT(OUT) :: hpsi (npwx*npol,m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  ! Workspace 
  !
  INTEGER :: ibnd, ig
  COMPLEX(DP) :: za
  !
  CALL start_clock ('stern')
  !
  !   compute the product of the hamiltonian with the h vector
  !
  hpsi=(0.0_DP,0.0_DP)
  !  
  IF(l_kinetic_only) THEN
     CALL k_psi( npwx, npw, m, psi, hpsi )
  ELSE
     CALL h_psi( npwx, npw, m, psi, hpsi )
  ENDIF
  !
  !   then we compute the operator H-epsilon S
  !
  DO ibnd=1,m
     za = CMPLX( -e(ibnd), 0._DP, KIND=DP )
     CALL ZAXPY(npw,za,psi(1,ibnd),1,hpsi(1,ibnd),1)
  ENDDO
  IF(npol==2) THEN
     DO ibnd=1,m
        za = CMPLX( -e(ibnd), 0._DP, KIND=DP )
        CALL ZAXPY(npw,za,psi(npwx+1,ibnd),1,hpsi(npwx+1,ibnd),1)
     ENDDO
  ENDIF
  !
  CALL apply_alpha_pv_to_m_wfcs(nbndval,m,psi,hpsi,alpha)
  !
  CALL stop_clock ('stern')
  !
END SUBROUTINE 
