!
! Copyright (C) 2015-2024 M. Govoni
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
  USE noncollin_module,     ONLY : npol,noncolin
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
  COMPLEX(DP), INTENT(IN) :: psi (npwx*npol,m)
  COMPLEX(DP), INTENT(OUT) :: hpsi (npwx*npol,m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  ! Workspace
  !
  INTEGER :: ibnd,ig
  COMPLEX(DP) :: za
  !
#if defined(__CUDA)
  CALL start_clock_gpu('stern')
#else
  CALL start_clock('stern')
#endif
  !
  ! compute the product of the hamiltonian with the h vector
  !
  IF(l_kinetic_only) THEN
     CALL k_psi( npwx, npw, m, psi, hpsi )
  ELSE
     !
     ! use h_psi_, i.e. h_psi without band parallelization, as wstat
     ! handles band parallelization separately in dfpt_module
     !
#if defined(__CUDA)
     !$acc host_data use_device(psi,hpsi)
     CALL h_psi__gpu( npwx, npw, m, psi, hpsi )
     !$acc end host_data
#else
     CALL h_psi_( npwx, npw, m, psi, hpsi )
#endif
  ENDIF
  !
  ! then we compute the operator H-epsilon S
  !
#if defined(__CUDA)
  !$acc parallel loop collapse(2) present(hpsi,e,psi)
  DO ibnd = 1,m
     DO ig = 1,npw
        hpsi(ig,ibnd) = hpsi(ig,ibnd)-e(ibnd)*psi(ig,ibnd)
        IF(noncolin) THEN
           hpsi(npwx+ig,ibnd) = hpsi(npwx+ig,ibnd)-e(ibnd)*psi(npwx+ig,ibnd)
        ENDIF
     ENDDO
  ENDDO
  !$acc end parallel
#else
  DO ibnd = 1,m
     za = -e(ibnd)
     CALL ZAXPY(npw,za,psi(1,ibnd),1,hpsi(1,ibnd),1)
     IF(noncolin) THEN
        CALL ZAXPY(npw,za,psi(npwx+1,ibnd),1,hpsi(npwx+1,ibnd),1)
     ENDIF
  ENDDO
#endif
  !
  CALL apply_alpha_pv_to_m_wfcs(nbndval,m,psi,hpsi,alpha)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('stern')
#else
  CALL stop_clock('stern')
#endif
  !
END SUBROUTINE
