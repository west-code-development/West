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
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE k_psi(lda,n,m,psi,hpsi)
  !-----------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian (kinetic-only)
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  USE kinds,            ONLY : DP
  USE gvect,            ONLY : gstart
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : npol,noncolin
#if defined(__CUDA)
  USE wvfct_gpum,       ONLY : g2kin=>g2kin_d
#else
  USE wvfct,            ONLY : g2kin
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda,n,m
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)
  !
  INTEGER :: ibnd,ig
  !
#if defined(__CUDA)
  CALL start_clock_gpu('k_psi')
#else
  CALL start_clock('k_psi')
#endif
  !
  ! ... Here we apply the kinetic energy (k+G)^2 psi
  !
#if defined(__CUDA)
  !$acc parallel loop collapse(2) present(hpsi,psi)
#else
  !$OMP PARALLEL DEFAULT(NONE) SHARED(m,n,hpsi,g2kin,psi,lda,noncolin) PRIVATE(ibnd,ig)
  !$OMP DO COLLAPSE(2)
#endif
  DO ibnd = 1,m
     DO ig = 1,lda
        IF(ig <= n) THEN
           hpsi(ig,ibnd) = g2kin(ig)*psi(ig,ibnd)
           IF(noncolin) THEN
              hpsi(lda+ig,ibnd) = g2kin(ig)*psi(lda+ig,ibnd)
           ENDIF
        ELSE
           hpsi(ig,ibnd) = 0._DP
           IF(noncolin) THEN
              hpsi(lda+ig,ibnd) = 0._DP
           ENDIF
        ENDIF
     ENDDO
  ENDDO
#if defined(__CUDA)
  !$acc end parallel
#else
  !$OMP ENDDO
  !$OMP END PARALLEL
#endif
  !
  ! ... Gamma-only trick: set to zero the imaginary part of hpsi at G=0
  !
  IF(gamma_only .AND. gstart == 2) THEN
#if defined(__CUDA)
     !$acc parallel loop present(hpsi)
#else
     !$OMP PARALLEL DO
#endif
     DO ibnd = 1,m
        hpsi(1,ibnd) = CMPLX(REAL(hpsi(1,ibnd),KIND=DP),0._DP,KIND=DP)
     ENDDO
#if defined(__CUDA)
     !$acc end parallel
#else
     !$OMP END PARALLEL DO
#endif
  ENDIF
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('k_psi')
#else
  CALL stop_clock('k_psi')
#endif
  !
END SUBROUTINE
