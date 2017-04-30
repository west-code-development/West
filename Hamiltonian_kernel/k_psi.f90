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
SUBROUTINE k_psi( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
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
  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : g2kin
  USE gvect,    ONLY : gstart
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY: npol, noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda*npol,m) 
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)   
  !
  INTEGER     :: ipol, ibnd, incr, ig
  !
  CALL start_clock( 'k_psi' )
  !  
  ! ... Here we apply the kinetic energy (k+G)^2 psi
  !
  IF (noncolin) THEN
!$OMP PARALLEL DEFAULT(none) SHARED(m,n,hpsi,g2kin,psi,lda) PRIVATE(ibnd,ig)
!$OMP DO 
     DO ibnd = 1,m
        DO ig = 1, n
           hpsi (ig, ibnd) = g2kin (ig) * psi (ig, ibnd)
        ENDDO
        DO ig = n+1, lda
           hpsi (ig, ibnd) = 0._DP
        ENDDO
        DO ig = 1, n
           hpsi (lda+ig, ibnd) = g2kin (ig) * psi (lda+ig, ibnd)
        ENDDO
        DO ig = n+1, lda
           hpsi (lda+ig, ibnd) = 0._DP
        ENDDO
     ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
  ELSE 
!$OMP PARALLEL DEFAULT(none) SHARED(m,n,hpsi,g2kin,psi,lda) PRIVATE(ibnd,ig)
!$OMP DO 
     DO ibnd = 1,m
        DO ig = 1, n
           hpsi (ig, ibnd) = g2kin (ig) * psi (ig, ibnd)
        ENDDO
        DO ig = n+1, lda
           hpsi (ig, ibnd) = 0._DP
        ENDDO
     ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
  ENDIF
  !
  ! ... Gamma-only trick: set to zero the imaginary part of hpsi at G=0
  !
  IF ( gamma_only .AND. gstart == 2 ) &
      hpsi(1,1:m) = CMPLX( DBLE( hpsi(1,1:m) ), 0.D0 ,kind=DP)
  !
  CALL stop_clock( 'k_psi' )
  !
  RETURN
  !
END SUBROUTINE k_psi
