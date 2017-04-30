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
MODULE fft_at_gamma
  !-----------------------------------------------------------------------
  !
  ! Everything is done following dffts
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nl,nlm
  USE gvecs,                ONLY : nls,nlsm
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE fft_types,            ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: double_invfft_gamma, double_fwfft_gamma 
  PUBLIC :: single_invfft_gamma, single_fwfft_gamma 
  !
  CONTAINS
  !
  !
  SUBROUTINE double_invfft_gamma(dfft,n,nx,a1,a2,b,cdriver)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          a1    = COMPLEX array containing ONE COMPLEX function in G space
    !          a2    = COMPLEX array containing ONE COMPLEX function in G space
    !          lda   = leading dimension of a1 or a2
    !          ldb   = leading dimension of b
    ! OUTPUT : b     = ONE COMPLEX array containing TWO REAL functions in R space
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n, nx
    COMPLEX(DP),INTENT(IN) :: a1(nx)
    COMPLEX(DP),INTENT(IN) :: a2(nx)
    COMPLEX(DP),INTENT(OUT) :: b(dfft%nnr)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
!$OMP PARALLEL private(ig)
!$OMP DO
    DO ig=1,dfft%nnr
       b(ig)= (0.0_DP,0.0_DP) 
    ENDDO
!$OMP ENDDO
!$OMP DO
    DO ig=1,n 
       b(nls (ig))=        a1(ig) + (0.0_DP,1.0_DP) * a2(ig)
       b(nlsm(ig))= DCONJG(a1(ig) - (0.0_DP,1.0_DP) * a2(ig))
    ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
    !
    CALL invfft(cdriver, b, dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE double_fwfft_gamma(dfft,n,nx,a,b1,b2,cdriver)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          a     = ONE COMPLEX array containing TWO REAL functions in R space
    !          lda   = leading dimension of a
    !          ldb   = leading dimension of b1 or b2 
    ! OUTPUT : b1     = ONE COMPLEX array containing ONE COMPLEX function in G space
    !          b2     = ONE COMPLEX array containing ONE COMPLEX function in G space
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n, nx
    COMPLEX(DP),INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP),INTENT(OUT) :: b1(nx)
    COMPLEX(DP),INTENT(OUT) :: b2(nx)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    COMPLEX(DP) :: fm, fp 
    !
    ! FFT call
    CALL fwfft(cdriver, a, dfft)
    ! Keep only G>=0
    !
!$OMP PARALLEL private(ig,fp,fm)
!$OMP DO
    DO ig = 1, n
       fp = ( a(nls (ig)) + a(nlsm(ig)) )*0.5_DP
       fm = ( a(nls (ig)) - a(nlsm(ig)) )*0.5_DP
       b1(ig) = CMPLX( DBLE(fp), DIMAG(fm), KIND=DP)
       b2(ig) = CMPLX(DIMAG(fp), -DBLE(fm), KIND=DP)
    ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
    !
    DO ig = (n+1), nx
       b1(ig) = (0.0_DP,0.0_DP)
       b2(ig) = (0.0_DP,0.0_DP)
    ENDDO
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_invfft_gamma(dfft,n,nx,a1,b,cdriver)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          a1     = ONE COMPLEX arrays containing ONE COMPLEX functions in G space
    !          lda   = leading dimension of a1
    !          ldb   = leading dimension of b          
    ! OUTPUT : b     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n, nx
    COMPLEX(DP),INTENT(IN) :: a1(nx)
    COMPLEX(DP),INTENT(OUT) :: b(dfft%nnr)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !
!$OMP PARALLEL private(ig)
!$OMP DO
    DO ig=1,dfft%nnr
       b(ig)= (0.0_DP,0.0_DP) 
    ENDDO
!$OMP ENDDO
!$OMP DO
    DO ig=1,n
       b(nls (ig))=        a1(ig)
       b(nlsm(ig))= DCONJG(a1(ig))
    ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
    !
    CALL invfft(cdriver, b, dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_fwfft_gamma(dfft,n,nx,a,b1,cdriver)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          a     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !          lda   = leading dimension of a
    !          ldb   = leading dimension of b1 
    ! OUTPUT : b1     = ONE COMPLEX array containing ONE COMPLEX functions in G space
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n, nx
    COMPLEX(DP),INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP),INTENT(OUT) :: b1(nx)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !
    CALL fwfft(cdriver, a, dfft)
    ! Keep only G>=0
    !
!$OMP PARALLEL private(ig)
!$OMP DO
    DO ig=1,n 
       b1(ig) = a(nls(ig))
    ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
    !
    DO ig = (n+1), nx
       b1(ig) = (0.0_DP,0.0_DP)
    ENDDO
    !
  END SUBROUTINE
  !
  !
END MODULE
