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
MODULE fft_at_k
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
  PUBLIC :: single_invfft_k, single_fwfft_k 
  !
  CONTAINS
  !
  !
  SUBROUTINE single_invfft_k(dfft,n,nx,a1,b,cdriver,igk)
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
    INTEGER,INTENT(IN),OPTIONAL :: igk(n)
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    b=0.0_DP
    IF(PRESENT(igk)) THEN
       DO ig=1,n
          b(nls(igk(ig)))=a1(ig)
       ENDDO
    ELSE
       DO ig=1,n
          b(nls(ig))=a1(ig)
       ENDDO
    ENDIF
    !
    CALL invfft(cdriver,b,dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_fwfft_k(dfft,n,nx,a,b1,cdriver,igk)
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
    INTEGER,INTENT(IN) :: n,nx
    COMPLEX(DP),INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP),INTENT(OUT) :: b1(nx)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    INTEGER,INTENT(IN),OPTIONAL :: igk(n)
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    CALL fwfft(cdriver, a, dfft)
    !
    b1=0.0_DP
    !
    IF(PRESENT(igk)) THEN 
       DO ig=1,n 
          b1(ig) = a(nls(igk(ig)))
       ENDDO
    ELSE
       DO ig=1,n 
          b1(ig) = a(nls(ig))
       ENDDO
    ENDIF
    !
  END SUBROUTINE
  !
END MODULE
