!
! Copyright (C) 2015-2025 M. Govoni
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
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE fft_types,            ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),PARAMETER :: z_0 = (0._DP,0._DP)
  COMPLEX(DP),PARAMETER :: z_i = (0._DP,1._DP)
  !
  CONTAINS
  !
  SUBROUTINE double_invfft_gamma(dfft,n,nx,a1,a2,b,cdriver)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          nx    = maximum number of PW
    !          a1    = COMPLEX array containing ONE COMPLEX function in G space
    !          a2    = COMPLEX array containing ONE COMPLEX function in G space
    ! OUTPUT : b     = ONE COMPLEX array containing TWO REAL functions in R space
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor),INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n,nx
    COMPLEX(DP),INTENT(IN) :: a1(nx)
    COMPLEX(DP),INTENT(IN) :: a2(nx)
    COMPLEX(DP),INTENT(OUT) :: b(dfft%nnr)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !$acc kernels present(b)
    b(:) = z_0
    !$acc end kernels
    !
    !$acc parallel loop present(b,dfft,dfft%nl,a1,a2,dfft%nlm)
    DO ig = 1,n
       b(dfft%nl(ig)) = a1(ig)+z_i*a2(ig)
       b(dfft%nlm(ig)) = CONJG(a1(ig)-z_i*a2(ig))
    ENDDO
    !$acc end parallel
    !
    !$acc host_data use_device(b)
    CALL invfft(cdriver,b,dfft)
    !$acc end host_data
    !
  END SUBROUTINE
  !
  SUBROUTINE double_fwfft_gamma(dfft,n,nx,a,b1,b2,cdriver)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          nx    = maximum number of PW
    !          a     = ONE COMPLEX array containing TWO REAL functions in R space
    ! OUTPUT : b1    = ONE COMPLEX array containing ONE COMPLEX function in G space
    !          b2    = ONE COMPLEX array containing ONE COMPLEX function in G space
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor),INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n,nx
    COMPLEX(DP),INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP),INTENT(OUT) :: b1(nx)
    COMPLEX(DP),INTENT(OUT) :: b2(nx)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    COMPLEX(DP) :: fm,fp
    !
    !$acc host_data use_device(a)
    CALL fwfft(cdriver,a,dfft)
    !$acc end host_data
    !
    ! Keep only G>=0
    !
    !$acc parallel loop present(a,dfft,dfft%nl,dfft%nlm,b1,b2)
    DO ig = 1,n
       fp = (a(dfft%nl(ig))+a(dfft%nlm(ig)))*0.5_DP
       fm = (a(dfft%nl(ig))-a(dfft%nlm(ig)))*0.5_DP
       b1(ig) = CMPLX(REAL(fp,KIND=DP),AIMAG(fm),KIND=DP)
       b2(ig) = CMPLX(AIMAG(fp),-REAL(fm,KIND=DP),KIND=DP)
    ENDDO
    !$acc end parallel
    !
    IF(nx > n) THEN
       !$acc kernels present(b1,b2)
       b1(n+1:nx) = z_0
       b2(n+1:nx) = z_0
       !$acc end kernels
    ENDIF
    !
  END SUBROUTINE
  !
  SUBROUTINE single_invfft_gamma(dfft,n,nx,a1,b,cdriver)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          nx    = maximum number of PW
    !          a1    = ONE COMPLEX arrays containing ONE COMPLEX functions in G space
    ! OUTPUT : b     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor),INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n,nx
    COMPLEX(DP),INTENT(IN) :: a1(nx)
    COMPLEX(DP),INTENT(OUT) :: b(dfft%nnr)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !$acc kernels present(b)
    b(:) = z_0
    !$acc end kernels
    !
    !$acc parallel loop present(b,dfft,dfft%nl,a1,dfft%nlm)
    DO ig = 1,n
       b(dfft%nl(ig)) = a1(ig)
       b(dfft%nlm(ig)) = CONJG(a1(ig))
    ENDDO
    !$acc end parallel
    !
    !$acc host_data use_device(b)
    CALL invfft(cdriver,b,dfft)
    !$acc end host_data
    !
  END SUBROUTINE
  !
  SUBROUTINE single_fwfft_gamma(dfft,n,nx,a,b1,cdriver)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          nx    = maximum number of PW
    !          a     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    ! OUTPUT : b1    = ONE COMPLEX array containing ONE COMPLEX functions in G space
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(fft_type_descriptor),INTENT(IN) :: dfft
    INTEGER,INTENT(IN) :: n,nx
    COMPLEX(DP),INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP),INTENT(OUT) :: b1(nx)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !$acc host_data use_device(a)
    CALL fwfft(cdriver,a,dfft)
    !$acc end host_data
    !
    ! Keep only G>=0
    !
    !$acc parallel loop present(b1,a,dfft,dfft%nl)
    DO ig = 1,n
       b1(ig) = a(dfft%nl(ig))
    ENDDO
    !$acc end parallel
    !
    IF(nx > n) THEN
       !$acc kernels present(b1)
       b1(n+1:nx) = z_0
       !$acc end kernels
    ENDIF
    !
  END SUBROUTINE
  !
END MODULE
