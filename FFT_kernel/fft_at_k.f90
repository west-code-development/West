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
MODULE fft_at_k
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
  !
  CONTAINS
  !
  SUBROUTINE single_invfft_k(dfft,n,nx,a1,b,cdriver,igk)
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
    INTEGER,INTENT(IN),OPTIONAL :: igk(n)
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !$acc kernels present(b)
    b(:) = z_0
    !$acc end kernels
    !
    IF(PRESENT(igk)) THEN
       !$acc parallel loop present(b,dfft,dfft%nl,igk,a1)
       DO ig = 1,n
          b(dfft%nl(igk(ig))) = a1(ig)
       ENDDO
       !$acc end parallel
    ELSE
       !$acc parallel loop present(b,dfft,dfft%nl,a1)
       DO ig = 1,n
          b(dfft%nl(ig)) = a1(ig)
       ENDDO
       !$acc end parallel
    ENDIF
    !
    !$acc host_data use_device(b)
    CALL invfft(cdriver,b,dfft)
    !$acc end host_data
    !
  END SUBROUTINE
  !
  SUBROUTINE single_fwfft_k(dfft,n,nx,a,b1,cdriver,igk)
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
    INTEGER,INTENT(IN),OPTIONAL :: igk(n)
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !$acc host_data use_device(a)
    CALL fwfft(cdriver,a,dfft)
    !$acc end host_data
    !
    IF(PRESENT(igk)) THEN
       !$acc parallel loop present(b1,a,dfft,dfft%nl,igk)
       DO ig = 1,n
          b1(ig) = a(dfft%nl(igk(ig)))
       ENDDO
       !$acc end parallel
    ELSE
       !$acc parallel loop present(b1,a,dfft,dfft%nl)
       DO ig = 1,n
          b1(ig) = a(dfft%nl(ig))
       ENDDO
       !$acc end parallel
    ENDIF
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
