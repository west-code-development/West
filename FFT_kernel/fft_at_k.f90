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
MODULE fft_at_k
  !-----------------------------------------------------------------------
  !
  ! Everything is done following dffts
  !
  USE kinds,                ONLY : DP
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE fft_types,            ONLY : fft_type_descriptor
#if defined(__CUDA)
  USE west_gpu,             ONLY : dfft_nl_d,dfft_nlm_d
#endif
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
#if defined(__CUDA)
    ATTRIBUTES(DEVICE) :: a1,b
#endif
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    INTEGER,INTENT(IN),OPTIONAL :: igk(n)
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !$acc kernels
    b(:) = z_0
    !$acc end kernels
    !
    IF(PRESENT(igk)) THEN
       !$acc parallel loop present(igk)
       DO ig = 1,n
#if defined(__CUDA)
          b(dfft_nl_d(igk(ig))) = a1(ig)
#else
          b(dfft%nl(igk(ig))) = a1(ig)
#endif
       ENDDO
       !$acc end parallel
    ELSE
       !$acc parallel loop
       DO ig = 1,n
#if defined(__CUDA)
          b(dfft_nl_d(ig)) = a1(ig)
#else
          b(dfft%nl(ig)) = a1(ig)
#endif
       ENDDO
       !$acc end parallel
    ENDIF
    !
    CALL invfft(cdriver,b,dfft)
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
#if defined(__CUDA)
    ATTRIBUTES(DEVICE) :: a,b1
#endif
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    INTEGER,INTENT(IN),OPTIONAL :: igk(n)
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    CALL fwfft(cdriver,a,dfft)
    !
    IF(PRESENT(igk)) THEN
       !$acc parallel loop present(igk)
       DO ig = 1,n
#if defined(__CUDA)
          b1(ig) = a(dfft_nl_d(igk(ig)))
#else
          b1(ig) = a(dfft%nl(igk(ig)))
#endif
       ENDDO
       !$acc end parallel
    ELSE
       !$acc parallel loop
       DO ig = 1,n
#if defined(__CUDA)
          b1(ig) = a(dfft_nl_d(ig))
#else
          b1(ig) = a(dfft%nl(ig))
#endif
       ENDDO
       !$acc end parallel
    ENDIF
    !
    IF(nx > n) THEN
       !$acc kernels
       b1(n+1:nx) = z_0
       !$acc end kernels
    ENDIF
    !
  END SUBROUTINE
  !
END MODULE
