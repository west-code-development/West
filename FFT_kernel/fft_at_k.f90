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
  PRIVATE
  !
  INTERFACE single_invfft_k
     MODULE PROCEDURE single_invfft_k_cpu
#if defined(__CUDA)
     MODULE PROCEDURE single_invfft_k_gpu
#endif
  END INTERFACE
  !
  INTERFACE single_fwfft_k
     MODULE PROCEDURE single_fwfft_k_cpu
#if defined(__CUDA)
     MODULE PROCEDURE single_fwfft_k_gpu
#endif
  END INTERFACE
  !
  PUBLIC :: single_invfft_k
  PUBLIC :: single_fwfft_k
  !
  COMPLEX(DP), PARAMETER :: z_0 = (0._DP,0._DP)
  !
  CONTAINS
  !
  !
  SUBROUTINE single_invfft_k_cpu(dfft,n,nx,a1,b,cdriver,igk)
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
    b = z_0
    IF(PRESENT(igk)) THEN
       DO ig=1,n
          b(dfft%nl(igk(ig)))=a1(ig)
       ENDDO
    ELSE
       DO ig=1,n
          b(dfft%nl(ig))=a1(ig)
       ENDDO
    ENDIF
    !
    CALL invfft(cdriver,b,dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_fwfft_k_cpu(dfft,n,nx,a,b1,cdriver,igk)
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
    CALL fwfft(cdriver,a,dfft)
    !
    b1 = z_0
    !
    IF(PRESENT(igk)) THEN
       DO ig=1,n
          b1(ig) = a(dfft%nl(igk(ig)))
       ENDDO
    ELSE
       DO ig=1,n
          b1(ig) = a(dfft%nl(ig))
       ENDDO
    ENDIF
    !
  END SUBROUTINE
  !
#if defined(__CUDA)
  SUBROUTINE single_invfft_k_gpu(dfft,n,nx,a1,b,cdriver,igk)
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
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: n, nx
    COMPLEX(DP), DEVICE, INTENT(IN) :: a1(nx)
    COMPLEX(DP), DEVICE, INTENT(OUT) :: b(dfft%nnr)
    CHARACTER(LEN=*), INTENT(IN) :: cdriver
    INTEGER, INTENT(IN), OPTIONAL :: igk(n)
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
          b(dfft_nl_d(igk(ig))) = a1(ig)
       ENDDO
       !$acc end parallel
    ELSE
       !$acc parallel loop
       DO ig = 1,n
          b(dfft_nl_d(ig)) = a1(ig)
       ENDDO
       !$acc end parallel
    ENDIF
    !
    CALL invfft(cdriver,b,dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_fwfft_k_gpu(dfft,n,nx,a,b1,cdriver,igk)
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
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: n,nx
    COMPLEX(DP), DEVICE, INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP), DEVICE, INTENT(OUT) :: b1(nx)
    CHARACTER(LEN=*), INTENT(IN) :: cdriver
    INTEGER, INTENT(IN), OPTIONAL :: igk(n)
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
          b1(ig) = a(dfft_nl_d(igk(ig)))
       ENDDO
       !$acc end parallel
    ELSE
       !$acc parallel loop
       DO ig = 1,n
          b1(ig) = a(dfft_nl_d(ig))
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
  !
#endif
END MODULE
