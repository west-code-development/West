!
! Copyright (C) 2015-2021 M. Govoni
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
#if defined(__CUDA)
  USE west_gpu,             ONLY : dfft_nl_d,dfft_nlm_d
#endif
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  INTERFACE double_invfft_gamma
     MODULE PROCEDURE double_invfft_gamma_cpu
#if defined(__CUDA)
     MODULE PROCEDURE double_invfft_gamma_gpu
#endif
  END INTERFACE
  !
  INTERFACE double_fwfft_gamma
     MODULE PROCEDURE double_fwfft_gamma_cpu
#if defined(__CUDA)
     MODULE PROCEDURE double_fwfft_gamma_gpu
#endif
  END INTERFACE
  !
  INTERFACE single_invfft_gamma
     MODULE PROCEDURE single_invfft_gamma_cpu
#if defined(__CUDA)
     MODULE PROCEDURE single_invfft_gamma_gpu
#endif
  END INTERFACE
  !
  INTERFACE single_fwfft_gamma
     MODULE PROCEDURE single_fwfft_gamma_cpu
#if defined(__CUDA)
     MODULE PROCEDURE single_fwfft_gamma_gpu
#endif
  END INTERFACE
  !
  PUBLIC :: double_invfft_gamma
  PUBLIC :: double_fwfft_gamma
  PUBLIC :: single_invfft_gamma
  PUBLIC :: single_fwfft_gamma
  !
  COMPLEX(DP), PARAMETER :: z_0 = (0._DP,0._DP)
  COMPLEX(DP), PARAMETER :: z_i = (0._DP,1._DP)
  !
  CONTAINS
  !
  !
  SUBROUTINE double_invfft_gamma_cpu(dfft,n,nx,a1,a2,b,cdriver)
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
    !$OMP PARALLEL PRIVATE(ig)
    !$OMP DO
    DO ig=1,dfft%nnr
       b(ig)= z_0
    ENDDO
    !$OMP ENDDO
    !$OMP DO
    DO ig=1,n
       b(dfft%nl (ig))=       a1(ig) + z_i * a2(ig)
       b(dfft%nlm(ig))= CONJG(a1(ig) - z_i * a2(ig))
    ENDDO
    !$OMP ENDDO
    !$OMP END PARALLEL
    !
    CALL invfft(cdriver, b, dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE double_fwfft_gamma_cpu(dfft,n,nx,a,b1,b2,cdriver)
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
    CALL fwfft(cdriver, a, dfft)
    !
    ! Keep only G>=0
    !
    !$OMP PARALLEL PRIVATE(ig,fp,fm)
    !$OMP DO
    DO ig = 1, n
       fp = ( a(dfft%nl (ig)) + a(dfft%nlm(ig)) )*0.5_DP
       fm = ( a(dfft%nl (ig)) - a(dfft%nlm(ig)) )*0.5_DP
       b1(ig) = CMPLX(REAL(fp,KIND=DP), AIMAG(fm), KIND=DP)
       b2(ig) = CMPLX(AIMAG(fp), -REAL(fm,KIND=DP), KIND=DP)
    ENDDO
    !$OMP ENDDO
    !$OMP END PARALLEL
    !
    DO ig = (n+1), nx
       b1(ig) = z_0
       b2(ig) = z_0
    ENDDO
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_invfft_gamma_cpu(dfft,n,nx,a1,b,cdriver)
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
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    !$OMP PARALLEL PRIVATE(ig)
    !$OMP DO
    DO ig=1,dfft%nnr
       b(ig)= z_0
    ENDDO
    !$OMP ENDDO
    !$OMP DO
    DO ig=1,n
       b(dfft%nl (ig))=       a1(ig)
       b(dfft%nlm(ig))= CONJG(a1(ig))
    ENDDO
    !$OMP ENDDO
    !$OMP END PARALLEL
    !
    CALL invfft(cdriver, b, dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_fwfft_gamma_cpu(dfft,n,nx,a,b1,cdriver)
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
    INTEGER,INTENT(IN) :: n, nx
    COMPLEX(DP),INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP),INTENT(OUT) :: b1(nx)
    CHARACTER(LEN=*),INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    CALL fwfft(cdriver, a, dfft)
    !
    ! Keep only G>=0
    !
    !$OMP PARALLEL PRIVATE(ig)
    !$OMP DO
    DO ig=1,n
       b1(ig) = a(dfft%nl(ig))
    ENDDO
    !$OMP ENDDO
    !$OMP END PARALLEL
    !
    DO ig = (n+1), nx
       b1(ig) = z_0
    ENDDO
    !
  END SUBROUTINE
  !
  !
#if defined(__CUDA)
  SUBROUTINE double_invfft_gamma_gpu(dfft,n,nx,a1,a2,b,cdriver)
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
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: n, nx
    COMPLEX(DP), DEVICE, INTENT(IN) :: a1(nx)
    COMPLEX(DP), DEVICE, INTENT(IN) :: a2(nx)
    COMPLEX(DP), DEVICE, INTENT(OUT) :: b(dfft%nnr)
    CHARACTER(LEN=*), INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig,dfft_nnr
    !
    dfft_nnr = dfft%nnr
    !
    !$acc parallel loop
    DO ig = 1,dfft_nnr
       b(ig) = z_0
    ENDDO
    !$acc end parallel
    !
    !$acc parallel loop
    DO ig = 1,n
       b(dfft_nl_d(ig)) = a1(ig)+z_i*a2(ig)
       b(dfft_nlm_d(ig)) = CONJG(a1(ig)-z_i*a2(ig))
    ENDDO
    !$acc end parallel
    !
    CALL invfft(cdriver,b,dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE double_fwfft_gamma_gpu(dfft,n,nx,a,b1,b2,cdriver)
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
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: n, nx
    COMPLEX(DP), DEVICE, INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP), DEVICE, INTENT(OUT) :: b1(nx)
    COMPLEX(DP), DEVICE, INTENT(OUT) :: b2(nx)
    CHARACTER(LEN=*), INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    COMPLEX(DP) :: fm,fp
    !
    CALL fwfft(cdriver,a,dfft)
    !
    ! Keep only G>=0
    !
    !$acc parallel loop
    DO ig = 1,n
       fp = (a(dfft_nl_d(ig))+a(dfft_nlm_d(ig)))*0.5_DP
       fm = (a(dfft_nl_d(ig))-a(dfft_nlm_d(ig)))*0.5_DP
       b1(ig) = CMPLX(REAL(fp,KIND=DP),AIMAG(fm),KIND=DP)
       b2(ig) = CMPLX(AIMAG(fp),-REAL(fm,KIND=DP),KIND=DP)
    ENDDO
    !$acc end parallel
    !
    IF(nx > n) THEN
       !$acc parallel loop
       DO ig = n+1,nx
          b1(ig) = z_0
          b2(ig) = z_0
       ENDDO
       !$acc end parallel
    ENDIF
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_invfft_gamma_gpu(dfft,n,nx,a1,b,cdriver)
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
    !
    ! Workspace
    !
    INTEGER :: ig,dfft_nnr
    !
    dfft_nnr = dfft%nnr
    !
    !$acc parallel loop
    DO ig = 1,dfft_nnr
       b(ig) = z_0
    ENDDO
    !$acc end parallel
    !
    !$acc parallel loop
    DO ig = 1,n
       b(dfft_nl_d(ig)) = a1(ig)
       b(dfft_nlm_d(ig)) = CONJG(a1(ig))
    ENDDO
    !$acc end parallel
    !
    CALL invfft(cdriver,b,dfft)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE single_fwfft_gamma_gpu(dfft,n,nx,a,b1,cdriver)
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
    INTEGER, INTENT(IN) :: n, nx
    COMPLEX(DP), DEVICE, INTENT(INOUT) :: a(dfft%nnr)
    COMPLEX(DP), DEVICE, INTENT(OUT) :: b1(nx)
    CHARACTER(LEN=*), INTENT(IN) :: cdriver
    !
    ! Workspace
    !
    INTEGER :: ig
    !
    CALL fwfft(cdriver,a,dfft)
    !
    ! Keep only G>=0
    !
    !$acc parallel loop
    DO ig = 1,n
       b1(ig) = a(dfft_nl_d(ig))
    ENDDO
    !$acc end parallel
    !
    IF(nx > n) THEN
       !$acc parallel loop
       DO ig = n+1,nx
          b1(ig) = z_0
       ENDDO
       !$acc end parallel
    ENDIF
    !
  END SUBROUTINE
  !
  !
#endif
END MODULE
