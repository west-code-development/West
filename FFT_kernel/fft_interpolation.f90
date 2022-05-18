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
! He Ma, Marco Govoni, Han Yang
!
MODULE fourier_interpolation
 ! -----------------------------------------------------------------
 !
 IMPLICIT NONE
 !
 CONTAINS
 !
 SUBROUTINE set_nl (dfft, ng, ngx, ndim, nl, igk)
    !
    ! INPUT  : dfft       = FFT descriptor
    !          ng, ngx    = actual number of PW, max
    !          ndim       = dimensions of nl, 1: (n), 2:(n,-n)
    ! OUTPUT : nl         = computed mapping from G to R space (i,e. from [1,n] to [1, n1*n2*n3] )
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : g
    USE cell_base,     ONLY : at
    USE fft_types,     ONLY : fft_type_descriptor
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN)  :: ng, ngx, ndim
    INTEGER, INTENT(OUT) :: nl(ndim,ngx) ! mapping from [1, npw] to [1, n1*n2*n3]
    INTEGER,INTENT(IN),OPTIONAL :: igk(ng)
    !
    ! Workspace
    !
    INTEGER  :: ig, i_dim
    INTEGER  :: h, k, l, hmax, kmax, lmax, hidx, kidx, lidx
    !
    nl = 0
    !
    ! Maximum Miller indexes associated to the R grid
    !
    hmax = (dfft%nr1 - 1) / 2
    kmax = (dfft%nr2 - 1) / 2
    lmax = (dfft%nr3 - 1) / 2
    !
    DO i_dim = 1, ndim
       DO ig = 1, ng
          !
          ! Get the current Miller index
          !
          IF(PRESENT(igk)) THEN
             SELECT CASE (i_dim)
             CASE(1)
                h = NINT(SUM(g(:,igk(ig)) * at(:,1)))
                k = NINT(SUM(g(:,igk(ig)) * at(:,2)))
                l = NINT(SUM(g(:,igk(ig)) * at(:,3)))
             CASE(2)
                h = NINT(SUM(-g(:,igk(ig)) * at(:,1)))
                k = NINT(SUM(-g(:,igk(ig)) * at(:,2)))
                l = NINT(SUM(-g(:,igk(ig)) * at(:,3)))
             CASE DEFAULT
                CALL errore("","Bad mapping given: only 1 .OR. 2",1)
             END SELECT
          ELSE
             SELECT CASE (i_dim)
             CASE(1)
                h = NINT(SUM(g(:,ig) * at(:,1)))
                k = NINT(SUM(g(:,ig) * at(:,2)))
                l = NINT(SUM(g(:,ig) * at(:,3)))
             CASE(2)
                h = NINT(SUM(-g(:,ig) * at(:,1)))
                k = NINT(SUM(-g(:,ig) * at(:,2)))
                l = NINT(SUM(-g(:,ig) * at(:,3)))
             CASE DEFAULT
                CALL errore("","Bad mapping given: only 1 .OR. 2",1)
             END SELECT
          ENDIF
          !
          IF( ABS(h) <= hmax .AND. ABS(k) <= kmax .AND. ABS(l) <= lmax ) THEN
             !
             ! derive position of G vector in (n1, n2, n3) grid
             !
             IF ( h >= 0 ) THEN
                hidx = h + 1
             ELSE
                hidx = h + dfft%nr1 + 1
             ENDIF
             !
             IF ( k >= 0 ) THEN
                kidx = k + 1
             ELSE
                kidx = k + dfft%nr2 + 1
             ENDIF
             !
             IF ( l >= 0 ) THEN
                lidx = l + 1
             ELSE
                lidx = l + dfft%nr3 + 1
             ENDIF
             !
             IF (hidx>dfft%nr1 .OR. kidx>dfft%nr2 .OR. lidx>dfft%nr3) CALL errore('setnl','Mesh too small?',ig)
             !
#if defined(__MPI) && ! defined(__USE_3D_FFT)
             nl(i_dim,ig) = lidx + ( dfft%isind (hidx + (kidx - 1) * dfft%nr1x) - 1) * dfft%nr3x
#else
             nl(i_dim,ig) = hidx + (kidx - 1) * dfft%nr1x + (lidx - 1) * dfft%nr1x * dfft%nr2x
#endif
          ENDIF
          !
       ENDDO
    ENDDO
    !
 END SUBROUTINE
 !
 SUBROUTINE single_interp_invfft_k(dfft,n,nx,a1,b,cdriver,nl)
   !
   USE kinds,                 ONLY : DP
   USE fft_interfaces,        ONLY : fwfft,invfft
   USE fft_types,             ONLY : fft_type_descriptor
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
   INTEGER,INTENT(IN) :: nl(1,nx)
   !
   ! Workspace
   !
   INTEGER :: ig
   !
   b=0.0_DP
   !
   DO ig=1,n
      b(nl(1,ig))=a1(ig)
   ENDDO
   !
   CALL invfft(cdriver,b,dfft)
   !
 END SUBROUTINE
 !
 SUBROUTINE single_interp_fwfft_k(dfft,n,nx,a,b1,cdriver,nl)
   !
   USE kinds,                 ONLY : DP
   USE fft_interfaces,        ONLY : fwfft,invfft
   USE fft_types,             ONLY : fft_type_descriptor
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
   INTEGER,INTENT(IN) :: nl(1,nx)
   !
   ! Workspace
   !
   INTEGER :: ig
   !
   CALL fwfft(cdriver, a, dfft)
   !
   b1=0.0_DP
   !
   DO ig=1,n
      b1(ig) = a(nl(1,ig))
   ENDDO
   !
 END SUBROUTINE
 !
 SUBROUTINE single_interp_invfft_gamma(dfft,n,nx,a1,b,cdriver,nl)
   !
   USE kinds,                 ONLY : DP
   USE fft_interfaces,        ONLY : fwfft,invfft
   USE fft_types,             ONLY : fft_type_descriptor
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
   INTEGER,INTENT(IN) :: nl(2,nx)
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
      b(nl(1,ig))=       a1(ig)
      b(nl(2,ig))= CONJG(a1(ig))
   ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
   !
   CALL invfft(cdriver, b, dfft)
   !
 END SUBROUTINE
 !
 SUBROUTINE single_interp_fwfft_gamma(dfft,n,nx,a,b1,cdriver,nl)
   !
   USE kinds,                 ONLY : DP
   USE fft_interfaces,        ONLY : fwfft,invfft
   USE fft_types,             ONLY : fft_type_descriptor
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
   INTEGER,INTENT(IN) :: nl(2,nx)
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
      b1(ig) = a(nl(1,ig))
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
END MODULE
