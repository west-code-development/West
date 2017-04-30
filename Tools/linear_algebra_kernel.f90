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
MODULE linear_algebra_kernel
  !-----------------------------------------------------------------------
  !
  ! Serial Algebra toolkit
  !
  USE kinds,                ONLY : DP 
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE D_SU2SP(n,np,a,ap)
    !-----------------------------------------------------------------------
      !
      ! Real symmetric matrix : from Symm up to packed storage 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER,INTENT(IN) :: n, np
      REAL(DP),INTENT(IN) :: a(n,n)
      REAL(DP),INTENT(OUT) :: ap(np) 
      !
      ! Workspace
      !
      INTEGER :: i1, i2, i3 
      !
      ap = 0._DP
      !
!$OMP PARALLEL SHARED(n,a,ap) PRIVATE(i2,i1,i3)
!$OMP DO
      DO i2 = 1, n
         DO i1 = 1, i2
            i3 = ( (i2-1) * i2 ) / 2 
            ap(i3+i1) = a(i1,i2)
         ENDDO
      ENDDO 
!$OMP END DO
!$OMP END PARALLEL
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE D_SP2SU(n,np,a,ap)
    !-----------------------------------------------------------------------
      !
      ! Real symmetric matrix : to Symm up from packed storage 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER,INTENT(IN) :: n, np
      REAL(DP),INTENT(OUT) :: a(n,n)
      REAL(DP),INTENT(IN) :: ap(np) 
      !
      ! Workspace
      !
      INTEGER :: i1, i2, i3 
      !
      a = 0._DP 
      !
!$OMP PARALLEL SHARED(n,a,ap) PRIVATE(i2,i1,i3)
!$OMP DO
      DO i2 = 1, n
         DO i1 = 1, i2
            i3 = ( (i2-1) * i2 ) / 2 
            a(i1,i2) = ap(i3+i1)
         ENDDO
      ENDDO 
!$OMP END DO
!$OMP END PARALLEL
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matdiago_dsy(n,a,e,l_just_ev)
    !-----------------------------------------------------------------------
      !
      ! Diago of a REAL(DP) SYMMETRIC MATRIX. 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER,INTENT(IN) :: n 
      REAL(DP),INTENT(INOUT) :: a(n,n)
      REAL(DP),INTENT(OUT) :: e(n) 
      LOGICAL,INTENT(IN) :: l_just_ev
      !
      ! Workspace
      !
      INTEGER :: lwork, info
      REAL(DP),ALLOCATABLE :: work(:)
      REAL(DP),ALLOCATABLE :: m(:,:)
      !
      IF(l_just_ev) THEN
         !
         ALLOCATE(m(n,n))
         m=a
         !
         ALLOCATE(work(1))
         !
         CALL DSYEV('N','U',n,m,n,e,work,-1,info)
         !
         lwork=INT(REAL(work(1),DP))
         DEALLOCATE(work)
         !
         ALLOCATE(work(lwork))
         !
         CALL DSYEV('N','U',n,m,n,e,work,lwork,info)
         !
         DEALLOCATE(work)
         DEALLOCATE(m)
         !
      ELSE
         !
         ALLOCATE(work(1))
         !
         CALL DSYEV('V','U',n,a,n,e,work,-1,info)
         !
         lwork=INT(REAL(work(1),DP))
         DEALLOCATE(work)
         !
         ALLOCATE(work(lwork))
         !
         CALL DSYEV('V','U',n,a,n,e,work,lwork,info)
         !
         DEALLOCATE(work)
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matdiago_zhe(n,a,e,l_just_ev)
    !-----------------------------------------------------------------------
      !
      ! Diagonalization of a COMPLEX(DP) HERMITIAN MATRIX. 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER, INTENT(IN) :: n 
      COMPLEX(DP), INTENT(INOUT) :: a(n,n)
      REAL(DP), INTENT(OUT) :: e(n) 
      LOGICAL, INTENT(IN) :: l_just_ev
      !
      ! Workspace
      !
      INTEGER :: lwork, info
      REAL(DP), ALLOCATABLE :: rwork(:), w(:)
      COMPLEX(DP), ALLOCATABLE :: m(:,:), work(:)
      !
      IF(l_just_ev) THEN
         !
         ALLOCATE(m(n,n))
         m=a
         !
         ALLOCATE(work(1))
         ALLOCATE(rwork(3*n-2))
         !
         CALL ZHEEV('N','U',n,m,n,e,work,-1,rwork,info)
         !
         lwork=INT(REAL(work(1),DP))
         DEALLOCATE(work)
         !
         ALLOCATE(work(lwork))
         !
         CALL ZHEEV('N','U',n,m,n,e,work,lwork,rwork,info)
         !
         DEALLOCATE(rwork)
         DEALLOCATE(work)
         DEALLOCATE(m)
         !
      ELSE
         !
         ALLOCATE(work(1))
         ALLOCATE(rwork(3*n-2))
         !
         CALL ZHEEV('V','U',n,a,n,e,work,-1,rwork,info)
         !
         lwork=INT(REAL(work(1),DP))
         DEALLOCATE(work)
         !
         ALLOCATE(work(lwork))
         !
         CALL ZHEEV('V','U',n,a,n,e,work,lwork,rwork,info)
         !
         DEALLOCATE(work)
         DEALLOCATE(rwork)
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_dsy(n,a)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a REAL(DP) SYMMETRIC MATRIX. 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER,INTENT(IN) :: n 
      REAL(DP),INTENT(INOUT) :: a(n,n)
      !
      ! Workspace
      !
      INTEGER,EXTERNAL :: ilaenv
      INTEGER :: i, j, nb, lwork, info
      INTEGER,ALLOCATABLE :: ipiv(:)
      REAL(DP),ALLOCATABLE :: work(:)
      !
      ! Calculate the optimal size of the workspace array.
      !
      nb = ilaenv(1, "DSYTRI", "U", n, -1, -1, -1)
      lwork = n * nb
      ALLOCATE( work(lwork) )
      ALLOCATE( ipiv(n) )
      !
      ! Invert the matrix.
      !
      CALL DSYTRF("U", n, a, n, ipiv, work, lwork, info)
      !if (info /= 0) stop "error in call to dsytrf"
      CALL DSYTRI("U", n, a, n, ipiv, work, info)
      !if (info /= 0) stop "error in call to dsytri"
      DEALLOCATE (work)
      DEALLOCATE (ipiv)
      ! Copy the upper triangular part of A to the lower.
      DO j = 1, n - 1
         DO i = j + 1, n
            a(i, j) = a(j, i)
         ENDDO
      ENDDO
      !
    END SUBROUTINE
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_zge(n,a)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a COMPLEX(DP) GENERIC MATRIX. 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER,INTENT(IN) :: n 
      COMPLEX(DP),INTENT(INOUT) :: a(n,n)
      !
      ! Workspace
      !
      INTEGER,EXTERNAL :: ilaenv
      INTEGER :: i, j, nb, lwork, info
      INTEGER,ALLOCATABLE :: ipiv(:)
      COMPLEX(DP),ALLOCATABLE :: work(:)
      !
      ! Calculate the optimal size of the workspace array.
      !
      nb = ilaenv(1, "ZGETRI", "N", n, -1, -1, -1)
      lwork = n * nb
      ALLOCATE( work(lwork) )
      ALLOCATE( ipiv(n) )
      !
      ! Invert the matrix.
      !
      CALL ZGETRF( n, n, a, n, ipiv, info)
      !if (info /= 0) stop "error in call to dsytrf"
      CALL ZGETRI( n, a, n, ipiv, work, lwork, info )
      !if (info /= 0) stop "error in call to dsytri"
      DEALLOCATE( work )
      DEALLOCATE( ipiv )
      !
    END SUBROUTINE
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_dge(n,a)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a REAL(DP) GENERIC MATRIX. 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER,INTENT(IN) :: n 
      REAL(DP),INTENT(INOUT) :: a(n,n)
      !
      ! Workspace
      !
      INTEGER,EXTERNAL :: ilaenv
      INTEGER :: i, j, nb, lwork, info
      INTEGER,ALLOCATABLE :: ipiv(:)
      REAL(DP),ALLOCATABLE :: work(:)
      !
      ! Calculate the optimal size of the workspace array.
      !
      nb = ilaenv(1, "DGETRI", "N", n, -1, -1, -1)
      lwork = n * nb
      ALLOCATE( work(lwork) )
      ALLOCATE( ipiv(n) )
      !
      ! Invert the matrix.
      !
      CALL DGETRF( n, n, a, n, ipiv, info)
      !if (info /= 0) stop "error in call to dsytrf"
      CALL DGETRI( n, a, n, ipiv, work, lwork, info )
      !if (info /= 0) stop "error in call to dsytri"
      DEALLOCATE( work )
      DEALLOCATE( ipiv )
      !
    END SUBROUTINE
    !
END MODULE
