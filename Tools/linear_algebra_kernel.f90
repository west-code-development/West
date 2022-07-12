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
MODULE linear_algebra_kernel
  !-----------------------------------------------------------------------
  !
  ! Serial linear algebra toolkit
  !
  USE kinds, ONLY : DP
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
      DO i2 = 1,n
         DO i1 = 1,i2
            i3 = ((i2-1)*i2)/2
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
      DO i2 = 1,n
         DO i1 = 1,i2
            i3 = ((i2-1)*i2)/2
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
      ! Diago of a REAL(DP) SYMMETRIC MATRIX
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
         m = a
         !
         ALLOCATE(work(1))
         !
         CALL DSYEV('N','U',n,m,n,e,work,-1,info)
         !
         lwork = CEILING(work(1))
         !
         DEALLOCATE(work)
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
         lwork = CEILING(work(1))
         !
         DEALLOCATE(work)
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
      ! Diagonalization of a COMPLEX(DP) HERMITIAN MATRIX
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: n
      COMPLEX(DP),INTENT(INOUT) :: a(n,n)
      REAL(DP),INTENT(OUT) :: e(n)
      LOGICAL,INTENT(IN) :: l_just_ev
      !
      ! Workspace
      !
      INTEGER :: lwork, info
      REAL(DP),ALLOCATABLE :: rwork(:)
      COMPLEX(DP),ALLOCATABLE :: m(:,:), work(:)
      !
      IF(l_just_ev) THEN
         !
         ALLOCATE(m(n,n))
         m = a
         !
         ALLOCATE(work(1))
         ALLOCATE(rwork(3*n-2))
         !
         CALL ZHEEV('N','U',n,m,n,e,work,-1,rwork,info)
         !
         lwork = CEILING(REAL(work(1),DP))
         !
         DEALLOCATE(work)
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
         lwork = CEILING(REAL(work(1),DP))
         !
         DEALLOCATE(work)
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
    SUBROUTINE matinvrs_zge(n,a)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a COMPLEX(DP) GENERIC MATRIX
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
      INTEGER :: nb, lwork, info
      INTEGER,ALLOCATABLE :: ipiv(:)
      COMPLEX(DP),ALLOCATABLE :: work(:)
      !
      ! Calculate optimal size of workspace
      !
      nb = ilaenv(1,'ZGETRI','N',n,-1,-1,-1)
      lwork = n*nb
      ALLOCATE(work(lwork))
      ALLOCATE(ipiv(n))
      !
      ! Invert matrix
      !
      CALL ZGETRF(n,n,a,n,ipiv,info)
      CALL ZGETRI(n,a,n,ipiv,work,lwork,info)
      !
      DEALLOCATE(work)
      DEALLOCATE(ipiv)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_dge(n,a)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a REAL(DP) GENERIC MATRIX
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
      INTEGER :: nb, lwork, info
      INTEGER,ALLOCATABLE :: ipiv(:)
      REAL(DP),ALLOCATABLE :: work(:)
      !
      ! Calculate optimal size of workspace
      !
      nb = ilaenv(1,'DGETRI','N',n,-1,-1,-1)
      lwork = n*nb
      ALLOCATE(work(lwork))
      ALLOCATE(ipiv(n))
      !
      ! Invert matrix
      !
      CALL DGETRF(n,n,a,n,ipiv,info)
      CALL DGETRI(n,a,n,ipiv,work,lwork,info)
      !
      DEALLOCATE(work)
      DEALLOCATE(ipiv)
      !
    END SUBROUTINE
    !
#if defined(__CUDA)
    !-----------------------------------------------------------------------
    SUBROUTINE matdiago_dsy_gpu(n,a,e,l_just_ev)
    !-----------------------------------------------------------------------
      !
      ! Diago of a REAL(DP) SYMMETRIC MATRIX
      !
      USE west_gpu, ONLY : cusolver_h
      USE cublas
      USE cusolverdn
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
      INTEGER :: lwork
      INTEGER :: info
      INTEGER :: info_d
      REAL(DP),ALLOCATABLE :: work(:)
      !$acc declare device_resident(work)
      !
      IF(l_just_ev) THEN
         !$acc data copyin(a) copyout(e,info_d)
         !$acc host_data use_device(a,e,info_d)
         info = cusolverDnDsyevd_bufferSize(cusolver_h,&
         & CUSOLVER_EIG_MODE_NOVECTOR,CUBLAS_FILL_MODE_UPPER,n,a,n,e,lwork)
         !
         ALLOCATE(work(lwork))
         !
         info = cusolverDnDsyevd(cusolver_h,CUSOLVER_EIG_MODE_NOVECTOR,&
         & CUBLAS_FILL_MODE_UPPER,n,a,n,e,work,lwork,info_d)
         !$acc end host_data
         !$acc end data
         !
         DEALLOCATE(work)
      ELSE
         !$acc data copy(a) copyout(e,info_d)
         !$acc host_data use_device(a,e,info_d)
         info = cusolverDnDsyevd_bufferSize(cusolver_h,&
         & CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,n,a,n,e,lwork)
         !
         ALLOCATE(work(lwork))
         !
         info = cusolverDnDsyevd(cusolver_h,CUSOLVER_EIG_MODE_VECTOR,&
         & CUBLAS_FILL_MODE_UPPER,n,a,n,e,work,lwork,info_d)
         !$acc end host_data
         !$acc end data
         !
         DEALLOCATE(work)
      ENDIF
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matdiago_zhe_gpu(n,a,e,l_just_ev)
    !-----------------------------------------------------------------------
      !
      ! Diagonalization of a COMPLEX(DP) HERMITIAN MATRIX
      !
      USE west_gpu, ONLY : cusolver_h
      USE cublas
      USE cusolverdn
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: n
      COMPLEX(DP),INTENT(INOUT) :: a(n,n)
      REAL(DP),INTENT(OUT) :: e(n)
      LOGICAL,INTENT(IN) :: l_just_ev
      !
      ! Workspace
      !
      INTEGER :: lwork
      INTEGER :: info
      INTEGER :: info_d
      COMPLEX(DP),ALLOCATABLE :: work(:)
      !$acc declare device_resident(work)
      !
      IF(l_just_ev) THEN
         !$acc data copyin(a) copyout(e,info_d)
         !$acc host_data use_device(a,e,info_d)
         info = cusolverDnZheevd_bufferSize(cusolver_h,&
         & CUSOLVER_EIG_MODE_NOVECTOR,CUBLAS_FILL_MODE_UPPER,n,a,n,e,lwork)
         !
         ALLOCATE(work(lwork))
         !
         info = cusolverDnZheevd(cusolver_h,CUSOLVER_EIG_MODE_NOVECTOR,&
         & CUBLAS_FILL_MODE_UPPER,n,a,n,e,work,lwork,info_d)
         !$acc end host_data
         !$acc end data
         !
         DEALLOCATE(work)
      ELSE
         !$acc data copy(a) copyout(e,info_d)
         !$acc host_data use_device(a,e,info_d)
         info = cusolverDnZheevd_bufferSize(cusolver_h,&
         & CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,n,a,n,e,lwork)
         !
         ALLOCATE(work(lwork))
         !
         info = cusolverDnZheevd(cusolver_h,CUSOLVER_EIG_MODE_VECTOR,&
         & CUBLAS_FILL_MODE_UPPER,n,a,n,e,work,lwork,info_d)
         !$acc end host_data
         !$acc end data
         !
         DEALLOCATE(work)
      ENDIF
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_zge_gpu(n,a_d)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a COMPLEX(DP) GENERIC MATRIX
      !
      USE west_gpu, ONLY : cusolver_h,work_c,piv_d
      USE cublas
      USE cusolverdn
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: n
      COMPLEX(DP),DEVICE,INTENT(INOUT) :: a_d(n,n)
      !
      ! Workspace
      !
      INTEGER :: info
      INTEGER :: i1
      INTEGER :: i2
      INTEGER :: info_d
      !
      ! Note piv_d => NULL
      !
      !$acc data copyout(info_d)
      !$acc host_data use_device(work_c,info_d)
      info = cusolverDnZgetrf(cusolver_h,n,n,a_d,n,work_c,piv_d,info_d)
      info = cusolverDnZtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n,a_d,&
      & n,work_c,SIZE(work_c),info_d)
      !$acc end host_data
      !$acc end data
      !
      !$acc parallel loop collapse(2) present(work_c)
      DO i2 = 1,n-1
         DO i1 = 2,n
            IF(i1 > i2) THEN
               work_c(i1+(i2-1)*n) = a_d(i1,i2)
               a_d(i1,i2) = (0._DP,0._DP)
            ENDIF
         ENDDO
      ENDDO
      !$acc end parallel
      !
      !$acc host_data use_device(work_c)
      CALL ZTRSM('R','L','N','U',n,n,(1._DP,0._DP),work_c,n,a_d,n)
      !$acc end host_data
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_dge_gpu(n,a_d)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a REAL(DP) GENERIC MATRIX
      !
      USE west_gpu, ONLY : cusolver_h,work_r,piv_d
      USE cublas
      USE cusolverdn
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: n
      REAL(DP),DEVICE,INTENT(INOUT) :: a_d(n,n)
      !
      ! Workspace
      !
      INTEGER :: info
      INTEGER :: i1
      INTEGER :: i2
      INTEGER :: info_d
      !
      ! Note piv_d => NULL
      !
      !$acc data copyout(info_d)
      !$acc host_data use_device(work_r,info_d)
      info = cusolverDnDgetrf(cusolver_h,n,n,a_d,n,work_r,piv_d,info_d)
      info = cusolverDnDtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n,a_d,&
      & n,work_r,SIZE(work_r),info_d)
      !$acc end host_data
      !$acc end data
      !
      !$acc parallel loop collapse(2) present(work_r)
      DO i2 = 1,n-1
         DO i1 = 2,n
            IF(i1 > i2) THEN
               work_r(i1+(i2-1)*n) = a_d(i1,i2)
               a_d(i1,i2) = 0._DP
            ENDIF
         ENDDO
      ENDDO
      !$acc end parallel
      !
      !$acc host_data use_device(work_r)
      CALL DTRSM('R','L','N','U',n,n,1._DP,work_r,n,a_d,n)
      !$acc end host_data
      !
    END SUBROUTINE
#endif
    !
END MODULE
