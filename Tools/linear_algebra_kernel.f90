!
! Copyright (C) 2015-2022 M. Govoni
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
  USE kinds, ONLY : DP,i8b
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
      ! Diagonalization of a REAL(DP) SYMMETRIC MATRIX
      !
#if defined(__CUDA)
      USE west_gpu, ONLY : cusolv_h,cusolv_p
      USE cublas
      USE cusolverdn
#endif
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
      INTEGER :: lwork, info, info_d
#if defined(__CUDA)
      INTEGER(i8b) :: n8, lh8, ld8
      REAL(DP),ALLOCATABLE :: work_h(:)
      INTEGER :: ev_mode
#else
      CHARACTER :: ev_mode
#endif
      REAL(DP),ALLOCATABLE :: work(:)
      !$acc declare device_resident(work)
      REAL(DP),ALLOCATABLE :: m(:,:)
      !
#if defined(__CUDA)
      IF(l_just_ev) THEN
         ev_mode = CUSOLVER_EIG_MODE_NOVECTOR
      ELSE
         ev_mode = CUSOLVER_EIG_MODE_VECTOR
      ENDIF
      !
      !$acc data copyin(a) copyout(e,info_d)
      !
#if defined(__PGI) && (__PGIC__ > 22 || (__PGIC__ == 22 && __PGIC_MINOR__ > 7))
      n8 = INT(n,KIND=i8b)
      !
      !$acc host_data use_device(a,e,info_d)
      info = cusolverDnXsyevd_buffersize(cusolv_h,cusolv_p,ev_mode,CUBLAS_FILL_MODE_UPPER,n8,&
           & cudaDataType(CUDA_R_64F),a,n8,cudaDataType(CUDA_R_64F),e,cudaDataType(CUDA_R_64F),&
           & ld8,lh8)
      !
      ALLOCATE(work(ld8/8))
      ALLOCATE(work_h(lh8/8))
      !
      info = cusolverDnXsyevd(cusolv_h,cusolv_p,ev_mode,CUBLAS_FILL_MODE_UPPER,n8,&
           & cudaDataType(CUDA_R_64F),a,n8,cudaDataType(CUDA_R_64F),e,cudaDataType(CUDA_R_64F),&
           & work,ld8,work_h,lh8,info_d)
      !$acc end host_data
      !
      DEALLOCATE(work)
      DEALLOCATE(work_h)
#else
      !$acc host_data use_device(a,e,info_d)
      info = cusolverDnDsyevd_bufferSize(cusolv_h,ev_mode,CUBLAS_FILL_MODE_UPPER,n,a,n,e,lwork)
      !$acc end host_data
      !
      ALLOCATE(work(lwork))
      !
      !$acc host_data use_device(a,e,work,info_d)
      info = cusolverDnDsyevd(cusolv_h,ev_mode,CUBLAS_FILL_MODE_UPPER,n,a,n,e,work,lwork,info_d)
      !$acc end host_data
      !
      DEALLOCATE(work)
#endif
      !
      IF(.NOT. l_just_ev) THEN
         !$acc update host(a)
      ENDIF
      !$acc end data
#else
      IF(l_just_ev) THEN
         ALLOCATE(m(n,n))
         m(:,:) = a
         !
         ev_mode = 'N'
      ELSE
         ev_mode = 'V'
      ENDIF
      !
      ALLOCATE(work(1))
      !
      CALL DSYEV(ev_mode,'U',n,a,n,e,work,-1,info)
      !
      lwork = CEILING(work(1))
      !
      DEALLOCATE(work)
      ALLOCATE(work(lwork))
      !
      CALL DSYEV(ev_mode,'U',n,a,n,e,work,lwork,info)
      !
      DEALLOCATE(work)
      !
      IF(l_just_ev) THEN
         a(:,:) = m
         DEALLOCATE(m)
      ENDIF
#endif
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matdiago_zhe(n,a,e,l_just_ev)
    !-----------------------------------------------------------------------
      !
      ! Diagonalization of a COMPLEX(DP) HERMITIAN MATRIX
      !
#if defined(__CUDA)
      USE west_gpu, ONLY : cusolv_h,cusolv_p
      USE cublas
      USE cusolverdn
#endif
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
      INTEGER :: lwork, info, info_d
#if defined(__CUDA)
      INTEGER(i8b) :: n8, lh8, ld8
      COMPLEX(DP),ALLOCATABLE :: work_h(:)
      INTEGER :: ev_mode
#else
      CHARACTER :: ev_mode
#endif
      REAL(DP),ALLOCATABLE :: rwork(:)
      COMPLEX(DP),ALLOCATABLE :: work(:)
      !$acc declare device_resident(work)
      COMPLEX(DP),ALLOCATABLE :: m(:,:)
      !
#if defined(__CUDA)
      IF(l_just_ev) THEN
         ev_mode = CUSOLVER_EIG_MODE_NOVECTOR
      ELSE
         ev_mode = CUSOLVER_EIG_MODE_VECTOR
      ENDIF
      !
      !$acc data copyin(a) copyout(e,info_d)
      !
#if defined(__PGI) && (__PGIC__ > 22 || (__PGIC__ == 22 && __PGIC_MINOR__ > 7))
      n8 = INT(n,KIND=i8b)
      !
      !$acc host_data use_device(a,e,info_d)
      info = cusolverDnXsyevd_buffersize(cusolv_h,cusolv_p,ev_mode,CUBLAS_FILL_MODE_UPPER,n8,&
           & cudaDataType(CUDA_C_64F),a,n8,cudaDataType(CUDA_R_64F),e,cudaDataType(CUDA_C_64F),&
           & ld8,lh8)
      !
      ALLOCATE(work(ld8/16))
      ALLOCATE(work_h(lh8/16))
      !
      info = cusolverDnXsyevd(cusolv_h,cusolv_p,ev_mode,CUBLAS_FILL_MODE_UPPER,n8,&
           & cudaDataType(CUDA_C_64F),a,n8,cudaDataType(CUDA_R_64F),e,cudaDataType(CUDA_C_64F),&
           & work,ld8,work_h,lh8,info_d)
      !$acc end host_data
      !
      DEALLOCATE(work)
      DEALLOCATE(work_h)
#else
      !$acc host_data use_device(a,e,info_d)
      info = cusolverDnZheevd_bufferSize(cusolv_h,ev_mode,CUBLAS_FILL_MODE_UPPER,n,a,n,e,lwork)
      !$acc end host_data
      !
      ALLOCATE(work(lwork))
      !
      !$acc host_data use_device(a,e,work,info_d)
      info = cusolverDnZheevd(cusolv_h,ev_mode,CUBLAS_FILL_MODE_UPPER,n,a,n,e,work,lwork,info_d)
      !$acc end host_data
      !
      DEALLOCATE(work)
#endif
      !
      IF(.NOT. l_just_ev) THEN
         !$acc update host(a)
      ENDIF
      !$acc end data
#else
      IF(l_just_ev) THEN
         ALLOCATE(m(n,n))
         m(:,:) = a
         !
         ev_mode = 'N'
      ELSE
         ev_mode = 'V'
      ENDIF
      !
      ALLOCATE(work(1))
      ALLOCATE(rwork(3*n-2))
      !
      CALL ZHEEV(ev_mode,'U',n,a,n,e,work,-1,rwork,info)
      !
      lwork = CEILING(REAL(work(1),KIND=DP))
      !
      DEALLOCATE(work)
      ALLOCATE(work(lwork))
      !
      CALL ZHEEV(ev_mode,'U',n,a,n,e,work,lwork,rwork,info)
      !
      DEALLOCATE(rwork)
      DEALLOCATE(work)
      !
      IF(l_just_ev) THEN
         a(:,:) = m
         DEALLOCATE(m)
      ENDIF
#endif
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_zge(n,a)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a COMPLEX(DP) GENERIC MATRIX
      !
#if defined(__CUDA)
      USE west_gpu, ONLY : cusolv_h,cusolv_p,l_inv,l_inv_h,piv,work_c,work_c_h
      USE cublas
      USE cusolverdn
#endif
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
      INTEGER(i8b) :: n8
      INTEGER,EXTERNAL :: ilaenv
      INTEGER :: nb, lwork, info, info_d, i1, i2, i3
      INTEGER,ALLOCATABLE :: ipiv(:)
      COMPLEX(DP) :: tmp
      COMPLEX(DP),ALLOCATABLE :: work(:)
      !
#if defined(__CUDA)
      n8 = INT(n,KIND=i8b)
      !
      !$acc data create(info_d)
      !$acc host_data use_device(a,piv,work_c,info_d)
#if defined(__PGI) && (__PGIC__ > 22 || (__PGIC__ == 22 && __PGIC_MINOR__ > 7))
      info = cusolverDnXgetrf(cusolv_h,cusolv_p,n8,n8,cudaDataType(CUDA_C_64F),a,n8,piv,&
           & cudaDataType(CUDA_C_64F),work_c,l_inv,work_c_h,l_inv_h,info_d)
#else
      info = cusolverDnZgetrf(cusolv_h,n,n,a,n,work_c,piv,info_d)
#endif
      !$acc end host_data
      !
      !$acc kernels present(work_c)
      work_c(1:n**2) = (0._DP,0._DP)
      !$acc end kernels
      !
      !$acc parallel loop present(work_c)
      DO i1 = 1,n
         work_c(i1+(i1-1)*n) = (1._DP,0._DP)
      ENDDO
      !$acc end parallel
      !
      !$acc host_data use_device(a,piv,work_c,info_d)
#if defined(__PGI) && (__PGIC__ > 22 || (__PGIC__ == 22 && __PGIC_MINOR__ > 7))
      info = cusolverDnXgetrs(cusolv_h,cusolv_p,CUBLAS_OP_N,n8,n8,cudaDataType(CUDA_C_64F),a,&
           & n8,piv,cudaDataType(CUDA_C_64F),work_c(1:n**2),n8,info_d)
#else
      info = cusolverDnZgetrs(cusolv_h,CUBLAS_OP_N,n,n,a,n,piv,work_c(1:n**2),n,info_d)
#endif
      !$acc end host_data
      !$acc end data
      !
      !$acc parallel loop collapse(2) present(a,work_c)
      DO i2 = 1,n
         DO i1 = 1,n
            a(i1,i2) = work_c(i1+(i2-1)*n)
         ENDDO
      ENDDO
      !$acc end parallel
#else
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
#endif
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE matinvrs_dge(n,a)
    !-----------------------------------------------------------------------
      !
      ! Inversion of a REAL(DP) GENERIC MATRIX
      !
#if defined(__CUDA)
      USE west_gpu, ONLY : cusolv_h,cusolv_p,l_inv,l_inv_h,piv,work_r,work_r_h
      USE cublas
      USE cusolverdn
#endif
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
      INTEGER(i8b) :: n8
      INTEGER,EXTERNAL :: ilaenv
      INTEGER :: nb, lwork, info, info_d, i1, i2, i3
      INTEGER,ALLOCATABLE :: ipiv(:)
      REAL(DP) :: tmp
      REAL(DP),ALLOCATABLE :: work(:)
      !
#if defined(__CUDA)
      n8 = INT(n,KIND=i8b)
      !
      !$acc data copyout(info_d)
      !$acc host_data use_device(a,piv,work_r,info_d)
#if defined(__PGI) && (__PGIC__ > 22 || (__PGIC__ == 22 && __PGIC_MINOR__ > 7))
      info = cusolverDnXgetrf(cusolv_h,cusolv_p,n8,n8,cudaDataType(CUDA_R_64F),a,n8,piv,&
           & cudaDataType(CUDA_R_64F),work_r,l_inv,work_r_h,l_inv_h,info_d)
#else
      info = cusolverDnDgetrf(cusolv_h,n,n,a,n,work_r,piv,info_d)
#endif
      !$acc end host_data
      !
      !$acc kernels present(work_r)
      work_r(1:n**2) = 0._DP
      !$acc end kernels
      !
      !$acc parallel loop present(work_r)
      DO i1 = 1,n
         work_r(i1+(i1-1)*n) = 1._DP
      ENDDO
      !$acc end parallel
      !
      !$acc host_data use_device(a,piv,work_r,info_d)
#if defined(__PGI) && (__PGIC__ > 22 || (__PGIC__ == 22 && __PGIC_MINOR__ > 7))
      info = cusolverDnXgetrs(cusolv_h,cusolv_p,CUBLAS_OP_N,n8,n8,cudaDataType(CUDA_R_64F),a,&
           & n8,piv,cudaDataType(CUDA_R_64F),work_r(1:n**2),n8,info_d)
#else
      info = cusolverDnDgetrs(cusolv_h,CUBLAS_OP_N,n,n,a,n,piv,work_r(1:n**2),n,info_d)
#endif
      !$acc end host_data
      !$acc end data
      !
      !$acc parallel loop collapse(2) present(a,work_r)
      DO i2 = 1,n
         DO i1 = 1,n
            a(i1,i2) = work_r(i1+(i2-1)*n)
         ENDDO
      ENDDO
      !$acc end parallel
#else
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
#endif
      !
    END SUBROUTINE
    !
END MODULE
