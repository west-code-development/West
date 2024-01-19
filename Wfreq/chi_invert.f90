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
MODULE chi_invert
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   CONTAINS
   !
   !-----------------------------------------------------------------------
   SUBROUTINE chi_invert_real(matilda,head,lambda,nma)
      !-----------------------------------------------------------------------
      !
      ! For each frequency I calculate X, ky, head and lambda
      !
      ! X = (1-B)^{-1}
      ! temph = X * wh
      ! templ = wl * X
      ! tempt = wl * temph
      !
      ! ky = 1 - f - Tr ( tempt ) / 3
      !
      ! Head = (1-ky) / ky
      !
      ! Lamba = X * B + temph * templ / (3*ky)
      !
      USE kinds,                 ONLY : DP
      USE linear_algebra_kernel, ONLY : matinvrs_dge
      USE westcom,               ONLY : n_pdep_eigen_to_use,l_macropol
#if defined(__CUDA)
      USE west_gpu,              ONLY : x_r,wh_r,wl_r,temph_r,templ_r,tempt_r
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      REAL(DP),INTENT(IN) :: matilda(nma,nma)
      REAL(DP),INTENT(OUT) :: head,lambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use)
      INTEGER,INTENT(IN) :: nma
      !
      ! Workspace
      !
      INTEGER :: i1,i2
      REAL(DP) :: f
      REAL(DP) :: ky
#if !defined(__CUDA)
      REAL(DP) :: tempt_r(3,3)
      REAL(DP),ALLOCATABLE :: wh_r(:,:)
      REAL(DP),ALLOCATABLE :: wl_r(:,:)
      REAL(DP),ALLOCATABLE :: x_r(:,:)
      REAL(DP),ALLOCATABLE :: temph_r(:,:)
      REAL(DP),ALLOCATABLE :: templ_r(:,:)
      !
      ALLOCATE(x_r(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      IF(l_macropol) THEN
         ALLOCATE(wh_r(n_pdep_eigen_to_use,3))
         ALLOCATE(wl_r(3,n_pdep_eigen_to_use))
         ALLOCATE(temph_r(n_pdep_eigen_to_use,3))
         ALLOCATE(templ_r(3,n_pdep_eigen_to_use))
      ENDIF
#endif
      !
      !$acc update device(matilda)
      !
      IF(l_macropol) THEN
         !
         f = 0._DP
         DO i1 = 1,3
            f = f+matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1)/3._DP
         ENDDO
         !
         !$acc parallel loop collapse(2) present(wh_r,matilda)
         DO i2 = 1,3
            DO i1 = 1,n_pdep_eigen_to_use
               wh_r(i1,i2) = matilda(i1,n_pdep_eigen_to_use+i2)
            ENDDO
         ENDDO
         !$acc end parallel
         !
         !$acc parallel loop collapse(2) present(wl_r,matilda)
         DO i2 = 1,n_pdep_eigen_to_use
            DO i1 = 1,3
               wl_r(i1,i2) = matilda(n_pdep_eigen_to_use+i1,i2)
            ENDDO
         ENDDO
         !$acc end parallel
         !
      ENDIF
      !
      ! X = (1-B)^{-1}
      !
      !$acc parallel loop collapse(2) present(x_r,matilda)
      DO i2 = 1,n_pdep_eigen_to_use
         DO i1 = 1,n_pdep_eigen_to_use
            IF(i1 == i2) THEN
               x_r(i1,i2) = 1._DP-matilda(i1,i2)
            ELSE
               x_r(i1,i2) = -matilda(i1,i2)
            ENDIF
         ENDDO
      ENDDO
      !$acc end parallel
      !
      CALL matinvrs_dge(n_pdep_eigen_to_use,x_r)
      !
      IF(l_macropol) THEN
         !
         ! temph = X * wh
         !
         !$acc host_data use_device(x_r,wh_r,wl_r,temph_r,templ_r,tempt_r)
         CALL DGEMM('N','N',n_pdep_eigen_to_use,3,n_pdep_eigen_to_use,1._DP,x_r,n_pdep_eigen_to_use,&
         & wh_r,n_pdep_eigen_to_use,0._DP,temph_r,n_pdep_eigen_to_use)
         !
         ! templ = wl * X
         !
         CALL DGEMM('N','N',3,n_pdep_eigen_to_use,n_pdep_eigen_to_use,1._DP,wl_r,3,x_r,&
         & n_pdep_eigen_to_use,0._DP,templ_r,3)
         !
         ! tempt = wl * temph
         !
         CALL DGEMM('N','N',3,3,n_pdep_eigen_to_use,1._DP,wl_r,3,temph_r,n_pdep_eigen_to_use,0._DP,&
         & tempt_r,3)
         !$acc end host_data
         !
         !$acc update host(tempt_r)
         !
         ! ky = 1 - f - Tr ( tempt ) / 3
         !
         ky = 1._DP-f
         DO i1 = 1,3
            ky = ky-tempt_r(i1,i1)/3._DP
         ENDDO
         !
         head = (1._DP-ky)/ky
         !
      ELSE
         head = 0._DP
      ENDIF
      !
      !$acc host_data use_device(x_r,matilda,temph_r,templ_r,lambda)
      CALL DGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,n_pdep_eigen_to_use,1._DP,x_r,&
      & n_pdep_eigen_to_use,matilda,nma,0._DP,lambda,n_pdep_eigen_to_use)
      !
      IF(l_macropol) THEN
         f = 1._DP/(3._DP*ky)
         CALL DGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,3,f,temph_r,n_pdep_eigen_to_use,&
         & templ_r,3,1._DP,lambda,n_pdep_eigen_to_use)
      ENDIF
      !$acc end host_data
      !
      !$acc update host(lambda)
      !
#if !defined(__CUDA)
      DEALLOCATE(x_r)
      IF(l_macropol) THEN
         DEALLOCATE(wh_r)
         DEALLOCATE(wl_r)
         DEALLOCATE(temph_r)
         DEALLOCATE(templ_r)
      ENDIF
#endif
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE chi_invert_complex(matilda,head,lambda,nma,l_gammaq)
      !-----------------------------------------------------------------------
      !
      ! For each frequency and q-point I calculate X, ky, head and lambda
      !
      ! X = (1-B)^{-1}
      ! temph = X * wh
      ! templ = wl * X
      ! tempt = wl * temph
      !
      ! ky = 1 - f - Tr ( tempt ) / 3
      !
      ! Head = (1-ky) / ky
      !
      ! Lamba = X * B + temph * templ / (3*ky)
      !
      USE kinds,                 ONLY : DP
      USE linear_algebra_kernel, ONLY : matinvrs_zge
      USE westcom,               ONLY : n_pdep_eigen_to_use,l_macropol
#if defined(__CUDA)
      USE west_gpu,              ONLY : x_c,wh_c,wl_c,temph_c,templ_c,tempt_c
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(IN) :: matilda(nma,nma)
      COMPLEX(DP),INTENT(OUT) :: head,lambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use)
      INTEGER,INTENT(IN) :: nma
      LOGICAL,INTENT(IN),OPTIONAL :: l_gammaq
      !
      ! Workspace
      !
      LOGICAL :: l_dohead
      INTEGER :: i1,i2
      COMPLEX(DP) :: f
      COMPLEX(DP) :: ky
#if !defined(__CUDA)
      COMPLEX(DP) :: tempt_c(3,3)
      COMPLEX(DP),ALLOCATABLE :: wh_c(:,:)
      COMPLEX(DP),ALLOCATABLE :: wl_c(:,:)
      COMPLEX(DP),ALLOCATABLE :: x_c(:,:)
      COMPLEX(DP),ALLOCATABLE :: temph_c(:,:)
      COMPLEX(DP),ALLOCATABLE :: templ_c(:,:)
#endif
      !
      COMPLEX(DP),PARAMETER :: one = (1._DP,0._DP)
      COMPLEX(DP),PARAMETER :: zero = (0._DP,0._DP)
      !
      IF(PRESENT(l_gammaq)) THEN
         l_dohead = l_macropol .AND. l_gammaq
      ELSE
         l_dohead = l_macropol
      ENDIF
      !
#if !defined(__CUDA)
      ALLOCATE(x_c(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      IF(l_dohead) THEN
         ALLOCATE(wh_c(n_pdep_eigen_to_use,3))
         ALLOCATE(wl_c(3,n_pdep_eigen_to_use))
         ALLOCATE(temph_c(n_pdep_eigen_to_use,3))
         ALLOCATE(templ_c(3,n_pdep_eigen_to_use))
      ENDIF
#endif
      !
      !$acc update device(matilda)
      !
      IF(l_dohead) THEN
         !
         f = zero
         DO i1 = 1,3
            f = f+matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1)/3._DP
         ENDDO
         !
         !$acc parallel loop collapse(2) present(wh_c,matilda)
         DO i2 = 1,3
            DO i1 = 1,n_pdep_eigen_to_use
               wh_c(i1,i2) = matilda(i1,n_pdep_eigen_to_use+i2)
            ENDDO
         ENDDO
         !$acc end parallel
         !
         !$acc parallel loop collapse(2) present(wl_c,matilda)
         DO i2 = 1,n_pdep_eigen_to_use
            DO i1 = 1,3
               wl_c(i1,i2) = matilda(n_pdep_eigen_to_use+i1,i2)
            ENDDO
         ENDDO
         !$acc end parallel
         !
      ENDIF
      !
      ! X = (1-B)^{-1}
      !
      !$acc parallel loop collapse(2) present(x_c,matilda)
      DO i2 = 1,n_pdep_eigen_to_use
         DO i1 = 1,n_pdep_eigen_to_use
            IF(i1 == i2) THEN
               x_c(i1,i2) = 1._DP-matilda(i1,i2)
            ELSE
               x_c(i1,i2) = -matilda(i1,i2)
            ENDIF
         ENDDO
      ENDDO
      !$acc end parallel
      !
      CALL matinvrs_zge(n_pdep_eigen_to_use,x_c)
      !
      IF(l_dohead) THEN
         !
         ! temph = X * wh
         !
         !$acc host_data use_device(x_c,wh_c,wl_c,temph_c,templ_c,tempt_c)
         CALL ZGEMM('N','N',n_pdep_eigen_to_use,3,n_pdep_eigen_to_use,one,x_c,n_pdep_eigen_to_use,&
         & wh_c,n_pdep_eigen_to_use,zero,temph_c,n_pdep_eigen_to_use)
         !
         ! templ = wl * X
         !
         CALL ZGEMM('N','N',3,n_pdep_eigen_to_use,n_pdep_eigen_to_use,one,wl_c,3,x_c,&
         & n_pdep_eigen_to_use,zero,templ_c,3)
         !
         ! tempt = wl * temph
         !
         CALL ZGEMM('N','N',3,3,n_pdep_eigen_to_use,one,wl_c,3,temph_c,n_pdep_eigen_to_use,zero,&
         & tempt_c,3)
         !$acc end host_data
         !
         !$acc update host(tempt_c)
         !
         ! ky = 1 - f - Tr ( tempt ) / 3
         !
         ky = one-f
         DO i1 = 1,3
            ky = ky-tempt_c(i1,i1)/3._DP
         ENDDO
         !
         head = (one-ky)/ky
         !
      ELSE
         head = zero
      ENDIF
      !
      !$acc host_data use_device(x_c,matilda,temph_c,templ_c,lambda)
      CALL ZGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,n_pdep_eigen_to_use,one,x_c,&
      & n_pdep_eigen_to_use,matilda,nma,zero,lambda,n_pdep_eigen_to_use)
      !
      IF(l_dohead) THEN
         f = one/(3._DP*ky)
         CALL ZGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,3,f,temph_c,n_pdep_eigen_to_use,&
         & templ_c,3,one,lambda,n_pdep_eigen_to_use)
      ENDIF
      !$acc end host_data
      !
      !$acc update host(lambda)
      !
#if !defined(__CUDA)
      DEALLOCATE(x_c)
      IF(l_dohead) THEN
         DEALLOCATE(wh_c)
         DEALLOCATE(wl_c)
         DEALLOCATE(temph_c)
         DEALLOCATE(templ_c)
      ENDIF
#endif
      !
   END SUBROUTINE
   !
END MODULE
