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
MODULE chi_invert
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   PUBLIC :: chi_invert_real, chi_invert_complex
#if defined(__CUDA)
   PUBLIC :: chi_invert_real_gpu, chi_invert_complex_gpu
#endif
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
      REAL(DP) :: tempt(3,3)
      REAL(DP),ALLOCATABLE :: body(:,:)
      REAL(DP),ALLOCATABLE :: wh(:,:)
      REAL(DP),ALLOCATABLE :: wl(:,:)
      REAL(DP),ALLOCATABLE :: x(:,:)
      REAL(DP),ALLOCATABLE :: temph(:,:)
      REAL(DP),ALLOCATABLE :: templ(:,:)
      !
      ALLOCATE(body(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      ALLOCATE(x(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      !
      DO i2 = 1,n_pdep_eigen_to_use
         DO i1 = 1,n_pdep_eigen_to_use
            body(i1,i2) = matilda(i1,i2)
         ENDDO
      ENDDO
      !
      IF(l_macropol) THEN
         !
         f = 0._DP
         DO i1 = 1,3
            f = f+matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1)/3._DP
         ENDDO
         !
         ALLOCATE(wh(n_pdep_eigen_to_use,3))
         DO i2 = 1,3
            DO i1 = 1,n_pdep_eigen_to_use
               wh(i1,i2) = matilda(i1,n_pdep_eigen_to_use+i2)
            ENDDO
         ENDDO
         !
         ALLOCATE(wl(3,n_pdep_eigen_to_use))
         DO i2 = 1,n_pdep_eigen_to_use
            DO i1 = 1,3
               wl(i1,i2) = matilda(n_pdep_eigen_to_use+i1,i2)
            ENDDO
         ENDDO
         !
      ENDIF
      !
      ! X = (1-B)^{-1}
      x(:,:) = -body
      DO i1 = 1,n_pdep_eigen_to_use
         x(i1,i1) = x(i1,i1)+1._DP
      ENDDO
      !
      CALL matinvrs_dge(n_pdep_eigen_to_use,x)
      !
      IF(l_macropol) THEN
         !
         ! temph = X * wh
         ALLOCATE(temph(n_pdep_eigen_to_use,3))
         CALL DGEMM('N','N',n_pdep_eigen_to_use,3,n_pdep_eigen_to_use,1._DP,x,n_pdep_eigen_to_use,&
         & wh,n_pdep_eigen_to_use,0._DP,temph,n_pdep_eigen_to_use)
         DEALLOCATE(wh)
         !
         ! templ = wl * X
         ALLOCATE(templ(3,n_pdep_eigen_to_use))
         CALL DGEMM('N','N',3,n_pdep_eigen_to_use,n_pdep_eigen_to_use,1._DP,wl,3,x,&
         & n_pdep_eigen_to_use,0._DP,templ,3)
         !
         ! tempt = wl * temph
         CALL DGEMM('N','N',3,3,n_pdep_eigen_to_use,1._DP,wl,3,temph,n_pdep_eigen_to_use,0._DP,tempt,3)
         DEALLOCATE(wl)
         !
         ! ky = 1 - f - Tr ( tempt ) / 3
         ky = 1._DP-f
         DO i1 = 1,3
            ky = ky-tempt(i1,i1)/3._DP
         ENDDO
         !
         head = (1._DP-ky)/ky
         !
      ELSE
         head = 0._DP
      ENDIF
      !
      CALL DGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,n_pdep_eigen_to_use,1._DP,x,n_pdep_eigen_to_use,&
      & body,n_pdep_eigen_to_use,0._DP,lambda,n_pdep_eigen_to_use)
      !
      IF(l_macropol) THEN
         CALL DGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,3,1._DP/(3._DP*ky),temph,&
         & n_pdep_eigen_to_use,templ,3,1._DP,lambda,n_pdep_eigen_to_use)
         DEALLOCATE(temph)
         DEALLOCATE(templ)
      ENDIF
      !
      DEALLOCATE(body,x)
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
      COMPLEX(DP) :: tempt(3,3)
      COMPLEX(DP),ALLOCATABLE :: body(:,:)
      COMPLEX(DP),ALLOCATABLE :: wh(:,:)
      COMPLEX(DP),ALLOCATABLE :: wl(:,:)
      COMPLEX(DP),ALLOCATABLE :: x(:,:)
      COMPLEX(DP),ALLOCATABLE :: temph(:,:)
      COMPLEX(DP),ALLOCATABLE :: templ(:,:)
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
      ALLOCATE(body(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      ALLOCATE(x(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      !
      DO i2 = 1,n_pdep_eigen_to_use
         DO i1 = 1,n_pdep_eigen_to_use
            body(i1,i2) = matilda(i1,i2)
         ENDDO
      ENDDO
      !
      IF(l_dohead) THEN
         !
         f = zero
         DO i1 = 1,3
            f = f+matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1)/3._DP
         ENDDO
         !
         ALLOCATE(wh(n_pdep_eigen_to_use,3))
         DO i2 = 1,3
            DO i1 = 1,n_pdep_eigen_to_use
               wh(i1,i2) = matilda(i1,n_pdep_eigen_to_use+i2)
            ENDDO
         ENDDO
         !
         ALLOCATE(wl(3,n_pdep_eigen_to_use))
         DO i2 = 1,n_pdep_eigen_to_use
            DO i1 = 1,3
               wl(i1,i2) = matilda(n_pdep_eigen_to_use+i1,i2)
            ENDDO
         ENDDO
         !
      ENDIF
      !
      ! X = (1-B)^{-1}
      x(:,:) = -body
      DO i1 = 1,n_pdep_eigen_to_use
         x(i1,i1) = x(i1,i1)+one
      ENDDO
      !
      CALL matinvrs_zge(n_pdep_eigen_to_use,x)
      !
      IF(l_dohead) THEN
         !
         ! temph = X * wh
         ALLOCATE(temph(n_pdep_eigen_to_use,3))
         CALL ZGEMM('N','N',n_pdep_eigen_to_use,3,n_pdep_eigen_to_use,one,x,n_pdep_eigen_to_use,&
         & wh,n_pdep_eigen_to_use,zero,temph,n_pdep_eigen_to_use)
         DEALLOCATE(wh)
         !
         ! templ = wl * X
         ALLOCATE(templ(3,n_pdep_eigen_to_use))
         CALL ZGEMM('N','N',3,n_pdep_eigen_to_use,n_pdep_eigen_to_use,one,wl,3,x,&
         & n_pdep_eigen_to_use,zero,templ,3)
         !
         ! tempt = wl * temph
         CALL ZGEMM('N','N',3,3,n_pdep_eigen_to_use,one,wl,3,temph,n_pdep_eigen_to_use,zero,tempt,3)
         DEALLOCATE(wl)
         !
         ! ky = 1 - f - Tr ( tempt ) / 3
         ky = one-f
         DO i1 = 1,3
            ky = ky-tempt(i1,i1)/3._DP
         ENDDO
         !
         head = (one-ky)/ky
         !
      ELSE
         head = zero
      ENDIF
      !
      CALL ZGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,n_pdep_eigen_to_use,one,x,n_pdep_eigen_to_use,&
      & body,n_pdep_eigen_to_use,zero,lambda,n_pdep_eigen_to_use)
      !
      IF(l_dohead) THEN
         CALL ZGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,3,one/(3._DP*ky),temph,&
         & n_pdep_eigen_to_use,templ,3,one,lambda,n_pdep_eigen_to_use)
         DEALLOCATE(temph)
         DEALLOCATE(templ)
      ENDIF
      !
      DEALLOCATE(body,x)
      !
   END SUBROUTINE
   !
#if defined(__CUDA)
   !-----------------------------------------------------------------------
   SUBROUTINE chi_invert_real_gpu(matilda,head,lambda,nma)
      !-----------------------------------------------------------------------
      !
      ! For each frequency I calculate X, ky, head and lambda
      !
      ! X = (1-B)^{-1}
      ! tmph = X * wh
      ! tmpl = wl * X
      ! tmpt = wl * tmph
      !
      ! ky = 1 - f - Tr ( tmpt ) / 3
      !
      ! Head = (1-ky) / ky
      !
      ! Lamba = X * B + tmph * tmpl / (3*ky)
      !
      USE kinds,                 ONLY : DP
      USE linear_algebra_kernel, ONLY : matinvrs_dge_gpu
      USE westcom,               ONLY : n_pdep_eigen_to_use,l_macropol
      USE west_gpu,              ONLY : body_r_d,x_r_d,wh_r_d,wl_r_d,tmph_r_d,tmpl_r_d,tmpt_r_d,lambda_r_d
      USE cublas
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
      REAL(DP) :: tmpt(3,3)
      !
      body_r_d = matilda(1:n_pdep_eigen_to_use,1:n_pdep_eigen_to_use)
      !
      IF(l_macropol) THEN
         !
         f = 0._DP
         DO i1 = 1,3
            f = f+matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1)/3._DP
         ENDDO
         !
         wh_r_d = matilda(1:n_pdep_eigen_to_use,n_pdep_eigen_to_use+1:n_pdep_eigen_to_use+3)
         wl_r_d = matilda(n_pdep_eigen_to_use+1:n_pdep_eigen_to_use+3,1:n_pdep_eigen_to_use)
         !
      ENDIF
      !
      ! X = (1-B)^{-1}
      !
      !$acc parallel loop collapse(2)
      DO i2 = 1,n_pdep_eigen_to_use
         DO i1 = 1,n_pdep_eigen_to_use
            IF(i1 == i2) THEN
               x_r_d(i1,i2) = 1._DP-body_r_d(i1,i2)
            ELSE
               x_r_d(i1,i2) = -body_r_d(i1,i2)
            ENDIF
         ENDDO
      ENDDO
      !$acc end parallel
      !
      CALL matinvrs_dge_gpu(n_pdep_eigen_to_use,x_r_d)
      !
      IF(l_macropol) THEN
         !
         ! tmph = X * wh
         CALL DGEMM('N','N',n_pdep_eigen_to_use,3,n_pdep_eigen_to_use,1._DP,x_r_d,&
         & n_pdep_eigen_to_use,wh_r_d,n_pdep_eigen_to_use,0._DP,tmph_r_d,n_pdep_eigen_to_use)
         !
         ! tmpl = wl * X
         CALL DGEMM('N','N',3,n_pdep_eigen_to_use,n_pdep_eigen_to_use,1._DP,wl_r_d,3,x_r_d,&
         & n_pdep_eigen_to_use,0._DP,tmpl_r_d,3)
         !
         ! tmpt = wl * tmph
         CALL DGEMM('N','N',3,3,n_pdep_eigen_to_use,1._DP,wl_r_d,3,tmph_r_d,n_pdep_eigen_to_use,&
         & 0._DP,tmpt_r_d,3)
         !
         tmpt = tmpt_r_d
         !
         ! ky = 1 - f - Tr ( tmpt ) / 3
         ky = 1._DP-f
         DO i1 = 1,3
            ky = ky-tmpt(i1,i1)/3._DP
         ENDDO
         !
         head = (1._DP-ky)/ky
         !
      ELSE
         head = 0._DP
      ENDIF
      !
      CALL DGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,n_pdep_eigen_to_use,1._DP,x_r_d,&
      & n_pdep_eigen_to_use,body_r_d,n_pdep_eigen_to_use,0._DP,lambda_r_d,n_pdep_eigen_to_use)
      !
      IF(l_macropol) THEN
         f = 1._DP/(3._DP*ky)
         CALL DGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,3,f,tmph_r_d,&
         & n_pdep_eigen_to_use,tmpl_r_d,3,1._DP,lambda_r_d,n_pdep_eigen_to_use)
      ENDIF
      !
      lambda = lambda_r_d
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE chi_invert_complex_gpu(matilda,head,lambda,nma,l_gammaq)
      !-----------------------------------------------------------------------
      !
      ! For each frequency and q-point I calculate X, ky, head and lambda
      !
      ! X = (1-B)^{-1}
      ! tmph = X * wh
      ! tmpl = wl * X
      ! tmpt = wl * tmph
      !
      ! ky = 1 - f - Tr ( tmpt ) / 3
      !
      ! Head = (1-ky) / ky
      !
      ! Lamba = X * B + tmph * tmpl / (3*ky)
      !
      USE kinds,                 ONLY : DP
      USE linear_algebra_kernel, ONLY : matinvrs_zge_gpu
      USE westcom,               ONLY : n_pdep_eigen_to_use,l_macropol
      USE west_gpu,              ONLY : body_c_d,x_c_d,wh_c_d,wl_c_d,tmph_c_d,tmpl_c_d,tmpt_c_d,lambda_c_d
      USE cublas
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
      COMPLEX(DP) :: tmpt(3,3)
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
      body_c_d = matilda(1:n_pdep_eigen_to_use,1:n_pdep_eigen_to_use)
      !
      IF(l_dohead) THEN
         !
         f = zero
         DO i1 = 1,3
            f = f+matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1)/3._DP
         ENDDO
         !
         wh_c_d = matilda(1:n_pdep_eigen_to_use,n_pdep_eigen_to_use+1:n_pdep_eigen_to_use+3)
         wl_c_d = matilda(n_pdep_eigen_to_use+1:n_pdep_eigen_to_use+3,1:n_pdep_eigen_to_use)
         !
      ENDIF
      !
      ! X = (1-B)^{-1}
      !
      !$acc parallel loop collapse(2)
      DO i2 = 1,n_pdep_eigen_to_use
         DO i1 = 1,n_pdep_eigen_to_use
            IF(i1 == i2) THEN
               x_c_d(i1,i2) = one-body_c_d(i1,i2)
            ELSE
               x_c_d(i1,i2) = -body_c_d(i1,i2)
            ENDIF
         ENDDO
      ENDDO
      !$acc end parallel
      !
      CALL matinvrs_zge_gpu(n_pdep_eigen_to_use,x_c_d)
      !
      IF(l_dohead) THEN
         !
         ! tmph = X * wh
         CALL ZGEMM('N','N',n_pdep_eigen_to_use,3,n_pdep_eigen_to_use,one,x_c_d,&
         & n_pdep_eigen_to_use,wh_c_d,n_pdep_eigen_to_use,zero,tmph_c_d,n_pdep_eigen_to_use)
         !
         ! tmpl = wl * X
         CALL ZGEMM('N','N',3,n_pdep_eigen_to_use,n_pdep_eigen_to_use,one,wl_c_d,3,x_c_d,&
         & n_pdep_eigen_to_use,zero,tmpl_c_d,3)
         !
         ! tmpt = wl * tmph
         CALL ZGEMM('N','N',3,3,n_pdep_eigen_to_use,one,wl_c_d,3,tmph_c_d,n_pdep_eigen_to_use,zero,&
         & tmpt_c_d,3)
         !
         tmpt = tmpt_c_d
         !
         ! ky = 1 - f - Tr ( tmpt ) / 3
         ky = one-f
         DO i1 = 1,3
            ky = ky-tmpt(i1,i1)/3._DP
         ENDDO
         !
         head = (one-ky)/ky
         !
      ELSE
         head = zero
      ENDIF
      !
      CALL ZGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,n_pdep_eigen_to_use,one,x_c_d,&
      & n_pdep_eigen_to_use,body_c_d,n_pdep_eigen_to_use,zero,lambda_c_d,n_pdep_eigen_to_use)
      !
      IF(l_dohead) THEN
         f = one/(3._DP*ky)
         CALL ZGEMM('N','N',n_pdep_eigen_to_use,n_pdep_eigen_to_use,3,f,tmph_c_d,&
         & n_pdep_eigen_to_use,tmpl_c_d,3,one,lambda_c_d,n_pdep_eigen_to_use)
      ENDIF
      !
      lambda = lambda_c_d
      !
   END SUBROUTINE
#endif
   !
END MODULE
