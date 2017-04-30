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
  USE westcom,               ONLY : west_prefix,n_pdep_eigen_to_use,l_macropol
  USE io_files,              ONLY : tmp_dir
  !
  ! I/O
  !
  REAL(DP),INTENT(IN) :: matilda(nma,nma)
  REAL(DP),INTENT(OUT) :: head,lambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use) 
  INTEGER,INTENT(IN) :: nma
  !
  ! Workspace
  ! 
  REAL(DP) :: f
  REAL(DP), ALLOCATABLE :: body(:,:)
  REAL(DP),ALLOCATABLE :: wh(:,:)
  REAL(DP),ALLOCATABLE :: wl(:,:)
  INTEGER :: i1,i2
  REAL(DP),ALLOCATABLE :: x(:,:)
  REAL(DP),ALLOCATABLE :: temph(:,:)
  REAL(DP),ALLOCATABLE :: templ(:,:)
  REAL(DP) :: tempt(3,3)
  REAL(DP) :: ky
  !
  ALLOCATE( body(n_pdep_eigen_to_use,n_pdep_eigen_to_use) )
  ALLOCATE( x(n_pdep_eigen_to_use,n_pdep_eigen_to_use) )
  !
  DO i2 = 1, n_pdep_eigen_to_use
     DO i1 = 1, n_pdep_eigen_to_use
        body(i1,i2) = matilda(i1,i2) 
     ENDDO
  ENDDO
  !
  IF(l_macropol) THEN
     !
     f = 0._DP
     DO i1 = 1, 3
        f = f + matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1) / 3._DP
     ENDDO
     !
     ALLOCATE( wh( n_pdep_eigen_to_use, 3)) 
     DO i2 = 1, 3
        DO i1 = 1, n_pdep_eigen_to_use
           wh(i1,i2) = matilda(i1,n_pdep_eigen_to_use+i2)
        ENDDO
     ENDDO
     !
     ALLOCATE( wl( 3, n_pdep_eigen_to_use)) 
     DO i2 = 1, n_pdep_eigen_to_use
        DO i1 = 1, 3
           wl(i1,i2) = matilda(n_pdep_eigen_to_use+i1,i2)
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! X = (1-B)^{-1}
  x(:,:) = -body(:,:)
  DO i1=1,n_pdep_eigen_to_use
     x(i1,i1) = x(i1,i1) + 1._DP
  ENDDO
  !
  CALL matinvrs_dge(n_pdep_eigen_to_use,x)
  !
  IF( l_macropol) THEN
     !
     ! temph = X * wh
     ALLOCATE( temph(n_pdep_eigen_to_use,3) )
     CALL DGEMM( 'N', 'N', n_pdep_eigen_to_use, 3, n_pdep_eigen_to_use, 1._DP, x, n_pdep_eigen_to_use, &
     & wh, n_pdep_eigen_to_use, 0._DP, temph, n_pdep_eigen_to_use )
     DEALLOCATE(wh)
     !
     ! templ = wl * X
     ALLOCATE( templ(3,n_pdep_eigen_to_use) )
     CALL DGEMM( 'N', 'N', 3, n_pdep_eigen_to_use, n_pdep_eigen_to_use, 1._DP, wl, 3, x, &
     & n_pdep_eigen_to_use, 0._DP, templ, 3 )
     !
     ! tempt = wl * temph
     CALL DGEMM( 'N', 'N', 3, 3, n_pdep_eigen_to_use, 1._DP, wl, 3, temph, n_pdep_eigen_to_use, 0._DP, tempt, 3 )
     DEALLOCATE(wl)
     !
     ! ky = 1 - f - Tr ( tempt ) / 3
     ky = 1._DP - f
     DO i1=1,3
        ky = ky - tempt(i1,i1) / 3._DP
     ENDDO
     !
     head = (1._DP - ky) / ky
     !
  ELSE
     head = 0._DP
  ENDIF
  !
  CALL DGEMM( 'N', 'N', n_pdep_eigen_to_use, n_pdep_eigen_to_use, n_pdep_eigen_to_use, 1._DP, x, n_pdep_eigen_to_use, &
  & body, n_pdep_eigen_to_use, 0._DP, lambda, n_pdep_eigen_to_use )
  !
  IF( l_macropol ) THEN
     CALL DGEMM( 'N', 'N', n_pdep_eigen_to_use, n_pdep_eigen_to_use, 3, 1._DP/(3._DP*ky), temph, &
     & n_pdep_eigen_to_use, templ, 3, 1._DP, lambda, n_pdep_eigen_to_use )
     DEALLOCATE( temph )
     DEALLOCATE( templ )
  ENDIF
  !
  DEALLOCATE( body, x) 
  !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE chi_invert_complex(matilda,head,lambda,nma)
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
  USE linear_algebra_kernel, ONLY : matinvrs_zge 
  USE westcom,               ONLY : west_prefix,n_pdep_eigen_to_use,l_macropol
  USE io_files,              ONLY : tmp_dir
  !
  ! I/O
  !
  COMPLEX(DP),INTENT(IN) :: matilda(nma,nma)
  COMPLEX(DP),INTENT(OUT) :: head,lambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use) 
  INTEGER,INTENT(IN) :: nma
  !
  ! Workspace
  ! 
  COMPLEX(DP) :: f
  COMPLEX(DP),ALLOCATABLE :: body(:,:)
  COMPLEX(DP),ALLOCATABLE :: wh(:,:)
  COMPLEX(DP),ALLOCATABLE :: wl(:,:)
  INTEGER :: i1,i2
  COMPLEX(DP),ALLOCATABLE :: x(:,:)
  COMPLEX(DP),ALLOCATABLE :: temph(:,:)
  COMPLEX(DP),ALLOCATABLE :: templ(:,:)
  COMPLEX(DP) :: tempt(3,3)
  COMPLEX(DP) :: ky,Zone,Zzero
  !
  ALLOCATE( body(n_pdep_eigen_to_use,n_pdep_eigen_to_use) )
  ALLOCATE( x(n_pdep_eigen_to_use,n_pdep_eigen_to_use) )
  !
  Zone = CMPLX(1._DP,0._DP,KIND=DP)
  Zzero = CMPLX(0._DP,0._DP,KIND=DP)
  !
  DO i2 = 1, n_pdep_eigen_to_use
     DO i1 = 1, n_pdep_eigen_to_use
        body(i1,i2) = matilda(i1,i2) 
     ENDDO
  ENDDO
  !
  IF(l_macropol) THEN
     !
     f = Zzero
     DO i1 = 1, 3
        f = f + matilda(n_pdep_eigen_to_use+i1,n_pdep_eigen_to_use+i1) / 3._DP
     ENDDO
     !
     ALLOCATE( wh( n_pdep_eigen_to_use, 3)) 
     DO i2 = 1, 3
        DO i1 = 1, n_pdep_eigen_to_use
           wh(i1,i2) = matilda(i1,n_pdep_eigen_to_use+i2)
        ENDDO
     ENDDO
     !
     ALLOCATE( wl( 3, n_pdep_eigen_to_use)) 
     DO i2 = 1, n_pdep_eigen_to_use
        DO i1 = 1, 3
           wl(i1,i2) = matilda(n_pdep_eigen_to_use+i1,i2)
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! X = (1-B)^{-1}
  x(:,:) = -body(:,:)
  DO i1=1,n_pdep_eigen_to_use
     x(i1,i1) = x(i1,i1) + Zone
  ENDDO
  !
  CALL matinvrs_zge(n_pdep_eigen_to_use,x)
  !
  IF( l_macropol) THEN
     !
     ! temph = X * wh
     ALLOCATE( temph(n_pdep_eigen_to_use,3) )
     CALL ZGEMM( 'N', 'N', n_pdep_eigen_to_use, 3, n_pdep_eigen_to_use, Zone, x, n_pdep_eigen_to_use, &
     & wh, n_pdep_eigen_to_use, Zzero, temph, n_pdep_eigen_to_use )
     DEALLOCATE(wh)
     !
     ! templ = wl * X
     ALLOCATE( templ(3,n_pdep_eigen_to_use) )
     CALL ZGEMM( 'N', 'N', 3, n_pdep_eigen_to_use, n_pdep_eigen_to_use, Zone, wl, 3, x, &
     & n_pdep_eigen_to_use, Zzero, templ, 3 )
     !
     ! tempt = wl * temph
     CALL ZGEMM( 'N', 'N', 3, 3, n_pdep_eigen_to_use, Zone, wl, 3, temph, n_pdep_eigen_to_use, Zzero, tempt, 3 )
     DEALLOCATE(wl)
     !
     ! ky = 1 - f - Tr ( tempt ) / 3
     ky = Zone - f
     DO i1=1,3
        ky = ky - tempt(i1,i1) / 3._DP
     ENDDO
     !
     head = (Zone - ky) / ky
     !
  ELSE
     head = Zzero
  ENDIF
  !
  CALL ZGEMM( 'N', 'N', n_pdep_eigen_to_use, n_pdep_eigen_to_use, n_pdep_eigen_to_use, Zone, x, n_pdep_eigen_to_use, &
  & body, n_pdep_eigen_to_use, Zzero, lambda, n_pdep_eigen_to_use )
  !
  IF( l_macropol ) THEN
     CALL ZGEMM( 'N', 'N', n_pdep_eigen_to_use, n_pdep_eigen_to_use, 3, Zone/(3._DP*ky), temph, &
     & n_pdep_eigen_to_use, templ, 3, Zone, lambda, n_pdep_eigen_to_use )
     DEALLOCATE( temph )
     DEALLOCATE( templ )
  ENDIF
  !
  DEALLOCATE( body, x )
  !
END SUBROUTINE
