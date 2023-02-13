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
MODULE wann_loc_wfc
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE wann_calc_proj(proj)
      !------------------------------------------------------------------------
      !
      USE kinds,                 ONLY : DP
      USE constants,             ONLY : tpi
      USE fft_base,              ONLY : dffts
      USE scatter_mod,           ONLY : scatter_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      REAL(DP),INTENT(OUT) :: proj(dffts%nnr,6)
      !
      ! Workspace
      !
      INTEGER :: il,ir,ix,iy,iz
      REAL(DP) :: nx,ny,nz,wcx,wcy,wcz,wsx,wsy,wsz
      !
      REAL(DP),ALLOCATABLE :: prod_gat(:)
      REAL(DP),ALLOCATABLE :: prod_distr(:)
      !
      proj = 0._DP
      !
      ALLOCATE(prod_gat(dffts%nr1x*dffts%nr2x*dffts%nr3x))
      ALLOCATE(prod_distr(dffts%nnr))
      !
      nx = REAL(dffts%nr1,KIND=DP)
      ny = REAL(dffts%nr2,KIND=DP)
      nz = REAL(dffts%nr3,KIND=DP)
      !
      DO il = 1,6
         !
         prod_gat = 0._DP
         prod_distr = 0._DP
         !
         ir = 0
         DO ix = 1,dffts%nr1
            !
            wcx = COS(tpi*REAL(ix-1,KIND=DP)/nx)
            wsx = SIN(tpi*REAL(ix-1,KIND=DP)/nx)
            !
            DO iy = 1,dffts%nr2
               !
               wcy = COS(tpi*REAL(iy-1,KIND=DP)/ny)
               wsy = SIN(tpi*REAL(iy-1,KIND=DP)/ny)
               !
               DO iz = 1,dffts%nr3
                  !
                  wcz = COS(tpi*REAL(iz-1,KIND=DP)/nz)
                  wsz = SIN(tpi*REAL(iz-1,KIND=DP)/nz)
                  !
                  ir = (iz-1)*(dffts%nr1x*dffts%nr2x) + (iy-1)*dffts%nr1x + ix
                  IF(il == 1) prod_gat(ir) = wcx
                  IF(il == 2) prod_gat(ir) = wsx
                  IF(il == 3) prod_gat(ir) = wcy
                  IF(il == 4) prod_gat(ir) = wsy
                  IF(il == 5) prod_gat(ir) = wcz
                  IF(il == 6) prod_gat(ir) = wsz
               ENDDO
               !
            ENDDO
            !
         ENDDO
         !
         CALL scatter_grid(dffts,prod_gat,prod_distr)
         !
         DO ir = 1,dffts%nnr
            proj(ir,il) = prod_distr(ir) / (nx*ny*nz)
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE(prod_gat)
      DEALLOCATE(prod_distr)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wann_joint_d(m,a,na,u)
      !------------------------------------------------------------------------
      !
      USE kinds,                 ONLY : DP
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: m
      INTEGER,INTENT(IN) :: na
      REAL(DP),INTENT(INOUT) :: a(m,m,na)
      COMPLEX(DP),INTENT(OUT) :: u(m,m)
      !
      ! Workspace
      !
      LOGICAL :: conv
      INTEGER :: ia,i,j,k,iter
      REAL(DP) :: sigma,sigma_old,tmpi,tmpj
      REAL(DP) :: c,s,h1,h2,e1,e2,g11,g22,g12,x,y,t,tau
      !
      INTEGER,PARAMETER :: itermax = 5000
      REAL(DP),PARAMETER :: spread_thr = 1.E-9_DP
      !
      u = 0._DP
      DO i = 1,m
         u(i,i) = 1._DP
      ENDDO
      !
      sigma_old = 0._DP
      DO ia = 1,na
         DO j = 1,m
            sigma_old = sigma_old+a(j,j,ia)*a(j,j,ia)
         ENDDO
      ENDDO
      !
      conv = .FALSE.
      !
      DO iter = 1,itermax
         !
         DO i = 1,m
            DO j = i+1,m
               !
               g11 = 0._DP
               g12 = 0._DP
               g22 = 0._DP
               !
               DO ia = 1,na
                  h1 = a(i,i,ia)-a(j,j,ia)
                  h2 = 2._DP*a(i,j,ia)
                  g11 = g11+h1*h1
                  g12 = g12+h1*h2
                  g22 = g22+h2*h2
               ENDDO
               !
               c = 1._DP
               s = 0._DP
               e1 = g11
               e2 = g22
               !
               IF(g12*g12 > 1.E-16_DP*ABS(g11*g22)) THEN
                  tau = 0.5_DP * (g22-g11) / g12
                  t = 1.0_DP / (ABS(tau) + SQRT(1._DP+tau**2))
                  IF(tau < 0._DP) t = -t
                  c = 1._DP / SQRT(1._DP+t**2)
                  s = t * c
                  e1 = e1 - t*g12
                  e2 = e2 + t*g12
               ENDIF
               !
               IF(e1 > e2) THEN
                  x = c
                  y = -s
               ELSE
                  x = s
                  y = c
               ENDIF
               !
               IF(x < 0._DP) THEN
                  x = -x
                  y = -y
               ENDIF
               !
               c = SQRT(0.5_DP*(x+1._DP))
               s = y / SQRT(2._DP*(x+1._DP))
               !
               DO ia = 1,na
                  DO k = 1,m
                     tmpi =  a(k,i,ia)*c + a(k,j,ia)*s
                     tmpj = -a(k,i,ia)*s + a(k,j,ia)*c
                     a(k,i,ia) = tmpi
                     a(k,j,ia) = tmpj
                  ENDDO
               ENDDO
               !
               DO ia = 1,na
                  DO k = 1,m
                     tmpi =  c*a(i,k,ia) + s*a(j,k,ia)
                     tmpj = -s*a(i,k,ia) + c*a(j,k,ia)
                     a(i,k,ia) = tmpi
                     a(j,k,ia) = tmpj
                  ENDDO
               ENDDO
               !
               DO k = 1,m
                  tmpi =  u(k,i)*c + u(k,j)*s
                  tmpj = -u(k,i)*s + u(k,j)*c
                  u(k,i) = CMPLX(tmpi,KIND=DP)
                  u(k,j) = CMPLX(tmpj,KIND=DP)
               ENDDO
               !
            ENDDO
         ENDDO
         !
         sigma = 0._DP
         DO ia = 1,na
            DO j = 1,m
               sigma = sigma + a(j,j,ia)*a(j,j,ia)
            ENDDO
         ENDDO
         !
         IF(ABS(sigma-sigma_old) < spread_thr) THEN
            conv = .TRUE.
            EXIT
         ELSE
            sigma_old = sigma
         ENDIF
         !
      ENDDO
      !
      IF(.NOT. conv) CALL errore('wann','convergence not achieved',itermax)
      !
    END SUBROUTINE
    !
END MODULE
