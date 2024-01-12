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
    SUBROUTINE wann_jade(m,a,na,u)
      !------------------------------------------------------------------------
      !
      ! Joint approximate diagonalization of eigen-matrices
      !
      ! Gygi et al., Computer Physics Communications 155, 1-6 (2003)
      !
      USE kinds,                 ONLY : DP
      USE io_global,             ONLY : stdout
      USE linear_algebra_kernel, ONLY : matdiago_dsy
      USE io_push,               ONLY : io_push_title
      USE westcom,               ONLY : wannier_tr_rel
#if defined(__CUDA)
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: m
      INTEGER,INTENT(IN) :: na
      REAL(DP),INTENT(INOUT) :: a(m,m,na)
      REAL(DP),INTENT(OUT) :: u(m,m)
      !
      ! Workspace
      !
      LOGICAL :: conv
      INTEGER :: ia,i,k,mwork,iter,sweep,p,q
      REAL(DP) :: sigma,sigma_old
      REAL(DP) :: c,s,h1,h2,e1,e2,g11,g22,g12,x,y,t,tau
      !
      INTEGER,ALLOCATABLE :: top(:),bot(:)
      REAL(DP),ALLOCATABLE :: ev(:)
      REAL(DP),ALLOCATABLE :: rot(:,:),aux(:,:)
      !$acc declare device_resident(rot,aux)
      !
      INTEGER,PARAMETER :: itermax = 100
      !
      REAL(DP) :: time_spent(2)
      REAL(DP), EXTERNAL :: get_clock
      CHARACTER(20), EXTERNAL :: human_readable_time
      !
#if defined(__CUDA)
      CALL start_clock_gpu('jade')
#else
      CALL start_clock('jade')
#endif
      !
      CALL io_push_title('Wannier (JADE)')
      !
      ALLOCATE(rot(m,m))
      ALLOCATE(aux(m,m))
      !
      ! Handle odd m
      !
      IF(MOD(m,2) == 0) THEN
         mwork = m
      ELSE
         mwork = m+1
      ENDIF
      !
      ALLOCATE(top(mwork/2))
      ALLOCATE(bot(mwork/2))
      !
      DO k = 1,mwork/2
         top(k) = k*2 - 1
         bot(k) = k*2
      ENDDO
      !
      u(:,:) = a(:,:,1)
      !
      ALLOCATE(ev(m))
      !
      CALL matdiago_dsy(m,u,ev,.FALSE.)
      !
      DEALLOCATE(ev)
      !
      !$acc enter data copyin(a,u,top,bot)
      !
      !$acc host_data use_device(u,a,aux)
      DO ia = 1,na
         CALL DGEMM('T','N',m,m,m,1._DP,u,m,a(1,1,ia),m,0._DP,aux,m)
         CALL DGEMM('N','N',m,m,m,1._DP,aux,m,u,m,0._DP,a(1,1,ia),m)
      ENDDO
      !$acc end host_data
      !
      ! Compute initial spread
      !
      sigma_old = 0._DP
      !$acc parallel loop collapse(2) reduction(+:sigma_old) present(a) copy(sigma_old)
      DO ia = 1,na
         DO i = 1,m
            sigma_old = sigma_old + a(i,i,ia)**2
         ENDDO
      ENDDO
      !$acc end parallel
      !
      conv = .FALSE.
      !
      DO iter = 1,itermax
         !
         time_spent(1) = get_clock('jade')
         !
         DO sweep = 1,m-1
            !
            !$acc kernels present(rot)
            rot(:,:) = 0._DP
            !$acc end kernels
            !
            !$acc parallel loop present(top,bot,rot,a)
            DO k = 1,mwork/2
               !
               p = MIN(top(k),bot(k))
               q = MAX(top(k),bot(k))
               !
               ! Handle odd m
               !
               IF(q <= m) THEN
                  !
                  ! Compute 2x2 matrix G
                  !
                  g11 = 0._DP
                  g12 = 0._DP
                  g22 = 0._DP
                  !
                  DO ia = 1,na
                     h1 = a(p,p,ia)-a(q,q,ia)
                     h2 = 2._DP*a(p,q,ia)
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
                  ! Compute eigenvalues and eigenvectors of G
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
                  ! Use the eigenvector with the largest eigenvalue
                  !
                  IF(e1 > e2) THEN
                     x = c
                     y = -s
                  ELSE
                     x = s
                     y = c
                  ENDIF
                  !
                  ! Choose x >= 0 to ensure small rotation angle
                  !
                  IF(x < 0._DP) THEN
                     x = -x
                     y = -y
                  ENDIF
                  !
                  ! Compute 2x2 rotation matrix R
                  !
                  c = SQRT(0.5_DP*(x+1._DP))
                  s = y / SQRT(2._DP*(x+1._DP))
                  !
                  rot(p,p) = c
                  rot(q,p) = s
                  rot(p,q) = -s
                  rot(q,q) = c
                  !
               ELSE
                  !
                  rot(p,p) = 1._DP
                  !
               ENDIF
               !
            ENDDO
            !$acc end parallel
            !
            ! Apply rotation R to rows and columns of A
            !
            !$acc host_data use_device(rot,a,aux)
            DO ia = 1,na
               CALL DGEMM('T','N',m,m,m,1._DP,rot,m,a(1,1,ia),m,0._DP,aux,m)
               CALL DGEMM('N','N',m,m,m,1._DP,aux,m,rot,m,0._DP,a(1,1,ia),m)
            ENDDO
            !$acc end host_data
            !
            ! Accumulate unitary transformation matrix U
            !
            !$acc host_data use_device(u,rot,aux)
            CALL DGEMM('N','N',m,m,m,1._DP,u,m,rot,m,0._DP,aux,m)
            !$acc end host_data
            !
            !$acc kernels present(u,aux)
            u(:,:) = aux
            !$acc end kernels
            !
            ! Go to next round of tournament
            !
            CALL wann_tournament(top,bot,mwork)
            !
            !$acc update device(top,bot)
            !
         ENDDO
         !
         ! Compute new spread
         !
         sigma = 0._DP
         !$acc parallel loop collapse(2) reduction(+:sigma) present(a) copy(sigma)
         DO ia = 1,na
            DO i = 1,m
               sigma = sigma + a(i,i,ia)**2
            ENDDO
         ENDDO
         !$acc end parallel
         !
         WRITE(stdout,"(/,5X,'                  *----------*            *-----------------*')")
         WRITE(stdout,"(  5X,'#     Iteration = | ', I8,' |','   ','Spread = | ', ES15.8,' |')") &
         & iter, sigma
         WRITE(stdout,"(  5X,'                  *----------*            *-----------------*')")
         !
         time_spent(2) = get_clock('jade')
         !
         WRITE(stdout,"(5X,'Time spent in last iteration ',A)") &
         & TRIM(human_readable_time(time_spent(2)-time_spent(1)))
         !
         ! Check convergence
         !
         IF(ABS((sigma-sigma_old)/sigma_old) < wannier_tr_rel) THEN
            conv = .TRUE.
            EXIT
         ELSE
            sigma_old = sigma
         ENDIF
         !
      ENDDO
      !
      !$acc exit data delete(top,bot) copyout(a,u)
      DEALLOCATE(rot)
      DEALLOCATE(aux)
      DEALLOCATE(top)
      DEALLOCATE(bot)
      !
      IF(.NOT. conv) WRITE(stdout,'(7X,"** WARNING : JADE not converged in ",I5," steps")') itermax
      !
#if defined(__CUDA)
      CALL stop_clock_gpu('jade')
#else
      CALL stop_clock('jade')
#endif
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wann_tournament(top,bot,m)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: m
      INTEGER,INTENT(INOUT) :: top(m/2)
      INTEGER,INTENT(INOUT) :: bot(m/2)
      !
      INTEGER,ALLOCATABLE :: new_top(:)
      INTEGER,ALLOCATABLE :: new_bot(:)
      !
      ALLOCATE(new_top(m/2))
      ALLOCATE(new_bot(m/2))
      !
      new_top(1) = top(1)
      new_top(3:m/2) = top(2:m/2-1)
      new_top(2) = bot(1)
      new_bot(1:m/2-1) = bot(2:m/2)
      new_bot(m/2) = top(m/2)
      !
      top(:) = new_top
      bot(:) = new_bot
      !
      DEALLOCATE(new_top)
      DEALLOCATE(new_bot)
      !
    END SUBROUTINE
    !
END MODULE
