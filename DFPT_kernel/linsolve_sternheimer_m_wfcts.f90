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
SUBROUTINE linsolve_sternheimer_m_wfcts ( nbndval, m, b, x, e, eprec, tr2, ierr )
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                 A * | x_i > = | b_i >       i=(1:m)
  !
  !     with :      A = ( H - E_i + alpha * P_v )
  !
  !                 where H is a complex hermitean matrix (H_{SCF}), E_v is a real scalar (energy of band v)
  !                 x and b are complex vectors
  !
  !     on input:
  !
  !                 nbndval  integer  number of valence bands.
  !
  !                 m        integer  number of wfcts to process.
  !
  !                 b        contains the right hand side vector
  !                          of the system.
  !
  !                 x        contains an estimate of the solution
  !                          vector.
  !
  !                 e        real     unperturbed eigenvalues.
  !
  !                 eprec    real     preconditioned energies.
  !
  !                 tr2      real     threshold
  !
  !     on output:  x        contains the refined estimate of the
  !                          solution vector.
  !
  !                 ierr     integer  error (if /=0 something went wrong)
  !
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE gvect,                ONLY : gstart
  USE control_flags,        ONLY : gamma_only
  USE westcom,              ONLY : n_dfpt_maxiter,alphapv_dfpt
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol,noncolin
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN)  :: nbndval, m
  REAL(DP),INTENT(IN)  :: e(m), eprec(m)
  COMPLEX(DP),INTENT(IN)  :: b (npwx*npol, m)
  COMPLEX(DP),INTENT(INOUT) :: x (npwx*npol, m)
  REAL(DP),INTENT(IN) :: tr2
  INTEGER, INTENT(OUT) :: ierr
  !
  ! Workspace
  !
  INTEGER :: iter, ibnd, lbnd
  COMPLEX(DP),ALLOCATABLE :: g(:,:), t(:,:), h(:,:), hold(:,:)
  COMPLEX(DP) ::  dcgamma, dclambda
  REAL(DP),ALLOCATABLE :: rho(:), rhoold(:), eu(:), a(:), c(:)
  REAL(DP) :: anorm
  LOGICAL,ALLOCATABLE :: is_conv(:)
  REAL(KIND=DP),EXTERNAL :: DDOT
  COMPLEX(KIND=DP),EXTERNAL :: ZDOTC
  !
  CALL start_clock ('linstern')
  !
  ! Zeros
  !
  ALLOCATE( g(npwx*npol,m), t(npwx*npol,m), h(npwx*npol,m), hold(npwx*npol,m) )
  ALLOCATE( rho(m), rhoold(m), eu(m), a(m), c(m), is_conv(m) )
  !
  g = 0.0_DP
  t = 0.0_DP
  h = 0.0_DP
  hold = 0.0_DP
  is_conv=.FALSE.
  ierr=0
  !
  ! Step 1, initialization of the loop
  !
  CALL apply_sternheimerop_to_m_wfcs(nbndval, x, g, e, alphapv_dfpt, m)
  !
  DO ibnd=1,m
     CALL ZAXPY(npw,(-1._DP,0._DP),b(1,ibnd),1,g(1,ibnd),1)
  ENDDO
  IF(noncolin) THEN
     DO ibnd=1,m
        CALL ZAXPY(npw,(-1._DP,0._DP),b(npwx+1,ibnd),1,g(npwx+1,ibnd),1)
     ENDDO
  ENDIF
  !
  ! Loop
  !
  DO iter = 1, n_dfpt_maxiter
     !
     lbnd = 0
     DO ibnd = 1, m
        IF (is_conv (ibnd) ) CYCLE
        lbnd = lbnd+1
        !
        CALL precondition_m_wfcts( 1, g(1,ibnd), h(1,ibnd), eprec(ibnd) )
        !
        IF (gamma_only) THEN
           rho(lbnd)=2.0_DP * DDOT(2*npw,h(1,ibnd),1,g(1,ibnd),1)
           IF(gstart==2) THEN
              rho(lbnd)=rho(lbnd)-REAL(h(1,ibnd),KIND=DP)*REAL(g(1,ibnd),KIND=DP)
           ENDIF
        ELSE
           rho(lbnd) = ZDOTC (npw, h(1,ibnd), 1, g(1,ibnd), 1)
           IF(noncolin) rho(lbnd) = rho(lbnd) + ZDOTC (npw, h(1+npwx,ibnd), 1, g(1+npwx,ibnd), 1)
        ENDIF
        !
     ENDDO
     !
     CALL mp_sum( rho(1:lbnd) , intra_bgrp_comm )
     !
     DO ibnd = m, 1, -1
        IF (is_conv(ibnd) ) CYCLE
        rho(ibnd)=rho(lbnd)
        lbnd = lbnd -1
        anorm = SQRT (rho (ibnd) )
        IF (anorm .LT. tr2) is_conv (ibnd) = .TRUE.
     ENDDO
     !
     !
     IF (ALL(is_conv(:))) EXIT
     !
     !
     lbnd = 0
     DO ibnd = 1, m
        IF (is_conv (ibnd) ) CYCLE
        !
        !
        CALL DSCAL (2*npwx*npol, - 1._DP, h (1, ibnd), 1)
        IF (iter .NE. 1) THEN
           dcgamma = CMPLX( rho (ibnd) / rhoold (ibnd), 0.0_DP, KIND=DP)
           CALL ZAXPY (npwx*npol, dcgamma, hold (1, ibnd), 1, h (1, ibnd), 1)
        ENDIF
        !
        lbnd = lbnd+1
        CALL ZCOPY(npwx*npol, h (1, ibnd), 1, hold (1, lbnd), 1)
        eu ( lbnd ) = e (ibnd )
        !
     ENDDO
     !
     !
     CALL apply_sternheimerop_to_m_wfcs(nbndval, hold, t, eu, alphapv_dfpt, lbnd)
     !
     lbnd=0
     DO ibnd = 1, m
        IF ( is_conv (ibnd) ) CYCLE
        lbnd=lbnd+1
        IF (gamma_only) THEN
           a(lbnd) = 2.0_DP * DDOT(2*npw,h(1,ibnd),1,g(1,ibnd),1)
           c(lbnd) = 2.0_DP * DDOT(2*npw,h(1,ibnd),1,t(1,lbnd),1)
           IF (gstart == 2) THEN
              a(lbnd)=a(lbnd)-REAL(h(1,ibnd),KIND=DP)*REAL(g(1,ibnd),KIND=DP)
              c(lbnd)=c(lbnd)-REAL(h(1,ibnd),KIND=DP)*REAL(t(1,lbnd),KIND=DP)
           ENDIF
        ELSE
           a(lbnd) = ZDOTC (npw, h(1,ibnd), 1, g(1,ibnd), 1)
           c(lbnd) = ZDOTC (npw, h(1,ibnd), 1, t(1,lbnd), 1)
           IF(noncolin) THEN
              a(lbnd) = a(lbnd) + ZDOTC (npw, h(1+npwx,ibnd), 1, g(1+npwx,ibnd), 1)
              c(lbnd) = c(lbnd) + ZDOTC (npw, h(1+npwx,ibnd), 1, t(1+npwx,lbnd), 1)
           ENDIF
        ENDIF
     ENDDO
     !
     CALL mp_sum(  a(1:lbnd), intra_bgrp_comm )
     CALL mp_sum(  c(1:lbnd), intra_bgrp_comm )
     !
     lbnd=0
     DO ibnd = 1, m
        IF ( is_conv (ibnd) ) CYCLE
        lbnd=lbnd+1
        dclambda = CMPLX( - a(lbnd) / c(lbnd), 0.0_DP, KIND=DP)
        !
        CALL ZAXPY(npwx*npol, dclambda, h(1,ibnd), 1, x(1,ibnd), 1)
        !
        CALL ZAXPY(npwx*npol, dclambda, t(1,lbnd), 1, g(1,ibnd), 1)
        !
        CALL ZCOPY(npwx*npol, h(1,ibnd), 1, hold(1,ibnd), 1)
        rhoold (ibnd) = rho (ibnd)
        !
     ENDDO ! on bands
     !
  ENDDO ! on iterations
  !
  DEALLOCATE( g, t, h, hold )
  DEALLOCATE( rho, rhoold, eu, a, c)
  !
  ierr=m-COUNT( is_conv(:) )
  !
  DEALLOCATE( is_conv )
  !
  CALL stop_clock ('linstern')
  !
END SUBROUTINE
!
#if defined(__CUDA)
!-----------------------------------------------------------------------
SUBROUTINE linsolve_sternheimer_m_wfcts_gpu(nbndval,m,b_d,x_d,e_d,eprec_d,tr2,ierr)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                 A * | x_i > = | b_i >       i=(1:m)
  !
  !     with :      A = ( H - E_i + alpha * P_v )
  !
  !                 where H is a complex hermitean matrix (H_{SCF}), E_v is a real scalar (energy of band v)
  !                 x and b are complex vectors
  !
  !     on input:
  !
  !                 nbndval  integer  number of valence bands.
  !
  !                 m        integer  number of wfcts to process.
  !
  !                 b        contains the right hand side vector
  !                          of the system.
  !
  !                 x        contains an estimate of the solution
  !                          vector.
  !
  !                 e        real     unperturbed eigenvalues.
  !
  !                 eprec    real     preconditioned energies.
  !
  !                 tr2      real     threshold
  !
  !     on output:  x        contains the refined estimate of the
  !                          solution vector.
  !
  !                 ierr     integer  error (if /=0 something went wrong)
  !
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : intra_bgrp_comm,nproc_bgrp
  USE mp,                   ONLY : mp_sum
  USE gvect,                ONLY : gstart
  USE control_flags,        ONLY : gamma_only
  USE westcom,              ONLY : n_dfpt_maxiter,alphapv_dfpt
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol,noncolin
  USE wvfct_gpum,           ONLY : g2kin_d
  USE west_cuda,            ONLY : g_d,t_d,h_d,eu_d,a_d,c_d,rho_d,rhoold_d,ibnd_d
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nbndval
  INTEGER, INTENT(IN) :: m
  REAL(DP), DEVICE, INTENT(IN) :: e_d(m)
  REAL(DP), DEVICE, INTENT(IN) :: eprec_d(m)
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: b_d(npwx*npol,m)
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: x_d(npwx*npol,m)
  REAL(DP), INTENT(IN) :: tr2
  INTEGER, INTENT(OUT) :: ierr
  !
  ! Workspace
  !
  LOGICAL, ALLOCATABLE :: is_conv(:)
  INTEGER :: ig, iter, ibnd, lbnd, nbnd_todo
  REAL(DP), ALLOCATABLE :: rho(:)
  REAL(DP) :: anorm
  REAL(DP) :: tmp_r, tmp2_r
  !
  CALL start_clock_gpu('linstern')
  !
  ALLOCATE(rho(m))
  ALLOCATE(is_conv(m))
  !
  is_conv = .FALSE.
  ierr = 0
  !
  ! Step 1, initialization of the loop
  !
  CALL apply_sternheimerop_to_m_wfcs(nbndval,x_d,g_d,e_d,alphapv_dfpt,m)
  !
  !$acc parallel loop collapse(2)
  DO ibnd = 1,m
     DO ig = 1,npw
        g_d(ig,ibnd) = g_d(ig,ibnd)-b_d(ig,ibnd)
        IF(noncolin) THEN
           g_d(npwx+ig,ibnd) = g_d(npwx+ig,ibnd)-b_d(npwx+ig,ibnd)
        ENDIF
     ENDDO
  ENDDO
  !$acc end parallel
  !
  ! From now on b_d is no longer needed, use it to store hold
  !
  CALL precompute_lbnd(m,is_conv,nbnd_todo,ibnd_d)
  !
  ! Loop
  !
  DO iter = 1,n_dfpt_maxiter
     !
     ! Preconditioning
     !
     h_d = 0.0_DP
     !
     !$acc parallel loop collapse(2)
     DO lbnd = 1,nbnd_todo
        DO ig = 1,npw
           ibnd = ibnd_d(lbnd)
           h_d(ig,ibnd) = g_d(ig,ibnd)/MAX(1.0_DP,g2kin_d(ig)/eprec_d(ibnd))
           IF(noncolin) THEN
              h_d(npwx+ig,ibnd) = g_d(npwx+ig,ibnd)/MAX(1.0_DP,g2kin_d(ig)/eprec_d(ibnd))
           ENDIF
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gamma_only) THEN
        !$acc parallel
        !$acc loop
        DO lbnd = 1,nbnd_todo
           ibnd = ibnd_d(lbnd)
           tmp_r = 0.0_DP
           !
           !$acc loop reduction(+:tmp_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h_d(ig,ibnd),KIND=DP)*REAL(g_d(ig,ibnd),KIND=DP) &
              & +AIMAG(h_d(ig,ibnd))*AIMAG(g_d(ig,ibnd))
           ENDDO
           !
           IF(gstart == 2) THEN
              rho_d(lbnd) = 2.0_DP*tmp_r-REAL(h_d(1,ibnd),KIND=DP)*REAL(g_d(1,ibnd),KIND=DP)
           ELSE
              rho_d(lbnd) = 2.0_DP*tmp_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel
        !$acc loop
        DO lbnd = 1,nbnd_todo
           ibnd = ibnd_d(lbnd)
           tmp_r = 0.0_DP
           !
           !$acc loop reduction(+:tmp_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h_d(ig,ibnd),KIND=DP)*REAL(g_d(ig,ibnd),KIND=DP) &
              & +AIMAG(h_d(ig,ibnd))*AIMAG(g_d(ig,ibnd))
              IF(noncolin) THEN
                 tmp_r = tmp_r+REAL(h_d(npwx+ig,ibnd),KIND=DP)*REAL(g_d(npwx+ig,ibnd),KIND=DP) &
                 & +AIMAG(h_d(npwx+ig,ibnd))*AIMAG(g_d(npwx+ig,ibnd))
              ENDIF
           ENDDO
           !
           rho_d(lbnd) = tmp_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     rho(1:nbnd_todo) = rho_d(1:nbnd_todo)
     !
     CALL mp_sum(rho(1:nbnd_todo),intra_bgrp_comm)
     !
     lbnd = nbnd_todo
     DO ibnd = m,1,-1
        IF(is_conv(ibnd)) CYCLE
        rho(ibnd) = rho(lbnd)
        lbnd = lbnd-1
        anorm = SQRT(rho(ibnd))
        IF(anorm < tr2) is_conv(ibnd) = .TRUE.
     ENDDO
     !
     IF(ALL(is_conv(:))) EXIT
     !
     rho_d(1:m) = rho(1:m)
     !
     CALL precompute_lbnd(m,is_conv,nbnd_todo,ibnd_d)
     !
     ! b_d used to store hold_d
     !
     DO lbnd = 1,nbnd_todo
        !$acc parallel loop
        DO ig = 1,npwx*npol
           ibnd = ibnd_d(lbnd)
           !
           IF(iter /= 1) THEN
              h_d(ig,ibnd) = -h_d(ig,ibnd)+rho_d(ibnd)/rhoold_d(ibnd)*b_d(ig,ibnd)
           ELSE
              h_d(ig,ibnd) = -h_d(ig,ibnd)
           ENDIF
           !
           b_d(ig,lbnd) = h_d(ig,ibnd)
        ENDDO
        !$acc end parallel
     ENDDO
     !
     !$acc parallel loop
     DO lbnd = 1,nbnd_todo
        ibnd = ibnd_d(lbnd)
        eu_d(lbnd) = e_d(ibnd)
     ENDDO
     !$acc end parallel
     !
     ! b_d used to store hold_d
     !
     CALL apply_sternheimerop_to_m_wfcs(nbndval,b_d,t_d,eu_d,alphapv_dfpt,nbnd_todo)
     !
     IF(gamma_only) THEN
        !$acc parallel
        !$acc loop
        DO lbnd = 1,nbnd_todo
           ibnd = ibnd_d(lbnd)
           tmp_r = 0.0_DP
           tmp2_r = 0.0_DP
           !
           !$acc loop reduction(+:tmp_r,tmp2_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h_d(ig,ibnd),KIND=DP)*REAL(g_d(ig,ibnd),KIND=DP) &
              & +AIMAG(h_d(ig,ibnd))*AIMAG(g_d(ig,ibnd))
              tmp2_r = tmp2_r+REAL(h_d(ig,ibnd),KIND=DP)*REAL(t_d(ig,lbnd),KIND=DP) &
              & +AIMAG(h_d(ig,ibnd))*AIMAG(t_d(ig,lbnd))
           ENDDO
           !
           IF(gstart == 2) THEN
              a_d(lbnd) = 2.0_DP*tmp_r-REAL(h_d(1,ibnd),KIND=DP)*REAL(g_d(1,ibnd),KIND=DP)
              c_d(lbnd) = 2.0_DP*tmp2_r-REAL(h_d(1,ibnd),KIND=DP)*REAL(t_d(1,lbnd),KIND=DP)
           ELSE
              a_d(lbnd) = 2.0_DP*tmp_r
              c_d(lbnd) = 2.0_DP*tmp2_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel
        !$acc loop
        DO lbnd = 1,nbnd_todo
           ibnd = ibnd_d(lbnd)
           tmp_r = 0.0_DP
           tmp2_r = 0.0_DP
           !
           !$acc loop reduction(+:tmp_r,tmp2_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h_d(ig,ibnd),KIND=DP)*REAL(g_d(ig,ibnd),KIND=DP) &
              & +AIMAG(h_d(ig,ibnd))*AIMAG(g_d(ig,ibnd))
              tmp2_r = tmp2_r+REAL(h_d(ig,ibnd),KIND=DP)*REAL(t_d(ig,lbnd),KIND=DP) &
              & +AIMAG(h_d(ig,ibnd))*AIMAG(t_d(ig,lbnd))
              IF(noncolin) THEN
                 tmp_r = tmp_r+REAL(h_d(npwx+ig,ibnd),KIND=DP)*REAL(g_d(npwx+ig,ibnd),KIND=DP) &
                 & +AIMAG(h_d(npwx+ig,ibnd))*AIMAG(g_d(npwx+ig,ibnd))
                 tmp2_r = tmp2_r+REAL(h_d(npwx+ig,ibnd),KIND=DP)*REAL(t_d(npwx+ig,lbnd),KIND=DP) &
                 & +AIMAG(h_d(npwx+ig,ibnd))*AIMAG(t_d(npwx+ig,lbnd))
              ENDIF
           ENDDO
           !
           a_d(lbnd) = tmp_r
           c_d(lbnd) = tmp2_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     IF(nproc_bgrp > 1) THEN
        CALL mp_sum(a_d(1:nbnd_todo),intra_bgrp_comm)
        CALL mp_sum(c_d(1:nbnd_todo),intra_bgrp_comm)
     ENDIF
     !
     ! b_d used to store hold_d
     !
     !$acc parallel loop collapse(2)
     DO lbnd = 1,nbnd_todo
        DO ig = 1,npwx*npol
           ibnd = ibnd_d(lbnd)
           x_d(ig,ibnd) = x_d(ig,ibnd)-a_d(lbnd)/c_d(lbnd)*h_d(ig,ibnd)
           g_d(ig,ibnd) = g_d(ig,ibnd)-a_d(lbnd)/c_d(lbnd)*t_d(ig,lbnd)
           b_d(ig,ibnd) = h_d(ig,ibnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc parallel loop
     DO lbnd = 1,nbnd_todo
        ibnd = ibnd_d(lbnd)
        rhoold_d(ibnd) = rho_d(ibnd)
     ENDDO
     !$acc end parallel
  ENDDO ! on iterations
  !
  ierr = m-COUNT(is_conv(:))
  !
  DEALLOCATE(rho)
  DEALLOCATE(is_conv)
  !
  CALL stop_clock_gpu('linstern')
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE precompute_lbnd(m,is_conv,nbnd_todo,ibnd_d)
  !----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m
  LOGICAL, INTENT(IN) :: is_conv(m)
  INTEGER, INTENT(OUT) :: nbnd_todo
  INTEGER, DEVICE, INTENT(OUT) :: ibnd_d(m)
  !
  ! Workspace
  !
  INTEGER :: ibnd
  INTEGER, ALLOCATABLE :: ibnd_h(:)
  !
  ALLOCATE(ibnd_h(m))
  !
  ibnd_h = 0
  nbnd_todo = 0
  !
  DO ibnd = 1, m
     IF(is_conv(ibnd)) CYCLE
     nbnd_todo = nbnd_todo+1
     ibnd_h(nbnd_todo) = ibnd
  ENDDO
  !
  ibnd_d = ibnd_h
  !
  DEALLOCATE(ibnd_h)
  !
END SUBROUTINE
#endif
