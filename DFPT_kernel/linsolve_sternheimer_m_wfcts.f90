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
  USE wvfct,                ONLY : g2kin
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
  INTEGER :: ig
  COMPLEX(DP),ALLOCATABLE :: g(:,:), t(:,:), h(:,:), hold(:,:)
  COMPLEX(DP) ::  dcgamma, dclambda
  REAL(DP),ALLOCATABLE :: rho(:), rhoold(:), eu(:), a(:), c(:)
  REAL(DP) :: energy, norma, anorm
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
  IF(npol==2) THEN
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
