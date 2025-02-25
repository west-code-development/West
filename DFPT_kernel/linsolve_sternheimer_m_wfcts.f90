!
! Copyright (C) 2015-2025 M. Govoni
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
SUBROUTINE linsolve_sternheimer_m_wfcts(nbndval,m,b,x,e,eprec,tr2,ierr)
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
  USE wvfct,                ONLY : g2kin
#if defined(__CUDA)
  USE west_gpu,             ONLY : is_conv,l2i_map,eu,a,c,rho,rhoold,g,t,h
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: nbndval
  INTEGER,INTENT(IN) :: m
  REAL(DP),INTENT(IN) :: e(m)
  REAL(DP),INTENT(IN) :: eprec(m)
  COMPLEX(DP),INTENT(INOUT) :: b(npwx*npol,m)
  COMPLEX(DP),INTENT(INOUT) :: x(npwx*npol,m)
  REAL(DP),INTENT(IN) :: tr2
  INTEGER, INTENT(OUT) :: ierr
  !
  ! Workspace
  !
  INTEGER :: ig, iter, ibnd, lbnd, not_conv
  REAL(DP) :: anorm
  REAL(DP) :: tmp_r, tmp2_r
#if !defined(__CUDA)
  COMPLEX(DP),ALLOCATABLE :: g(:,:),t(:,:),h(:,:)
  REAL(DP),ALLOCATABLE :: rho(:),rhoold(:),eu(:),a(:),c(:)
  INTEGER,ALLOCATABLE :: l2i_map(:)
  LOGICAL,ALLOCATABLE :: is_conv(:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('linstern')
#else
  CALL start_clock('linstern')
#endif
  !
#if !defined(__CUDA)
  ALLOCATE(g(npwx*npol,m))
  ALLOCATE(t(npwx*npol,m))
  ALLOCATE(h(npwx*npol,m))
  ALLOCATE(rho(m))
  ALLOCATE(rhoold(m))
  ALLOCATE(eu(m))
  ALLOCATE(a(m))
  ALLOCATE(c(m))
  ALLOCATE(l2i_map(m))
  ALLOCATE(is_conv(m))
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('linstern')
#else
  CALL start_clock('linstern')
#endif
  !
  ierr = 0
  !
  !$acc kernels present(g)
  g(:,:) = 0._DP
  !$acc end kernels
  !
  !$acc kernels present(t)
  t(:,:) = 0._DP
  !$acc end kernels
  !
  !$acc kernels present(h)
  h(:,:) = 0._DP
  !$acc end kernels
  !
  !$acc kernels present(is_conv)
  is_conv(:) = .FALSE.
  !$acc end kernels
  !
  ! Step 1, initialization of the loop
  !
  CALL apply_sternheimerop_to_m_wfcs(nbndval,x,g,e,alphapv_dfpt,m)
  !
  !$acc parallel loop collapse(2) present(g,b)
  DO ibnd = 1,m
     DO ig = 1,npw
        g(ig,ibnd) = g(ig,ibnd)-b(ig,ibnd)
        IF(noncolin) THEN
           g(npwx+ig,ibnd) = g(npwx+ig,ibnd)-b(npwx+ig,ibnd)
        ENDIF
     ENDDO
  ENDDO
  !$acc end parallel
  !
  ! From now on b is no longer needed, use it to store hold
  !
  CALL map_lbnd2ibnd(m,is_conv,not_conv,l2i_map)
  !
  ! Loop
  !
  DO iter = 1,n_dfpt_maxiter
     !
     ! Preconditioning
     !
     !$acc parallel loop collapse(2) present(l2i_map,h,g,g2kin,eprec)
     DO lbnd = 1,not_conv
        DO ig = 1,npw
           ibnd = l2i_map(lbnd)
           h(ig,ibnd) = g(ig,ibnd)/MAX(1._DP,g2kin(ig)/eprec(ibnd))
           IF(noncolin) THEN
              h(npwx+ig,ibnd) = g(npwx+ig,ibnd)/MAX(1._DP,g2kin(ig)/eprec(ibnd))
           ENDIF
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gamma_only) THEN
        !$acc parallel present(l2i_map,h,g,rho)
        !$acc loop
        DO lbnd = 1,not_conv
           ibnd = l2i_map(lbnd)
           tmp_r = 0._DP
           !
           !$acc loop reduction(+:tmp_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h(ig,ibnd),KIND=DP)*REAL(g(ig,ibnd),KIND=DP) &
              & +AIMAG(h(ig,ibnd))*AIMAG(g(ig,ibnd))
           ENDDO
           !
           IF(gstart == 2) THEN
              rho(lbnd) = 2._DP*tmp_r-REAL(h(1,ibnd),KIND=DP)*REAL(g(1,ibnd),KIND=DP)
           ELSE
              rho(lbnd) = 2._DP*tmp_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel present(l2i_map,h,g,rho)
        !$acc loop
        DO lbnd = 1,not_conv
           ibnd = l2i_map(lbnd)
           tmp_r = 0._DP
           !
           !$acc loop reduction(+:tmp_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h(ig,ibnd),KIND=DP)*REAL(g(ig,ibnd),KIND=DP) &
              & +AIMAG(h(ig,ibnd))*AIMAG(g(ig,ibnd))
              IF(noncolin) THEN
                 tmp_r = tmp_r+REAL(h(npwx+ig,ibnd),KIND=DP)*REAL(g(npwx+ig,ibnd),KIND=DP) &
                 & +AIMAG(h(npwx+ig,ibnd))*AIMAG(g(npwx+ig,ibnd))
              ENDIF
           ENDDO
           !
           rho(lbnd) = tmp_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc host_data use_device(rho)
     CALL mp_sum(rho(1:not_conv),intra_bgrp_comm)
     !$acc end host_data
     !
     lbnd = not_conv
     !
     !$acc serial present(is_conv,rho) copy(lbnd,not_conv)
     DO ibnd = m,1,-1
        IF(is_conv(ibnd)) CYCLE
        rho(ibnd) = rho(lbnd)
        lbnd = lbnd-1
        anorm = rho(ibnd)
        IF(anorm < tr2) is_conv(ibnd) = .TRUE.
     ENDDO
     !$acc end serial
     !
     CALL map_lbnd2ibnd(m,is_conv,not_conv,l2i_map)
     !
     IF(not_conv == 0) EXIT
     !
     ! hold stored in b
     !
     DO lbnd = 1,not_conv
        !$acc parallel loop present(l2i_map,h,rho,rhoold,b)
        DO ig = 1,npwx*npol
           ibnd = l2i_map(lbnd)
           !
           IF(iter /= 1) THEN
              h(ig,ibnd) = -h(ig,ibnd)+rho(ibnd)/rhoold(ibnd)*b(ig,ibnd)
           ELSE
              h(ig,ibnd) = -h(ig,ibnd)
           ENDIF
           !
           b(ig,lbnd) = h(ig,ibnd)
        ENDDO
        !$acc end parallel
     ENDDO
     !
     !$acc parallel loop present(l2i_map,eu,e)
     DO lbnd = 1,not_conv
        ibnd = l2i_map(lbnd)
        eu(lbnd) = e(ibnd)
     ENDDO
     !$acc end parallel
     !
     ! hold stored in b
     !
     CALL apply_sternheimerop_to_m_wfcs(nbndval,b,t,eu,alphapv_dfpt,not_conv)
     !
     IF(gamma_only) THEN
        !$acc parallel present(l2i_map,h,g,t,a,c)
        !$acc loop
        DO lbnd = 1,not_conv
           ibnd = l2i_map(lbnd)
           tmp_r = 0._DP
           tmp2_r = 0._DP
           !
           !$acc loop reduction(+:tmp_r,tmp2_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h(ig,ibnd),KIND=DP)*REAL(g(ig,ibnd),KIND=DP) &
              & +AIMAG(h(ig,ibnd))*AIMAG(g(ig,ibnd))
              tmp2_r = tmp2_r+REAL(h(ig,ibnd),KIND=DP)*REAL(t(ig,lbnd),KIND=DP) &
              & +AIMAG(h(ig,ibnd))*AIMAG(t(ig,lbnd))
           ENDDO
           !
           IF(gstart == 2) THEN
              a(lbnd) = 2._DP*tmp_r-REAL(h(1,ibnd),KIND=DP)*REAL(g(1,ibnd),KIND=DP)
              c(lbnd) = 2._DP*tmp2_r-REAL(h(1,ibnd),KIND=DP)*REAL(t(1,lbnd),KIND=DP)
           ELSE
              a(lbnd) = 2._DP*tmp_r
              c(lbnd) = 2._DP*tmp2_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel present(l2i_map,h,g,t,a,c)
        !$acc loop
        DO lbnd = 1,not_conv
           ibnd = l2i_map(lbnd)
           tmp_r = 0._DP
           tmp2_r = 0._DP
           !
           !$acc loop reduction(+:tmp_r,tmp2_r)
           DO ig = 1,npw
              tmp_r = tmp_r+REAL(h(ig,ibnd),KIND=DP)*REAL(g(ig,ibnd),KIND=DP) &
              & +AIMAG(h(ig,ibnd))*AIMAG(g(ig,ibnd))
              tmp2_r = tmp2_r+REAL(h(ig,ibnd),KIND=DP)*REAL(t(ig,lbnd),KIND=DP) &
              & +AIMAG(h(ig,ibnd))*AIMAG(t(ig,lbnd))
              IF(noncolin) THEN
                 tmp_r = tmp_r+REAL(h(npwx+ig,ibnd),KIND=DP)*REAL(g(npwx+ig,ibnd),KIND=DP) &
                 & +AIMAG(h(npwx+ig,ibnd))*AIMAG(g(npwx+ig,ibnd))
                 tmp2_r = tmp2_r+REAL(h(npwx+ig,ibnd),KIND=DP)*REAL(t(npwx+ig,lbnd),KIND=DP) &
                 & +AIMAG(h(npwx+ig,ibnd))*AIMAG(t(npwx+ig,lbnd))
              ENDIF
           ENDDO
           !
           a(lbnd) = tmp_r
           c(lbnd) = tmp2_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc host_data use_device(a,c)
     CALL mp_sum(a(1:not_conv),intra_bgrp_comm)
     CALL mp_sum(c(1:not_conv),intra_bgrp_comm)
     !$acc end host_data
     !
     ! hold stored in b
     !
     !$acc parallel loop collapse(2) present(l2i_map,x,a,c,h,g,t,b)
     DO lbnd = 1,not_conv
        DO ig = 1,npwx*npol
           ibnd = l2i_map(lbnd)
           x(ig,ibnd) = x(ig,ibnd)-a(lbnd)/c(lbnd)*h(ig,ibnd)
           g(ig,ibnd) = g(ig,ibnd)-a(lbnd)/c(lbnd)*t(ig,lbnd)
           b(ig,ibnd) = h(ig,ibnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc parallel loop present(l2i_map,rhoold,rho)
     DO lbnd = 1,not_conv
        ibnd = l2i_map(lbnd)
        rhoold(ibnd) = rho(ibnd)
     ENDDO
     !$acc end parallel
  ENDDO ! on iterations
  !
  ierr = not_conv
  !
#if !defined(__CUDA)
  DEALLOCATE(g)
  DEALLOCATE(t)
  DEALLOCATE(h)
  DEALLOCATE(rho)
  DEALLOCATE(rhoold)
  DEALLOCATE(eu)
  DEALLOCATE(a)
  DEALLOCATE(c)
  DEALLOCATE(l2i_map)
  DEALLOCATE(is_conv)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('linstern')
#else
  CALL stop_clock('linstern')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE map_lbnd2ibnd(m,is_conv,not_conv,l2i_map)
  !----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m
  LOGICAL, INTENT(IN) :: is_conv(m)
  INTEGER, INTENT(OUT) :: not_conv
  INTEGER, INTENT(OUT) :: l2i_map(m)
  !
  ! Workspace
  !
  INTEGER :: ibnd
  !
  not_conv = 0
  !
  !$acc serial present(is_conv,l2i_map) copy(not_conv)
  DO ibnd = 1, m
     IF(is_conv(ibnd)) CYCLE
     not_conv = not_conv+1
     l2i_map(not_conv) = ibnd
  ENDDO
  !$acc end serial
  !
END SUBROUTINE
