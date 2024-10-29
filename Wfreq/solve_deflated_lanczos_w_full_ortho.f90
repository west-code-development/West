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
SUBROUTINE solve_deflated_lanczos_w_full_ortho(nbnd_to_deflate, NRHS, NLSTEPS, b, alpha_s, beta_s, q_s, bnorm)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE gvect,                ONLY : gstart
  USE pwcom,                ONLY : npw,npwx
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
#if defined(__CUDA)
  USE west_gpu,             ONLY : beta,alpha,r,tmp_r,tmp_c
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: nbnd_to_deflate
  INTEGER,INTENT(IN) :: NRHS
  INTEGER,INTENT(IN) :: NLSTEPS
  COMPLEX(DP),INTENT(IN) :: b(npwx*npol,NRHS)
  REAL(DP),INTENT(OUT) :: alpha_s(NLSTEPS,NRHS)
  REAL(DP),INTENT(OUT) :: beta_s(NLSTEPS-1,NRHS)
  COMPLEX(DP),INTENT(OUT) :: q_s(npwx*npol,NRHS,NLSTEPS)
  REAL(DP),INTENT(OUT) :: bnorm(NRHS)
  !
  ! Workspace
  !
  INTEGER :: ia,ig,ip,il
  REAL(DP) :: reduce_r
  COMPLEX(DP) :: reduce_c
  !
#if !defined(__CUDA)
  REAL(DP),ALLOCATABLE :: beta(:)
  REAL(DP),ALLOCATABLE :: alpha(:)
  REAL(DP),ALLOCATABLE :: tmp_r(:)
  COMPLEX(DP),ALLOCATABLE :: tmp_c(:)
  COMPLEX(DP),ALLOCATABLE :: r(:,:)
#endif
  !
#if !defined(__CUDA)
  ALLOCATE(beta(NRHS))
  ALLOCATE(alpha(NRHS))
  ALLOCATE(tmp_r(NRHS))
  ALLOCATE(r(npwx*npol,NRHS))
  IF(.NOT. gamma_only) ALLOCATE(tmp_c(NRHS))
#endif
  !
  ! INIT
  !
  !$acc kernels present(r,b)
  r(:,:) = b
  !$acc end kernels
  !
#if defined(__CUDA)
  CALL start_clock_gpu("lan_H")
#else
  CALL start_clock( "lan_H" )
#endif
  !
  ! FROM R TO Q & BETA
  !
  IF(gamma_only) THEN
     !$acc parallel vector_length(1024) present(r,beta)
     !$acc loop
     DO ip = 1,NRHS
        reduce_r = 0.0_DP
        !$acc loop reduction(+:reduce_r)
        DO ig = 1,npw
           reduce_r = reduce_r+REAL(r(ig,ip),KIND=DP)**2+AIMAG(r(ig,ip))**2
        ENDDO
        IF(gstart == 2) THEN
           beta(ip) = 2.0_DP*reduce_r-REAL(r(1,ip),KIND=DP)**2
        ELSE
           beta(ip) = 2.0_DP*reduce_r
        ENDIF
     ENDDO
     !$acc end parallel
  ELSE
     !$acc parallel vector_length(1024) present(r,beta)
     !$acc loop
     DO ip = 1,NRHS
        reduce_r = 0.0_DP
        !$acc loop reduction(+:reduce_r)
        DO ig = 1,npw
           reduce_r = reduce_r+REAL(r(ig,ip),KIND=DP)**2+AIMAG(r(ig,ip))**2
           IF(noncolin) THEN
              reduce_r = reduce_r+REAL(r(ig+npwx,ip),KIND=DP)**2+AIMAG(r(ig+npwx,ip))**2
           ENDIF
        ENDDO
        beta(ip) = reduce_r
     ENDDO
     !$acc end parallel
  ENDIF
  !
  !$acc update host(beta)
  !
  CALL mp_sum(beta,intra_bgrp_comm)
  !
  DO ip = 1,NRHS
     beta(ip) = SQRT(beta(ip))
     bnorm(ip) = beta(ip)
  ENDDO
  !
  ! --- Lanczos iterations
  !
  DO il = 1,NLSTEPS
     !
     !$acc update device(beta)
     !
     !$acc parallel loop collapse(2) present(q_s,r,beta)
     DO ip = 1,NRHS
        DO ig = 1,npwx*npol
           q_s(ig,ip,il) = r(ig,ip)/beta(ip)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! NEW Q WAS FOUND !
     !
     ! APPLY H  q --> r
     !
     ! use h_psi_, i.e. h_psi without band parallelization, as west
     ! handles band parallelization by itself
     !
#if defined(__CUDA)
     !$acc host_data use_device(q_s,r)
     CALL h_psi__gpu(npwx,npw,NRHS,q_s(:,:,il),r)
     !$acc end host_data
#else
     CALL h_psi_(npwx,npw,NRHS,q_s(:,:,il),r)
#endif
     CALL apply_alpha_pc_to_m_wfcs(nbnd_to_deflate,NRHS,r,(1.0_DP,0.0_DP))
     !
     ! use beta
     !
     IF(il > 1) THEN
        !$acc parallel loop collapse(2) present(q_s,r,beta)
        DO ip = 1,NRHS
           DO ig = 1,npw
              r(ig,ip) = r(ig,ip)-q_s(ig,ip,il-1)*beta(ip)
              IF(noncolin) THEN
                 r(ig+npwx,ip) = r(ig+npwx,ip)-q_s(ig+npwx,ip,il-1)*beta(ip)
              ENDIF
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     ! get alpha
     !
     IF(gamma_only) THEN
        !$acc parallel vector_length(1024) present(q_s,r,alpha)
        !$acc loop
        DO ip = 1,NRHS
           reduce_r = 0.0_DP
           !$acc loop reduction(+:reduce_r)
           DO ig = 1,npw
              reduce_r = reduce_r+REAL(q_s(ig,ip,il),KIND=DP)*REAL(r(ig,ip),KIND=DP) &
              & +AIMAG(q_s(ig,ip,il))*AIMAG(r(ig,ip))
           ENDDO
           IF(gstart == 2) THEN
              alpha(ip) = 2.0_DP*reduce_r-REAL(q_s(1,ip,il),KIND=DP)*REAL(r(1,ip),KIND=DP)
           ELSE
              alpha(ip) = 2.0_DP*reduce_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel vector_length(1024) present(q_s,r,alpha)
        !$acc loop
        DO ip = 1,NRHS
           reduce_r = 0.0_DP
           !$acc loop reduction(+:reduce_r)
           DO ig = 1,npw
              reduce_r = reduce_r+REAL(q_s(ig,ip,il),KIND=DP)*REAL(r(ig,ip),KIND=DP) &
              & +AIMAG(q_s(ig,ip,il))*AIMAG(r(ig,ip))
              IF(noncolin) THEN
                 reduce_r = reduce_r+REAL(q_s(ig+npwx,ip,il),KIND=DP)*REAL(r(ig+npwx,ip),KIND=DP) &
                 & +AIMAG(q_s(ig+npwx,ip,il))*AIMAG(r(ig+npwx,ip))
              ENDIF
           ENDDO
           alpha(ip) = reduce_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc update host(alpha)
     !
     CALL mp_sum(alpha,intra_bgrp_comm)
     !
     DO ip = 1,NRHS
        alpha_s(il,ip) = alpha(ip)
     ENDDO
     !
     ! use alpha
     !
     !$acc update device(alpha)
     !
     !$acc parallel loop collapse(2) present(q_s,r,alpha)
     DO ip = 1,NRHS
        DO ig = 1,npw
           r(ig,ip) = r(ig,ip)-q_s(ig,ip,il)*alpha(ip)
           IF(noncolin) THEN
              r(ig+npwx,ip) = r(ig+npwx,ip)-q_s(ig+npwx,ip,il)*alpha(ip)
           ENDIF
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! Enforce full ortho
     !
     DO ia = 1,il-2
        IF(gamma_only) THEN
           !$acc parallel vector_length(1024) present(q_s,r,tmp_r)
           !$acc loop
           DO ip = 1,NRHS
              reduce_r = 0.0_DP
              !$acc loop reduction(+:reduce_r)
              DO ig = 1,npw
                 reduce_r = reduce_r+REAL(q_s(ig,ip,ia),KIND=DP)*REAL(r(ig,ip),KIND=DP) &
                 & +AIMAG(q_s(ig,ip,ia))*AIMAG(r(ig,ip))
              ENDDO
              IF(gstart == 2) THEN
                 tmp_r(ip) = 2.0_DP*reduce_r-REAL(q_s(1,ip,ia),KIND=DP)*REAL(r(1,ip),KIND=DP)
              ELSE
                 tmp_r(ip) = 2.0_DP*reduce_r
              ENDIF
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(tmp_r)
           CALL mp_sum(tmp_r,intra_bgrp_comm)
           !$acc end host_data
           !
           !$acc parallel loop collapse(2) present(q_s,r,tmp_r)
           DO ip = 1,NRHS
              DO ig = 1,npw
                 r(ig,ip) = r(ig,ip)-q_s(ig,ip,ia)*tmp_r(ip)
                 IF(noncolin) THEN
                    r(ig+npwx,ip) = r(ig+npwx,ip)-q_s(ig+npwx,ip,ia)*tmp_r(ip)
                 ENDIF
              ENDDO
           ENDDO
           !$acc end parallel
        ELSE
           !$acc parallel vector_length(1024) present(q_s,r,tmp_c)
           !$acc loop
           DO ip = 1,NRHS
              reduce_c = 0.0_DP
              !$acc loop reduction(+:reduce_c)
              DO ig = 1,npw
                 reduce_c = reduce_c+CONJG(q_s(ig,ip,ia))*r(ig,ip)
                 IF(noncolin) THEN
                    reduce_c = reduce_c+CONJG(q_s(ig+npwx,ip,ia))*r(ig+npwx,ip)
                 ENDIF
              ENDDO
              tmp_c(ip) = reduce_c
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(tmp_c)
           CALL mp_sum(tmp_c,intra_bgrp_comm)
           !$acc end host_data
           !
           !$acc parallel loop collapse(2) present(q_s,r,tmp_c)
           DO ip = 1,NRHS
              DO ig = 1,npw
                 r(ig,ip) = r(ig,ip)-q_s(ig,ip,ia)*tmp_c(ip)
                 IF(noncolin) THEN
                    r(ig+npwx,ip) = r(ig+npwx,ip)-q_s(ig+npwx,ip,ia)*tmp_c(ip)
                 ENDIF
              ENDDO
           ENDDO
           !$acc end parallel
        ENDIF
     ENDDO
     !
     ! get beta
     !
     IF(gamma_only) THEN
        !$acc parallel vector_length(1024) present(r,beta)
        !$acc loop
        DO ip = 1,NRHS
           reduce_r = 0.0_DP
           !$acc loop reduction(+:reduce_r)
           DO ig = 1,npw
              reduce_r = reduce_r+REAL(r(ig,ip),KIND=DP)**2+AIMAG(r(ig,ip))**2
           ENDDO
           IF(gstart == 2) THEN
              beta(ip) = 2.0_DP*reduce_r-REAL(r(1,ip),KIND=DP)**2
           ELSE
              beta(ip) = 2.0_DP*reduce_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel vector_length(1024) present(r,beta)
        !$acc loop
        DO ip = 1,NRHS
           reduce_r = 0.0_DP
           !$acc loop reduction(+:reduce_r)
           DO ig = 1,npw
              reduce_r = reduce_r+REAL(r(ig,ip),KIND=DP)**2+AIMAG(r(ig,ip))**2
              IF(noncolin) THEN
                 reduce_r = reduce_r+REAL(r(ig+npwx,ip),KIND=DP)**2+AIMAG(r(ig+npwx,ip))**2
              ENDIF
           ENDDO
           beta(ip) = reduce_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc update host(beta)
     !
     CALL mp_sum(beta,intra_bgrp_comm)
     !
     DO ip = 1,NRHS
        beta(ip) = SQRT(beta(ip))
     ENDDO
     !
     IF(il < NLSTEPS) THEN
        DO ip = 1,NRHS
           beta_s(il,ip) = beta(ip)
        ENDDO
     ENDIF
     !
  ENDDO
  !
#if !defined(__CUDA)
  DEALLOCATE(beta)
  DEALLOCATE(alpha)
  DEALLOCATE(tmp_r)
  DEALLOCATE(r)
  IF(.NOT. gamma_only) DEALLOCATE(tmp_c)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu("lan_H")
#else
  CALL stop_clock("lan_H")
#endif
  !
END SUBROUTINE
