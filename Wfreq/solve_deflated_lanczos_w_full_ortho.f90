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
  INTEGER :: ia,ip,il
  !
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  REAL(DP),ALLOCATABLE :: beta(:),alpha(:),raux(:)
  COMPLEX(DP),ALLOCATABLE :: ZA(:)
  COMPLEX(DP),ALLOCATABLE :: r(:,:)
  !
  ALLOCATE( beta(NRHS),alpha(NRHS),raux(NRHS),ZA(NRHS),r(npwx*npol,NRHS) )
  !
  ! INIT
  r=b
  q_s=0._DP
  !
  CALL start_clock( "lan_H" )
  !
  ! FROM R TO Q & BETA
  !
  IF(gamma_only) THEN
     DO ip = 1, NRHS
        beta(ip) = 2._DP * DDOT(2*npw, r(1,ip), 1, r(1,ip), 1)
        IF(gstart==2) beta(ip) = beta(ip) - REAL(r(1,ip),KIND=DP) * REAL(r(1,ip),KIND=DP)
     ENDDO
  ELSE
     DO ip = 1, NRHS
        beta(ip) = DDOT(2*npw, r(1,ip), 1, r(1,ip), 1)
     ENDDO
     IF(noncolin) THEN
        DO ip = 1, NRHS
           beta(ip) = beta(ip) + DDOT(2*npw, r(1+npwx,ip), 1, r(1+npwx,ip), 1)
        ENDDO
     ENDIF
  ENDIF
  CALL mp_sum(beta,intra_bgrp_comm)
  beta(:) = SQRT( beta(:) )
  bnorm = beta
  !
  !
  ! --- Lanczos iterations
  !
  DO il=1,NLSTEPS
     !
!$OMP PARALLEL DO
     DO ip = 1, NRHS
        q_s(:,ip,il) = r(:,ip) / beta(ip)
     ENDDO
!$OMP END PARALLEL DO
     !
     ! NEW Q WAS FOUND !
     !
     ! APPLY H  q --> r
     !
     ! use h_psi_, i.e. h_psi without band parallelization, as wfreq
     ! handles band parallelization separately in solve_wfreq and solve_gfreq
     !
     CALL h_psi_( npwx, npw, NRHS, q_s(1,1,il), r )
     CALL apply_alpha_pc_to_m_wfcs(nbnd_to_deflate,NRHS,r,(1._DP,0._DP))
     !
     ! use beta
     !
     IF(il>1) THEN
        ZA(:)=CMPLX( -beta(:), 0._DP, KIND=DP )
        DO ip = 1, NRHS
           CALL ZAXPY(npw,ZA(ip),q_s(1,ip,il-1),1,r(1,ip),1)
        ENDDO
        IF(noncolin) THEN
           DO ip = 1, NRHS
              CALL ZAXPY(npw,ZA(ip),q_s(1+npwx,ip,il-1),1,r(1+npwx,ip),1)
           ENDDO
        ENDIF
     ENDIF
     !
     ! get alpha
     !
     IF(gamma_only) THEN
        DO ip = 1, NRHS
           alpha(ip) = 2._DP * DDOT(2*npw, q_s(1,ip,il), 1, r(1,ip), 1)
           IF(gstart==2) alpha(ip) = alpha(ip) - REAL(q_s(1,ip,il),KIND=DP) * REAL(r(1,ip),KIND=DP)
        ENDDO
     ELSE
        DO ip = 1, NRHS
           alpha(ip) = ZDOTC(npw, q_s(1,ip,il), 1, r(1,ip), 1)
        ENDDO
        IF(noncolin) THEN
           DO ip = 1, NRHS
              alpha(ip) = alpha(ip) + ZDOTC(npw, q_s(1+npwx,ip,il), 1, r(1+npwx,ip), 1)
           ENDDO
        ENDIF
     ENDIF
     !
     CALL mp_sum(alpha,intra_bgrp_comm)
     DO ip = 1,NRHS
        alpha_s(il,ip) = alpha(ip)
     ENDDO
     !
     ! use alpha
     !
     ZA(:)=CMPLX( -alpha(:), 0._DP, KIND=DP )
     DO ip = 1, NRHS
        CALL ZAXPY(npw,ZA(ip),q_s(1,ip,il),1,r(1,ip),1)
     ENDDO
     IF(noncolin) THEN
        DO ip = 1, NRHS
           CALL ZAXPY(npw,ZA(ip),q_s(1+npwx,ip,il),1,r(1+npwx,ip),1)
        ENDDO
     ENDIF
     !
     ! Enforce full ortho
     !
     DO ia=1,il-2
        IF(gamma_only) THEN
           DO ip=1,NRHS
              raux(ip) = 2._DP * DDOT(2*npw, q_s(1,ip,ia), 1, r(1,ip), 1)
              IF(gstart==2) raux(ip) = raux(ip) - REAL(q_s(1,ip,ia),KIND=DP) * REAL(r(1,ip),KIND=DP)
           ENDDO
           CALL mp_sum(raux,intra_bgrp_comm)
           ZA(:)=CMPLX( -raux(:), 0._DP, KIND=DP )
        ELSE
           DO ip=1,NRHS
              ZA(ip) = -ZDOTC(npw, q_s(1,ip,ia), 1, r(1,ip), 1)
           ENDDO
           IF(noncolin) THEN
              DO ip=1,NRHS
                 ZA(ip) = ZA(ip) -ZDOTC(npw, q_s(1+npwx,ip,ia), 1, r(1+npwx,ip), 1)
              ENDDO
           ENDIF
           CALL mp_sum(ZA,intra_bgrp_comm)
        ENDIF
        DO ip=1,NRHS
           CALL ZAXPY(npw,ZA(ip),q_s(1,ip,ia),1,r(1,ip),1)
        ENDDO
        IF(noncolin) THEN
           DO ip=1,NRHS
              CALL ZAXPY(npw,ZA(ip),q_s(1+npwx,ip,ia),1,r(1+npwx,ip),1)
           ENDDO
        ENDIF
     ENDDO
     !
     ! get beta
     !
     IF(gamma_only) THEN
        DO ip = 1, NRHS
           beta(ip) = 2._DP * DDOT(2*npw, r(1,ip), 1, r(1,ip), 1)
           IF(gstart==2) beta(ip) = beta(ip) - REAL(r(1,ip),KIND=DP) * REAL(r(1,ip),KIND=DP)
        ENDDO
     ELSE
        DO ip = 1, NRHS
           beta(ip) = DDOT(2*npw, r(1,ip), 1, r(1,ip), 1)
        ENDDO
        IF(noncolin) THEN
           DO ip = 1, NRHS
              beta(ip) = beta(ip) + DDOT(2*npw, r(1+npwx,ip), 1, r(1+npwx,ip), 1)
           ENDDO
        ENDIF
     ENDIF
     CALL mp_sum(beta,intra_bgrp_comm)
     beta(:) = SQRT( beta(:) )
     IF(il<NLSTEPS) THEN
        DO ip = 1, NRHS
           beta_s(il,ip) = beta(ip)
        ENDDO
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE( beta,alpha,raux,ZA,r )
  !
  CALL stop_clock( "lan_H" )
  !
END SUBROUTINE
!
#if defined(__CUDA)
!-----------------------------------------------------------------------
SUBROUTINE solve_deflated_lanczos_w_full_ortho_gpu(nbnd_to_deflate, NRHS, NLSTEPS, b, alpha_s, beta_s, q_s_d, bnorm)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : intra_bgrp_comm,nproc_bgrp
  USE mp,                   ONLY : mp_sum
  USE gvect,                ONLY : gstart
  USE pwcom,                ONLY : npw,npwx
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
  USE west_gpu,             ONLY : r,tmp_r,tmp_c
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nbnd_to_deflate
  INTEGER, INTENT(IN) :: NRHS
  INTEGER, INTENT(IN) :: NLSTEPS
  COMPLEX(DP), INTENT(IN) :: b(npwx*npol,NRHS)
  REAL(DP), INTENT(OUT) :: alpha_s(NLSTEPS,NRHS)
  REAL(DP), INTENT(OUT) :: beta_s(NLSTEPS-1,NRHS)
  COMPLEX(DP), DEVICE, INTENT(OUT) :: q_s_d(npwx*npol,NRHS,NLSTEPS)
  REAL(DP), INTENT(OUT) :: bnorm(NRHS)
  !
  ! Workspace
  !
  INTEGER :: ia,ig,ip,il
  REAL(DP) :: reduce_r
  COMPLEX(DP) :: reduce_c
  !
  ! INIT
  !
  q_s_d = (0.0_DP,0.0_DP)
  r = b
  !$acc update device(r)
  !
  CALL start_clock_gpu("lan_H")
  !
  ! FROM R TO Q & BETA
  !
  IF(gamma_only) THEN
     !$acc parallel vector_length(1024) present(r,tmp_r)
     !$acc loop
     DO ip = 1,NRHS
        reduce_r = 0.0_DP
        !$acc loop reduction(+:reduce_r)
        DO ig = 1,npw
           reduce_r = reduce_r+REAL(r(ig,ip),KIND=DP)**2+AIMAG(r(ig,ip))**2
        ENDDO
        IF(gstart == 2) THEN
           tmp_r(ip) = 2.0_DP*reduce_r-REAL(r(1,ip),KIND=DP)**2
        ELSE
           tmp_r(ip) = 2.0_DP*reduce_r
        ENDIF
     ENDDO
     !$acc end parallel
  ELSE
     !$acc parallel vector_length(1024) present(r,tmp_r)
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
        tmp_r(ip) = reduce_r
     ENDDO
     !$acc end parallel
  ENDIF
  !
  !$acc update host(tmp_r)
  !
  CALL mp_sum(tmp_r,intra_bgrp_comm)
  !
  DO ip = 1,NRHS
     tmp_r(ip) = SQRT(tmp_r(ip))
     bnorm(ip) = tmp_r(ip)
  ENDDO
  !
  ! --- Lanczos iterations
  !
  DO il = 1,NLSTEPS
     !
     !$acc update device(tmp_r)
     !
     !$acc parallel loop collapse(2) present(r,tmp_r)
     DO ip = 1,NRHS
        DO ig = 1,npwx*npol
           q_s_d(ig,ip,il) = r(ig,ip)/tmp_r(ip)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! NEW Q WAS FOUND !
     !
     ! APPLY H  q --> r
     !
     ! use h_psi__gpu, i.e. h_psi_gpu without band parallelization, as wfreq
     ! handles band parallelization separately in solve_wfreq and solve_gfreq
     !
     !$acc host_data use_device(r)
     CALL h_psi__gpu(npwx,npw,NRHS,q_s_d(1,1,il),r)
     CALL apply_alpha_pc_to_m_wfcs(nbnd_to_deflate,NRHS,r,(1.0_DP,0.0_DP))
     !$acc end host_data
     !
     ! use beta
     !
     IF(il > 1) THEN
        !$acc parallel loop collapse(2) present(r,tmp_r)
        DO ip = 1,NRHS
           DO ig = 1,npw
              r(ig,ip) = r(ig,ip)-q_s_d(ig,ip,il-1)*tmp_r(ip)
              IF(noncolin) THEN
                 r(ig+npwx,ip) = r(ig+npwx,ip)-q_s_d(ig+npwx,ip,il-1)*tmp_r(ip)
              ENDIF
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     ! get alpha
     !
     IF(gamma_only) THEN
        !$acc parallel vector_length(1024) present(r,tmp_r)
        !$acc loop
        DO ip = 1,NRHS
           reduce_r = 0.0_DP
           !$acc loop reduction(+:reduce_r)
           DO ig = 1,npw
              reduce_r = reduce_r+REAL(q_s_d(ig,ip,il),KIND=DP)*REAL(r(ig,ip),KIND=DP) &
              & +AIMAG(q_s_d(ig,ip,il))*AIMAG(r(ig,ip))
           ENDDO
           IF(gstart == 2) THEN
              tmp_r(ip) = 2.0_DP*reduce_r-REAL(q_s_d(1,ip,il),KIND=DP)*REAL(r(1,ip),KIND=DP)
           ELSE
              tmp_r(ip) = 2.0_DP*reduce_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel vector_length(1024) present(r,tmp_r)
        !$acc loop
        DO ip = 1,NRHS
           reduce_r = 0.0_DP
           !$acc loop reduction(+:reduce_r)
           DO ig = 1,npw
              reduce_r = reduce_r+REAL(q_s_d(ig,ip,il),KIND=DP)*REAL(r(ig,ip),KIND=DP) &
              & +AIMAG(q_s_d(ig,ip,il))*AIMAG(r(ig,ip))
              IF(noncolin) THEN
                 reduce_r = reduce_r+REAL(q_s_d(ig+npwx,ip,il),KIND=DP)*REAL(r(ig+npwx,ip),KIND=DP) &
                 & +AIMAG(q_s_d(ig+npwx,ip,il))*AIMAG(r(ig+npwx,ip))
              ENDIF
           ENDDO
           tmp_r(ip) = reduce_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc update host(tmp_r)
     !
     CALL mp_sum(tmp_r,intra_bgrp_comm)
     !
     DO ip = 1,NRHS
        alpha_s(il,ip) = tmp_r(ip)
     ENDDO
     !
     ! use alpha
     !
     !$acc update device(tmp_r)
     !
     !$acc parallel loop collapse(2) present(r,tmp_r)
     DO ip = 1,NRHS
        DO ig = 1,npw
           r(ig,ip) = r(ig,ip)-q_s_d(ig,ip,il)*tmp_r(ip)
           IF(noncolin) THEN
              r(ig+npwx,ip) = r(ig+npwx,ip)-q_s_d(ig+npwx,ip,il)*tmp_r(ip)
           ENDIF
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! Enforce full ortho
     !
     DO ia = 1,il-2
        IF(gamma_only) THEN
           !$acc parallel vector_length(1024) present(r,tmp_r)
           !$acc loop
           DO ip = 1,NRHS
              reduce_r = 0.0_DP
              !$acc loop reduction(+:reduce_r)
              DO ig = 1,npw
                 reduce_r = reduce_r+REAL(q_s_d(ig,ip,ia),KIND=DP)*REAL(r(ig,ip),KIND=DP) &
                 & +AIMAG(q_s_d(ig,ip,ia))*AIMAG(r(ig,ip))
              ENDDO
              IF(gstart == 2) THEN
                 tmp_r(ip) = 2.0_DP*reduce_r-REAL(q_s_d(1,ip,ia),KIND=DP)*REAL(r(1,ip),KIND=DP)
              ELSE
                 tmp_r(ip) = 2.0_DP*reduce_r
              ENDIF
           ENDDO
           !$acc end parallel
           !
           IF(nproc_bgrp > 1) THEN
              !$acc host_data use_device(tmp_r)
              CALL mp_sum(tmp_r,intra_bgrp_comm)
              !$acc end host_data
           ENDIF
           !
           !$acc parallel loop collapse(2) present(r,tmp_r)
           DO ip = 1,NRHS
              DO ig = 1,npw
                 r(ig,ip) = r(ig,ip)-q_s_d(ig,ip,ia)*tmp_r(ip)
                 IF(noncolin) THEN
                    r(ig+npwx,ip) = r(ig+npwx,ip)-q_s_d(ig+npwx,ip,ia)*tmp_r(ip)
                 ENDIF
              ENDDO
           ENDDO
           !$acc end parallel
        ELSE
           !$acc parallel vector_length(1024) present(r,tmp_c)
           !$acc loop
           DO ip = 1,NRHS
              reduce_c = 0.0_DP
              !$acc loop reduction(+:reduce_c)
              DO ig = 1,npw
                 reduce_c = reduce_c+CONJG(q_s_d(ig,ip,ia))*r(ig,ip)
                 IF(noncolin) THEN
                    reduce_c = reduce_c+CONJG(q_s_d(ig+npwx,ip,ia))*r(ig+npwx,ip)
                 ENDIF
              ENDDO
              tmp_c(ip) = reduce_c
           ENDDO
           !$acc end parallel
           !
           IF(nproc_bgrp > 1) THEN
              !$acc host_data use_device(tmp_c)
              CALL mp_sum(tmp_c,intra_bgrp_comm)
              !$acc end host_data
           ENDIF
           !
           !$acc parallel loop collapse(2) present(r,tmp_c)
           DO ip = 1,NRHS
              DO ig = 1,npw
                 r(ig,ip) = r(ig,ip)-q_s_d(ig,ip,ia)*tmp_c(ip)
                 IF(noncolin) THEN
                    r(ig+npwx,ip) = r(ig+npwx,ip)-q_s_d(ig+npwx,ip,ia)*tmp_c(ip)
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
        !$acc parallel vector_length(1024) present(r,tmp_r)
        !$acc loop
        DO ip = 1,NRHS
           reduce_r = 0.0_DP
           !$acc loop reduction(+:reduce_r)
           DO ig = 1,npw
              reduce_r = reduce_r+REAL(r(ig,ip),KIND=DP)**2+AIMAG(r(ig,ip))**2
           ENDDO
           IF(gstart == 2) THEN
              tmp_r(ip) = 2.0_DP*reduce_r-REAL(r(1,ip),KIND=DP)**2
           ELSE
              tmp_r(ip) = 2.0_DP*reduce_r
           ENDIF
        ENDDO
        !$acc end parallel
     ELSE
        !$acc parallel vector_length(1024) present(r,tmp_r)
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
           tmp_r(ip) = reduce_r
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc update host(tmp_r)
     !
     CALL mp_sum(tmp_r,intra_bgrp_comm)
     !
     DO ip = 1,NRHS
        tmp_r(ip) = SQRT(tmp_r(ip))
     ENDDO
     !
     IF(il < NLSTEPS) THEN
        DO ip = 1,NRHS
           beta_s(il,ip) = tmp_r(ip)
        ENDDO
     ENDIF
     !
  ENDDO
  !
  CALL stop_clock_gpu("lan_H")
  !
END SUBROUTINE
#endif
