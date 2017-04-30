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
SUBROUTINE solve_deflated_lanczos_w_full_ortho ( nbnd_to_deflate, NRHS, NLSTEPS, b, alpha_s, beta_s, q_s, bnorm)  
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
  INTEGER :: i1,i2,i3,ip,ig,il,is
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
  beta(:) = DSQRT( beta(:) )
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
     CALL h_psi( npwx, npw, NRHS, q_s(1,1,il), r )
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
     DO i1=1,il-2
        IF(gamma_only) THEN
           DO ip=1,NRHS
              raux(ip) = 2._DP * DDOT(2*npw, q_s(1,ip,i1), 1, r(1,ip), 1)
              IF(gstart==2) raux(ip) = raux(ip) - REAL(q_s(1,ip,i1),KIND=DP) * REAL(r(1,ip),KIND=DP)
           ENDDO
           CALL mp_sum(raux,intra_bgrp_comm)
           ZA(:)=CMPLX( -raux(:), 0._DP, KIND=DP )
        ELSE
           DO ip=1,NRHS
              ZA(ip) = -ZDOTC(npw, q_s(1,ip,i1), 1, r(1,ip), 1)
           ENDDO
           IF(noncolin) THEN
              DO ip=1,NRHS
                 ZA(ip) = ZA(ip) -ZDOTC(npw, q_s(1+npwx,ip,i1), 1, r(1+npwx,ip), 1)
              ENDDO
           ENDIF 
           CALL mp_sum(ZA,intra_bgrp_comm)
        ENDIF
        DO ip=1,NRHS
           CALL ZAXPY(npw,ZA(ip),q_s(1,ip,i1),1,r(1,ip),1) 
        ENDDO
        IF(noncolin) THEN
           DO ip=1,NRHS
              CALL ZAXPY(npw,ZA(ip),q_s(1+npwx,ip,i1),1,r(1+npwx,ip),1) 
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
     beta(:) = DSQRT( beta(:) )
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
