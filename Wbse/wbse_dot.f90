!
! Copyright (C) 2015-2023 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_dot(x,y,m,nks,dotp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE mp_global,            ONLY : inter_bgrp_comm,nbgrp,my_bgrp_id,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE pwcom,                ONLY : isk,wg,nspin,npw,npwx,ngk
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m,nks
  COMPLEX(DP), INTENT(IN) :: x(npwx,m,nks)
  COMPLEX(DP), INTENT(IN) :: y(npwx,m,nks)
  COMPLEX(DP), INTENT(OUT) :: dotp(nspin)
  !
  ! Workspace
  !
  INTEGER :: ig, lbnd, ibnd, iks, is, current_spin, nbndval
  REAL(DP) :: tmp_r
  COMPLEX(DP) :: tmp_c
  !
  DO is = 1, nspin
     !
     tmp_r = 0._DP
     tmp_c = (0._DP,0._DP)
     !
     DO iks = 1, nks
        !
        npw = ngk(iks)
        nbndval = nbnd_occ(iks)
        current_spin = isk(iks)
        IF(current_spin /= is) CYCLE
        !
        IF(gamma_only) THEN
           !
           !$acc parallel loop collapse(2) reduction(+:tmp_r) present(wg,x,y) copy(tmp_r)
           DO lbnd = 1, m
              DO ig = 1, npw
                 ibnd = nbgrp*(lbnd-1) + my_bgrp_id + 1 + n_trunc_bands
                 tmp_r = tmp_r + wg(ibnd,iks)*(REAL(x(ig,lbnd,iks),KIND=DP)*REAL(y(ig,lbnd,iks),KIND=DP) &
                 & + AIMAG(x(ig,lbnd,iks))*AIMAG(y(ig,lbnd,iks)))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           tmp_r = tmp_r*2._DP
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:tmp_r) present(wg,x,y) copy(tmp_r)
              DO lbnd = 1, m
                 ibnd = nbgrp*(lbnd-1) + my_bgrp_id + 1 + n_trunc_bands
                 tmp_r = tmp_r - wg(ibnd,iks)*REAL(x(1,lbnd,iks),KIND=DP)*REAL(y(1,lbnd,iks),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
        ELSE
           !
           !$acc parallel loop collapse(2) reduction(+:tmp_c) present(wg,x,y) copy(tmp_c)
           DO lbnd = 1, m
              DO ig = 1, npw
                 ibnd = nbgrp*(lbnd-1) + my_bgrp_id + 1 + n_trunc_bands
                 tmp_c = tmp_c + wg(ibnd,iks)*CONJG(x(ig,lbnd,iks))*y(ig,lbnd,iks)
              ENDDO
           ENDDO
           !$acc end parallel
           !
        ENDIF
        !
     ENDDO
     !
     IF(gamma_only) THEN
        CALL mp_sum(tmp_r,intra_bgrp_comm)
        CALL mp_sum(tmp_r,inter_bgrp_comm)
        dotp(is) = CMPLX(tmp_r*nspin/2._DP,KIND=DP)
     ELSE
        CALL mp_sum(tmp_c,intra_bgrp_comm)
        CALL mp_sum(tmp_c,inter_bgrp_comm)
        dotp(is) = tmp_c*nspin/2._DP
     ENDIF
     !
  ENDDO
  !
END SUBROUTINE
