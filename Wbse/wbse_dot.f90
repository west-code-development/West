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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_dot(x,y,m,dotp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE pwcom,                ONLY : wg,nspin,npw,npwx,ngk
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands
  USE distribution_center,  ONLY : kpt_pool,band_group
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m
  COMPLEX(DP), INTENT(IN) :: x(npwx,m,kpt_pool%nloc)
  COMPLEX(DP), INTENT(IN) :: y(npwx,m,kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: dotp(nspin)
  !
  ! Workspace
  !
  INTEGER :: ig, lbnd, ibnd, iks, iks_g, nbndval, band_group_myoffset
  REAL(DP) :: tmp_r
  !
  band_group_myoffset = band_group%myoffset
  !
  dotp(:) = (0._DP,0._DP)
  !
  DO iks = 1, kpt_pool%nloc
     !
     iks_g = kpt_pool%l2g(iks)
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     tmp_r = 0._DP
     !
     !$acc parallel loop collapse(2) reduction(+:tmp_r) present(wg,x,y) copy(tmp_r)
     DO lbnd = 1, m
        DO ig = 1, npw
           !
           ibnd = band_group_myoffset+lbnd+n_trunc_bands
           !
           tmp_r = tmp_r + wg(ibnd,iks)*2._DP*(REAL(x(ig,lbnd,iks),KIND=DP)*REAL(y(ig,lbnd,iks),KIND=DP) &
           & + AIMAG(x(ig,lbnd,iks))*AIMAG(y(ig,lbnd,iks)))
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gstart == 2) THEN
        !$acc parallel loop reduction(+:tmp_r) present(wg,x,y) copy(tmp_r)
        DO lbnd = 1, m
           !
           ibnd = band_group_myoffset+lbnd+n_trunc_bands
           !
           tmp_r = tmp_r - wg(ibnd,iks)*REAL(x(1,lbnd,iks),KIND=DP)*REAL(y(1,lbnd,iks),KIND=DP)
           !
        ENDDO
        !$acc end parallel
     ENDIF
     !
     dotp(iks_g) = CMPLX(tmp_r*nspin/2._DP,KIND=DP)
     !
  ENDDO
  !
  CALL mp_sum(dotp,intra_bgrp_comm)
  CALL mp_sum(dotp,inter_bgrp_comm)
  CALL mp_sum(dotp,inter_pool_comm)
  !
END SUBROUTINE
