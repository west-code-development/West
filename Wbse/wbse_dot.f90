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
SUBROUTINE wbse_dot (x,y,npwx,nbnd,nks,wbse_dot_out)
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE mp,                   ONLY : mp_sum
  USE pwcom,                ONLY : isk,wg
  USE westcom,              ONLY : nbnd_occ
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart
  USE klist,                ONLY : ngk
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (IN)     :: npwx,nbnd,nks
  COMPLEX(DP), INTENT (IN) :: x(npwx,nbnd,nks)
  COMPLEX(DP), INTENT (IN) :: y(npwx,nbnd,nks)
  COMPLEX(DP), INTENT (OUT):: wbse_dot_out(nspin)
  !
  INTEGER     :: ibnd, iks, is, current_spin, nbndval
  COMPLEX(DP) :: temp
  REAL(DP), EXTERNAL    :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  wbse_dot_out(:) = (0._DP,0._DP)
  !
  DO is = 1, nspin
     !
     temp = (0._DP, 0._DP)
     !
     DO iks = 1, nks
        !
        nbndval = nbnd_occ(iks)
        !
        current_spin = isk(iks)
        !
        IF (current_spin /= is) CYCLE
        !
        IF (gamma_only) THEN
           !
           DO ibnd=1, nbndval
              !
              temp = temp + 2._DP*wg(ibnd,iks)*DDOT(2*ngk(iks),x(:,ibnd,iks),1,y(:,ibnd,iks),1)
              !
              IF (gstart==2) temp = temp - wg(ibnd,iks)*REAL(x(1,ibnd,iks),KIND=DP)*REAL(y(1,ibnd,iks),KIND=DP)
              !
           ENDDO
           !
        ELSE
           !
           DO ibnd=1, nbndval
              !
              temp = temp + wg(ibnd,iks) * ZDOTC(ngk(iks),x(:,ibnd,iks),1,y(:,ibnd,iks),1)
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
     CALL mp_sum(temp, inter_pool_comm)
     CALL mp_sum(temp, intra_bgrp_comm)
     !
     wbse_dot_out(is) = temp*nspin/2._DP
     !
  ENDDO
  !
END SUBROUTINE
