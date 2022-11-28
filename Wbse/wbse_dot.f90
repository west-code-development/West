!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE wbse_dot(x,y,nbnd,nks,wbse_dot_out)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE pwcom,                ONLY : isk,wg,nspin,npw,npwx,ngk
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nbnd,nks
  COMPLEX(DP), INTENT(IN) :: x(npwx,nbnd,nks)
  COMPLEX(DP), INTENT(IN) :: y(npwx,nbnd,nks)
  COMPLEX(DP), INTENT(OUT) :: wbse_dot_out(nspin)
  !
  ! Workspace
  !
  INTEGER :: ig, ibnd, iks, is, current_spin, nbndval
  COMPLEX(DP) :: reduce
  !
  DO is = 1, nspin
     !
     reduce = (0._DP,0._DP)
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
           !$acc parallel loop collapse(2) reduction(+:reduce) present(wg,x,y) copy(reduce)
           DO ibnd = 1, nbndval
              DO ig = 1, npw
                 reduce = reduce + wg(ibnd,iks)*(REAL(x(ig,ibnd,iks),KIND=DP)*REAL(y(ig,ibnd,iks),KIND=DP) &
                 & + AIMAG(x(ig,ibnd,iks))*AIMAG(y(ig,ibnd,iks)))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:reduce) present(wg,x,y) copy(reduce)
              DO ibnd = 1, nbndval
                 reduce = reduce - 0.5_DP*wg(ibnd,iks)*REAL(x(1,ibnd,iks),KIND=DP)*REAL(y(1,ibnd,iks),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
           reduce = reduce*2._DP
           !
        ELSE
           !
           !$acc parallel loop collapse(2) reduction(+:reduce) present(wg,x,y) copy(reduce)
           DO ibnd = 1, nbndval
              DO ig = 1, npw
                 reduce = reduce + wg(ibnd,iks)*CONJG(x(ig,ibnd,iks))*y(ig,ibnd,iks)
              ENDDO
           ENDDO
           !$acc end parallel
           !
        ENDIF
        !
     ENDDO
     !
     CALL mp_sum(reduce,intra_bgrp_comm)
     !
     wbse_dot_out(is) = reduce*nspin/2._DP
     !
  ENDDO
  !
END SUBROUTINE
