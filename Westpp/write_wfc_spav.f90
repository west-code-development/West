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
! Marco Govoni
!
! -------------------------------------------------------------------
SUBROUTINE wfc_spav(wfc, r0, nr, rmax, spav)
  ! -----------------------------------------------------------------
  !
  ! r0 and rmax are in a.u
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : tpiba
  USE gvect,                 ONLY : g,gstart,ngm
  USE control_flags,         ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP),INTENT(IN) :: wfc(ngm)
  REAL(DP),INTENT(IN) :: r0(3)
  INTEGER,INTENT(IN) :: nr
  REAL(DP),INTENT(IN) :: rmax
  REAL(DP),INTENT(OUT) :: spav(nr+1)
  !
  ! Workspace
  !
  INTEGER :: ig, ir
  REAL(DP) :: gmod, arg, gr
  REAL(DP) :: rmod(nr+1)
  COMPLEX(DP) :: f_times_exp_igdotr0
  !
  DO ir = 1, nr+1
     rmod(ir) = rmax/nr*(ir-1)
  ENDDO
  !
  spav = 0._DP
  !
  IF(gamma_only) THEN
     !
     IF(gstart == 2) THEN
        spav(:) = spav(:) + REAL(wfc(1), KIND=DP)
     ENDIF
     !
     DO ig = gstart, ngm
        !
        arg = ( r0(1)*g(1,ig) + r0(2)*g(2,ig) + r0(3)*g(3,ig) ) * tpiba
        f_times_exp_igdotr0 = wfc(ig)*CMPLX(COS(arg), SIN(arg), KIND=DP)
        gmod = SQRT(g(1,ig)**2 + g(2,ig)**2 + g(3,ig)**2) * tpiba
        !
        spav(1) = spav(1) + 2._DP * REAL(f_times_exp_igdotr0, KIND=DP)
        !
        DO ir =2, nr+1
           !
           gr = gmod * rmod(ir)
           spav( ir ) = spav(ir) + 2._DP * REAL(f_times_exp_igdotr0, KIND=DP) * SIN(gr) / gr
           !
        ENDDO
        !
     ENDDO
     !
  ELSE
     !
     CALL errore('spav', 'spav for non gamma_only not yet implemented', '1')
     !
  ENDIF
  !
END SUBROUTINE
!
! -------------------------------------------------------------------
SUBROUTINE write_wfc_spav(fname, wfc, r0, nr, rmax)
  ! -----------------------------------------------------------------
  !
  ! r0 and rmax are in a.u
  !
  USE kinds,                 ONLY : DP
  USE mp_global,             ONLY : me_bgrp,root_bgrp,intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE gvect,                 ONLY : ngm
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  CHARACTER(LEN=*),INTENT(IN) :: fname
  COMPLEX(DP),INTENT(IN) :: wfc(ngm)
  REAL(DP),INTENT(IN) :: r0(3)
  INTEGER,INTENT(IN) :: nr
  REAL(DP),INTENT(IN) :: rmax
  !
  ! Workspace
  !
  INTEGER :: iu
  INTEGER :: ir
  REAL(DP) :: dr
  REAL(DP) :: spav(nr+1)
  !
  CALL wfc_spav(wfc, r0, nr, rmax, spav)
  !
  CALL mp_sum(spav, intra_bgrp_comm)
  !
  IF(me_bgrp == root_bgrp) THEN
     OPEN(NEWUNIT=iu, FILE=TRIM(ADJUSTL(fname)))
     WRITE(iu, '(a,3f16.6)') '# r        f(r) : R0 (a.u)', r0(:)
     dr = rmax/nr
     DO ir = 1, nr+1
        WRITE(iu, '(2es14.6)') dr*(ir-1), spav(ir)
     ENDDO
     CLOSE(iu)
  ENDIF
  !
END SUBROUTINE
