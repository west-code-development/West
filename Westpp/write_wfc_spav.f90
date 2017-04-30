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
! -------------------------------------------------------------------
SUBROUTINE write_wfc_spav ( iu, fname, wfc, r0, nr, rmax )
  ! -----------------------------------------------------------------
  !
  ! r0 and rmax are in a.u 
  !
  USE pwcom,                 ONLY : npw,npwx,tpiba,alat
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : celldm, at, bg, omega
  USE ions_base,             ONLY : nat, tau, atm, ityp
  USE mp_global,             ONLY : me_bgrp,root_bgrp,intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE constants,             ONLY : fpi, tpi
  USE gvect,                 ONLY : g, gstart, ngm 
  USE control_flags,         ONLY : gamma_only 
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: iu
  CHARACTER(LEN=*),INTENT(IN) :: fname
  COMPLEX(DP),INTENT(IN) :: wfc(ngm)
  REAL(DP),INTENT(IN) :: r0(3)
  INTEGER,INTENT(IN) :: nr
  REAL(DP),INTENT(IN) :: rmax
  ! 
  ! Workspace
  !
  REAL(DP) :: rmod(nr+1)
  REAL(DP) :: spav(nr+1)
  INTEGER :: ig, ir
  REAL(DP) :: gmod, arg, gr, sinx_by_x
  COMPLEX(DP) :: f_times_exp_igdotr0
  ! 
  DO ir = 1, nr+1
     rmod(ir) = DBLE(ir-1) * rmax / DBLE(nr) 
  ENDDO  
  !
  spav = 0._DP
  !
  IF ( gamma_only ) THEN 
     !
     IF(gstart == 2) THEN 
        spav(:) = spav(:) + REAL( wfc(ig) , KIND=DP )
     ENDIF 
     !
     DO ig = gstart, ngm
        !
        arg = ( r0(1)*g(1,ig) + r0(2)*g(2,ig) + r0(3)*g(3,ig) ) * tpiba
        f_times_exp_igdotr0 = wfc(ig)*CMPLX( COS(arg), SIN(arg), KIND=DP) 
        gmod = SQRT(g(1,ig)**2 + g(2,ig)**2 + g(3,ig)**2) * tpiba
        !
        spav(1) = spav(1) + 2._DP * REAL( f_times_exp_igdotr0, KIND=DP )
        !
        DO ir =2, nr+1 
           !
           gr = gmod * rmod(ir) 
           spav( ir ) = spav(ir) + 2._DP * REAL( f_times_exp_igdotr0, KIND=DP ) * SIN( gr ) / gr 
           !
        ENDDO
        !
     ENDDO
     !
  ELSE
     !
     CALL errore('spav','spav for non gamma_only not yet implemented','1')
     !
  ENDIF
  !
  CALL mp_sum( spav, intra_bgrp_comm )
  !
  IF( me_bgrp == root_bgrp ) THEN
     OPEN(UNIT=iu,FILE=TRIM(ADJUSTL(fname)))
     WRITE(iu,'(a,3f16.6)') "# r        f(r) : R0 (a.u)", r0(:)
     DO ir = 1, nr+1
        WRITE(iu,'(2es14.6)') rmod(ir), spav(ir)
     ENDDO
     CLOSE(iu)
  ENDIF
  !
END SUBROUTINE
