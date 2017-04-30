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
SUBROUTINE set_eprec(m,wfc,eprec)
  !-----------------------------------------------------------------------
  !
  ! Set eprec, for precondiconditioning
  !
  USE kinds,                 ONLY : DP
  USE wvfct,                 ONLY : g2kin
  USE noncollin_module,      ONLY : noncolin,npol
  USE pwcom,                 ONLY : npw,npwx
  USE mp,                    ONLY : mp_sum
  USE mp_global,             ONLY : intra_bgrp_comm
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m
  COMPLEX(DP),INTENT(IN) :: wfc(npwx*npol,m)
  REAL(DP),INTENT(OUT) :: eprec(m)
  !
  ! Workspace
  !
  INTEGER :: ibnd, ig 
  !
  eprec = 0._DP
  !
  DO ibnd = 1, m
     !
     DO ig = 1, npw
        eprec(ibnd) = eprec(ibnd) + CONJG( wfc(ig,ibnd) ) * wfc(ig,ibnd) * g2kin(ig)
     ENDDO
     !
     IF (noncolin) THEN 
        DO ig = 1, npw
           eprec(ibnd) = eprec(ibnd) + CONJG( wfc(ig+npwx,ibnd) ) * wfc(ig+npwx,ibnd) * g2kin(ig)
        ENDDO
     ENDIF 
     !
  ENDDO
  !
  CALL mp_sum(eprec,intra_bgrp_comm)
  !
  eprec = eprec * 1.35_DP
  !
END SUBROUTINE
