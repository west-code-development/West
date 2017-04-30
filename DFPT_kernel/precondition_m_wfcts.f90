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
SUBROUTINE precondition_m_wfcts (m,f,pf,eprec)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE wvfct,                 ONLY : g2kin
  USE noncollin_module,      ONLY : noncolin,npol
  USE pwcom,                 ONLY : npw,npwx
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m
  COMPLEX(DP),INTENT(IN) :: f(npwx*npol,m)
  COMPLEX(DP),INTENT(OUT) :: pf(npwx*npol,m)
  REAL(DP),INTENT(IN) :: eprec(m)
  !
  ! Workspace
  !
  INTEGER :: ibnd, ig
  !
  pf = 0._DP
  !
  DO ibnd=1,m
     !
!$OMP PARALLEL DO
     DO ig=1,npw 
        pf(ig,ibnd) = f(ig,ibnd) / MAX(1._DP,g2kin(ig)/eprec(ibnd))
     ENDDO
!$OMP END PARALLEL DO
     IF( noncolin ) THEN 
!$OMP PARALLEL DO
        DO ig=1,npw
           pf(npwx+ig,ibnd) = f(npwx+ig,ibnd) / MAX(1._DP,g2kin(ig)/eprec(ibnd))
        ENDDO
!$OMP END PARALLEL DO
     ENDIF
     !
  ENDDO
  !
END SUBROUTINE
