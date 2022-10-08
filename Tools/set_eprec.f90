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
SUBROUTINE set_eprec(m,wfc,eprec)
  !-----------------------------------------------------------------------
  !
  ! Set eprec, for precondiconditioning
  !
  USE kinds,                 ONLY : DP
#if defined(__CUDA)
  USE wvfct_gpum,            ONLY : g2kin=>g2kin_d
#else
  USE wvfct,                 ONLY : g2kin
#endif
  USE noncollin_module,      ONLY : noncolin,npol
  USE pwcom,                 ONLY : npw,npwx
  USE mp,                    ONLY : mp_sum
  USE mp_global,             ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m
  COMPLEX(DP),INTENT(IN) :: wfc(npwx*npol,m)
  REAL(DP),INTENT(OUT) :: eprec(m)
#if defined(__CUDA)
  ATTRIBUTES(DEVICE) :: wfc,eprec
#endif
  !
  ! Workspace
  !
  INTEGER :: ibnd,ig
  REAL(DP) :: reduce
  REAL(DP),PARAMETER :: factor = 1.35_DP
  !
  !$acc parallel vector_length(1024)
  !$acc loop
  DO ibnd = 1,m
     reduce = 0._DP
     !$acc loop reduction(+:reduce)
     DO ig = 1,npw
        reduce = reduce+CONJG(wfc(ig,ibnd))*wfc(ig,ibnd)*g2kin(ig)
        IF(noncolin) THEN
           reduce = reduce+CONJG(wfc(ig+npwx,ibnd))*wfc(ig+npwx,ibnd)*g2kin(ig)
        ENDIF
     ENDDO
     eprec(ibnd) = reduce*factor
  ENDDO
  !$acc end parallel
  !
  CALL mp_sum(eprec,intra_bgrp_comm)
  !
END SUBROUTINE
