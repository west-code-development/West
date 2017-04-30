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
SUBROUTINE set_isz( isz_mode, isz )
  !-----------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE constants,              ONLY : pi,fpi
  USE pwcom,                  ONLY : tpiba2,npw
  USE gvect,                  ONLY : g,gstart
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : intra_bgrp_comm
  USE io_global,              ONLY : stdout
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: isz_mode
  REAL(DP),INTENT(OUT) :: isz
  !
  ! Workspace
  ! 
  REAL(DP) :: gammafact, g2
  INTEGER :: ig, partial
  !
  WRITE(stdout,'(5x,"isz_mode = ",i6)') isz_mode
  !
  SELECT CASE ( isz_mode )
     !
     CASE(1) ! spherical region
        !
        isz = ( (6._DP * pi * pi / omega )**(1._DP/3._DP) ) / ( 2._DP * pi * pi )
        !
     CASE(2) ! gygi-baldereschi  
        !
        gammafact = 1._DP 
        !
        partial = 0._DP
        DO ig = gstart, npw
           g2 = ( g(1,ig)*g(1,ig) + g(2,ig)*g(2,ig) + g(3,ig)*g(3,ig) ) * tpiba2
           partial = partial + DEXP( - gammafact * g2 ) / g2
        ENDDO
        !
        CALL mp_sum( partial, intra_bgrp_comm )
        ! 
        isz = 1._DP / ( fpi * SQRT( pi * gammafact ) ) - 2._DP * partial / omega + gammafact / omega  
        !
     CASE DEFAULT
        !
        isz = 0._DP
        !
  END SELECT
  !
END SUBROUTINE
