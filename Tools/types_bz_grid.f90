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
! Marco Govoni, Matteo Gerosa
!
!-----------------------------------------------------------------------
MODULE types_bz_grid
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE class_bz_grid,   ONLY : bz_grid
  !
  IMPLICIT NONE
  !
  TYPE(bz_grid) :: k_grid
  TYPE(bz_grid) :: q_grid
  !
  CONTAINS
  !
  !
  FUNCTION findG(g0, unit_type) RESULT(ig0)
     !
     ! ... ig0 is the index of G (unit_type = [ "cryst", "cart"])
     ! ... if on exit ig0 == 0 --> G is not found
     !
     USE cell_base,        ONLY : bg
     USE gvect,            ONLY : g, ngm
     USE constants,        ONLY : eps8
     !
     IMPLICIT NONE
     !
     ! I/O
     !
     REAL(DP), INTENT(IN) :: g0(3)
     CHARACTER(LEN=*), INTENT(IN) :: unit_type
     !
     ! Workspace
     !
     REAL(DP) :: gtemp(3)
     INTEGER :: ig, ig0
     !
     SELECT CASE(unit_type)
     CASE("cryst","cart")
     CASE DEFAULT
        CALL errore( "types_bz_grid", "unit_type not supported, supported only cryst or cart", 1 )
     END SELECT
     !
     gtemp = g0
     IF( unit_type == "cryst" ) CALL cryst_to_cart( 1, gtemp, bg, 1)
     !
     ! gtemp is in cart
     !
     ig0 = 0
     DO ig = 1, ngm
        IF ( ALL ( ABS( g(:,ig) - gtemp(:) ) < eps8 ) ) THEN
           ig0 = ig
           EXIT
        ENDIF
     ENDDO
     !
  END FUNCTION
  !
  !
  SUBROUTINE compute_phase( g0, unit_type, phase )
     !
     ! ... phase(r) = exp(-iG_0*r)  (allocated externally)
     !
     USE fft_base,         ONLY : dffts
     USE fft_interfaces,   ONLY : invfft
     USE mp,               ONLY : mp_max
     USE mp_bands,         ONLY : intra_bgrp_comm
     !
     IMPLICIT NONE
     !
     ! I/O
     !
     REAL(DP), INTENT(IN) :: g0(3)
     CHARACTER(LEN=*), INTENT(IN) :: unit_type
     COMPLEX(DP), INTENT(OUT) :: phase(:)
     !
     ! Workspace
     !
     INTEGER :: ig0, glob_ig0
     !
     SELECT CASE(unit_type)
     CASE("cryst","cart")
     CASE DEFAULT
        CALL errore( "types_bz_grid", "unit_type not supported, supported only cryst or cart", 1 )
     END SELECT
     !
     ig0 = findG(g0,unit_type)
     glob_ig0 = ig0
     CALL mp_max( glob_ig0, intra_bgrp_comm )
     !
     IF( glob_ig0 == 0 ) CALL errore( "types_bz_grid", "G0 not found", 1 )
     !
     ! phase(r) = exp(-iG_0*r)
     !
     phase = (0._DP, 0._DP)
     !
     IF ( ig0 /= 0 ) THEN
        phase( dffts%nl(ig0) ) = (1._DP, 0._DP)
     ENDIF
     CALL invfft( 'Wave', phase, dffts )
     phase(1:dffts%nnr) = DCONJG( phase(1:dffts%nnr) )
     !
  END SUBROUTINE
  !
END MODULE
