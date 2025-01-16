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
MODULE check_ovl_wfc
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    SUBROUTINE read_bisection_loc(current_spin, numband, bisec_loc)
      !
      USE io_global,     ONLY : ionode, ionode_id
      USE mp,            ONLY : mp_bcast
      USE mp_world,      ONLY : world_comm
      USE westcom,       ONLY : bisection_info
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: current_spin, numband
      INTEGER, INTENT(INOUT) :: bisec_loc(numband)
      !
      INTEGER :: ibnd, num_localized_orb,iunit
      CHARACTER :: my_spin
      CHARACTER(LEN=256) :: file_bisection
      !
      ! read eigenvalues from file
      !
      bisec_loc(:) = 0
      !
      WRITE(my_spin,'(i1)') current_spin
      !
      file_bisection = TRIM(bisection_info)//'.'//my_spin
      !
      IF(ionode) THEN
         OPEN(NEWUNIT=iunit, FILE=TRIM(file_bisection), FORM='FORMATTED', STATUS='OLD')
         !
         READ(iunit,*) num_localized_orb
         !
         DO ibnd = 1, num_localized_orb
            READ(iunit,*) bisec_loc(ibnd)
         ENDDO
         !
         CLOSE(iunit)
      ENDIF
      !
      CALL mp_bcast(num_localized_orb,ionode_id,world_comm)
      CALL mp_bcast(bisec_loc,ionode_id,world_comm)
      !
    END SUBROUTINE
    !
    SUBROUTINE check_ovl_bisection(orbital_i, orbital_j, ovl_value)
      !
      ! input is the bitwise number of orb_i and orb_j
      !
      USE kinds,                ONLY : DP
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: orbital_i
      INTEGER, INTENT(IN) :: orbital_j
      REAL(DP), INTENT(OUT):: ovl_value
      !
      INTEGER :: loc_i, loc_j
      LOGICAL :: l_p_left_i, l_p_left_j, l_p_right_i, l_p_right_j
      INTEGER :: int_p_left_i, int_p_left_j, int_p_right_i, int_p_right_j
      !
      loc_i = orbital_i
      loc_j = orbital_j
      !
      DO WHILE(loc_i /= 0 .AND. loc_j /= 0)
         !
         ! get the weight of projections for each state
         !
         int_p_right_i = AND(loc_i, 1)
         int_p_left_i  = AND(loc_i, 2)
         int_p_right_j = AND(loc_j, 1)
         int_p_left_j  = AND(loc_j, 2)
         !
         IF(int_p_right_i == 0) THEN
            l_p_right_i = .FALSE.
         ELSE
            l_p_right_i = .TRUE.
         ENDIF
         !
         IF(int_p_left_i == 0) THEN
            l_p_left_i = .FALSE.
         ELSE
            l_p_left_i = .TRUE.
         ENDIF
         !
         IF(int_p_right_j == 0) THEN
            l_p_right_j = .FALSE.
         ELSE
            l_p_right_j = .TRUE.
         ENDIF
         !
         IF(int_p_left_j == 0) THEN
            l_p_left_j = .FALSE.
         ELSE
            l_p_left_j = .TRUE.
         ENDIF
         !
         ! return false as soon as the states are found to be separated
         !
         IF(.NOT. ((l_p_right_i .AND. l_p_right_j) .OR. (l_p_left_i .AND. l_p_left_j))) THEN
            !
            ovl_value = -1._DP
            !
            RETURN
            !
         ENDIF
         !
         loc_i = ISHFT(loc_i, -2)
         loc_j = ISHFT(loc_j, -2)
         !
      ENDDO
      !
      ovl_value = 1._DP
      !
    END SUBROUTINE
    !
    SUBROUTINE check_ovl_wannier(orb_i, orb_j, ovl_value)
      !
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dffts
      USE cell_base,            ONLY : omega
      USE mp_global,            ONLY : intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: orb_i(dffts%nnr)
      REAL(DP), INTENT(IN) :: orb_j(dffts%nnr)
      REAL(DP), INTENT(OUT) :: ovl_value
      !
      INTEGER :: dffts_nnr, ir
      REAL(DP) :: summ_ib, summ_jb, summ_ij
      !
      dffts_nnr = dffts%nnr
      !
      summ_ib = 0._DP
      summ_jb = 0._DP
      summ_ij = 0._DP
      !
      !$acc parallel loop reduction(+:summ_ib,summ_jb,summ_ij) present(orb_i,orb_j) copy(summ_ib,summ_jb,summ_ij)
      DO ir = 1, dffts_nnr
         summ_ib = summ_ib + orb_i(ir)**4
         summ_jb = summ_jb + orb_j(ir)**4
         summ_ij = summ_ij + orb_i(ir)**2 * orb_j(ir)**2
      ENDDO
      !$acc end parallel
      !
      summ_ib = summ_ib * omega / (dffts%nr1*dffts%nr2*dffts%nr3)
      summ_jb = summ_jb * omega / (dffts%nr1*dffts%nr2*dffts%nr3)
      summ_ij = summ_ij * omega / (dffts%nr1*dffts%nr2*dffts%nr3)
      !
      CALL mp_sum(summ_ib, intra_bgrp_comm)
      CALL mp_sum(summ_jb, intra_bgrp_comm)
      CALL mp_sum(summ_ij, intra_bgrp_comm)
      !
      ovl_value = summ_ij / SQRT(summ_ib*summ_jb)
      !
    END SUBROUTINE
    !
END MODULE
