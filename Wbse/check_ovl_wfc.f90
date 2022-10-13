!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE check_ovl_wfc
  !
  IMPLICIT NONE
  !
  INTERFACE check_ovl_wannier
     !
     MODULE PROCEDURE check_ovl_wannier_real, check_ovl_wannier_cmplx
     !
  ENDINTERFACE
  !
  PUBLIC :: check_ovl_bisection, read_bisection_loc
  !
  CONTAINS
    !
    SUBROUTINE read_bisection_loc(current_spin, numband, bisec_loc)
      !
      USE io_global,     ONLY : stdout, ionode, ionode_id
      USE mp,            ONLY : mp_bcast, mp_barrier
      USE mp_world,      ONLY : world_comm
      USE bse_module,    ONLY : et_qp
      USE wvfct,         ONLY : nbnd, et
      USE pwcom,         ONLY : nks
      USE lsda_mod,      ONLY : lsda
      USE westcom,       ONLY : bisection_info
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)    :: current_spin, numband
      INTEGER, INTENT(INOUT) :: bisec_loc(numband)
      !
      INTEGER :: ibnd, num_localized_orb,iunit
      CHARACTER(LEN=3) :: my_spin
      CHARACTER(LEN=256) :: file_bisection
      !
      ! read eigenvalues from file
      !
      CALL mp_barrier(world_comm)
      !
      bisec_loc(:) = 0
      !
      WRITE(my_spin,'(i1)') current_spin
      !
      file_bisection =  TRIM(bisection_info)//'.'//TRIM(my_spin)
      !file_bisection = 'bisection_localization_'//TRIM(my_spin)//'.dat'
      !
      if (ionode) then
         !
         open(NEWUNIT=iunit, file =TRIM(file_bisection),form = 'formatted',status = 'old')
         !
         read(iunit,*) num_localized_orb
         !
         do ibnd = 1, num_localized_orb
            !
            read (iunit, * ) bisec_loc(ibnd)
            !
         enddo
         !
         close (iunit)
         !
      endif
      !
      call mp_bcast (num_localized_orb,ionode_id,world_comm )
      call mp_bcast (bisec_loc,ionode_id,world_comm)
      !
      return
      !
    END SUBROUTINE
    !
    SUBROUTINE check_ovl_bisection (orbital_i, orbital_j, ovl_value)
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
      !
      LOGICAL :: l_p_left_i, l_p_left_j, l_p_right_i, l_p_right_j
      INTEGER :: int_p_left_i, int_p_left_j, int_p_right_i, int_p_right_j
      !
      loc_i = orbital_i
      loc_j = orbital_j
      !
      DO WHILE ( (loc_i /= 0) .AND. (loc_j /= 0) )
         !
         ! get the weight of projections for each state
         !
         int_p_right_i = AND(loc_i, 1)
         int_p_left_i  = AND(loc_i, 2)
         int_p_right_j = AND(loc_j, 1)
         int_p_left_j  = AND(loc_j, 2)
         !
         IF (int_p_right_i == 0) THEN
            l_p_right_i = .FALSE.
         ELSE
            l_p_right_i = .TRUE.
         ENDIF
         !
         IF (int_p_left_i == 0) THEN
            l_p_left_i = .FALSE.
         ELSE
            l_p_left_i = .TRUE.
         ENDIF
         !
         IF (int_p_right_j == 0) THEN
            l_p_right_j = .FALSE.
         ELSE
            l_p_right_j = .TRUE.
         ENDIF
         !
         IF (int_p_left_j == 0) THEN
            l_p_left_j = .FALSE.
         ELSE
            l_p_left_j = .TRUE.
         ENDIF
         !
         ! return false as soon as the states are found to be separated
         !
         IF ( .NOT.( ( l_p_right_i .AND. l_p_right_j ) .OR. ( l_p_left_i .AND. l_p_left_j ) ) ) THEN
            !
            ovl_value = -1.0_DP
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
      ovl_value = 1.0_DP
      !
      RETURN
      !
    END SUBROUTINE
    !
    SUBROUTINE check_ovl_wannier_real (orb_real_i, orb_real_j, ovl_value)
      !
      USE kinds,                ONLY : DP
      USE control_flags,        ONLY : gamma_only
      USE fft_base,             ONLY : dfftp,dffts
      USE pwcom,                ONLY : omega
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
      USE mp_global,            ONLY : inter_image_comm
      USE bse_module,           ONLY : ovl_thr
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN)  :: orb_real_i(dfftp%nnr)
      REAL(DP), INTENT(IN)  :: orb_real_j(dfftp%nnr)
      REAL(DP), INTENT(OUT) :: ovl_value
      !
      INTEGER   :: ir
      !
      REAL (DP) :: summ_ib, summ_jb, summ_ij
      !
      ovl_value = 0.0_DP
      summ_ib   = 0.0_DP
      summ_jb   = 0.0_DP
      summ_ij   = 0.0_DP
      !
      DO ir = 1, dfftp%nnr
         !
         summ_ib = summ_ib +  orb_real_i(ir)**4
         summ_jb = summ_jb +  orb_real_j(ir)**4
         summ_ij = summ_ij +  orb_real_i(ir)**2 * orb_real_j(ir)**2
         !
      ENDDO
      !
      summ_ib = summ_ib * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      summ_jb = summ_jb * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      summ_ij = summ_ij * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      !
      CALL mp_sum(summ_ib, intra_bgrp_comm)
      CALL mp_sum(summ_jb, intra_bgrp_comm)
      CALL mp_sum(summ_ij, intra_bgrp_comm)
      !
      ovl_value = summ_ij/sqrt(summ_ib*summ_jb)
      !
      RETURN
      !
    END SUBROUTINE
    !
    SUBROUTINE check_ovl_wannier_cmplx (orb_cmpl_i, orb_cmpl_j, ovl_value)
      !
      USE kinds,                ONLY : DP
      USE control_flags,        ONLY : gamma_only
      USE fft_base,             ONLY : dfftp,dffts
      USE pwcom,                ONLY : omega
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
      USE mp_global,            ONLY : inter_image_comm
      USE bse_module,           ONLY : ovl_thr
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(IN)  :: orb_cmpl_i(dfftp%nnr)
      COMPLEX(DP), INTENT(IN)  :: orb_cmpl_j(dfftp%nnr)
      REAL(DP),    INTENT(OUT) :: ovl_value
      !
      INTEGER   :: ir
      !
      REAL (DP) :: summ_ib, summ_jb, summ_ij
      !
      ovl_value = 0.0_DP
      summ_ib = 0.0_DP
      summ_jb = 0.0_DP
      summ_ij = 0.0_DP
      !
      DO ir = 1, dfftp%nnr
         !
         summ_ib = summ_ib + ((DBLE(orb_cmpl_i(ir)))**4 + (AIMAG(orb_cmpl_i(ir)))**4)
         summ_jb = summ_jb + ((DBLE(orb_cmpl_j(ir)))**4 + (AIMAG(orb_cmpl_j(ir)))**4)
         summ_ij = summ_ij + (((DBLE(orb_cmpl_i(ir)))**2 + (AIMAG(orb_cmpl_i(ir)))**2) &
                           *  ((DBLE(orb_cmpl_j(ir)))**2 + (AIMAG(orb_cmpl_j(ir)))**2))
         !
      ENDDO
      !
      summ_ib = summ_ib * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      summ_jb = summ_jb * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      summ_ij = summ_ij * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      !
      CALL mp_sum(summ_ib, intra_bgrp_comm)
      CALL mp_sum(summ_jb, intra_bgrp_comm)
      CALL mp_sum(summ_ij, intra_bgrp_comm)
      !
      ovl_value = summ_ij/sqrt(summ_ib*summ_jb)
      !
      RETURN
      !
    END SUBROUTINE
    !
END MODULE
