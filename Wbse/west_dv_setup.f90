!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE west_dv_setup (l_bse_calc)
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares some variables which are needed for derivatives
  !  1) Set the nonlinear core correction 
  !  2) Computes dmuxc (derivative of the XC potential)
  !  3) Set gradient correction (GGA) if needed
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : ntyp => nsp
  USE fft_base,              ONLY : dfftp
  USE uspp_param,            ONLY : upf
  USE spin_orb,              ONLY : domag
  USE uspp,                  ONLY : nlcc_any
  USE noncollin_module,      ONLY : noncolin
  USE eqv,                   ONLY : dmuxc
  USE funct,                 ONLY : dft_is_gradient, exx_is_active
  USE wavefunctions_module,  ONLY : psic
  USE lsda_mod,              ONLY : nspin
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: l_bse_calc 
  !
  CALL start_clock ('wbse_dv_setup')
  !
  IF (l_bse_calc) THEN ! BSE
     !
     dmuxc = (0.0_DP, 0.0_DP)
     !
     RETURN
     ! 
  ENDIF 
  ! 
  ! 0) allocate dmuxc
  ! 
  ALLOCATE(dmuxc (dfftp%nnr, nspin, nspin))
  !
  ! 1) Set the nonlinear core correction
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  ! 2) Compute the derivative of the XC potential
  !
  CALL setup_dmuxc()
  !
  ! 3) Setup gradient correction
  !
  IF (dft_is_gradient()) THEN
     !
     IF (noncolin .AND. domag) THEN
        IF (.NOT.ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
        psic(:) = (0.0_dp, 0.0_dp)
     ENDIF
     !
     CALL setup_dgc()
     !
     IF (ALLOCATED(psic)) DEALLOCATE(psic)
     !
  ENDIF
  !
  CALL stop_clock ('wbse_dv_setup')
  !
  RETURN
  !
END SUBROUTINE west_dv_setup
