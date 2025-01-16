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
! Marco Govoni
!
!-----------------------------------------------------------------------
PROGRAM westpp
  !-----------------------------------------------------------------------
  !
  ! This is the main program that generates post-processing data for WEST.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  USE westcom,              ONLY : westpp_calculation
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WESTPP'
  INTEGER :: i
  LOGICAL :: lgate(11)
  !
  ! *** START ***
  !
  CALL check_stop_init( )
  !
  ! Initialize MPI, clocks, print initial messages
  !
#if defined(__MPI)
  CALL mp_startup( start_images = .TRUE. )
#endif
  !
  CALL west_environment_start( code )
  !
  CALL west_readin( code )
  !
  CALL westpp_setup( )
  !
  lgate = .FALSE.
  DO i = 1, 11
     !
     ! Rho --> density
     IF( westpp_calculation(i:i) == 'r' .OR. westpp_calculation(i:i) == 'R' ) lgate(1) = .TRUE.
     !
     ! Wavefunction
     IF( westpp_calculation(i:i) == 'w' .OR. westpp_calculation(i:i) == 'W' ) lgate(2) = .TRUE.
     !
     ! Eigenpotentials
     IF( westpp_calculation(i:i) == 'e' .OR. westpp_calculation(i:i) == 'E' ) lgate(3) = .TRUE.
     !
     ! Screened exact exchange
     IF( westpp_calculation(i:i) == 's' .OR. westpp_calculation(i:i) == 'S' ) lgate(4) = .TRUE.
     !
     ! Dipole matrix
     IF( westpp_calculation(i:i) == 'd' .OR. westpp_calculation(i:i) == 'D' ) lgate(5) = .TRUE.
     !
     ! Localization factor and inverse participation ratio
     IF( westpp_calculation(i:i) == 'l' .OR. westpp_calculation(i:i) == 'L' ) lgate(6) = .TRUE.
     !
     ! Exciton
     IF( westpp_calculation(i:i) == 'x' .OR. westpp_calculation(i:i) == 'X' ) lgate(7) = .TRUE.
     !
     ! Density response to exciton
     IF( westpp_calculation(i:i) == 'p' .OR. westpp_calculation(i:i) == 'P' ) lgate(8) = .TRUE.
     !
     ! Boys / Wannier localization
     IF( westpp_calculation(i:i) == 'b' .OR. westpp_calculation(i:i) == 'B' ) lgate(9) = .TRUE.
     !
     ! Decomposition of BSE/TDDFT excited state and calculation of dipole moments
     IF( westpp_calculation(i:i) == 'c' .OR. westpp_calculation(i:i) == 'C' ) lgate(10) = .TRUE.
     !
     ! Spin multiplicity of BSE/TDDFT excited state (<S^2>, for nspin==2 only)
     IF( westpp_calculation(i:i) == 'm' .OR. westpp_calculation(i:i) == 'M' ) lgate(11) = .TRUE.
     !
  ENDDO
  !
  IF( lgate(1) ) CALL do_rho( )
  IF( lgate(2) ) CALL do_wfc2( )
  IF( lgate(3) ) CALL do_eigenpot2( )
  IF( lgate(4) ) CALL do_sxx( )
  IF( lgate(5) ) CALL do_dip( )
  IF( lgate(6) ) CALL do_loc( )
  IF( lgate(7) ) CALL do_exc( )
  IF( lgate(8) ) CALL do_resp( )
  IF( lgate(9) ) CALL do_wann( )
  IF( lgate(10) ) CALL do_exc_comp( )
  IF( lgate(11) ) CALL do_exc_spin( )
  !
  CALL exx_ungo( )
  !
  CALL clean_scratchfiles( )
  !
  CALL west_print_clocks( )
  !
  CALL west_environment_end( code )
  !
  CALL mp_global_end( )
  !
END PROGRAM
