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
PROGRAM wbsepp
  !-----------------------------------------------------------------------
  !
  ! This is the main program that generates post-processing data for Wbse.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  USE westcom,              ONLY : wbsepp_calculation
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WBSEPP'
  INTEGER :: i
  LOGICAL :: lgate(4)
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
  CALL wbsepp_readin( )
  !
  CALL wbsepp_setup( )
  !
  lgate = .FALSE.
  DO i = 1, 4
     IF( wbsepp_calculation(i:i) == 's' .OR. wbsepp_calculation(i:i) == 'S' ) lgate(1) = .TRUE. ! Spectrum
     IF( wbsepp_calculation(i:i) == 'p' .OR. wbsepp_calculation(i:i) == 'P' ) lgate(2) = .TRUE. ! Eig decomposition
     IF( wbsepp_calculation(i:i) == 'e' .OR. wbsepp_calculation(i:i) == 'E' ) lgate(3) = .TRUE. ! Exciton states
     IF( wbsepp_calculation(i:i) == 'd' .OR. wbsepp_calculation(i:i) == 'D' ) lgate(4) = .TRUE. ! Density response
  ENDDO
  !
  IF( lgate(1) ) CALL do_spectrum( )
  IF( lgate(2) ) CALL do_eig_decomp( )
  IF( lgate(3) ) CALL do_exc( )
  IF( lgate(4) ) CALL do_density_resp( )
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
