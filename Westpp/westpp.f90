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
PROGRAM westpp
  !-----------------------------------------------------------------------
  ! 
  ! This is the main program that generates post-processing data for WEST.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end 
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE westcom,              ONLY : westpp_calculation
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WESTPP'
  INTEGER :: i
  LOGICAL :: lgate(4)
  !
  ! *** START *** 
  !
  CALL check_stop_init ()
  !
  ! Initialize MPI, clocks, print initial messages
  !
#if defined(__MPI)
  CALL mp_startup ( start_images = .TRUE. )
#endif
  !
  CALL west_environment_start ( code )
  !
  CALL westpp_readin( )
  !
  CALL westpp_setup( )
  !
  lgate = .FALSE. 
  DO i = 1, 8
     IF( westpp_calculation(i:i) == 'r' .OR. westpp_calculation(i:i) == 'R' ) lgate(1) = .TRUE. ! Rho --> density 
     IF( westpp_calculation(i:i) == 'w' .OR. westpp_calculation(i:i) == 'W' ) lgate(2) = .TRUE. ! Wavefunction 
     IF( westpp_calculation(i:i) == 'e' .OR. westpp_calculation(i:i) == 'E' ) lgate(3) = .TRUE. ! Eigenpotentials  
     IF( westpp_calculation(i:i) == 's' .OR. westpp_calculation(i:i) == 'S' ) lgate(4) = .TRUE. ! SXX  
  ENDDO
  !
  IF( lgate(1) ) CALL do_rho( )
  IF( lgate(2) ) CALL do_wfc2( )
  IF( lgate(3) ) CALL do_eigenpot2( )
  IF( lgate(4) ) CALL do_sxx( )
  !
  CALL clean_scratchfiles( )
  !
  CALL print_clock(' ')
  !
  CALL west_environment_end( code )
  !
  CALL mp_global_end()
  !
  STOP
  !
9000 FORMAT (/5x,'Program ',a12,' starts ...',/5x,'Today is ',a9,' at ',a9)
  !
END PROGRAM
