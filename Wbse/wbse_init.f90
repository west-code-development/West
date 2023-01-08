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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
PROGRAM wbse_init
  !-----------------------------------------------------------------------
  !
  ! This is the main program that calculates the static screening.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup,mp_global_end
  USE west_environment,     ONLY : west_environment_start,west_environment_end
  USE qbox_interface,       ONLY : init_qbox,finalize_qbox
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WBSE_INIT'
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
  CALL wbse_init_readin( )
  !
  CALL wbse_init_setup( )
  !
  CALL init_qbox( )
  !
  CALL wbse_init_methods( )
  !
  CALL finalize_qbox( )
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
