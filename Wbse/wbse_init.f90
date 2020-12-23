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
PROGRAM wbse_init
  !-----------------------------------------------------------------------
  !
  ! This is the main program that calculates the static screening.
  !
  !USE io_global,            ONLY : stdout
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  !USE mp,                   ONLY : mp_sum,mp_barrier
  !USE qbox_interface,       ONLY : init_qbox_interface, finalize_qbox_interface
  USE qbox_interface,       ONLY : init_qbox,finalize_qbox
  USE westcom,              ONLY : l_test_ovl, use_wstat_pdep,use_qbox
  !wbsecom combined into westcom
  !USE wbsecom,              ONLY : l_test_ovl, use_wstat_pdep
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WBSE_INIT'
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
  CALL wbse_init_readin ( )
  !
  IF (l_test_ovl) THEN
     !
     use_qbox = .False.
     !
  ENDIF
  !
  CALL wbse_setup ( )
  !
  IF ( use_qbox ) THEN
     !
     CALL init_qbox()
     !CALL init_qbox_interface()
     !
     CALL wbse_init_methods (.TRUE.)
     !
     CALL finalize_qbox()
     !CALL finalize_qbox_interface()
     !
  ELSEIF (use_wstat_pdep) THEN
     !
     CALL wbse_init_methods (.False.)
     !
  ELSE
     !
     CALL wbse_init_fock_energy()
     !
  ENDIF
  !
  CALL exx_ungo ( )
  !
  CALL wbse_clear ( )
  !
  CALL clean_scratchfiles( )
  !
  CALL west_print_clocks( )
  !CALL print_clock(' ')
  !
  CALL west_environment_end( code )
  !
  CALL mp_global_end()
  !
!9000 FORMAT (/5x,'Program ',a12,' starts ...',/5x,'Today is ',a9,' at ',a9)
  !
END PROGRAM
