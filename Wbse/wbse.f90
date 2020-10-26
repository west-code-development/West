!
! Copyright (C) 2015-2017 M. Govoni 
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
PROGRAM wbse
  !-----------------------------------------------------------------------
  ! 
  ! This is the main program that calculates the static screening.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE wbsecom,              ONLY : l_davidson, l_lanzcos 
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WBSE'
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
  CALL wbse_readin ( )
  !
  CALL wbse_setup ( )
  !
  IF (l_davidson) THEN
     CALL wbse_davidson_diago ( )
  ENDIF
  !
  IF (l_lanzcos) THEN
     CALL wbse_lanzcos_diago ( )
  ENDIF
  !
  CALL exx_ungo ( )
  !
  CALL clean_scratchfiles( )
  !
  CALL west_print_clocks( )
  !
  CALL west_environment_end( code )
  !
  CALL mp_global_end()
  !
END PROGRAM
