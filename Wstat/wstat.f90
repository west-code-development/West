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
PROGRAM wstat
  !-----------------------------------------------------------------------
  ! 
  ! This is the main program that calculates the static screening.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE wavefunctions_module, ONLY : evc
  USE function3d
  USE pwcom,                ONLY : npw,npwx
  USE westcom,              ONLY : n_pdep_eigen
  ! 
  IMPLICIT NONE
  !
  INTEGER :: nx, ny, nz, i
  CHARACTER(LEN=9) :: code = 'WSTAT'
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
  CALL wstat_readin ( )
  !
  CALL wstat_setup ( )
  !
  CALL mp_barrier( world_comm )
  !
  nx = n_pdep_eigen
  ny = n_pdep_eigen
  nz = n_pdep_eigen
  !
  PRINT*, 'nx, ny, nz = ', nx, ny, nz
  DO i = 1, 100
      !
      PRINT*, 'Writing ', i, '%'
      CALL write_function3d( 'wfcl.f3d', nx, ny, nz, npw, npwx, evc(:, 3))
      !
  ENDDO
  !
  CALL mp_barrier( world_comm )
  !
  DO i = 1, 100
      !
      PRINT*, 'Reading ', i, '%'
      CALL read_function3d ( 'wfcl.f3d', nx, ny, nz, npw, npwx, evc(:, 3))
      !
  ENDDO
  !
  CALL mp_barrier( world_comm )
  !
  PRINT*, 'Reading finished'
  !
  CALL west_print_clocks( )
  !
  PRINT*, 'STOP'
  !
  STOP
  !
  CALL davidson_diago ( )
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
