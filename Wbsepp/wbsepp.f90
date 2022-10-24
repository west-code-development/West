!
! Copyright (C) 2015-2021 M. Govoni
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
  ! This is the main program that calculates the static screening.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  USE westcom,              ONLY : l_eig_decomp, l_lz_spec, l_exc_plot, l_exc_rho_res_plot
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'Wbsepp'
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
  CALL wbsepp_readin ( )
  !
  CALL wbsepp_setup ( )
  !
  IF (l_eig_decomp) CALL wbsepp_decompose_eig_contributions ( )
  IF (l_exc_plot)   CALL wbsepp_plot_exc ()
  IF (l_exc_rho_res_plot) CALL wbsepp_plot_charged_density_res_exc ()
  !TODO: several version of wbsepp_meg()
  !IF (l_meg) CALL wbsepp_meg ()
  IF (l_lz_spec) CALL wbsepp_ads_spectrum ()
  !
  CALL exx_ungo ( )
  !
  CALL clean_scratchfiles( )
  !
  CALL print_clock(' ')
  !
  CALL west_environment_end( code )
  !
  CALL mp_global_end()
  !
END PROGRAM
