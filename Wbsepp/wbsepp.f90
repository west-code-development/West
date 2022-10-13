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
PROGRAM wbse_pp
  !-----------------------------------------------------------------------
  !
  ! This is the main program that calculates the static screening.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE westcom,            ONLY : l_meg, l_eig_decomp, l_lz_spec, l_exc_plot, &
                                   l_exc_rho_res_plot
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
9000 FORMAT (/5x,'Program ',a12,' starts ...',/5x,'Today is ',a9,' at ',a9)
  !
  !
END PROGRAM
    !
    !
    SUBROUTINE wbsepp_setup
      !
      USE io_global,              ONLY : stdout
      USE types_coulomb,          ONLY : pot3D
      USE westcom,                ONLY : alphapv_dfpt,npwq,wstat_save_dir,west_prefix,nbnd_occ,&
                                       & n_pdep_basis,n_pdep_eigen,n_pdep_times,l_use_ecutrho,nbndval0x
      USE mp,                     ONLY : mp_max
      USE mp_global,              ONLY : intra_bgrp_comm
      USE pwcom,                  ONLY : npw,npwx
      USE kinds,                  ONLY : DP
      USE wavefunctions,          ONLY : evc
      USE gvect,                  ONLY : gstart,g,ig_l2g,ngm,ngmx
      USE constants,              ONLY : e2,fpi
      USE cell_base,              ONLY : tpiba2, alat
      USE io_files,               ONLY : tmp_dir
      USE westcom,                ONLY : l_meg, l_eig_decomp, l_lz_spec, l_exc_plot, &
                                         l_exc_rho_res_plot, wbse_save_dir
      !
      IMPLICIT NONE
      !
      REAL(DP) :: q(3)
      REAL(DP) :: qq
      COMPLEX(DP),EXTERNAL :: get_alpha_pv
      INTEGER :: ig
      !
      CALL do_setup ( )
      !
      l_use_ecutrho = .FALSE.
      !
      CALL set_npwq()
      !
      CALL pot3D%init('Dense',.FALSE.,'gb')
      !
      CALL set_nbndocc()
      !
      nbndval0x = nbnd_occ(1)
      !
      wbse_save_dir = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.save'
      !
    END SUBROUTINE

SUBROUTINE wbsepp_readin()
  !
  USE pwcom
  USE westcom
  USE ions_base,        ONLY : nat
  USE uspp,             ONLY : okvan
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : nproc,mpime,root
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: iunit =5, ios
  !
  CALL start_clock('wbsepp_readin')
  !
  CALL fetch_input_yml(2,(/1,8/),.TRUE.,.FALSE.)
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  CALL read_pwout( )
  !
  ! PW checks
  !
  IF (domag) CALL errore('wbse_readin','domag version not available',1)
  IF (okvan) CALL errore('wbse_readin','ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) CALL errore('wbse_readin', 'double grid not implemented',1)
  !
  CALL stop_clock('wbsepp_readin')
  !
END SUBROUTINE
