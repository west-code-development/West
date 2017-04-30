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
SUBROUTINE fetch_namelist(num_namelists,driver)
  !-----------------------------------------------------------------------
  !
  USE pwcom
  USE westcom
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : mpime,root,world_comm
  USE mp_global,        ONLY : nimage
  USE io_push,          ONLY : io_push_title,io_push_value,io_push_bar,io_push_es0,io_push_c256 
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: num_namelists
  INTEGER, INTENT(IN) :: driver(num_namelists)
  !
  ! Workspace
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: iunit=5
  INTEGER :: i
  INTEGER :: numsp
  !
  ! NAMELISTS
  !
  ! 1
  NAMELIST /input_west/ &
      & qe_prefix, &
      & west_prefix, &
      & outdir
  ! 2
  NAMELIST /wstat_control/ &
      & wstat_calculation, &
      & n_pdep_eigen, &
      & n_pdep_times, &
      & n_pdep_maxiter, &
      & n_dfpt_maxiter, &
      & n_pdep_read_from_file, &
      & trev_pdep, &
      & trev_pdep_rel, &
      & tr2_dfpt, &
      & l_kinetic_only, &
      & l_minimize_exx_if_active, &
      & l_use_ecutrho
  ! 3
  NAMELIST /wfreq_control/ &
      & wfreq_calculation, &
      & n_pdep_eigen_to_use, &
      & qp_bandrange, &
      & macropol_calculation, &
      & n_lanczos, &
      & n_imfreq, &
      & n_refreq, &
      & ecut_imfreq, &
      & ecut_refreq, &
      & wfreq_eta, &
      & n_secant_maxiter, &
      & trev_secant, & 
      & l_enable_lanczos, &
      & l_enable_gwetot, &
      & div_kind_hf, &
      & o_restart_time, &
      & ecut_spectralf, & 
      & n_spectralf
  ! 4
  NAMELIST /westpp_control/ &
      & westpp_calculation, &
      & westpp_range, &
      & westpp_format, &
      & westpp_sign, &
      & westpp_n_pdep_eigen_to_use, &
      & westpp_r0, &
      & westpp_nr, &
      & westpp_rmax, &
      & westpp_epsinfty
  ! 
  !
  CALL start_clock('fetch_nml')
  !
  ! Connect The INPUT FILE TO unit 5, then skip the title line
  !
  IF ( mpime==root ) THEN 
     CALL input_from_file ( )
  ENDIF
  !
  ! NAMELIST 1 : INPUT_WEST
  !
  IF(ANY(driver(:)==1)) THEN
     !
     ! DEFAULTS
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     qe_prefix = 'pwscf'
     west_prefix = 'west'
     !
     ! READ
     !
     IF ( mpime == root) READ(iunit,input_west)
     tmp_dir = trimcheck (outdir)
     !
     ! BCAST
     !
     CALL mp_bcast(qe_prefix,root,world_comm)
     prefix=qe_prefix
     CALL mp_bcast(west_prefix,root,world_comm)
     CALL mp_bcast(tmp_dir,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title("I/O Summary : input_west")
     !
     numsp = 14
     CALL io_push_c256('qe_prefix',qe_prefix,numsp)
     CALL io_push_c256('west_prefix',qe_prefix,numsp)
     CALL io_push_c256('outdir',outdir,numsp)
     !
     CALL io_push_bar()
     !
  ENDIF
  !
  ! NAMELIST 2 : WSTAT_CONTROL
  !
  IF(ANY(driver(:)==2)) THEN
     !
     ! DEFAULTS
     !
     wstat_calculation        = 'S'
     n_pdep_eigen             = 1
     n_pdep_times             = 4
     n_pdep_maxiter           = 100
     n_dfpt_maxiter           = 250
     n_pdep_read_from_file    = 0
     trev_pdep                = 1.d-3
     trev_pdep_rel            = 1.d-1
     tr2_dfpt                 = 1.d-12
     l_kinetic_only           = .false.
     l_minimize_exx_if_active = .false. 
     l_use_ecutrho            = .false.
     !
     ! READ
     !
     IF ( mpime == root) READ(iunit,wstat_control)
     !
     ! BCAST
     !
     CALL mp_bcast(wstat_calculation,root,world_comm)
     CALL mp_bcast(n_pdep_eigen,root,world_comm)
     CALL mp_bcast(n_pdep_times,root,world_comm)
     CALL mp_bcast(n_pdep_maxiter,root,world_comm)
     CALL mp_bcast(n_dfpt_maxiter,root,world_comm)
     CALL mp_bcast(n_pdep_read_from_file,root,world_comm)
     CALL mp_bcast(trev_pdep,root,world_comm)
     CALL mp_bcast(trev_pdep_rel,root,world_comm)
     CALL mp_bcast(tr2_dfpt,root,world_comm)
     CALL mp_bcast(l_kinetic_only,root,world_comm)
     CALL mp_bcast(l_minimize_exx_if_active,root,world_comm)
     CALL mp_bcast(l_use_ecutrho,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title('I/O Summary : wstat_control')
     !
     numsp=30
     CALL io_push_value('wstat_calculation',wstat_calculation,numsp)
     CALL io_push_value('n_pdep_eigen',n_pdep_eigen,numsp)
     CALL io_push_value('n_pdep_times',n_pdep_times,numsp)
     CALL io_push_value('n_pdep_maxiter',n_pdep_maxiter,numsp)
     CALL io_push_value('n_dfpt_maxiter',n_dfpt_maxiter,numsp)
     CALL io_push_value('n_pdep_read_from_file',n_pdep_read_from_file,numsp)
     CALL io_push_es0('trev_pdep',trev_pdep,numsp)
     CALL io_push_es0('trev_pdep_rel',trev_pdep_rel,numsp)
     CALL io_push_es0('tr2_dfpt',tr2_dfpt,numsp)
     CALL io_push_value('l_kinetic_only',l_kinetic_only,numsp)
     CALL io_push_value('l_minimize_exx_if_active',l_minimize_exx_if_active,numsp)
     CALL io_push_value('l_use_ecutrho',l_use_ecutrho,numsp)
     !
     CALL io_push_bar()
     !
     ! CHECKS 
     !
     SELECT CASE(wstat_calculation) 
     CASE('r','R','s','S')
     CASE DEFAULT
        CALL errore('fetch_nml','Err: wstat_calculation /= S or R',1)
     END SELECT
     !
     IF( n_pdep_times < 2 ) CALL errore('fetch_nml','Err: n_pdep_times<2',1) 
     IF( n_pdep_eigen < 1 ) CALL errore('fetch_nml','Err: n_pdep_eigen<1',1)
     IF( n_pdep_eigen*n_pdep_times < nimage ) CALL errore('fetch_nml','Err: n_pdep_eigen*n_pdep_times<nimage',1) 
     IF( n_pdep_maxiter < 1 ) CALL errore('fetch_nml','Err: n_pdep_maxiter<1',1) 
     IF( n_dfpt_maxiter < 1 ) CALL errore('fetch_nml','Err: n_dfpt_maxiter<1',1) 
     IF( n_pdep_read_from_file < 0 ) CALL errore('fetch_nml','Err: n_pdep_read_from_file<0',1) 
     IF( n_pdep_read_from_file > n_pdep_eigen ) CALL errore('fetch_nml','Err: n_pdep_read_from_file>n_pdep_eigen',1) 
     IF(tr2_dfpt<=0._DP) CALL errore('fetch_nml','Err: tr2_dfpt<0.',1)
     IF(trev_pdep<=0._DP) CALL errore('fetch_nml','Err: trev_pdep<0.',1)
     IF(trev_pdep_rel<=0._DP) CALL errore('fetch_nml','Err: trev_pdep_rel<0.',1)
     !
  ENDIF
  !
  ! NAMELIST 3 : WFREQ_CONTROL
  !
  IF(ANY(driver(:)==3)) THEN
     !
     ! DEFAULTS
     !
     wfreq_calculation       = 'XWGQ'
     n_pdep_eigen_to_use     = 2
     qp_bandrange            = (/ 1, 2 /)
     macropol_calculation    = 'N'
     n_lanczos               = 20
     n_imfreq                = 10
     n_refreq                = 10
     ecut_imfreq             = 1._DP
     ecut_refreq             = 2._DP
     wfreq_eta               = 0.003675_DP
     n_secant_maxiter        = 1
     trev_secant             = 0.003675_DP
     l_enable_lanczos        = .TRUE.
     l_enable_gwetot         = .FALSE.
     div_kind_hf             = 2 
     o_restart_time          = 0._DP
     ecut_spectralf          = (/ -2._DP, 2._DP /)
     n_spectralf             = 10
     !
     ! READ
     !
     IF ( mpime == root) READ(iunit,wfreq_control)
     !
     ! BCAST
     !
     CALL mp_bcast(wfreq_calculation,root,world_comm)
     CALL mp_bcast(n_pdep_eigen_to_use,root,world_comm)
     CALL mp_bcast(qp_bandrange,root,world_comm)
     CALL mp_bcast(macropol_calculation,root,world_comm)
     CALL mp_bcast(n_lanczos,root,world_comm)
     CALL mp_bcast(n_imfreq,root,world_comm)
     CALL mp_bcast(n_refreq,root,world_comm)
     CALL mp_bcast(ecut_imfreq,root,world_comm)
     CALL mp_bcast(ecut_refreq,root,world_comm)
     CALL mp_bcast(wfreq_eta,root,world_comm)
     CALL mp_bcast(n_secant_maxiter,root,world_comm)
     CALL mp_bcast(trev_secant,root,world_comm)
     CALL mp_bcast(l_enable_lanczos,root,world_comm)
     CALL mp_bcast(l_enable_gwetot,root,world_comm)
     CALL mp_bcast(div_kind_hf,root,world_comm)
     CALL mp_bcast(o_restart_time,root,world_comm)
     CALL mp_bcast(ecut_spectralf,root,world_comm)
     CALL mp_bcast(n_spectralf,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title('I/O Summary : wfreq_control')
     !
     numsp=40
     CALL io_push_value('wfreq_calculation',wfreq_calculation,numsp)
     CALL io_push_value('n_pdep_eigen_to_use',n_pdep_eigen_to_use,numsp)
     CALL io_push_value('qp_bandrange(1)',qp_bandrange(1),numsp)
     CALL io_push_value('qp_bandrange(2)',qp_bandrange(2),numsp)
     CALL io_push_value('macropol_calculation',macropol_calculation,numsp)
     CALL io_push_value('n_lanczos',n_lanczos,numsp)
     CALL io_push_value('n_imfreq',n_imfreq,numsp)
     CALL io_push_value('n_refreq',n_refreq,numsp)
     CALL io_push_value('ecut_imfreq [Ry]',ecut_imfreq,numsp)
     CALL io_push_value('ecut_refreq [Ry]',ecut_refreq,numsp)
     CALL io_push_value('wfreq_eta [Ry]',wfreq_eta,numsp)
     CALL io_push_value('n_secant_maxiter',n_secant_maxiter,numsp)
     CALL io_push_value('trev_secant [Ry]',trev_secant,numsp)
     CALL io_push_value('l_enable_lanczos',l_enable_lanczos,numsp)
     CALL io_push_value('l_enable_gwetot',l_enable_gwetot,numsp)
     CALL io_push_value('div_kind_hf',div_kind_hf,numsp)
     CALL io_push_value('o_restart_time [min]',o_restart_time,numsp)
     CALL io_push_value('ecut_spectralf(1) [Ry]',ecut_spectralf(1),numsp)
     CALL io_push_value('ecut_spectralf(2) [Ry]',ecut_spectralf(2),numsp)
     CALL io_push_value('n_spectralf',n_spectralf,numsp)
     !
     CALL io_push_bar()
     !
     ! CHECKS 
     !
     IF( n_lanczos < 2 ) CALL errore('fetch_nml','Err: n_lanczos<2',1) 
     IF( n_pdep_eigen_to_use < 1 ) CALL errore('fetch_nml','Err: n_pdep_eigen_to_use<1',1) 
     IF( n_pdep_eigen_to_use > n_pdep_eigen ) CALL errore('fetch_nml','Err: n_pdep_eigen_to_use>n_pdep_eigen',1) 
     IF( n_imfreq < 1 ) CALL errore('fetch_nml','Err: n_imfreq<1',1) 
     IF( n_refreq < 1 ) CALL errore('fetch_nml','Err: n_refreq<1',1) 
     IF( n_spectralf < 2 ) CALL errore('fetch_nml','Err: n_spectralf<1',1) 
     IF( qp_bandrange(1) < 1 ) CALL errore('fetch_nml','Err: qp_bandrange(1)<1',1) 
     IF( qp_bandrange(2) < 1 ) CALL errore('fetch_nml','Err: qp_bandrange(2)<1',1) 
     IF( qp_bandrange(2) < qp_bandrange(1) ) CALL errore('fetch_nml','Err: qp_bandrange(2)<qp_bandrange(1)',1) 
     IF( ecut_imfreq<=0._DP) CALL errore('fetch_nml','Err: ecut_imfreq<0.',1)
     IF( ecut_refreq<=0._DP) CALL errore('fetch_nml','Err: ecut_imfreq<0.',1)
     IF( ecut_spectralf(2)<ecut_spectralf(1)) CALL errore('fetch_nml','Err: ecut_spectralf(2)<ecut_spectralf(1)',1)
     IF( wfreq_eta<=0._DP) CALL errore('fetch_nml','Err: wfreq_eta<0.',1)
     IF( n_secant_maxiter < 1 ) CALL errore('fetch_nml','Err: n_secant_maxiter<1',1) 
     IF( trev_secant<=0._DP) CALL errore('fetch_nml','Err: trev_secant<0.',1)
     SELECT CASE(macropol_calculation) 
     CASE('N','n','C','c')
     CASE DEFAULT
        CALL errore('fetch_nml','Err: macropol_calculation /= (N,C)',1)
     END SELECT
     !
  ENDIF
  !
  ! NAMELIST 4 : WESTPP_CONTROL
  !
  IF(ANY(driver(:)==4)) THEN
     !
     ! DEFAULTS
     !
     westpp_calculation         = 'r'
     westpp_range               = (/ 1, 2 /)
     westpp_format              = 'C'
     westpp_sign                = .FALSE.
     westpp_n_pdep_eigen_to_use = 1
     westpp_r0                  = (/ 0.d0, 0.d0, 0.d0 /)
     westpp_nr                  = 100
     westpp_rmax                = 1.d0
     westpp_epsinfty            = 1.d0
     !
     ! READ
     !
     IF ( mpime == root) READ(iunit,westpp_control)
     !
     ! BCAST
     !
     CALL mp_bcast(westpp_calculation,root,world_comm)
     CALL mp_bcast(westpp_range,root,world_comm)
     CALL mp_bcast(westpp_format,root,world_comm)
     CALL mp_bcast(westpp_sign,root,world_comm)
     CALL mp_bcast(westpp_n_pdep_eigen_to_use,root,world_comm)
     CALL mp_bcast(westpp_r0,root,world_comm)
     CALL mp_bcast(westpp_nr,root,world_comm)
     CALL mp_bcast(westpp_rmax,root,world_comm)
     CALL mp_bcast(westpp_epsinfty,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title('I/O Summary : westpp_control')
     !
     numsp=40
     CALL io_push_value('westpp_calculation',westpp_calculation,numsp)
     CALL io_push_value('westpp_range(1)',westpp_range(1),numsp)
     CALL io_push_value('westpp_range(2)',westpp_range(2),numsp)
     CALL io_push_value('westpp_format',westpp_format,numsp)
     CALL io_push_value('westpp_sign',westpp_sign,numsp)
     CALL io_push_value('westpp_n_pdep_eigen_to_use',westpp_n_pdep_eigen_to_use,numsp)
     CALL io_push_value('westpp_r0(1)',westpp_r0(1),numsp)
     CALL io_push_value('westpp_r0(2)',westpp_r0(2),numsp)
     CALL io_push_value('westpp_r0(3)',westpp_r0(3),numsp)
     CALL io_push_value('westpp_nr',westpp_nr,numsp)
     CALL io_push_value('westpp_rmax',westpp_rmax,numsp)
     CALL io_push_value('westpp_epsinfty',westpp_epsinfty,numsp)
     !
     CALL io_push_bar()
     !
     ! CHECKS 
     !
     IF( westpp_range(1) < 1 ) CALL errore('fetch_nml','Err: westpp_range(1)<1',1)
     IF( westpp_range(2) < 1 ) CALL errore('fetch_nml','Err: westpp_range(2)<1',1)
     IF( westpp_range(2) < westpp_range(1) ) CALL errore('fetch_nml','Err: westpp_range(2)<westpp_range(1)',1)
     IF( westpp_nr < 1 ) CALL errore('fetch_nml','Err: westpp_nr<1',1)
     IF( westpp_n_pdep_eigen_to_use < 1 ) CALL errore('fetch_nml','Err: westpp_n_pdep_eigen_to_use<1',1)
     IF( westpp_rmax < 0.d0 ) CALL errore('fetch_nml','Err: westpp_rmax<0',1)
     IF( westpp_epsinfty < 1.d0 ) CALL errore('fetch_nml','Err: westpp_epsinfty<1',1)
     !
  ENDIF
  !
  CALL stop_clock('fetch_nml')
  !
END SUBROUTINE
