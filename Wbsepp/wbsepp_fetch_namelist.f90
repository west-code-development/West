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
SUBROUTINE wbsepp_fetch_namelist(num_namelists,driver)
  !-----------------------------------------------------------------------
  !
  USE pwcom
  USE westcom
  USE wbseppcom
  USE wbsecom,          ONLY : n_plep_read_from_file, macropol_dfpt 
  USE qbox_interface
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
      & outdir, &
      & l_load_qbox_wfc, &
      & qbox_ks_wfc_filename
  ! 2
  NAMELIST /wbsepp_input/ &
      & wbsepp_type, & ! 1:eig_decomposion, 2:meg, 3:ads_spect
      !
      ! for eig_decomposion
      !
      & n_plep_read_from_file, &
      & macropol_dfpt, & 
      !
      ! for exc plot
      !
      & r0_input,  &
      & iexc_plot, & 
      !
      ! for lanzcos
      !
      & itermax, &
      & itermax0, &
      & ipol, &
      & sym_op, &
      & units, &
      & verbosity, &
      & extrapolation, & 
      & start, &
      & end, & 
      & increment, &
      & epsil, &
      & spin_channel
      ! 
      ! for meg
      !
  ! 
  CALL start_clock('wbsepp_fetch_nml')
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
     l_load_qbox_wfc = .false.
     qbox_ks_wfc_filename = 'qb.xml'
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
     CALL mp_bcast(l_load_qbox_wfc,root,world_comm)
     CALL mp_bcast(qbox_ks_wfc_filename,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title("I/O Summary : input_west")
     !
     numsp = 23
     CALL io_push_value('qe_prefix',qe_prefix,numsp)
     CALL io_push_value('west_prefix',west_prefix,numsp)
     CALL io_push_value('outdir',outdir,numsp)
     CALL io_push_value('l_load_qbox_wfc',l_load_qbox_wfc,numsp)
     CALL io_push_value('qbox_ks_wfc_filename',qbox_ks_wfc_filename,numsp)
     !
     CALL io_push_bar()
     !
  ENDIF
  !
  IF(ANY(driver(:)==2)) THEN
     !
     ! DEFAULTS
     !
     wbsepp_type = 0
     !
     n_plep_read_from_file = 0
     macropol_dfpt = .FALSE.
     !
     r0_input = (/ 0.d0, 0.d0, 0.d0 /)
     iexc_plot = 1
     !
     itermax  = 1000
     itermax0 = 1000
     extrapolation = "no"
     end = 2.50d0
     increment = 0.001d0
     start = 0.0d0
     epsil = 0.02
     ipol = 1
     sym_op = 0
     verbosity = 0
     units = 0
     spin_channel = 1
     !
     ! READ
     !
     IF ( mpime == root) READ(iunit,wbsepp_input)
     !
     ! BCAST
     !
     CALL mp_bcast(wbsepp_type,root,world_comm)
     ! 
     CALL mp_bcast(n_plep_read_from_file,root,world_comm)
     CALL mp_bcast(macropol_dfpt,root,world_comm)
     !
     CALL mp_bcast(r0_input,root,world_comm)
     CALL mp_bcast(iexc_plot,root,world_comm)
     !
     CALL mp_bcast(itermax,root,world_comm)
     CALL mp_bcast(itermax0,root,world_comm)
     CALL mp_bcast(ipol,root,world_comm)
     CALL mp_bcast(sym_op,root,world_comm)
     CALL mp_bcast(units,root,world_comm)
     CALL mp_bcast(verbosity,root,world_comm)
     CALL mp_bcast(extrapolation,root,world_comm)
     CALL mp_bcast(start,root,world_comm)
     CALL mp_bcast(end,root,world_comm)
     CALL mp_bcast(increment,root,world_comm)
     CALL mp_bcast(epsil,root,world_comm)
     CALL mp_bcast(spin_channel,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title("I/O Summary : wbsepp_input")
     !
     l_meg        = .FALSE.
     l_eig_decomp = .FALSE.
     l_lz_spec    = .FALSE.
     l_exc_plot   = .FALSE.
     l_exc_rho_res_plot = .FALSE.
     !
     SELECT CASE(wbsepp_type)
     CASE( 1 )
         l_eig_decomp = .TRUE.
     CASE( 2 )
         l_lz_spec  = .TRUE.
     CASE( 3 )
         l_meg      = .TRUE.
     CASE( 4 )
         l_exc_plot = .TRUE.
     CASE( 5 )
         l_exc_rho_res_plot = .TRUE.
     CASE DEFAULT
        CALL errore('wbsepp_fetch_nml','Err: wbsepp_type /= 1,2,3', wbsepp_type)
     END SELECT
     !
     IF ((spin_channel < 1) .OR. (spin_channel > 2)) THEN
        CALL errore('wbsepp_fetch_nml','Err: spin_channel/= 1,2', spin_channel)
     ENDIF
     ! 
     ! DISPLAY
     !
     numsp=30
     !
     IF (l_eig_decomp) THEN
        !  
        CALL io_push_value('wbsepp_type',wbsepp_type,numsp)
        CALL io_push_value('n_plep_read_from_file',n_plep_read_from_file,numsp)
        CALL io_push_value('macropol_dfpt',macropol_dfpt,numsp)
        !
     ENDIF
     !
     IF (l_exc_rho_res_plot) THEN
        !  
        CALL io_push_value('wbsepp_type',wbsepp_type,numsp)
        CALL io_push_value('n_plep_read_from_file',n_plep_read_from_file,numsp)
        CALL io_push_value('iexc_plot',iexc_plot,numsp)
        !
     ENDIF
     !
     IF (l_exc_plot ) THEN
        !  
        CALL io_push_value('wbsepp_type',wbsepp_type,numsp)
        CALL io_push_value('n_plep_read_from_file',n_plep_read_from_file,numsp)
        CALL io_push_value('r0_input(1) [alat]',r0_input(1),numsp)
        CALL io_push_value('r0_input(2) [alat]',r0_input(2),numsp)
        CALL io_push_value('r0_input(3) [alat]',r0_input(3),numsp)
        CALL io_push_value('iexc_plot',iexc_plot,numsp)
        !
     ENDIF
     !
     IF (l_lz_spec) THEN
        !
        CALL io_push_value('wbsepp_type',wbsepp_type,numsp)
        CALL io_push_value('lziter_max', itermax,numsp)
        CALL io_push_value('lziter_cal', itermax0,numsp)
        CALL io_push_value('l_extrapolation', extrapolation,numsp)
        CALL io_push_value('w_start', start,numsp)
        CALL io_push_value('w_stop', end,numsp)
        CALL io_push_value('w_delta', increment,numsp)
        CALL io_push_value('w_unit', units,numsp)
        CALL io_push_value('ipol_index',ipol,numsp)
        CALL io_push_value('sym_opt',sym_op,numsp)
        CALL io_push_value('spin_channel',spin_channel,numsp)
        !
     ENDIF
     ! 
  ENDIF
  !
  CALL stop_clock('wbsepp_fetch_nml')
  !
END SUBROUTINE
