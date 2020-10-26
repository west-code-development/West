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
SUBROUTINE wbse_fetch_namelist( num_drivers, driver, verbose )
  !-----------------------------------------------------------------------
  !
  USE json_module,      ONLY : json_file
  USE westcom
  USE wbsecom
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : mpime,root,world_comm
  USE mp_global,        ONLY : nimage
  USE io_push,          ONLY : io_push_title,io_push_value,io_push_bar,io_push_es0,io_push_c512 
  USE gvect,            ONLY : ecutrho
  USE start_k,          ONLY : nk1, nk2, nk3
  USE control_flags,    ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: num_drivers
  INTEGER, INTENT(IN) :: driver(num_drivers)
  LOGICAL, INTENT(IN) :: verbose
  !
  ! Workspace
  !
  TYPE(json_file) :: json
  INTEGER :: i
  INTEGER :: iiarg, nargs
  INTEGER :: numsp
  LOGICAL :: found
  CHARACTER (LEN=512) :: input_file
  CHARACTER(LEN=512), EXTERNAL :: trimcheck
  CHARACTER(LEN=:),ALLOCATABLE :: cval
  REAL(DP) :: rval 
  INTEGER :: ival 
  INTEGER,ALLOCATABLE :: ivec(:)
  REAL(DP),ALLOCATABLE :: rvec(:)
  LOGICAL :: lval 
  INTEGER :: iunit
  ! 
  CALL start_clock('fetch_input')
  !
  ! PRESETS
  !
  ! ** input_west **
  IF ( ANY(driver(:)==1) ) THEN
     qe_prefix = 'pwscf'
     west_prefix = 'west'
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     !l_load_qbox_wfc      = .false.
     !qbox_ks_wfc_filename = 'qb.xml' 
  ENDIF
  !
  ! ** wbse_init **
  IF ( ANY(driver(:)==2) ) THEN
     ! 
     ! ** WARNING ** : In order to properly initialize these variables, this driver 
     !                 can be called only after: 
     !                 - fetch_input( driver 1 ) 
     !                 - read_pwout() 
     !
     !wbse_which_screening    = 'PDEP'
     !wbse_init_calculation   = 'S'
     !n_pdep_eigen            =  1
     !chi_kernel              = 'CHI'
     !which_spin_channel      = 'up'
     !which_overlap_method    = 'N'!B!W
     !overlap_thr             =  0.0_DP
     !bisec_ovl_thr_filename  = 'bisection.dat'
     !qbox_bisec_wf_filename  = 'qb.xml'
     which_bse_method         = 'PDEP'
     wbse_init_calculation    = 'S'
     n_pdep_eigen             =  1
     l_use_bisection_thr      = .false.
     l_use_localise_repr      = .false.
     overlap_thr              =  0.0
     which_spin_channel       =  0
     !
     chi_kernel               = 'CHI'!, "CHI_RPA"
!     qbox_bisec_wfc_filename  = 'qb.bi.xml'
     !
  ENDIF
  !
  ! ** wbse_control  **
  IF ( ANY(driver(:)==3) ) THEN
     ! 
     ! ** WARNING ** : In order to properly initialize these variables, this driver 
     !                 can be called only after: 
     !                 - fetch_input( driver 1 ) 
     !                 - read_pwout() 
     !
     !wbse_calculation        = 'S'
     !wbse_approx_method      = 'IPA' ! 'BSE' or 'TDDFT'
     !wbse_diag_method        = 'david'
     !n_plep_eigen            =  1
     !n_plep_times            =  4
     !n_plep_maxiter          =  10
     !n_plep_read_from_file   =  0
     !l_qp_correction         = .TRUE.
     !trev_plep               =  0.01
     !trev_plep_rel           =  0.1 
     !scissor_ope             =  0.0
     !eps_macro               =  1.0
     !ipol_input              = 'X' ! Y Z or XYZ  
     !n_lanczos_steps         =  1
     !macropol_dfpt           = .TRUE.
     !
     wbse_calculation         = 'S'
     n_plep_eigen             = 1
     n_plep_times             = 4
     n_plep_maxiter           = 100
     n_plep_read_from_file    = 0
     spin_excitation          = 'singlet'
     ! 
     l_diag_term_only         = .false.
     l_bse_calculation        = .true.
     l_preconditioning        = .true.
     l_qp_correction          = .false.
     !
     trev_plep                =  1.D-3
     trev_plep_rel            =  1.D-1
     scissor_ope              =  0.0
     eps_macro                =  1.0
     !
     wbse_diag_method         = 'david'
     !
     ! LZ
     !
     ipol_input               = 'XX'
     n_lzstep                 = 0
     macropol_dfpt            = .false.
     !
  ENDIF
  ! 
  ! ** wbsepp_control **
  !IF ( ANY(driver(:)==4) ) THEN
  !   westpp_r0                  = (/ 0.d0, 0.d0, 0.d0 /)
  !   westpp_nr                  = 100
  !   westpp_rmax                = 1.d0
  !   westpp_epsinfty            = 1.d0
  !ENDIF
  !
  ! READ the input
  !
  IF ( mpime==root ) THEN
     !
     !
     CALL json%initialize()
     CALL json%load_file( filename = main_input_file )
     !
     IF ( ANY(driver(:)==1) ) THEN 
        CALL json%get('input_west.qe_prefix', cval, found)
        IF( found ) qe_prefix = cval  
        CALL json%get('input_west.west_prefix', cval, found) 
        IF( found )  west_prefix = cval  
        CALL json%get('input_west.outdir', cval, found) 
        IF( found ) outdir = cval
!        CALL json%get('input_west.l_load_qbox_wfc', lval, found)
!        IF( found )  l_load_qbox_wfc = lval 
!        CALL json%get('input_west.qbox_ks_wfc_filename', cval, found)
!        IF( found )  qbox_ks_wf_filename = cval  
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN
        CALL json%get('wbse_init.which_bse_method', cval, found) 
        IF( found ) which_bse_method  = cval  
        CALL json%get('wbse_init.wbse_init_calculation', cval, found) 
        IF( found ) wbse_init_calculation  = cval
        CALL json%get('wbse_init.n_pdep_eigen', ival, found)
        IF( found ) n_pdep_eigen  = ival 
        CALL json%get('wbse_init.l_use_bisection_thr', lval, found)
        IF( found ) l_use_bisection_thr  = lval  
        CALL json%get('wbse_init.l_use_bisection_thr', lval, found)
        IF( found ) l_use_bisection_thr  = lval
        CALL json%get('wbse_init.l_use_localise_repr', lval, found)
        IF( found ) l_use_localise_repr  = lval
        CALL json%get('wbse_init.overlap_thr', rval, found)
        IF( found ) overlap_thr  = rval
        CALL json%get('wbse_init.which_spin_channel', ival, found)
        IF( found ) which_spin_channel  = ival
        CALL json%get('wbse_init.chi_kernel', cval, found)
        IF( found ) chi_kernel  = cval
!        CALL json%get('wbse_init.qbox_bisec_wfc_filename', cval, found)
!        IF( found ) qbox_bisec_wfc_filename  = cval
        !CALL json%get('wbse_init.wbse_which_screening', cval, found) 
        !IF( found ) wbse_which_screening  = cval  
        !CALL json%get('wbse_init.wbse_init_calculation', cval, found) 
        !IF( found ) wbse_init_calculation = cval
        !CALL json%get('wbse_init.n_pdep_eigen', ival, found)
        !IF( found ) n_pdep_eigen = ival  
        !CALL json%get('wbse_init.chi_kernel', cval, found) 
        !IF( found ) chi_kernel = cval  
        !CALL json%get('wbse_init.which_spin_channel', cval, found)
        !IF( found ) which_spin_channel = cval
        !CALL json%get('wbse_init.which_overlap_method', cval, found)
        !IF( found ) which_overlap_method = cval
        !CALL json%get('wbse_init.overlap_thr', rval, found)
        !IF( found ) overlap_thr = rval
        !CALL json%get('wbse_init.bisec_ovl_thr_filename', cval, found)
        !IF( found ) bisec_ovl_thr_filename = cval
        !CALL json%get('wbse_init.qbox_bisec_wf_filename', cval, found)
        !IF( found ) qbox_bisec_wf_filename = cval
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN 
        CALL json%get('wbse_control.wbse_calculation', cval, found)
        IF( found ) wbse_calculation = cval
        CALL json%get('wbse_control.n_plep_eigen', ival, found)
        IF( found ) n_plep_eigen = ival
        CALL json%get('wbse_control.n_plep_times', ival, found)
        IF( found ) n_plep_times = ival
        CALL json%get('wbse_control.n_plep_maxiter', ival, found)
        IF( found ) n_plep_maxiter = ival
        CALL json%get('wbse_control.n_plep_read_from_file', ival, found)
        IF( found ) n_plep_read_from_file = ival
        CALL json%get('wbse_control.spin_excitation', cval, found)
        IF( found ) spin_excitation = cval
        CALL json%get('wbse_control.l_diag_term_only', lval, found)
        IF( found ) l_diag_term_only = lval
        CALL json%get('wbse_control.l_bse_calculation', lval, found)
        IF( found ) l_bse_calculation = lval
        CALL json%get('wbse_control.l_preconditioning', lval, found)
        IF( found ) l_preconditioning = lval
        CALL json%get('wbse_control.l_qp_correction', lval, found)
        IF( found ) l_qp_correction = lval
        CALL json%get('wbse_control.trev_plep', rval, found)
        IF( found ) trev_plep = rval
        CALL json%get('wbse_control.trev_plep_rel', rval, found)
        IF( found ) trev_plep_rel = rval
        CALL json%get('wbse_control.scissor_ope', rval, found)
        IF( found ) scissor_ope = rval
        CALL json%get('wbse_control.eps_macro', rval, found)
        IF( found ) eps_macro = rval
        CALL json%get('wbse_control.wbse_diag_method', cval, found)
        IF( found ) wbse_diag_method = cval
        CALL json%get('wbse_control.ipol_input', cval, found)
        IF( found ) ipol_input = cval
        CALL json%get('wbse_control.n_lzstep', ival, found)
        IF( found ) n_lzstep = ival
        CALL json%get('wbse_control.macropol_dfpt', lval, found)
        IF( found ) macropol_dfpt = lval
!        CALL json%get('wbse_control.wbse_calculation', cval, found)
!        IF( found ) wbse_calculation = cval
!        CALL json%get('wbse_control.which_approx_method', cval, found)
!        IF( found ) wbse_approx_method = cval
!        CALL json%get('wbse_control.wbse_diag_method', cval, found)
!        IF( found ) wbse_diag_method = cval
!        CALL json%get('wbse_control.n_plep_eigen', ival, found)
!        IF( found ) n_plep_eigen = ival
!        CALL json%get('wbse_control.n_plep_times', ival, found)
!        IF( found ) n_plep_times = ival
!        CALL json%get('wbse_control.n_plep_maxiter', ival, found)
!        IF( found ) n_plep_maxiter = ival
!        CALL json%get('wbse_control.n_plep_read_from_file', ival, found)
!        IF( found ) n_plep_read_from_file = ival
!        CALL json%get('wbse_control.l_qp_correction', lval, found)
!        IF( found ) l_qp_correction = lval
!        CALL json%get('wbse_control.trev_plep', rval, found)
!        IF( found ) trev_plep = rval
!        CALL json%get('wbse_control.trev_plep_rel', rval, found)
!        IF( found ) trev_plep_rel = rval
!        CALL json%get('wbse_control.scissor_ope', rval, found)
!        IF( found ) scissor_ope = rval
!        CALL json%get('wbse_control.eps_macro', rval, found)
!        IF( found ) eps_macro = rval
!        CALL json%get('wbse_control.ipol_input', cval, found)
!        IF( found ) ipol_input = cval
!        CALL json%get('wbse_control.n_lanczos_steps', ival, found)
!        IF( found ) n_lanczos_steps = ival
!        CALL json%get('wbse_control.macropol_dfpt', lval, found)
!        IF( found ) macropol_dfpt = lval
     ENDIF
     !
     CALL json%destroy()
     !
  ENDIF
  !
  ! BCAST & CHECKS
  !
  IF ( ANY(driver(:)==1) ) THEN
     !
     CALL mp_bcast(qe_prefix,root,world_comm)
     prefix=qe_prefix
     CALL mp_bcast(west_prefix,root,world_comm)
     tmp_dir = trimcheck (outdir)
     CALL mp_bcast(tmp_dir,root,world_comm)
     !
     !CALL mp_bcast(l_load_qbox_wfc,root,world_comm)
     !CALL mp_bcast(qbox_ks_wf_filename,root,world_comm)
!     CALL mp_bcast(l_load_qbox_wfc,root,world_comm)
!     CALL mp_bcast(qbox_ks_wfc_filename,root,world_comm)
     !
  ENDIF
  !
  IF ( ANY(driver(:)==2) ) THEN
     !
     !CALL mp_bcast(wbse_which_screening,root,world_comm)
     !CALL mp_bcast(wbse_init_calculation,root,world_comm)
     !CALL mp_bcast(n_pdep_eigen,root,world_comm)
     !CALL mp_bcast(chi_kernel,root,world_comm)
     !CALL mp_bcast(which_spin_channel,root,world_comm)
     !CALL mp_bcast(which_overlap_method,root,world_comm)
     !CALL mp_bcast(overlap_thr,root,world_comm)
     !CALL mp_bcast(bisec_ovl_thr_filename,root,world_comm)
     !CALL mp_bcast(qbox_bisec_wf_filename,root,world_comm)
     CALL mp_bcast(which_bse_method,root,world_comm)
     CALL mp_bcast(wbse_init_calculation,root,world_comm)
     CALL mp_bcast(n_pdep_eigen,root,world_comm)
     CALL mp_bcast(l_use_bisection_thr,root,world_comm)
     CALL mp_bcast(l_use_localise_repr,root,world_comm)
     CALL mp_bcast(overlap_thr,root,world_comm)
     CALL mp_bcast(which_spin_channel,root,world_comm)
     CALL mp_bcast(chi_kernel,root,world_comm)
!     CALL mp_bcast(qbox_bisec_wfc_filename,root,world_comm)
     !
     ! CHECKS 
     !
     l_test_ovl = .FALSE.
     !
     IF (wbse_init_calculation == 'I') l_test_ovl = .TRUE.
     ! 
!     use_qbox = .FALSE.
     use_wstat_pdep = .FALSE.
     SELECT CASE(which_bse_method)
     CASE('finite-field','FF')
!     use_qbox = .TRUE.
     CASE('PDEP','pdep')
     use_wstat_pdep = .TRUE.
     CASE DEFAULT
        CALL errore('wbse_fetch_nml','Err: which_bse_method != FF/PDEP', 1)
     END SELECT
     !
  ENDIF
  !
  IF ( ANY(driver(:)==3) ) THEN
     !
     !CALL mp_bcast(wbse_calculation,root,world_comm)
     !CALL mp_bcast(wbse_approx_method,root,world_comm)
     !CALL mp_bcast(wbse_diag_method,root,world_comm)
     !CALL mp_bcast(n_plep_eigen,root,world_comm)
     !CALL mp_bcast(n_plep_times,root,world_comm)
     !CALL mp_bcast(n_plep_maxiter,root,world_comm)
     !CALL mp_bcast(n_plep_read_from_file,root,world_comm)
     !CALL mp_bcast(l_qp_correction,root,world_comm)
     !CALL mp_bcast(trev_plep,root,world_comm)
     !CALL mp_bcast(trev_plep_rel,root,world_comm)
     !CALL mp_bcast(scissor_ope,root,world_comm)
     !CALL mp_bcast(eps_macro,root,world_comm)
     !CALL mp_bcast(eps_macro,root,world_comm)
     !CALL mp_bcast(ipol_input,root,world_comm)
     !CALL mp_bcast(n_lanczos_steps,root,world_comm)
     !CALL mp_bcast(macropol_dfpt,root,world_comm)
     !
     CALL mp_bcast(wbse_calculation,root,world_comm)
     CALL mp_bcast(n_plep_eigen,root,world_comm)
     CALL mp_bcast(n_plep_times,root,world_comm)
     CALL mp_bcast(n_plep_maxiter,root,world_comm)
     CALL mp_bcast(n_plep_read_from_file,root,world_comm)
     CALL mp_bcast(spin_excitation,root,world_comm)
     CALL mp_bcast(l_bse_calculation,root,world_comm)
     CALL mp_bcast(l_diag_term_only,root,world_comm)
     CALL mp_bcast(l_preconditioning,root,world_comm)
     CALL mp_bcast(wbse_diag_method,root,world_comm)
     CALL mp_bcast(trev_plep,root,world_comm)
     CALL mp_bcast(trev_plep_rel,root,world_comm)
     CALL mp_bcast(scissor_ope,root,world_comm)
     CALL mp_bcast(eps_macro,root,world_comm)
     CALL mp_bcast(l_qp_correction,root,world_comm)
     CALL mp_bcast(ipol_input,root,world_comm)
     CALL mp_bcast(n_lzstep,root,world_comm)
     CALL mp_bcast(macropol_dfpt,root,world_comm)
     !
     ! CHECKS 
     !
     SELECT CASE(wbse_calculation)
     CASE('r','R','s','S')
     CASE DEFAULT
        CALL errore('wbse_fetch_nml','Err: wbse_calculation /= S or R',1)
     END SELECT
     !
     wlz_calculation = wbse_calculation
     l_davidson = .FALSE.
     l_lanzcos  = .FALSE.
     !
     SELECT CASE(wbse_diag_method)
     CASE('david', 'davidson')
         l_davidson = .TRUE.
     CASE('lanzcos', 'Lanzcos')
         l_lanzcos  = .TRUE.
     CASE DEFAULT
        CALL errore('wbse_fetch_nml','Err: wbse_diag_method /= david or lanzcos',1)
     END SELECT
     !
     SELECT CASE(spin_excitation)
     CASE('s', 'S', 'singlet')
         l_bse_triplet = .false.
     CASE('t', 'T', 'triplet')
         l_bse_triplet = .true.
     CASE DEFAULT
        CALL errore('wbse_fetch_nml','Err: spin_excitation /= s or t',1)
     END SELECT
     !
     IF (l_davidson) THEN
        IF( n_plep_times < 2 ) CALL errore('wbse_fetch_nml','Err: n_plep_times<2',1)
        IF( n_plep_eigen < 1 ) CALL errore('wbse_fetch_nml','Err: n_plep_eigen<1',1)
        IF( n_plep_eigen*n_plep_times < nimage ) CALL errore('wbse_fetch_nml','Err: n_plep_eigen*n_plep_times<nimage',1)
        IF( n_plep_maxiter < 1 ) CALL errore('wbse_fetch_nml','Err: n_plep_maxiter<1',1)
        IF( n_plep_read_from_file < 0 ) CALL errore('wbse_fetch_nml','Err: n_plep_read_from_file<0',1)
        IF( n_plep_read_from_file > n_plep_eigen ) CALL errore('wbse_fetch_nml','Err: n_plep_read_from_file>n_plep_eigen',1)
     ENDIF
     !
  ENDIF
  !
  ! REPORT
  !
  IF ( verbose ) THEN
     !
     IF ( ANY(driver(:)==1) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title("I/O Summary : input_west")
        !
        numsp = 14
        CALL io_push_c512('qe_prefix',qe_prefix,numsp)
        CALL io_push_c512('west_prefix',qe_prefix,numsp)
        CALL io_push_c512('outdir',outdir,numsp)
!        CALL io_push_value('l_load_qbox_wfc',l_load_qbox_wfc,numsp)
!        CALL io_push_c512('qbox_ks_wfc_filename',qbox_ks_wfc_filename,numsp)
        !
        CALL io_push_bar()
        !
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title('I/O Summary : wbse_init')
        !
        numsp = 23
        CALL io_push_value('wbse_init_calculation',wbse_init_calculation,numsp)
        CALL io_push_value('n_pdep_eigen',n_pdep_eigen,numsp)
        CALL io_push_value('chi_kernel',chi_kernel,numsp)
        CALL io_push_value('which_spin_channel',which_spin_channel,numsp)
        CALL io_push_value('l_use_localise_repr',l_use_localise_repr,numsp)
        CALL io_push_value('overlap_thr',overlap_thr,numsp)
        CALL io_push_value('l_use_bisection_thr',l_use_bisection_thr,numsp)
!        CALL io_push_value('qbox_bisec_wfc_filename',qbox_bisec_wfc_filename,numsp)
        !
        CALL io_push_bar()
        !
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title('I/O Summary : wbse_control')
        !
        numsp=23
        !
        IF (l_davidson) THEN
           !  
           CALL io_push_value('wbse_diag_method',wbse_diag_method,numsp)
           CALL io_push_value('wbse_calculation',wbse_calculation,numsp)
           CALL io_push_value('n_plep_eigen',n_plep_eigen,numsp)
           CALL io_push_value('n_plep_times',n_plep_times,numsp)
           CALL io_push_value('n_plep_maxiter',n_plep_maxiter,numsp)
           CALL io_push_value('n_plep_read_from_file',n_plep_read_from_file,numsp)
           CALL io_push_value('l_bse_calculation',l_bse_calculation,numsp)
           CALL io_push_value('l_diag_term_only',l_diag_term_only,numsp)
           CALL io_push_value('l_preconditioning',l_preconditioning,numsp)
           CALL io_push_value('l_qp_correction',l_qp_correction,numsp)
           CALL io_push_value('spin_excitation', spin_excitation,numsp)
           CALL io_push_value('trev_plep',trev_plep,numsp)
           CALL io_push_value('trev_plep_rel',trev_plep_rel,numsp)
           CALL io_push_value('scissor_ope',scissor_ope,numsp)
           CALL io_push_value('eps_macro',eps_macro,numsp)
           !
        ENDIF
        !
        IF (l_lanzcos) THEN
           !
           CALL io_push_value('wbse_diag_method',wbse_diag_method,numsp)
           CALL io_push_value('wbse_calculation',wbse_calculation,numsp)
           CALL io_push_value('ipol_input',ipol_input,numsp)
           CALL io_push_value('n_lzstep',n_lzstep,numsp)
           CALL io_push_value('l_bse_calculation',l_bse_calculation,numsp)
           CALL io_push_value('l_diag_term_only',l_diag_term_only,numsp)
           CALL io_push_value('l_qp_correction',l_qp_correction,numsp)
           CALL io_push_value('spin_excitation', spin_excitation,numsp)
           CALL io_push_value('scissor_ope',scissor_ope,numsp)
           CALL io_push_value('eps_macro',eps_macro,numsp)
           CALL io_push_value('macropol_dfpt',macropol_dfpt,numsp)
           ! 
        ENDIF
        !
        CALL io_push_bar()
        !
     ENDIF
     !
!     IF ( ANY(driver(:)==4) ) THEN
        !
        ! REPORT
        !
!        CALL io_push_title('I/O Summary : westpp_control')
        !
        !
!        CALL io_push_bar()
        !
!     ENDIF
     !
     IF( mpime == root ) THEN
        !
        CALL json%initialize()
        CALL json%load_file(filename=TRIM(logfile))
        !
        CALL add_wbse_intput_parameters_to_json_file( num_drivers, driver, json )
        ! 
        OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
        CALL json%print_file( iunit )
        CLOSE( iunit )
        CALL json%destroy()
        !
     ENDIF
     !
  ENDIF
  !
  CALL stop_clock('fetch_input')
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE add_wbse_intput_parameters_to_json_file( num_drivers, driver, json )
  !-----------------------------------------------------------------------
  !
  USE json_module,      ONLY : json_file
  USE pwcom
  USE westcom
  USE wbsecom
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : mpime,root,world_comm
  USE mp_global,        ONLY : nimage
  USE io_push,          ONLY : io_push_title,io_push_value,io_push_bar,io_push_es0,io_push_c512 
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: num_drivers
  INTEGER, INTENT(IN) :: driver(num_drivers)
  TYPE(json_file), INTENT(INOUT) :: json 
  !
  IF ( mpime == root ) THEN
     !
     IF ( ANY(driver(:)==1) ) THEN
        !
        CALL json%add('input.input_west.qe_prefix',TRIM(qe_prefix))
        CALL json%add('input.input_west.west_prefix',TRIM(west_prefix))
        CALL json%add('input.input_west.outdir',TRIM(outdir))
!        CALL json%add('input.input_west.l_load_qbox_wfc',l_load_qbox_wfc)
!        CALL json%add('input.input_west.qbox_ks_wfc_filename',TRIM(qbox_ks_wfc_filename))
        !
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN
        !
        CALL json%add('input.wbse_init.wbse_init_calculation',TRIM(wbse_init_calculation))
        CALL json%add('input.wbse_ini.n_pdep_eigen',n_pdep_eigen)
        CALL json%add('input.wbse_ini.chi_kernel',TRIM(chi_kernel))
        CALL json%add('input.wbse_ini.which_spin_channel',which_spin_channel)
        CALL json%add('input.wbse_ini.l_use_localise_repr',l_use_localise_repr)
        CALL json%add('input.wbse_ini.overlap_thr',overlap_thr)
        CALL json%add('input.wbse_ini.l_use_bisection_thr',l_use_bisection_thr)
!        CALL json%add('input.wbse_ini.qbox_bisec_wfc_filename',TRIM(qbox_bisec_wfc_filename))
        !
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN
        !
        CALL json%add('input.wbse_control.wbse_calculation',TRIM(wbse_calculation))
        CALL json%add('input.wbse_control.n_plep_eigen',n_plep_eigen)
        CALL json%add('input.wbse_control.n_plep_times',n_plep_times)
        CALL json%add('input.wbse_control.n_plep_maxiter',n_plep_maxiter)
        CALL json%add('input.wbse_control.n_plep_read_from_file',n_plep_read_from_file)
        CALL json%add('input.wbse_control.spin_excitation',TRIM(spin_excitation))
        CALL json%add('input.wbse_control.l_diag_term_only',l_diag_term_only)
        CALL json%add('input.wbse_control.l_bse_calculation',l_bse_calculation)
        CALL json%add('input.wbse_control.l_qp_correction',l_qp_correction)
        CALL json%add('input.wbse_control.trev_plep',trev_plep)
        CALL json%add('input.wbse_control.trev_plep_rel',trev_plep_rel)
        CALL json%add('input.wbse_control.scissor_ope',scissor_ope)
        CALL json%add('input.wbse_control.eps_macro',eps_macro)
        CALL json%add('input.wbse_control.wbse_diag_method',TRIM(wbse_diag_method))
        CALL json%add('input.wbse_control.ipol_input',TRIM(ipol_input))
        CALL json%add('input.wbse_control.n_lzstep',n_lzstep)
        CALL json%add('input.wbse_control.macropol_dfpt',macropol_dfpt)
        !
     ENDIF
     !
     IF ( ANY(driver(:)==4) ) THEN
        !
        !
     ENDIF
     !
  ENDIF
  !
END SUBROUTINE
