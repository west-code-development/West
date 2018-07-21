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
SUBROUTINE wbse_fetch_namelist(num_namelists,driver)
  !-----------------------------------------------------------------------
  !
  USE pwcom
  USE westcom
  USE wbsecom
!  USE qbox_interface
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
      & qe_prefix,      &
      & west_prefix,    &
      & outdir,         &
      & l_load_qbox_wfc,&
      & qbox_ks_wfc_filename
  ! 2
  NAMELIST /wbse_init/ &
      & which_bse_method, &
      & wbse_init_calculation, &
      & n_pdep_eigen, &
      & chi_kernel,  &
      & which_spin_channel, &
      & overlap_thr, &
      & l_use_localise_repr, &
      & l_use_bisection_thr, &
      & qbox_bisec_wfc_filename
  ! 3
  NAMELIST /wbse_control/ &
      & wbse_calculation, &
      & n_plep_eigen, &
      & n_plep_times, &
      & n_plep_maxiter, &
      & n_plep_read_from_file, &
      & spin_excitation,&
      & l_bse_calculation, &
      & l_qp_correction, &
      & l_diag_term_only,  &
      & l_preconditioning, &
      & trev_plep, &
      & trev_plep_rel, &
      & scissor_ope, &
      & eps_macro, &
      & wbse_diag_method, &
      & ipol_input, &
      & n_lzstep, &
      & macropol_dfpt
  ! 4
  NAMELIST /qbox_control/ &
      & nrowmax, &
      & xml_file, &
      & xc, &
      & alpha_pbe0, &
      & amplitude, &
      & wf_dyn, &
      & btHF, &
      & blHF, &
      & nitscf, &
      & nite, &
      & io
  ! 
  CALL start_clock('wbse_fetch_nml')
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
     which_bse_method = 'finite_field'
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
     CALL mp_bcast(which_bse_method,root,world_comm)
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
     CALL io_push_value('which_bse_method',which_bse_method,numsp)
     !
     CALL io_push_bar()
     !
     ! CHECK 
     !
     use_qbox = .FALSE. 
     use_wstat_pdep = .FALSE. 
     SELECT CASE(which_bse_method)
     CASE('finite-field','FF')
     use_qbox = .TRUE. 
     CASE('PDEP','pdep')
     use_wstat_pdep = .TRUE. 
     CASE DEFAULT
        CALL errore('wbse_fetch_nml','Err: which_bse_method != FF/PDEP', 1)
     END SELECT
     !
  ENDIF
  !
  IF(ANY(driver(:)==2)) THEN
     !
     ! DEFAULTS
     !
     wbse_init_calculation   = 'S'
     n_pdep_eigen            = 1
     l_use_bisection_thr     = .false.
     l_use_localise_repr     = .false.
     overlap_thr             = 0.0
     which_spin_channel      = 0
     !
     chi_kernel              = 'CHI'!, "CHI_RPA"
     qbox_bisec_wfc_filename = 'qb.bi.xml'
     !
     ! READ
     !
     IF ( mpime == root) READ(iunit,wbse_init)
     !
     ! BCAST
     !
     CALL mp_bcast(wbse_init_calculation,root,world_comm)
     CALL mp_bcast(l_use_bisection_thr,root,world_comm)
     CALL mp_bcast(overlap_thr,root,world_comm)
     CALL mp_bcast(chi_kernel,root,world_comm)
     CALL mp_bcast(n_pdep_eigen,root,world_comm)
     CALL mp_bcast(which_spin_channel,root,world_comm)
     CALL mp_bcast(l_use_localise_repr,root,world_comm)
     CALL mp_bcast(qbox_bisec_wfc_filename,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title("I/O Summary : wbse_init")
     !
     numsp = 23
     CALL io_push_value('wbse_init_calculation',wbse_init_calculation,numsp)
     CALL io_push_value('n_pdep_eigen',n_pdep_eigen,numsp)
     CALL io_push_value('chi_kernel',chi_kernel,numsp)
     CALL io_push_value('which_spin_channel',which_spin_channel,numsp)
     CALL io_push_value('l_use_localise_repr',l_use_localise_repr,numsp)
     CALL io_push_value('overlap_thr',overlap_thr,numsp)
     CALL io_push_value('l_use_bisection_thr',l_use_bisection_thr,numsp)
     CALL io_push_value('qbox_bisec_wfc_filename',qbox_bisec_wfc_filename,numsp)
     !
     CALL io_push_bar()
     !
     l_test_ovl = .FALSE.
     !
     IF (wbse_init_calculation == 'I') l_test_ovl = .TRUE.
     !  
  ENDIF
  !
  ! NAMELIST 3 : WBSE_CONTROL
  !
  IF(ANY(driver(:)==3)) THEN
     !
     ! DEFAULTS
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
     l_preconditioning        = .false.
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
     ! READ
     !
     IF ( mpime == root) READ(iunit,wbse_control)
     !
     ! BCAST
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
     ! DISPLAY
     !
     CALL io_push_title('I/O Summary : WBSE_control')
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
  ! NAMELIST 4 : QBOX_CONTROL
  !
  IF(ANY(driver(:)==4) .AND. use_qbox) THEN
     !
     ! DEFAULTS
     !
     nrowmax    = 0             ! band parallelization control parameter, see Qbox manual
     xml_file   = 'qb.init.xml' ! xml file generated by Qbox calculation, used to initialize Qbox
     xc         = 'PBE'
     alpha_pbe0 = 0.25          ! alpha for PBE0 calculation
     wf_dyn     = 'PSDA'        ! wavefunction update algorithm, see Qbox manual
     blHF       = '2 2 2'       ! bisection level HF calculation, see Qbox manual
     btHF       = 0             ! bisection thres HF calculation, see Qbox manual
     nitscf     = 10            ! number of self-consistent iterations
     nite       = 0             ! number of non-self-consistent iterations between each SCF iteration
     io         = 'cube'
     !
     ! READ
     !
     IF ( mpime == root) READ(iunit,qbox_control)
     !
     ! BCAST
     !
     CALL mp_bcast(nrowmax,root,world_comm)
     CALL mp_bcast(xml_file,root,world_comm)
     CALL mp_bcast(xc,root,world_comm)
     CALL mp_bcast(alpha_pbe0,root,world_comm)
     CALL mp_bcast(wf_dyn,root,world_comm)
     CALL mp_bcast(blHF,root,world_comm)
     CALL mp_bcast(btHF,root,world_comm)
     CALL mp_bcast(nitscf,root,world_comm)
     CALL mp_bcast(nite,root,world_comm)
     CALL mp_bcast(io,root,world_comm)
     !
     ! DISPLAY
     !
     CALL io_push_title('I/O Summary : qbox_control')
     !
     numsp=23
     IF ( nrowmax /=0 ) CALL io_push_value('nrowmax',nrowmax,numsp)
     CALL io_push_value('xml_file',xml_file,numsp)
     CALL io_push_value('xc',xc,numsp)
     IF ( TRIM(xc) == 'PBE0' ) CALL io_push_value('alpha_pbe0',alpha_pbe0,numsp)
     CALL io_push_value('wf_dyn',wf_dyn,numsp)
     CALL io_push_value('blHF',blHF,numsp)
     CALL io_push_value('btHF',btHF,numsp)
     CALL io_push_value('nitscf',nitscf,numsp)
     CALL io_push_value('nite',nite,numsp)
     CALL io_push_value('io',io,numsp)
     !
  ENDIF
  !
  CALL stop_clock('wbse_fetch_nml')
  !
END SUBROUTINE
