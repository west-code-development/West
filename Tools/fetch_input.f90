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
SUBROUTINE fetch_input( num_drivers, driver, verbose )
  !-----------------------------------------------------------------------
  !
  USE json_module,      ONLY : json_file
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
  INTEGER, INTENT(IN) :: num_drivers
  INTEGER, INTENT(IN) :: driver(num_drivers)
  LOGICAL, INTENT(IN) :: verbose
  !
  ! Workspace
  !
  type(json_file) :: json
  INTEGER :: i
  INTEGER :: iiarg, nargs
  INTEGER :: numsp
  LOGICAL :: found
  CHARACTER (LEN=256) :: input_file
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=:),ALLOCATABLE :: cval
  REAL(DP) :: rval 
  INTEGER :: ival 
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
  ENDIF
  !
  ! ** wstat_control **
  IF ( ANY(driver(:)==2) ) THEN
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
  ENDIF
  !
  ! ** wfreq_control **
  IF ( ANY(driver(:)==3) ) THEN
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
  ENDIF
  ! 
  ! ** westpp_control **
  IF ( ANY(driver(:)==4) ) THEN
     westpp_calculation         = 'r'
     westpp_range               = (/ 1, 2 /)
     westpp_format              = 'C'
     westpp_sign                = .FALSE.
     westpp_n_pdep_eigen_to_use = 1
     westpp_r0                  = (/ 0.d0, 0.d0, 0.d0 /)
     westpp_nr                  = 100
     westpp_rmax                = 1.d0
     westpp_epsinfty            = 1.d0
  ENDIF
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
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN 
        CALL json%get('wstat_control.wstat_calculation', cval, found) 
        IF( found ) wstat_calculation = cval  
        CALL json%get('wstat_control.n_pdep_eigen', ival, found) 
        IF( found ) n_pdep_eigen = ival  
        CALL json%get('wstat_control.n_pdep_times', ival, found) 
        IF( found ) n_pdep_times = ival  
        CALL json%get('wstat_control.n_pdep_maxiter', ival, found) 
        IF( found ) n_pdep_maxiter = ival  
        CALL json%get('wstat_control.n_dfpt_maxiter', ival, found) 
        IF( found ) n_dfpt_maxiter = ival  
        CALL json%get('wstat_control.n_pdep_read_from_file', ival, found) 
        IF( found ) n_pdep_read_from_file = ival  
        CALL json%get('wstat_control.trev_pdep', rval, found) 
        IF( found ) trev_pdep = rval  
        CALL json%get('wstat_control.trev_pdep_rel', rval, found) 
        IF( found ) trev_pdep_rel = rval  
        CALL json%get('wstat_control.tr2_dfpt', rval, found) 
        IF( found ) tr2_dfpt = rval  
        CALL json%get('wstat_control.l_kinetic_only', lval, found) 
        IF( found ) l_kinetic_only = lval  
        CALL json%get('wstat_control.l_minimize_exx_if_active', lval, found) 
        IF( found ) l_minimize_exx_if_active = lval  
        CALL json%get('wstat_control.l_use_ecutrho', lval, found) 
        IF( found ) l_use_ecutrho = lval  
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN 
        CALL json%get('wfreq_control.wfreq_calculation', cval, found)
        IF( found ) wfreq_calculation = cval  
        CALL json%get('wfreq_control.n_pdep_eigen_to_use', ival, found) 
        IF( found ) n_pdep_eigen_to_use = ival  
        CALL json%get('wfreq_control.qp_bandrange(1)', rval, found) 
        IF( found ) qp_bandrange(1) = rval  
        CALL json%get('wfreq_control.qp_bandrange(2)', rval, found) 
        IF( found ) qp_bandrange(2) = rval  
        CALL json%get('wfreq_control.macropol_calculation', cval, found) 
        IF( found ) macropol_calculation = cval  
        CALL json%get('wfreq_control.n_lanczos', ival, found) 
        IF( found ) n_lanczos = ival  
        CALL json%get('wfreq_control.n_imfreq', ival, found) 
        IF( found ) n_imfreq = ival  
        CALL json%get('wfreq_control.n_refreq', ival, found) 
        IF( found ) n_refreq = ival  
        CALL json%get('wfreq_control.ecut_imfreq', rval, found) 
        IF( found ) ecut_imfreq = rval  
        CALL json%get('wfreq_control.ecut_refreq', rval, found) 
        IF( found ) ecut_refreq = rval  
        CALL json%get('wfreq_control.wfreq_eta', rval, found) 
        IF( found ) wfreq_eta = rval  
        CALL json%get('wfreq_control.n_secant_maxiter', ival, found) 
        IF( found ) n_secant_maxiter = ival  
        CALL json%get('wfreq_control.trev_secant', rval, found) 
        IF( found ) trev_secant = rval  
        CALL json%get('wfreq_control.l_enable_lanczos', lval, found) 
        IF( found ) l_enable_lanczos = lval  
        CALL json%get('wfreq_control.l_enable_gwetot', lval, found) 
        IF( found ) l_enable_gwetot = lval  
        CALL json%get('wfreq_control.div_kind_hf', ival, found) 
        IF( found ) div_kind_hf = ival  
        CALL json%get('wfreq_control.o_restart_time', rval, found) 
        IF( found ) o_restart_time = rval  
        CALL json%get('wfreq_control.ecut_spectralf(1)', rval, found) 
        IF( found ) ecut_spectralf(1) = rval  
        CALL json%get('wfreq_control.ecut_spectralf(2)', rval, found) 
        IF( found ) ecut_spectralf(2) = rval  
        CALL json%get('wfreq_control.n_spectralf', ival, found) 
        IF( found ) n_spectralf = ival  
     ENDIF
     !
     IF ( ANY(driver(:)==4) ) THEN 
        CALL json%get('westpp_control.westpp_calculation', cval, found) 
        IF( found ) westpp_calculation = cval
        CALL json%get('westpp_control.westpp_range(1)', rval, found) 
        IF( found ) westpp_range(1) = rval  
        CALL json%get('westpp_control.westpp_range(2)', rval, found) 
        IF( found ) westpp_range(2) = rval  
        CALL json%get('westpp_control.westpp_format', cval, found) 
        IF( found ) westpp_format = cval
        CALL json%get('westpp_control.westpp_sign', lval, found) 
        IF( found ) westpp_sign = lval  
        CALL json%get('westpp_control.westpp_n_pdep_eigen_to_use', ival, found) 
        IF( found ) westpp_n_pdep_eigen_to_use = ival  
        CALL json%get('westpp_control.westpp_r0(1)', rval, found) 
        IF( found ) westpp_r0(1) = rval  
        CALL json%get('westpp_control.westpp_r0(2)', rval, found) 
        IF( found ) westpp_r0(2) = rval  
        CALL json%get('westpp_control.westpp_r0(3)', rval, found) 
        IF( found ) westpp_r0(3) = rval  
        CALL json%get('westpp_control.westpp_nr', ival, found) 
        IF( found ) westpp_nr = ival  
        CALL json%get('westpp_control.westpp_rmax', rval, found) 
        IF( found ) westpp_rmax = rval  
        CALL json%get('westpp_control.westpp_epsinfty', rval, found) 
        IF( found ) westpp_epsinfty = rval  
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
  ENDIF
  !
  IF ( ANY(driver(:)==2) ) THEN
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
     ! CHECKS 
     !
     SELECT CASE(wstat_calculation) 
     CASE('r','R','s','S')
     CASE DEFAULT
        CALL errore('fetch_input','Err: wstat_calculation /= S or R',1)
     END SELECT
     !
     IF( n_pdep_times < 2 ) CALL errore('fetch_input','Err: n_pdep_times<2',1) 
     IF( n_pdep_eigen < 1 ) CALL errore('fetch_input','Err: n_pdep_eigen<1',1)
     IF( n_pdep_eigen*n_pdep_times < nimage ) CALL errore('fetch_input','Err: n_pdep_eigen*n_pdep_times<nimage',1) 
     IF( n_pdep_maxiter < 1 ) CALL errore('fetch_input','Err: n_pdep_maxiter<1',1) 
     IF( n_dfpt_maxiter < 1 ) CALL errore('fetch_input','Err: n_dfpt_maxiter<1',1) 
     IF( n_pdep_read_from_file < 0 ) CALL errore('fetch_input','Err: n_pdep_read_from_file<0',1) 
     IF( n_pdep_read_from_file > n_pdep_eigen ) CALL errore('fetch_input','Err: n_pdep_read_from_file>n_pdep_eigen',1) 
     IF(tr2_dfpt<=0._DP) CALL errore('fetch_input','Err: tr2_dfpt<0.',1)
     IF(trev_pdep<=0._DP) CALL errore('fetch_input','Err: trev_pdep<0.',1)
     IF(trev_pdep_rel<=0._DP) CALL errore('fetch_input','Err: trev_pdep_rel<0.',1)
     !
  ENDIF
  !
  IF ( ANY(driver(:)==3) ) THEN
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
     ! CHECKS 
     !
     IF( n_lanczos < 2 ) CALL errore('fetch_input','Err: n_lanczos<2',1) 
     IF( n_pdep_eigen_to_use < 1 ) CALL errore('fetch_input','Err: n_pdep_eigen_to_use<1',1) 
     IF( n_pdep_eigen_to_use > n_pdep_eigen ) CALL errore('fetch_input','Err: n_pdep_eigen_to_use>n_pdep_eigen',1) 
     IF( n_imfreq < 1 ) CALL errore('fetch_input','Err: n_imfreq<1',1) 
     IF( n_refreq < 1 ) CALL errore('fetch_input','Err: n_refreq<1',1) 
     IF( n_spectralf < 2 ) CALL errore('fetch_input','Err: n_spectralf<1',1) 
     IF( qp_bandrange(1) < 1 ) CALL errore('fetch_input','Err: qp_bandrange(1)<1',1) 
     IF( qp_bandrange(2) < 1 ) CALL errore('fetch_input','Err: qp_bandrange(2)<1',1) 
     IF( qp_bandrange(2) < qp_bandrange(1) ) CALL errore('fetch_input','Err: qp_bandrange(2)<qp_bandrange(1)',1) 
     IF( ecut_imfreq<=0._DP) CALL errore('fetch_input','Err: ecut_imfreq<0.',1)
     IF( ecut_refreq<=0._DP) CALL errore('fetch_input','Err: ecut_imfreq<0.',1)
     IF( ecut_spectralf(2)<ecut_spectralf(1)) CALL errore('fetch_input','Err: ecut_spectralf(2)<ecut_spectralf(1)',1)
     IF( wfreq_eta<=0._DP) CALL errore('fetch_input','Err: wfreq_eta<0.',1)
     IF( n_secant_maxiter < 1 ) CALL errore('fetch_input','Err: n_secant_maxiter<1',1) 
     IF( trev_secant<=0._DP) CALL errore('fetch_input','Err: trev_secant<0.',1)
     SELECT CASE(macropol_calculation) 
     CASE('N','n','C','c')
     CASE DEFAULT
        CALL errore('fetch_input','Err: macropol_calculation /= (N,C)',1)
     END SELECT
     !
  ENDIF
  !
  IF ( ANY(driver(:)==4) ) THEN
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
     ! CHECKS 
     !
     IF( westpp_range(1) < 1 ) CALL errore('fetch_input','Err: westpp_range(1)<1',1)
     IF( westpp_range(2) < 1 ) CALL errore('fetch_input','Err: westpp_range(2)<1',1)
     IF( westpp_range(2) < westpp_range(1) ) CALL errore('fetch_input','Err: westpp_range(2)<westpp_range(1)',1)
     IF( westpp_nr < 1 ) CALL errore('fetch_input','Err: westpp_nr<1',1)
     IF( westpp_n_pdep_eigen_to_use < 1 ) CALL errore('fetch_input','Err: westpp_n_pdep_eigen_to_use<1',1)
     IF( westpp_rmax < 0.d0 ) CALL errore('fetch_input','Err: westpp_rmax<0',1)
     IF( westpp_epsinfty < 1.d0 ) CALL errore('fetch_input','Err: westpp_epsinfty<1',1)
     !
  ENDIF
  !
  ! REPORT
  !
  IF ( verbose ) THEN
     !
     CALL json%initialize()
     !
     CALL json%load_file(filename=TRIM(logfile))
     !
     IF ( ANY(driver(:)==1) ) THEN
        !
        ! REPORT
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
        CALL json%add('input.input_west.qe_prefix',TRIM(qe_prefix))
        CALL json%add('input.input_west.west_prefix',TRIM(west_prefix))
        CALL json%add('input.input_west.outdir',TRIM(outdir))
        !
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN
        !
        ! REPORT
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
        CALL json%add('input.wstat_control.wstat_calculation',TRIM(wstat_calculation))
        CALL json%add('input.wstat_control.n_pdep_eigen',n_pdep_eigen)
        CALL json%add('input.wstat_control.n_pdep_times',n_pdep_times)
        CALL json%add('input.wstat_control.n_pdep_maxiter',n_pdep_maxiter)
        CALL json%add('input.wstat_control.n_dfpt_maxiter',n_dfpt_maxiter)
        CALL json%add('input.wstat_control.n_pdep_read_from_file',n_pdep_read_from_file)
        CALL json%add('input.wstat_control.trev_pdep',trev_pdep)
        CALL json%add('input.wstat_control.trev_pdep_rel',trev_pdep_rel)
        CALL json%add('input.wstat_control.tr2_dfpt',tr2_dfpt)
        CALL json%add('input.wstat_control.l_kinetic_only',l_kinetic_only)
        CALL json%add('input.wstat_control.l_minimize_exx_if_active',l_minimize_exx_if_active)
        CALL json%add('input.wstat_control.l_use_ecutrho',l_use_ecutrho)
        !
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN
        !
        ! REPORT
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
        CALL json%add('input.wfreq_control.wfreq_calculation',TRIM(wfreq_calculation))
        CALL json%add('input.wfreq_control.n_pdep_eigen_to_use',n_pdep_eigen_to_use)
        CALL json%add('input.wfreq_control.qp_bandrange',qp_bandrange)
        CALL json%add('input.wfreq_control.macropol_calculation',macropol_calculation)
        CALL json%add('input.wfreq_control.n_lanczos',n_lanczos)
        CALL json%add('input.wfreq_control.n_imfreq',n_imfreq)
        CALL json%add('input.wfreq_control.n_refreq',n_refreq)
        CALL json%add('input.wfreq_control.ecut_imfreq',ecut_imfreq)
        CALL json%add('input.wfreq_control.ecut_refreq',ecut_refreq)
        CALL json%add('input.wfreq_control.wfreq_eta',wfreq_eta)
        CALL json%add('input.wfreq_control.n_secant_maxiter',n_secant_maxiter)
        CALL json%add('input.wfreq_control.trev_secant',trev_secant)
        CALL json%add('input.wfreq_control.l_enable_lanczos',l_enable_lanczos)
        CALL json%add('input.wfreq_control.l_enable_gwetot',l_enable_gwetot)
        CALL json%add('input.wfreq_control.div_kind_hf',div_kind_hf)
        CALL json%add('input.wfreq_control.o_restart_time',o_restart_time)
        CALL json%add('input.wfreq_control.ecut_spectralf',ecut_spectralf)
        CALL json%add('input.wfreq_control.n_spectralf',n_spectralf)
        !
     ENDIF
     !
     IF ( ANY(driver(:)==4) ) THEN
        !
        ! REPORT
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
        CALL json%add('input.westpp_control.westpp_calculation',TRIM(westpp_calculation))
        CALL json%add('input.westpp_control.westpp_range',westpp_range)
        CALL json%add('input.westpp_control.westpp_format',TRIM(westpp_format))
        CALL json%add('input.westpp_control.westpp_sign',westpp_sign)
        CALL json%add('input.westpp_control.westpp_n_pdep_eigen_to_use',westpp_n_pdep_eigen_to_use)
        CALL json%add('input.westpp_control.westpp_r0',westpp_r0)
        CALL json%add('input.westpp_control.westpp_nr',westpp_nr)
        CALL json%add('input.westpp_control.westpp_rmax',westpp_rmax)
        CALL json%add('input.westpp_control.westpp_epsinfty',westpp_epsinfty)
        !
     ENDIF
     !
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print_file( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
  CALL stop_clock('fetch_input')
  !
END SUBROUTINE
