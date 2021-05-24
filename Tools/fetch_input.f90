!
! Copyright (C) 2015-2019 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file License
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file: 
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE add_intput_parameters_to_json_file( num_drivers, driver, json )
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
        !
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN
        !
        CALL json%add('input.wstat_control.wstat_calculation',TRIM(wstat_calculation))
        CALL json%add('input.wstat_control.n_pdep_eigen',n_pdep_eigen)
        CALL json%add('input.wstat_control.n_pdep_times',n_pdep_times)
        CALL json%add('input.wstat_control.n_pdep_maxiter',n_pdep_maxiter)
        CALL json%add('input.wstat_control.n_dfpt_maxiter',n_dfpt_maxiter)
        CALL json%add('input.wstat_control.n_pdep_read_from_file',n_pdep_read_from_file)
        CALL json%add('input.wstat_control.n_steps_write_restart',n_steps_write_restart)
        CALL json%add('input.wstat_control.trev_pdep',trev_pdep)
        CALL json%add('input.wstat_control.trev_pdep_rel',trev_pdep_rel)
        CALL json%add('input.wstat_control.tr2_dfpt',tr2_dfpt)
        CALL json%add('input.wstat_control.l_kinetic_only',l_kinetic_only)
        CALL json%add('input.wstat_control.l_minimize_exx_if_active',l_minimize_exx_if_active)
        CALL json%add('input.wstat_control.l_use_ecutrho',l_use_ecutrho)
        CALL json%add('input.wstat_control.qlist', qlist) 
        !
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN
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
        CALL json%add('input.wfreq_control.o_restart_time',o_restart_time)
        CALL json%add('input.wfreq_control.ecut_spectralf',ecut_spectralf)
        CALL json%add('input.wfreq_control.n_spectralf',n_spectralf)
        !
     ENDIF
     !
     IF ( ANY(driver(:)==4) ) THEN
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
     IF ( ANY(driver(:)==5) ) THEN
        !
        CALL json%add('input.server_control.document',TRIM(document))
        !
     ENDIF
     !
  ENDIF
  !
END SUBROUTINE


SUBROUTINE fetch_input_yml( num_drivers, driver, verbose, debug )
  !
  USE west_version, ONLY : start_forpy, end_forpy
  USE io_push,      ONLY : io_push_title,io_push_value,io_push_bar,io_push_es0,io_push_c512 
  USE forpy_mod,    ONLY: call_py, call_py_noret, import_py, module_py
  USE forpy_mod,    ONLY: tuple, tuple_create 
  USE forpy_mod,    ONLY: dict, dict_create 
  USE forpy_mod,    ONLY: list, list_create 
  USE forpy_mod,    ONLY: object, cast
  USE forpy_mod,    ONLY : exception_matches, KeyError, err_clear, err_print 
  USE westcom
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_bcast, mp_barrier
  USE mp_world,         ONLY : mpime,root,world_comm
  USE mp_global,        ONLY : nimage
  USE gvect,            ONLY : ecutrho
  USE start_k,          ONLY : nk1, nk2, nk3
  USE control_flags,    ONLY : gamma_only
  USE json_module,      ONLY : json_file
  USE pwcom,            ONLY : nelec
  !
  IMPLICIT NONE 
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: num_drivers
  INTEGER, INTENT(IN) :: driver(num_drivers)
  LOGICAL, INTENT(IN) :: verbose
  LOGICAL, INTENT(IN) :: debug
  !
  INTEGER :: IERR
  TYPE(tuple) :: args
  TYPE(dict) :: kwargs
  TYPE(module_py) :: pymod
  TYPE(object) :: return_obj, tmp_obj
  TYPE(dict) :: return_dict
  TYPE(list) :: tmp_list
  INTEGER :: list_len
  INTEGER :: i
  INTEGER :: nq
  INTEGER :: numsp 
  CHARACTER(LEN=512), EXTERNAL :: trimcheck
  CHARACTER(LEN=:),ALLOCATABLE :: cvalue
  TYPE(json_file) :: json
  INTEGER :: iunit, lenc
  INTEGER, PARAMETER :: DUMMY_DEFAULT = -1210
  !
  CALL start_clock('fetch_input')
  !
  IF ( mpime==root ) THEN 
     ! 
     IERR = import_py(pymod, "west_fetch_input")
     !
     IF ( ANY(driver(:)==1) ) THEN 
        !  
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)) )
        IERR = args%setitem(1, "input_west" )
        IERR = args%setitem(2, verbose )
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, "read_keyword_from_file", args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, "qe_prefix"); qe_prefix = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, "west_prefix"); west_prefix = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, "outdir"); outdir = TRIM(ADJUSTL(cvalue))
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN
        !
        IF ( gamma_only ) THEN
           nq = 1 
        ELSE
           nq = nk1*nk2*nk3
        ENDIF
        !  
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)) )
        IERR = args%setitem(1, "wstat_control" )
        IERR = args%setitem(2, verbose )
        IERR = dict_create(kwargs)
        IERR = kwargs%setitem("nq",nq)
        IERR = kwargs%setitem("nelec",nelec)
        !
        IERR = call_py(return_obj, pymod, "read_keyword_from_file", args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, "wstat_calculation"); wstat_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(n_pdep_eigen, "n_pdep_eigen", DUMMY_DEFAULT)
        IERR = return_dict%get(n_pdep_times, "n_pdep_times", DUMMY_DEFAULT)
        IERR = return_dict%get(n_pdep_maxiter, "n_pdep_maxiter", DUMMY_DEFAULT)
        IERR = return_dict%get(n_dfpt_maxiter, "n_dfpt_maxiter", DUMMY_DEFAULT)
        IERR = return_dict%get(n_pdep_read_from_file, "n_pdep_read_from_file", DUMMY_DEFAULT)
        IERR = return_dict%get(n_steps_write_restart, "n_steps_write_restart", DUMMY_DEFAULT)
        IERR = return_dict%getitem(trev_pdep, "trev_pdep")
        IERR = return_dict%getitem(trev_pdep_rel, "trev_pdep_rel")
        IERR = return_dict%getitem(tr2_dfpt, "tr2_dfpt")
        IERR = return_dict%getitem(l_kinetic_only, "l_kinetic_only")
        IERR = return_dict%getitem(l_minimize_exx_if_active, "l_minimize_exx_if_active")
        IERR = return_dict%getitem(l_use_ecutrho, "l_use_ecutrho")
        IERR = return_dict%getitem(tmp_obj, "qlist")
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%len(list_len)
        IF( ALLOCATED(qlist) ) DEALLOCATE(qlist)
        ALLOCATE(qlist(list_len))
        DO i = 0, list_len-1 ! Python indices start at 0
           IERR = tmp_list%getitem(qlist(i+1), i) ! Fortran indices start at 1 
        ENDDO
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN
        !  
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)) )
        IERR = args%setitem(1, "wfreq_control" )
        IERR = args%setitem(2, verbose )
        IERR = dict_create(kwargs)
        IERR = kwargs%setitem("nelec",nelec)
        IERR = kwargs%setitem("ecutrho",ecutrho)
        !
        IERR = call_py(return_obj, pymod, "read_keyword_from_file", args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, "wfreq_calculation"); wfreq_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(n_pdep_eigen_to_use, "n_pdep_eigen_to_use", DUMMY_DEFAULT)
        IERR = return_dict%getitem(tmp_obj, "qp_bandrange")
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%len(list_len)
        IERR = tmp_list%getitem(qp_bandrange(1), 0) ! Fortran indices start at 1 
        IERR = tmp_list%getitem(qp_bandrange(2), 1) ! Fortran indices start at 1 
        CALL tmp_list%destroy 
        CALL tmp_obj%destroy
        IERR = return_dict%getitem(cvalue, "macropol_calculation"); macropol_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(n_lanczos, "n_lanczos", DUMMY_DEFAULT)
        IERR = return_dict%get(n_imfreq, "n_imfreq", DUMMY_DEFAULT)
        IERR = return_dict%get(n_refreq, "n_refreq", DUMMY_DEFAULT)
        IERR = return_dict%getitem(ecut_imfreq, "ecut_imfreq")
        IERR = return_dict%getitem(ecut_refreq, "ecut_refreq")
        IERR = return_dict%getitem(wfreq_eta, "wfreq_eta")
        IERR = return_dict%get(n_secant_maxiter, "n_secant_maxiter", DUMMY_DEFAULT)
        IERR = return_dict%getitem(trev_secant, "trev_secant")
        IERR = return_dict%getitem(l_enable_lanczos, "l_enable_lanczos")
        IERR = return_dict%getitem(l_enable_gwetot, "l_enable_gwetot")
        IERR = return_dict%getitem(o_restart_time, "o_restart_time")
        IERR = return_dict%getitem(tmp_obj, "ecut_spectralf")
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%getitem(ecut_spectralf(1), 0) ! Fortran indices start at 1 
        IERR = tmp_list%getitem(ecut_spectralf(2), 1) ! Fortran indices start at 1 
        CALL tmp_list%destroy 
        CALL tmp_obj%destroy 
        IERR = return_dict%get(n_spectralf, "n_spectralf", DUMMY_DEFAULT)
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF ( ANY(driver(:)==4) ) THEN
        !  
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)) )
        IERR = args%setitem(1, "westpp_control" )
        IERR = args%setitem(2, verbose )
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, "read_keyword_from_file", args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, "westpp_calculation"); westpp_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(tmp_obj, "westpp_range")
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%getitem(westpp_range(1), 0) ! Fortran indices start at 1 
        IERR = tmp_list%getitem(westpp_range(2), 1) ! Fortran indices start at 1 
        CALL tmp_list%destroy 
        CALL tmp_obj%destroy 
        IERR = return_dict%getitem(cvalue, "westpp_format"); westpp_format = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(westpp_sign, "westpp_sign")
        IERR = return_dict%get(westpp_n_pdep_eigen_to_use, "westpp_n_pdep_eigen_to_use", DUMMY_DEFAULT)
        IERR = return_dict%getitem(tmp_obj, "westpp_r0")
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%getitem(westpp_r0(1), 0) ! Fortran indices start at 1 
        IERR = tmp_list%getitem(westpp_r0(2), 1) ! Fortran indices start at 1 
        IERR = tmp_list%getitem(westpp_r0(3), 2) ! Fortran indices start at 1 
        CALL tmp_list%destroy 
        CALL tmp_obj%destroy 
        IERR = return_dict%get(westpp_nr, "westpp_nr", DUMMY_DEFAULT)
        IERR = return_dict%getitem(westpp_rmax, "westpp_rmax")
        IERR = return_dict%getitem(westpp_epsinfty, "westpp_epsinfty")
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF ( ANY(driver(:)==5) ) THEN
        !  
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)) )
        IERR = args%setitem(1, "server_control" )
        IERR = args%setitem(2, verbose )
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, "read_keyword_from_file", args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, "document"); document = TRIM(ADJUSTL(cvalue))
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     CALL pymod%destroy
     !
     !
  ENDIF
  !
  ! BCAST & CHECKS
  !
  IF ( ANY(driver(:)==1) ) THEN
     !
     CALL mp_bcast(qe_prefix,root,world_comm); prefix=qe_prefix
     CALL mp_bcast(west_prefix,root,world_comm)
     CALL mp_bcast(outdir,root,world_comm); tmp_dir = trimcheck (outdir)
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
     CALL mp_bcast(n_steps_write_restart,root,world_comm)
     CALL mp_bcast(trev_pdep,root,world_comm)
     CALL mp_bcast(trev_pdep_rel,root,world_comm)
     CALL mp_bcast(tr2_dfpt,root,world_comm)
     CALL mp_bcast(l_kinetic_only,root,world_comm)
     CALL mp_bcast(l_minimize_exx_if_active,root,world_comm)
     CALL mp_bcast(l_use_ecutrho,root,world_comm)
     IF(mpime == root) nq = SIZE(qlist)
     CALL mp_bcast(nq,root,world_comm)
     IF(mpime /= root) THEN 
        IF( ALLOCATED(qlist) ) DEALLOCATE(qlist)
        ALLOCATE(qlist(nq))
     ENDIF  
     CALL mp_bcast(qlist,root,world_comm)
     !
     ! CHECKS 
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
     IF( n_pdep_eigen == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_pdep_eigen')
     IF( n_pdep_times == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_pdep_times')
     IF( n_pdep_maxiter == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_pdep_maxiter')
     IF( n_dfpt_maxiter == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_dfpt_maxiter')
     IF( n_pdep_read_from_file == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_pdep_read_from_file')
     IF( n_steps_write_restart == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_steps_write_restart')
     IF(gamma_only) THEN
        IF (SIZE(qlist)/=1) CALL errore('fetch_input','Err: SIZE(qlist)/=1.',1)
     ELSE 
        IF (SIZE(qlist)>nk1*nk2*nk3) CALL errore('fetch_input','Err: SIZE(qlist)>nk1*nk2*nk3.',1)
     ENDIF
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
     IF( n_secant_maxiter < 0 ) CALL errore('fetch_input','Err: n_secant_maxiter<0',1) 
     IF( trev_secant<=0._DP) CALL errore('fetch_input','Err: trev_secant<0.',1)
     IF( n_pdep_eigen_to_use == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_pdep_eigen_to_use')
     IF( n_lanczos == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_lanczos')
     IF( n_imfreq == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_imfreq')
     IF( n_refreq == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_refreq')
     IF( n_secant_maxiter == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_secant_maxiter')
     IF( n_spectralf == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read n_spectralf')
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
     IF( westpp_rmax < 0._DP ) CALL errore('fetch_input','Err: westpp_rmax<0.',1)
     IF( westpp_epsinfty < 1._DP ) CALL errore('fetch_input','Err: westpp_epsinfty<1.',1)
     IF( westpp_n_pdep_eigen_to_use == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read westpp_n_pdep_eigen_to_use')
     IF( westpp_nr == DUMMY_DEFAULT ) CALL errore('fetch_input','Err: cannot read westpp_nr')
     !
  ENDIF
  !
  IF ( ANY(driver(:)==5) ) THEN
     !
     lenc = LEN(document)
     CALL mp_bcast(lenc,root,world_comm)
     IF(mpime/=root) ALLOCATE(CHARACTER(LEN=lenc) :: document)
     CALL mp_bcast(document,root,world_comm)
     !
  ENDIF
  !
  CALL mp_barrier(world_comm)
  !
  ! REPORT
  !
  IF ( debug .AND. mpime == root ) THEN
     !
     IF ( ANY(driver(:)==1) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title("I/O Summary : input_west")
        !
        numsp = 14
        CALL io_push_c512('qe_prefix',qe_prefix,numsp)
        CALL io_push_c512('west_prefix',west_prefix,numsp)
        CALL io_push_c512('outdir',outdir,numsp)
        !
        CALL io_push_bar()
        !
     ENDIF
     !
     IF ( ANY(driver(:)==2) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title('I/O Summary : wstat_control')
        !
        numsp = 30
        CALL io_push_value('wstat_calculation',wstat_calculation,numsp)
        CALL io_push_value('n_pdep_eigen',n_pdep_eigen,numsp)
        CALL io_push_value('n_pdep_times',n_pdep_times,numsp)
        CALL io_push_value('n_pdep_maxiter',n_pdep_maxiter,numsp)
        CALL io_push_value('n_dfpt_maxiter',n_dfpt_maxiter,numsp)
        CALL io_push_value('n_pdep_read_from_file',n_pdep_read_from_file,numsp)
        CALL io_push_value('n_steps_write_restart',n_steps_write_restart,numsp)
        CALL io_push_es0('trev_pdep',trev_pdep,numsp)
        CALL io_push_es0('trev_pdep_rel',trev_pdep_rel,numsp)
        CALL io_push_es0('tr2_dfpt',tr2_dfpt,numsp)
        CALL io_push_value('l_kinetic_only',l_kinetic_only,numsp)
        CALL io_push_value('l_minimize_exx_if_active',l_minimize_exx_if_active,numsp)
        CALL io_push_value('l_use_ecutrho',l_use_ecutrho,numsp)
        DO i = 1, SIZE(qlist) 
           CALL io_push_value('qlist',qlist(i),numsp)
        ENDDO
        !
        CALL io_push_bar()
        !
     ENDIF
     !
     IF ( ANY(driver(:)==3) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title('I/O Summary : wfreq_control')
        !
        numsp = 40
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
        CALL io_push_value('o_restart_time [min]',o_restart_time,numsp)
        CALL io_push_value('ecut_spectralf(1) [Ry]',ecut_spectralf(1),numsp)
        CALL io_push_value('ecut_spectralf(2) [Ry]',ecut_spectralf(2),numsp)
        CALL io_push_value('n_spectralf',n_spectralf,numsp)
        !
        CALL io_push_bar()
        !
     ENDIF
     !
     IF ( ANY(driver(:)==4) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title('I/O Summary : westpp_control')
        !
        numsp = 40
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
     ENDIF
     !
     IF ( ANY(driver(:)==5) ) THEN
        !
        ! REPORT
        !
        CALL io_push_title('I/O Summary : server_control')
        !
        numsp = 40
        CALL io_push_value('document',document,numsp)
        !
        CALL io_push_bar()
        !
     ENDIF
     !
  ENDIF
  !
  IF ( verbose .AND. mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     !
     CALL add_intput_parameters_to_json_file( num_drivers, driver, json )
     ! 
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     CALL json%destroy()
     !
     !
  ENDIF
  !
  CALL stop_clock('fetch_input')
  !
END SUBROUTINE
