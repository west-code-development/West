!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE add_intput_parameters_to_json_file(num_drivers, driver, json)
  !-----------------------------------------------------------------------
  !
  USE json_module,      ONLY : json_file
  USE westcom,          ONLY : qe_prefix,west_prefix,outdir,wstat_calculation,n_pdep_eigen,&
                             & n_pdep_times,n_pdep_maxiter,n_dfpt_maxiter,n_pdep_read_from_file,&
                             & n_steps_write_restart,trev_pdep,trev_pdep_rel,tr2_dfpt,&
                             & l_kinetic_only,l_minimize_exx_if_active,l_use_ecutrho,qlist,&
                             & wfreq_calculation,n_pdep_eigen_to_use,qp_bandrange,qp_bands,&
                             & macropol_calculation,n_lanczos,n_imfreq,n_refreq,ecut_imfreq,&
                             & ecut_refreq,wfreq_eta,n_secant_maxiter,trev_secant,l_enable_lanczos,&
                             & l_qdet_verbose,l_enable_off_diagonal,o_restart_time,ecut_spectralf,&
                             & n_spectralf,westpp_calculation,westpp_range,westpp_format,&
                             & westpp_sign,westpp_n_pdep_eigen_to_use,westpp_r0,westpp_nr,&
                             & westpp_rmax,westpp_epsinfty,westpp_box,document,&
                             & wbse_init_calculation,localization,wfc_from_qbox,bisection_info,&
                             & chi_kernel,overlap_thr,spin_channel,wbse_calculation,solver,&
                             & qp_correction,scissor_ope,n_liouville_eigen,n_liouville_times,&
                             & n_liouville_maxiter,n_liouville_read_from_file,trev_liouville,&
                             & trev_liouville_rel,ipol_input,wbse_epsinfty,spin_excitation,&
                             & l_preconditioning,l_reduce_io,wbsepp_calculation,wbsepp_r0,&
                             & iexc_plot,n_lanczos_to_use,n_extrapolation,range,broad
  USE mp_world,         ONLY : mpime,root
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: num_drivers
  INTEGER, INTENT(IN) :: driver(num_drivers)
  TYPE(json_file), INTENT(INOUT) :: json
  !
  IF(mpime == root) THEN
     !
     IF(ANY(driver(:)==1)) THEN
        !
        CALL json%add('input.input_west.qe_prefix',TRIM(qe_prefix))
        CALL json%add('input.input_west.west_prefix',TRIM(west_prefix))
        CALL json%add('input.input_west.outdir',TRIM(outdir))
        !
     ENDIF
     !
     IF(ANY(driver(:)==2)) THEN
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
        CALL json%add('input.wstat_control.qlist',qlist)
        !
     ENDIF
     !
     IF(ANY(driver(:)==3)) THEN
        !
        CALL json%add('input.wfreq_control.wfreq_calculation',TRIM(wfreq_calculation))
        CALL json%add('input.wfreq_control.n_pdep_eigen_to_use',n_pdep_eigen_to_use)
        CALL json%add('input.wfreq_control.qp_bandrange',qp_bandrange)
        CALL json%add('input.wfreq_control.qp_bands',qp_bands)
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
        CALL json%add('input.wfreq_control.l_enable_off_diagonal',l_enable_off_diagonal)
        CALL json%add('input.wfreq_control.o_restart_time',o_restart_time)
        CALL json%add('input.wfreq_control.ecut_spectralf',ecut_spectralf)
        CALL json%add('input.wfreq_control.n_spectralf',n_spectralf)
        !
     ENDIF
     !
     IF(ANY(driver(:)==4)) THEN
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
        CALL json%add('input.westpp_control.westpp_box',westpp_box)
        !
     ENDIF
     !
     IF(ANY(driver(:)==5)) THEN
        !
        CALL json%add('input.server_control.document',TRIM(document))
        !
     ENDIF
     !
     IF(ANY(driver(:)==6)) THEN
        !
        CALL json%add('input.wbse_init_control.wbse_init_calculation',TRIM(wbse_init_calculation))
        CALL json%add('input.wbse_init_control.localization',TRIM(localization))
        CALL json%add('input.wbse_init_control.wfc_from_qbox',TRIM(wfc_from_qbox))
        CALL json%add('input.wbse_init_control.bisection_info',TRIM(bisection_info))
        CALL json%add('input.wbse_init_control.chi_kernel',TRIM(chi_kernel))
        CALL json%add('input.wbse_init_control.overlap_thr',overlap_thr)
        CALL json%add('input.wbse_init_control.spin_channel',spin_channel)
        !
     ENDIF
     !
     IF(ANY(driver(:)==7)) THEN
        !
        CALL json%add('input.wbse_control.wbse_calculation',TRIM(wbse_calculation))
        CALL json%add('input.wbse_control.solver',TRIM(solver))
        CALL json%add('input.wbse_control.qp_correction',TRIM(qp_correction))
        CALL json%add('input.wbse_control.scissor_ope',scissor_ope)
        CALL json%add('input.wbse_control.n_liouville_eigen',n_liouville_eigen)
        CALL json%add('input.wbse_control.n_liouville_times',n_liouville_times)
        CALL json%add('input.wbse_control.n_liouville_maxiter',n_liouville_maxiter)
        CALL json%add('input.wbse_control.n_liouville_read_from_file',n_liouville_read_from_file)
        CALL json%add('input.wbse_control.trev_liouville',trev_liouville)
        CALL json%add('input.wbse_control.trev_liouville_rel',trev_liouville_rel)
        CALL json%add('input.wbse_control.n_lanczos',n_lanczos)
        CALL json%add('input.wbse_control.n_steps_write_restart',n_steps_write_restart)
        CALL json%add('input.wbse_control.ipol_input',TRIM(ipol_input))
        CALL json%add('input.wbse_control.macropol_calculation',TRIM(macropol_calculation))
        CALL json%add('input.wbse_control.wbse_epsinfty',wbse_epsinfty)
        CALL json%add('input.wbse_control.spin_excitation',TRIM(spin_excitation))
        CALL json%add('input.wbse_control.l_preconditioning',l_preconditioning)
        CALL json%add('input.wbse_control.l_reduce_io',l_reduce_io)
        !
     ENDIF
     !
     IF(ANY(driver(:)==8)) THEN
        !
        CALL json%add('input.wbsepp_control.wbsepp_calculation',TRIM(wbsepp_calculation))
        CALL json%add('input.wbsepp_control.n_liouville_read_from_file',n_liouville_read_from_file)
        CALL json%add('input.wbsepp_control.ipol_input',TRIM(ipol_input))
        CALL json%add('input.wbsepp_control.spin_channel',spin_channel)
        CALL json%add('input.wbsepp_control.wbsepp_r0',wbsepp_r0)
        CALL json%add('input.wbsepp_control.iexc_plot',iexc_plot)
        CALL json%add('input.wbsepp_control.n_lanczos_to_use',n_lanczos_to_use)
        CALL json%add('input.wbsepp_control.n_extrapolation',n_extrapolation)
        CALL json%add('input.wbsepp_control.range',range)
        CALL json%add('input.wbsepp_control.broad',broad)
        !
     ENDIF
     !
  ENDIF
  !
END SUBROUTINE
!
SUBROUTINE fetch_input_yml(num_drivers, driver, verbose)
  !
  USE forpy_mod,        ONLY : call_py,import_py,module_py,tuple,tuple_create,dict,dict_create,&
                             & list,object,cast
  USE westcom,          ONLY : qe_prefix,west_prefix,outdir,wstat_calculation,n_pdep_eigen,&
                             & n_pdep_times,n_pdep_maxiter,n_dfpt_maxiter,n_pdep_read_from_file,&
                             & n_steps_write_restart,trev_pdep,trev_pdep_rel,tr2_dfpt,&
                             & l_kinetic_only,l_minimize_exx_if_active,l_use_ecutrho,qlist,&
                             & wfreq_calculation,n_pdep_eigen_to_use,qp_bandrange,qp_bands,&
                             & macropol_calculation,n_lanczos,n_imfreq,n_refreq,ecut_imfreq,&
                             & ecut_refreq,wfreq_eta,n_secant_maxiter,trev_secant,l_enable_lanczos,&
                             & l_qdet_verbose,l_enable_off_diagonal,o_restart_time,ecut_spectralf,&
                             & n_spectralf,westpp_calculation,westpp_range,westpp_format,&
                             & westpp_sign,westpp_n_pdep_eigen_to_use,westpp_r0,westpp_nr,&
                             & westpp_rmax,westpp_epsinfty,westpp_box,document,&
                             & wbse_init_calculation,localization,wfc_from_qbox,bisection_info,&
                             & chi_kernel,overlap_thr,spin_channel,wbse_calculation,solver,&
                             & qp_correction,scissor_ope,n_liouville_eigen,n_liouville_times,&
                             & n_liouville_maxiter,n_liouville_read_from_file,trev_liouville,&
                             & trev_liouville_rel,ipol_input,wbse_epsinfty,spin_excitation,&
                             & l_preconditioning,l_reduce_io,wbsepp_calculation,wbsepp_r0,&
                             & iexc_plot,n_lanczos_to_use,n_extrapolation,range,broad,&
                             & main_input_file,logfile
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : tmp_dir,prefix
  USE mp,               ONLY : mp_bcast,mp_barrier
  USE mp_world,         ONLY : mpime,root,world_comm
  USE mp_global,        ONLY : nimage
  USE gvect,            ONLY : ecutrho
  USE start_k,          ONLY : nk1,nk2,nk3
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
  !
  INTEGER :: IERR
  TYPE(tuple) :: args
  TYPE(dict) :: kwargs
  TYPE(module_py) :: pymod
  TYPE(object) :: return_obj,tmp_obj
  TYPE(dict) :: return_dict
  TYPE(list) :: tmp_list
  INTEGER :: list_len
  INTEGER :: i
  INTEGER :: nq
  INTEGER :: n_qp_bands
  CHARACTER(LEN=512), EXTERNAL :: trimcheck
  CHARACTER(LEN=:),ALLOCATABLE :: cvalue
  TYPE(json_file) :: json
  INTEGER :: iunit, lenc
  INTEGER, PARAMETER :: DUMMY_DEFAULT = -1210
  !
  CALL start_clock('fetch_input')
  !
  IF(mpime == root) THEN
     !
     IERR = import_py(pymod, 'west_fetch_input')
     !
     IF(ANY(driver(:)==1)) THEN
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'input_west')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'qe_prefix'); qe_prefix = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'west_prefix'); west_prefix = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'outdir'); outdir = TRIM(ADJUSTL(cvalue))
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF(ANY(driver(:)==2)) THEN
        !
        IF(gamma_only) THEN
           nq = 1
        ELSE
           nq = nk1*nk2*nk3
        ENDIF
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'wstat_control')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        IERR = kwargs%setitem('nq',nq)
        IERR = kwargs%setitem('nelec',nelec)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'wstat_calculation'); wstat_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(n_pdep_eigen, 'n_pdep_eigen', DUMMY_DEFAULT)
        IERR = return_dict%get(n_pdep_times, 'n_pdep_times', DUMMY_DEFAULT)
        IERR = return_dict%get(n_pdep_maxiter, 'n_pdep_maxiter', DUMMY_DEFAULT)
        IERR = return_dict%get(n_dfpt_maxiter, 'n_dfpt_maxiter', DUMMY_DEFAULT)
        IERR = return_dict%get(n_pdep_read_from_file, 'n_pdep_read_from_file', DUMMY_DEFAULT)
        IERR = return_dict%get(n_steps_write_restart, 'n_steps_write_restart', DUMMY_DEFAULT)
        IERR = return_dict%getitem(trev_pdep, 'trev_pdep')
        IERR = return_dict%getitem(trev_pdep_rel, 'trev_pdep_rel')
        IERR = return_dict%getitem(tr2_dfpt, 'tr2_dfpt')
        IERR = return_dict%getitem(l_kinetic_only, 'l_kinetic_only')
        IERR = return_dict%getitem(l_minimize_exx_if_active, 'l_minimize_exx_if_active')
        IERR = return_dict%getitem(l_use_ecutrho, 'l_use_ecutrho')
        IERR = return_dict%getitem(tmp_obj, 'qlist')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%len(list_len)
        IF(ALLOCATED(qlist)) DEALLOCATE(qlist)
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
     IF(ANY(driver(:)==3)) THEN
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'wfreq_control')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        IERR = kwargs%setitem('nelec',nelec)
        IERR = kwargs%setitem('ecutrho',ecutrho)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'wfreq_calculation'); wfreq_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(n_pdep_eigen_to_use, 'n_pdep_eigen_to_use', DUMMY_DEFAULT)
        IERR = return_dict%getitem(tmp_obj, 'qp_bandrange')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%len(list_len)
        IERR = tmp_list%getitem(qp_bandrange(1), 0) ! Fortran indices start at 1
        IERR = tmp_list%getitem(qp_bandrange(2), 1) ! Fortran indices start at 1
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        IERR = return_dict%getitem(tmp_obj, 'qp_bands')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%len(list_len)
        IF(ALLOCATED(qp_bands)) DEALLOCATE(qp_bands)
        ALLOCATE(qp_bands(list_len))
        DO i = 0, list_len-1 ! Python indices start at 0
           IERR = tmp_list%getitem(qp_bands(i+1), i) ! Fortran indices start at 1
        ENDDO
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        IERR = return_dict%getitem(cvalue, 'macropol_calculation'); macropol_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(n_lanczos, 'n_lanczos', DUMMY_DEFAULT)
        IERR = return_dict%get(n_imfreq, 'n_imfreq', DUMMY_DEFAULT)
        IERR = return_dict%get(n_refreq, 'n_refreq', DUMMY_DEFAULT)
        IERR = return_dict%getitem(ecut_imfreq, 'ecut_imfreq')
        IERR = return_dict%getitem(ecut_refreq, 'ecut_refreq')
        IERR = return_dict%getitem(wfreq_eta, 'wfreq_eta')
        IERR = return_dict%get(n_secant_maxiter, 'n_secant_maxiter', DUMMY_DEFAULT)
        IERR = return_dict%getitem(trev_secant, 'trev_secant')
        IERR = return_dict%getitem(l_enable_lanczos, 'l_enable_lanczos')
        IERR = return_dict%getitem(l_qdet_verbose, 'l_qdet_verbose')
        IERR = return_dict%getitem(l_enable_off_diagonal, 'l_enable_off_diagonal')
        IERR = return_dict%getitem(o_restart_time, 'o_restart_time')
        IERR = return_dict%getitem(tmp_obj, 'ecut_spectralf')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%getitem(ecut_spectralf(1), 0) ! Fortran indices start at 1
        IERR = tmp_list%getitem(ecut_spectralf(2), 1) ! Fortran indices start at 1
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        IERR = return_dict%get(n_spectralf, 'n_spectralf', DUMMY_DEFAULT)
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF(ANY(driver(:)==4)) THEN
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'westpp_control')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'westpp_calculation'); westpp_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(tmp_obj, 'westpp_range')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%getitem(westpp_range(1), 0) ! Fortran indices start at 1
        IERR = tmp_list%getitem(westpp_range(2), 1) ! Fortran indices start at 1
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        IERR = return_dict%getitem(cvalue, 'westpp_format'); westpp_format = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(westpp_sign, 'westpp_sign')
        IERR = return_dict%get(westpp_n_pdep_eigen_to_use, 'westpp_n_pdep_eigen_to_use', DUMMY_DEFAULT)
        IERR = return_dict%getitem(tmp_obj, 'westpp_r0')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%getitem(westpp_r0(1), 0) ! Fortran indices start at 1
        IERR = tmp_list%getitem(westpp_r0(2), 1) ! Fortran indices start at 1
        IERR = tmp_list%getitem(westpp_r0(3), 2) ! Fortran indices start at 1
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        IERR = return_dict%get(westpp_nr, 'westpp_nr', DUMMY_DEFAULT)
        IERR = return_dict%getitem(westpp_rmax, 'westpp_rmax')
        IERR = return_dict%getitem(westpp_epsinfty, 'westpp_epsinfty')
        IERR = return_dict%getitem(tmp_obj, 'westpp_box')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%getitem(westpp_box(1), 0)
        IERR = tmp_list%getitem(westpp_box(2), 1)
        IERR = tmp_list%getitem(westpp_box(3), 2)
        IERR = tmp_list%getitem(westpp_box(4), 3)
        IERR = tmp_list%getitem(westpp_box(5), 4)
        IERR = tmp_list%getitem(westpp_box(6), 5)
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF(ANY(driver(:)==5)) THEN
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'server_control')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'document'); document = TRIM(ADJUSTL(cvalue))
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF(ANY(driver(:)==6)) THEN
        !
        IF(ALLOCATED(qlist)) DEALLOCATE(qlist)
        ALLOCATE(qlist(1))
        qlist = (/1/)
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'wbse_init_control')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'wbse_init_calculation'); wbse_init_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'localization'); localization = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'wfc_from_qbox'); wfc_from_qbox = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'bisection_info'); bisection_info = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'chi_kernel'); chi_kernel = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(overlap_thr, 'overlap_thr')
        IERR = return_dict%get(spin_channel, 'spin_channel', DUMMY_DEFAULT)
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF(ANY(driver(:)==7)) THEN
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'wbse_control')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'wbse_calculation'); wbse_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'solver'); solver = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'qp_correction'); qp_correction = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(scissor_ope, 'scissor_ope')
        IERR = return_dict%get(n_liouville_eigen, 'n_liouville_eigen', DUMMY_DEFAULT)
        IERR = return_dict%get(n_liouville_times, 'n_liouville_times', DUMMY_DEFAULT)
        IERR = return_dict%get(n_liouville_maxiter, 'n_liouville_maxiter', DUMMY_DEFAULT)
        IERR = return_dict%get(n_liouville_read_from_file, 'n_liouville_read_from_file', DUMMY_DEFAULT)
        IERR = return_dict%getitem(trev_liouville, 'trev_liouville')
        IERR = return_dict%getitem(trev_liouville_rel, 'trev_liouville_rel')
        IERR = return_dict%get(n_lanczos,'n_lanczos', DUMMY_DEFAULT)
        IERR = return_dict%get(n_steps_write_restart, 'n_steps_write_restart', DUMMY_DEFAULT)
        IERR = return_dict%getitem(cvalue, 'ipol_input'); ipol_input = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(cvalue, 'macropol_calculation'); macropol_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(wbse_epsinfty, 'wbse_epsinfty')
        IERR = return_dict%getitem(cvalue, 'spin_excitation'); spin_excitation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%getitem(l_preconditioning, 'l_preconditioning')
        IERR = return_dict%getitem(l_reduce_io, 'l_reduce_io')
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     IF(ANY(driver(:)==8)) THEN
        !
        IF(ALLOCATED(qlist)) DEALLOCATE(qlist)
        ALLOCATE(qlist(1))
        qlist = (/1/)
        !
        IERR = tuple_create(args, 3)
        IERR = args%setitem(0, TRIM(ADJUSTL(main_input_file)))
        IERR = args%setitem(1, 'wbsepp_control')
        IERR = args%setitem(2, verbose)
        IERR = dict_create(kwargs)
        !
        IERR = call_py(return_obj, pymod, 'read_keyword_from_file', args, kwargs)
        IERR = cast(return_dict, return_obj)
        !
        CALL args%destroy
        CALL kwargs%destroy
        CALL return_obj%destroy
        !
        IERR = return_dict%getitem(cvalue, 'wbsepp_calculation'); wbsepp_calculation = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(n_liouville_read_from_file, 'n_liouville_read_from_file', DUMMY_DEFAULT)
        IERR = return_dict%getitem(cvalue, 'ipol_input'); ipol_input = TRIM(ADJUSTL(cvalue))
        IERR = return_dict%get(spin_channel, 'spin_channel', DUMMY_DEFAULT)
        IERR = return_dict%getitem(tmp_obj, 'wbsepp_r0')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%len(list_len)
        IERR = tmp_list%getitem(wbsepp_r0(1), 0) ! Fortran indices start at 1
        IERR = tmp_list%getitem(wbsepp_r0(2), 1) ! Fortran indices start at 1
        IERR = tmp_list%getitem(wbsepp_r0(3), 2) ! Fortran indices start at 1
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        IERR = return_dict%get(iexc_plot, 'iexc_plot', DUMMY_DEFAULT)
        IERR = return_dict%get(n_lanczos_to_use, 'n_lanczos_to_use', DUMMY_DEFAULT)
        IERR = return_dict%get(n_extrapolation, 'n_extrapolation', DUMMY_DEFAULT)
        IERR = return_dict%getitem(tmp_obj, 'range')
        IERR = cast(tmp_list,tmp_obj)
        IERR = tmp_list%len(list_len)
        IERR = tmp_list%getitem(range(1), 0) ! Fortran indices start at 1
        IERR = tmp_list%getitem(range(2), 1) ! Fortran indices start at 1
        IERR = tmp_list%getitem(range(3), 2) ! Fortran indices start at 1
        CALL tmp_list%destroy
        CALL tmp_obj%destroy
        IERR = return_dict%getitem(broad, 'broad')
        !
        CALL return_dict%destroy
        !
     ENDIF
     !
     CALL pymod%destroy
     !
  ENDIF
  !
  ! BCAST & CHECKS
  !
  IF(ANY(driver(:)==1)) THEN
     !
     CALL mp_bcast(qe_prefix,root,world_comm); prefix = qe_prefix
     CALL mp_bcast(west_prefix,root,world_comm)
     CALL mp_bcast(outdir,root,world_comm); tmp_dir = trimcheck(outdir)
     !
  ENDIF
  !
  IF(ANY(driver(:)==2)) THEN
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
        IF(ALLOCATED(qlist)) DEALLOCATE(qlist)
        ALLOCATE(qlist(nq))
     ENDIF
     CALL mp_bcast(qlist,root,world_comm)
     !
     ! CHECKS
     !
     IF(n_pdep_times < 2) CALL errore('fetch_input','Err: n_pdep_times<2',1)
     IF(n_pdep_eigen < 1) CALL errore('fetch_input','Err: n_pdep_eigen<1',1)
     IF(n_pdep_eigen*n_pdep_times < nimage) CALL errore('fetch_input','Err: n_pdep_eigen*n_pdep_times<nimage',1)
     IF(n_pdep_maxiter < 1) CALL errore('fetch_input','Err: n_pdep_maxiter<1',1)
     IF(n_dfpt_maxiter < 1) CALL errore('fetch_input','Err: n_dfpt_maxiter<1',1)
     IF(n_pdep_read_from_file < 0) CALL errore('fetch_input','Err: n_pdep_read_from_file<0',1)
     IF(n_pdep_read_from_file > n_pdep_eigen) CALL errore('fetch_input','Err: n_pdep_read_from_file>n_pdep_eigen',1)
     IF(tr2_dfpt <= 0._DP) CALL errore('fetch_input','Err: tr2_dfpt<0.',1)
     IF(trev_pdep <= 0._DP) CALL errore('fetch_input','Err: trev_pdep<0.',1)
     IF(trev_pdep_rel <= 0._DP) CALL errore('fetch_input','Err: trev_pdep_rel<0.',1)
     IF(n_pdep_eigen == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_pdep_eigen',1)
     IF(n_pdep_times == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_pdep_times',1)
     IF(n_pdep_maxiter == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_pdep_maxiter',1)
     IF(n_dfpt_maxiter == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_dfpt_maxiter',1)
     IF(n_pdep_read_from_file == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_pdep_read_from_file',1)
     IF(n_steps_write_restart == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_steps_write_restart',1)
     IF(gamma_only) THEN
        IF(SIZE(qlist)/=1) CALL errore('fetch_input','Err: SIZE(qlist)/=1',1)
     ELSE
        IF(SIZE(qlist)>nk1*nk2*nk3) CALL errore('fetch_input','Err: SIZE(qlist)>nk1*nk2*nk3',1)
     ENDIF
     !
  ENDIF
  !
  IF(ANY(driver(:)==3)) THEN
     !
     CALL mp_bcast(wfreq_calculation,root,world_comm)
     CALL mp_bcast(n_pdep_eigen_to_use,root,world_comm)
     CALL mp_bcast(qp_bandrange,root,world_comm)
     IF(mpime == root) n_qp_bands = SIZE(qp_bands)
     CALL mp_bcast(n_qp_bands,root,world_comm)
     IF(mpime /= root) THEN
        IF(ALLOCATED(qp_bands)) DEALLOCATE(qp_bands)
        ALLOCATE(qp_bands(n_qp_bands))
     ENDIF
     CALL mp_bcast(qp_bands,root,world_comm)
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
     CALL mp_bcast(l_qdet_verbose,root,world_comm)
     CALL mp_bcast(l_enable_off_diagonal,root,world_comm)
     CALL mp_bcast(o_restart_time,root,world_comm)
     CALL mp_bcast(ecut_spectralf,root,world_comm)
     CALL mp_bcast(n_spectralf,root,world_comm)
     !
     ! CHECKS
     !
     IF(.NOT. gamma_only) THEN
        DO i = 1, 9
           IF(wfreq_calculation(i:i) == 'H') CALL errore('fetch_input','Err: QDET requires gamma_only',1)
        ENDDO
     ENDIF
     IF(.NOT. l_enable_off_diagonal) THEN
        DO i = 1, 9
           IF(wfreq_calculation(i:i) == 'H') CALL errore('fetch_input','Err: QDET requires l_enable_off_diagonal',1)
        ENDDO
     ENDIF
     IF(n_lanczos < 2) CALL errore('fetch_input','Err: n_lanczos<2',1)
     IF(n_pdep_eigen_to_use < 1) CALL errore('fetch_input','Err: n_pdep_eigen_to_use<1',1)
     IF(n_pdep_eigen_to_use > n_pdep_eigen) CALL errore('fetch_input','Err: n_pdep_eigen_to_use>n_pdep_eigen',1)
     IF(n_imfreq < 1) CALL errore('fetch_input','Err: n_imfreq<1',1)
     IF(n_refreq < 1) CALL errore('fetch_input','Err: n_refreq<1',1)
     IF(n_spectralf < 2) CALL errore('fetch_input','Err: n_spectralf<1',1)
     IF(qp_bandrange(1) < 1) CALL errore('fetch_input','Err: qp_bandrange(1)<1',1)
     IF(qp_bandrange(2) < 1) CALL errore('fetch_input','Err: qp_bandrange(2)<1',1)
     IF(qp_bandrange(2) < qp_bandrange(1)) CALL errore('fetch_input','Err: qp_bandrange(2)<qp_bandrange(1)',1)
     IF(qp_bands(1) /= 0) THEN
        DO i = 0, SIZE(qp_bands)-1 ! Python indices start at 0
           IF(qp_bands(i+1) < 1) CALL errore('fetch_input','Err: qp_bands<1',1)
           IF(i /= SIZE(qp_bands)-1) THEN
              IF(qp_bands(i+1) >= qp_bands(i+2)) CALL errore('fetch_input','Err: qp_bands must be sorted in ascending order',1)
           ENDIF
        ENDDO
     ENDIF
     IF(ecut_imfreq <= 0._DP) CALL errore('fetch_input','Err: ecut_imfreq<0.',1)
     IF(ecut_refreq <= 0._DP) CALL errore('fetch_input','Err: ecut_imfreq<0.',1)
     IF(ecut_spectralf(2) < ecut_spectralf(1)) CALL errore('fetch_input','Err: ecut_spectralf(2)<ecut_spectralf(1)',1)
     IF(wfreq_eta <= 0._DP) CALL errore('fetch_input','Err: wfreq_eta<0.',1)
     IF(n_secant_maxiter < 0) CALL errore('fetch_input','Err: n_secant_maxiter<0',1)
     IF(trev_secant <= 0._DP) CALL errore('fetch_input','Err: trev_secant<0.',1)
     IF(l_enable_off_diagonal .AND. .NOT. gamma_only) CALL errore('fetch_input','Err: off-diagonal implemented for gamma only',1)
     IF(n_pdep_eigen_to_use == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_pdep_eigen_to_use',1)
     IF(n_lanczos == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_lanczos',1)
     IF(n_imfreq == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_imfreq',1)
     IF(n_refreq == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_refreq',1)
     IF(n_secant_maxiter == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_secant_maxiter',1)
     IF(n_spectralf == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_spectralf',1)
     !
     SELECT CASE(macropol_calculation)
     CASE('N','n','C','c')
     CASE DEFAULT
        CALL errore('fetch_input','Err: macropol_calculation/=(N,C)',1)
     END SELECT
     !
  ENDIF
  !
  IF(ANY(driver(:)==4)) THEN
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
     CALL mp_bcast(westpp_box,root,world_comm)
     !
     ! CHECKS
     !
     IF(westpp_range(1) < 1) CALL errore('fetch_input','Err: westpp_range(1)<1',1)
     IF(westpp_range(2) < 1) CALL errore('fetch_input','Err: westpp_range(2)<1',1)
     IF(westpp_range(2) < westpp_range(1)) CALL errore('fetch_input','Err: westpp_range(2)<westpp_range(1)',1)
     IF(westpp_nr < 1) CALL errore('fetch_input','Err: westpp_nr<1',1)
     IF(westpp_n_pdep_eigen_to_use < 1) CALL errore('fetch_input','Err: westpp_n_pdep_eigen_to_use<1',1)
     IF(westpp_rmax < 0._DP) CALL errore('fetch_input','Err: westpp_rmax<0.',1)
     IF(westpp_epsinfty < 1._DP) CALL errore('fetch_input','Err: westpp_epsinfty<1.',1)
     IF(westpp_n_pdep_eigen_to_use == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch westpp_n_pdep_eigen_to_use',1)
     IF(westpp_nr == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch westpp_nr',1)
     IF(westpp_box(1) > westpp_box(2) .OR. westpp_box(3) > westpp_box(4) .OR. westpp_box(5) > westpp_box(6)) &
     & CALL errore('fetch_input','Err: invalid westpp_box',1)
     !
  ENDIF
  !
  IF(ANY(driver(:)==5)) THEN
     !
     lenc = LEN(document)
     CALL mp_bcast(lenc,root,world_comm)
     IF(mpime /= root) ALLOCATE(CHARACTER(LEN=lenc) :: document)
     CALL mp_bcast(document,root,world_comm)
     !
  ENDIF
  !
  IF(ANY(driver(:)==6)) THEN
     !
     IF(mpime == root) nq = SIZE(qlist)
     CALL mp_bcast(nq,root,world_comm)
     IF(mpime /= root) THEN
        IF(ALLOCATED(qlist)) DEALLOCATE(qlist)
        ALLOCATE(qlist(nq))
     ENDIF
     CALL mp_bcast(qlist,root,world_comm)
     !
     CALL mp_bcast(wbse_init_calculation,root,world_comm)
     CALL mp_bcast(localization,root,world_comm)
     CALL mp_bcast(wfc_from_qbox,root,world_comm)
     CALL mp_bcast(bisection_info,root,world_comm)
     CALL mp_bcast(chi_kernel,root,world_comm)
     CALL mp_bcast(overlap_thr,root,world_comm)
     CALL mp_bcast(spin_channel,root,world_comm)
     !
     IF(.NOT. gamma_only) CALL errore('fetch_input','Err: BSE requires gamma_only',1)
     IF(spin_channel == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch spin_channel',1)
     !
     SELECT CASE(TRIM(localization))
     CASE('N','n','B','b')
     CASE DEFAULT
        CALL errore('fetch_input','Err: localization/=(N,B)',1)
     END SELECT
     !
     IF(spin_channel < 0 .OR. spin_channel > 2) THEN
        CALL errore('fetch_input','Err: spin_channel/=0,1,2',spin_channel)
     ENDIF
     !
     !use_qbox = .FALSE.
     !use_wstat_pdep = .FALSE.
     !SELECT CASE(TRIM(which_bse_method))
     !CASE('finite-field','FF')
     !use_qbox = .TRUE.
     !CASE('PDEP','pdep')
     !use_wstat_pdep = .TRUE.
     !CASE DEFAULT
     !   CALL errore('fetch_input','Err: which_bse_method != FF/PDEP', 1)
     !END SELECT
     !
     !CALL mp_bcast(use_qbox,root,world_comm)
     !CALL mp_bcast(use_wstat_pdep,root,world_comm)
     !
  ENDIF
  !
  IF(ANY(driver(:)==7)) THEN
     !
     CALL mp_bcast(wbse_calculation,root,world_comm)
     CALL mp_bcast(solver,root,world_comm)
     CALL mp_bcast(qp_correction,root,world_comm)
     CALL mp_bcast(scissor_ope,root,world_comm)
     CALL mp_bcast(n_liouville_eigen,root,world_comm)
     CALL mp_bcast(n_liouville_times,root,world_comm)
     CALL mp_bcast(n_liouville_maxiter,root,world_comm)
     CALL mp_bcast(n_liouville_read_from_file,root,world_comm)
     CALL mp_bcast(trev_liouville,root,world_comm)
     CALL mp_bcast(trev_liouville_rel,root,world_comm)
     CALL mp_bcast(n_lanczos,root,world_comm)
     CALL mp_bcast(n_steps_write_restart,root,world_comm)
     CALL mp_bcast(ipol_input,root,world_comm)
     CALL mp_bcast(macropol_calculation,root,world_comm)
     CALL mp_bcast(wbse_epsinfty,root,world_comm)
     CALL mp_bcast(spin_excitation,root,world_comm)
     CALL mp_bcast(l_preconditioning,root,world_comm)
     CALL mp_bcast(l_reduce_io,root,world_comm)
     !
     ! CHECKS
     !
     IF(.NOT. gamma_only) CALL errore('fetch_input','Err: BSE requires gamma_only',1)
     IF(n_liouville_eigen == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_liouville_eigen',1)
     IF(n_liouville_times == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_liouville_times',1)
     IF(n_liouville_maxiter == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_liouville_maxiter',1)
     IF(n_liouville_read_from_file == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_liouville_read_from_file',1)
     IF(n_lanczos == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_lanczos',1)
     IF(n_steps_write_restart == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_steps_write_restart',1)
     !
     SELECT CASE(macropol_calculation)
     CASE('N','n','C','c')
     CASE DEFAULT
        CALL errore('fetch_input','Err: macropol_calculation/=(N,C)',1)
     END SELECT
     !
     SELECT CASE(TRIM(solver))
     CASE('BSE','bse','TDDFT','tddft')
     CASE DEFAULT
        CALL errore('fetch_input','Err: solver must be BSE or TDDFT',1)
     END SELECT
     !
     SELECT CASE(wbse_calculation)
     CASE('D','d','L','l')
     CASE DEFAULT
        CALL errore('fetch_input','Err: wbse_calculation/=(D,L)',1)
     END SELECT
     !
     SELECT CASE(spin_excitation)
     CASE('s','S','singlet','t','T','triplet')
     CASE DEFAULT
        CALL errore('fetch_input','Err: spin_excitation/=(S,T)',1)
     END SELECT
     !
     IF(wbse_calculation == 'D' .OR. wbse_calculation == 'd') THEN
        IF(n_liouville_times < 2) CALL errore('fetch_input','Err: n_liouville_times<2',1)
        IF(n_liouville_eigen < 1) CALL errore('fetch_input','Err: n_liouville_eigen<1',1)
        IF(n_liouville_eigen*n_liouville_times < nimage) &
        & CALL errore('fetch_input','Err: n_liouville_eigen*n_liouville_times<nimage',1)
        IF(n_liouville_maxiter < 1) CALL errore('fetch_input','Err: n_liouville_maxiter<1',1)
        IF(n_liouville_read_from_file < 0) CALL errore('fetch_input','Err: n_liouville_read_from_file<0',1)
        IF(n_liouville_read_from_file > n_liouville_eigen) &
        & CALL errore('fetch_input','Err: n_liouville_read_from_file>n_liouville_eigen',1)
     ENDIF
     !
  ENDIF
  !
  IF(ANY(driver(:)==8)) THEN
     !
     IF(mpime == root) nq = SIZE(qlist)
     CALL mp_bcast(nq,root,world_comm)
     IF(mpime /= root) THEN
        IF(ALLOCATED(qlist)) DEALLOCATE(qlist)
        ALLOCATE(qlist(nq))
     ENDIF
     CALL mp_bcast(qlist,root,world_comm)
     !
     CALL mp_bcast(wbsepp_calculation,root,world_comm)
     CALL mp_bcast(n_liouville_read_from_file,root,world_comm)
     CALL mp_bcast(ipol_input,root,world_comm)
     CALL mp_bcast(spin_channel,root,world_comm)
     CALL mp_bcast(wbsepp_r0,root,world_comm)
     CALL mp_bcast(iexc_plot,root,world_comm)
     CALL mp_bcast(n_lanczos_to_use,root,world_comm)
     CALL mp_bcast(n_extrapolation,root,world_comm)
     CALL mp_bcast(range,root,world_comm)
     CALL mp_bcast(broad,root,world_comm)
     !
     ! CHECKS
     !
     IF(.NOT. gamma_only) CALL errore('fetch_input','Err: BSE requires gamma_only',1)
     IF(n_liouville_read_from_file == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_liouville_read_from_file',1)
     IF(spin_channel == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch spin_channel',1)
     IF(iexc_plot == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch iexc_plot',1)
     IF(n_lanczos_to_use == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_lanczos_to_use',1)
     IF(n_extrapolation == DUMMY_DEFAULT) CALL errore('fetch_input','Err: cannot fetch n_extrapolation',1)
     !
     IF(spin_channel < 1 .OR. spin_channel > 2) CALL errore('fetch_input','Err: spin_channel/=1,2',spin_channel)
     !
  ENDIF
  !
  CALL mp_barrier(world_comm)
  !
  IF(verbose .AND. mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     !
     CALL add_intput_parameters_to_json_file(num_drivers, driver, json)
     !
     OPEN(NEWUNIT=iunit, FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
     !
  ENDIF
  !
  CALL stop_clock('fetch_input')
  !
END SUBROUTINE
