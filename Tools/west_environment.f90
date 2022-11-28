!
! Copyright (C) 2015-2022 M. Govoni
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
MODULE west_environment
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: west_environment_start
  PUBLIC :: west_environment_end
  !
  !-----------------------------------------------------------------------
CONTAINS
  !-----------------------------------------------------------------------
  !
  SUBROUTINE west_environment_start( code )
    !
    USE io_global,             ONLY : stdout, meta_ionode
    USE io_files,              ONLY : crash_file, nd_nmbr
    USE mp_images,             ONLY : me_image, my_image_id
    USE westcom,               ONLY : savedir, logfile, outdir, west_prefix
    USE base64_module,         ONLY : base64_init
    USE json_string_utilities, ONLY : lowercase_string
    USE west_version,          ONLY : start_forpy
#if defined(__HDF5)
    USE hdf5_qe,               ONLY : phdf5_start => initialize_hdf5
    USE qeh5_base_module,      ONLY : hdf5_start => initialize_hdf5
#endif
#if defined(__CUDA)
    USE west_gpu,              ONLY : west_gpu_start
#endif
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: code
    !
    LOGICAL :: exst, debug = .FALSE.
    CHARACTER(LEN=80) :: uname
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    INTEGER :: ios, crashunit, n_json
    INTEGER, PARAMETER :: n_json_max = 100
    LOGICAL, EXTERNAL :: matches
    !
    ! ... Intel compilers v .ge.8 allocate a lot of stack space
    ! ... Stack limit is often small, thus causing SIGSEGV and crash
    !
#if defined(__INTEL_COMPILER)
    CALL remove_stack_limit ( )
#endif
    CALL start_forpy()
    !
    ! Input from (-i), output from (-o)
    !
    CALL parse_command_arguments()
    CALL fetch_input_yml(1,(/1/),.FALSE.)
    !
    savedir = TRIM(outdir) // TRIM(west_prefix) // '.' // TRIM(lowercase_string(code)) // '.save/'
    logfile = TRIM(savedir) // TRIM(lowercase_string(code)) // '.json'
    !
    ! Do not overwrite existing JSON file
    !
    n_json = 1
    !
    INQUIRE(FILE=TRIM(logfile),EXIST=exst)
    !
    DO WHILE(exst)
       logfile = TRIM(savedir)//TRIM(lowercase_string(code))//'_'//TRIM(int_to_char(n_json))//'.json'
       INQUIRE(FILE=TRIM(logfile),EXIST=exst)
       n_json = n_json+1
       IF(n_json > n_json_max) CALL errore(TRIM(code),'Too many JSON files',1)
    ENDDO
    !
    CALL my_mkdir( TRIM(savedir) )
    !
    ! ... use .FALSE. to disable all clocks except the total cpu time clock
    ! ... use .TRUE. to enable clocks
    !
    CALL init_clocks( .TRUE. )
    CALL start_clock( TRIM(code) )
    !
    ! ... for compatibility with PWSCF
    !
#if defined(__MPI)
    nd_nmbr = TRIM( int_to_char( me_image+1 ))
#else
    nd_nmbr = ' '
#endif
    !
    IF( meta_ionode ) THEN
       !
       ! ... search for file CRASH and delete it
       !
       INQUIRE( FILE=TRIM(crash_file), EXIST=exst )
       IF( exst ) THEN
          OPEN( NEWUNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD',IOSTAT=ios )
          IF (ios==0) THEN
             CLOSE( UNIT=crashunit, STATUS='DELETE', IOSTAT=ios )
          ELSE
             WRITE(stdout,'(5x,"Remark: CRASH file could not be deleted")')
          ENDIF
       ENDIF
       !
    ELSE
       ! ... one processor per image (other than meta_ionode)
       ! ... or, for debugging purposes, all processors,
       ! ... open their own standard output file
#if defined(DEBUG) || defined(__DEBUG)
       debug = .TRUE.
#endif
       IF (debug ) THEN
          uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // &
               trim(int_to_char( me_image))
          OPEN( UNIT = stdout, FILE = TRIM(uname),STATUS='UNKNOWN')
       ELSE
#if defined(_WIN32)
          OPEN( UNIT = stdout, FILE='NUL:', STATUS='UNKNOWN' )
#else
          OPEN( UNIT = stdout, FILE='/dev/null', STATUS='UNKNOWN' )
#endif
       ENDIF
       !
    ENDIF
    !
    ! Initialize base64 tables
    CALL base64_init()
    !
    CALL west_opening_message( code )
#if defined(__MPI)
    CALL report_parallel_status( )
#else
    CALL errore(TRIM(code), 'West needs MPI to run', 1 )
#endif
    !
#if defined(__CUDA)
    CALL west_gpu_start( )
#endif
    !
#if defined(__HDF5)
    CALL hdf5_start( )
    CALL phdf5_start( )
#endif
    !
  END SUBROUTINE
  !
  SUBROUTINE west_environment_end( code )
    !
    USE io_global,        ONLY : stdout, meta_ionode
    USE json_module,      ONLY : json_file
    USE mp_world,         ONLY : mpime,root,world_comm
    USE mp,               ONLY : mp_barrier
    USE westcom,          ONLY : logfile
    USE west_version,     ONLY : end_forpy
#if defined(__HDF5)
    USE hdf5_qe,          ONLY : phdf5_end => finalize_hdf5
    USE qeh5_base_module, ONLY : hdf5_end => finalize_hdf5
#endif
#if defined(__CUDA)
    USE west_gpu,         ONLY : west_gpu_end
#endif
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: code
    INTEGER :: iunit
    TYPE(json_file) :: json
    CHARACTER(LEN=9) :: cdate, ctime
    CHARACTER(LEN=80) :: time_str
    LOGICAL :: found
    !
#if defined(_HDF5)
    CALL phdf5_end( )
    CALL hdf5_end( )
#endif
    !
#if defined(__CUDA)
    CALL west_gpu_end( )
#endif
    !
    IF ( meta_ionode ) WRITE( stdout, * )
    !
    CALL stop_clock(  TRIM(code) )
    CALL print_clock( TRIM(code) )
    !
    CALL date_and_tim( cdate, ctime )
    !
    time_str = 'This run was terminated on:  ' // ctime // ' ' // cdate
    !
    IF( meta_ionode ) THEN
       WRITE( stdout,*)
       WRITE( stdout,3334) time_str
       WRITE( stdout,3335)
    ENDIF
    !
    IF( meta_ionode ) THEN
       WRITE( stdout,'(A)')      '   JOB DONE.'
       WRITE( stdout,3335)
    ENDIF
3334 FORMAT(3X,A60,/)
3335 FORMAT('=',78('-'),'=')
    FLUSH(stdout)
    !
    IF( mpime == root ) THEN
      !
      CALL json%initialize()
      CALL json%load(filename=TRIM(logfile))
      !
      CALL json%update('runjob.completed', .TRUE., found)
      CALL json%add('runjob.endtime', TRIM(ctime) )
      CALL json%add('runjob.enddate', TRIM(cdate) )
      !
      OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
      CALL json%print( iunit )
      CLOSE( iunit )
      !
      CALL json%destroy()
      !
    ENDIF
    !
    CALL end_forpy()
    !
    CALL mp_barrier(world_comm)
    !
  END SUBROUTINE
  !
  SUBROUTINE west_opening_message( code )
    !
    USE json_module,     ONLY : json_file
    USE io_global,       ONLY : stdout
    USE global_version,  ONLY : version_number
    USE west_version,    ONLY : west_version_number, west_git_revision
    USE mp_world,        ONLY : mpime,root
    USE westcom,         ONLY : logfile
    USE base64_module,   ONLY : islittleendian
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=*), INTENT(IN) :: code
    !
    ! Workspace
    !
    TYPE(json_file) :: json
    INTEGER :: iunit
    CHARACTER(LEN=9) :: cdate, ctime
    !
    CALL date_and_tim( cdate, ctime )
    !
    IF ( TRIM(west_git_revision) /= 'unknown' ) THEN
       WRITE( stdout, '(/5X,"Program ",A," v. ",A," git rev. ",A," starts on ",A9," at ",A9)' ) &
       & TRIM(code), TRIM(west_version_number), TRIM(west_git_revision), cdate, ctime
    ELSE
       WRITE( stdout, '(/5X,"Program ",A," v. ",A," starts on ",A9," at ",A9)' ) &
       & TRIM(code), TRIM(west_version_number), cdate, ctime
    ENDIF
    !
    WRITE( stdout, '(/5X,"This program is part of the open-source West suite",&
    &/5X,"for massively parallel calculations of excited states in materials; please cite", &
    &/9X,"""M. Govoni et al., J. Chem. Theory Comput. 11, 2680 (2015);",&
    &/9X," URL http://www.west-code.org"", ", &
    &/5X,"in publications or presentations arising from this work.")' )
    !
    WRITE( stdout, '(/5X,"Based on the Quantum ESPRESSO v. ",A)') TRIM(version_number)
    !
    IF( islittleendian() ) THEN
       WRITE( stdout, '(/5X,"I/O is Little Endian",A)') ""
    ELSE
       WRITE( stdout, '(/5X,"I/O is Big Endian",A)') ""
    ENDIF
    !
    IF( mpime == root ) THEN
      !
      CALL json%initialize()
      !
      CALL json%add('runjob.startdate', TRIM(cdate) )
      CALL json%add('runjob.starttime', TRIM(ctime) )
      CALL json%add('runjob.completed', .FALSE. )
      CALL json%add('software.package', 'WEST' )
      CALL json%add('software.program', TRIM(code) )
      CALL json%add('software.version', TRIM(west_version_number) )
      IF( TRIM(west_git_revision) /= 'unknown' ) CALL json%add('software.westgit', TRIM(west_git_revision) )
      CALL json%add('software.website','http://www.west-code.org')
      CALL json%add('software.citation','M. Govoni et al., J. Chem. Theory Comput. 11, 2680 (2015).')
      CALL json%add('software.qeversion', TRIM(version_number) )
      CALL json%add('config.io.islittleendian', islittleendian() )
      !
      OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
      CALL json%print( iunit )
      CLOSE( iunit )
      !
      CALL json%destroy()
      !
    ENDIF
    !
  END SUBROUTINE
  !
  SUBROUTINE report_parallel_status()
     !-----------------------------------------------------------------------
     !
     ! ... Report the mpi/openmp status
     !
     USE json_module,      ONLY : json_file
     USE io_global,        ONLY : stdout
     USE mp_global,        ONLY : nimage,npool,nbgrp,nproc_bgrp
     USE mp_world,         ONLY : nproc,mpime,root
     USE io_push,          ONLY : io_push_title,io_push_bar
     USE westcom,          ONLY : logfile
     !
     IMPLICIT NONE
     !
#if defined(_OPENMP)
     INTEGER, EXTERNAL :: omp_get_max_threads
#endif
     !
     INTEGER :: nth, ncores
     TYPE(json_file) :: json
     INTEGER :: iunit
     !
#if defined(_OPENMP)
     nth = omp_get_max_threads()
#else
     nth = 1
#endif
     ncores = nproc * nth
     !
     CALL io_push_title('**MPI** Parallelization Status')
     WRITE(stdout, "(5x, '   ',i14,'      ',4i14)") nproc, nimage, npool, nbgrp, nproc_bgrp
     CALL io_push_bar()
     WRITE(stdout, "(5x, '                N         =         I      X      P      X      B      X      Z')")
     WRITE(stdout, "(5x, '                ^                   ^             ^             ^             ^')")
     WRITE(stdout, "(5x, '                |                   |             |             |             |')")
     WRITE(stdout, "(5x, '              #rnk                  |             |             |             |')")
     WRITE(stdout, "(5x, '                                 #image           |             |             |')")
     WRITE(stdout, "(5x, '                                                #pool           |             |')")
     WRITE(stdout, "(5x, '                                                              #bgrp           |')")
     WRITE(stdout, "(5x, '                                                                            #R&G')")
     CALL io_push_bar()
     !
#if defined(_OPENMP)
     WRITE(stdout, "(5x, '**OPENMP** Parallelization Status')")
     WRITE(stdout, "(5x, '#thr/rnk               = ',i12)") nth
     CALL io_push_bar()
     WRITE(stdout, "(5x, '#prc = (rnk) * (thr) = ',i12)") ncores
     CALL io_push_bar()
#else
     WRITE(stdout, "(5x, '#prc = ',i12)") ncores
     CALL io_push_bar()
#endif
    !
    IF( mpime == root ) THEN
       !
       CALL json%initialize()
       !
       CALL json%load(filename=TRIM(logfile))
       !
       CALL json%add('parallel.nranks', nproc )
       CALL json%add('parallel.nimage', nimage )
       CALL json%add('parallel.npool', npool )
       CALL json%add('parallel.nbgrp', nbgrp )
       CALL json%add('parallel.nrg', nproc_bgrp )
       CALL json%add('parallel.nproc', ncores )
#if defined(_OPENMP)
       CALL json%add('parallel.nthreads', nth )
#endif
       !
       OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
       CALL json%print( iunit )
       CLOSE( iunit )
       !
       CALL json%destroy()
       !
    ENDIF
    !
  END SUBROUTINE
  !
END MODULE
