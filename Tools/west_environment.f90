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
MODULE west_environment
  !-----------------------------------------------------------------------
  !
  USE io_global,       ONLY: stdout, meta_ionode
  !
  IMPLICIT NONE
  !
  SAVE
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
    USE kinds,           ONLY: DP
    USE io_files,        ONLY: crash_file, nd_nmbr
    USE mp_images,       ONLY: me_image, my_image_id, root_image, nimage
    !
    CHARACTER(LEN=*), INTENT(IN) :: code
    !
    LOGICAL           :: exst, debug = .false.
    CHARACTER(LEN=80) :: uname
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    INTEGER :: ios, crashunit
    INTEGER, EXTERNAL :: find_free_unit
    !
    ! ... Intel compilers v .ge.8 allocate a lot of stack space
    ! ... Stack limit is often small, thus causing SIGSEGV and crash
    !
#if defined(__INTEL_COMPILER)
    CALL remove_stack_limit ( )
#endif
    !
    ! ... use ".FALSE." to disable all clocks except the total cpu time clock
    ! ... use ".TRUE."  to enable clocks
    !
    CALL init_clocks( .TRUE. )
    CALL start_clock( TRIM(code) )
    !
    ! ... for compatibility with PWSCF
    !
#if defined(__MPI)
    nd_nmbr = TRIM ( int_to_char( me_image+1 ))
#else
    nd_nmbr = ' '
#endif
    !
    IF( meta_ionode ) THEN
       !
       ! ...  search for file CRASH and delete it
       !
       INQUIRE( FILE=TRIM(crash_file), EXIST=exst )
       IF( exst ) THEN
          crashunit = find_free_unit()
          OPEN( UNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD',IOSTAT=ios )
          IF (ios==0) THEN
             CLOSE( UNIT=crashunit, STATUS='DELETE', IOSTAT=ios )
          ELSE
             WRITE(stdout,'(5x,"Remark: CRASH file could not be deleted")')
          END IF
       END IF
       !
    ELSE
       ! ... one processor per image (other than meta_ionode)
       ! ... or, for debugging purposes, all processors,
       ! ... open their own standard output file
#if defined(DEBUG)
       debug = .true.
#endif
       IF (debug ) THEN
          uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // &
               trim(int_to_char( me_image))
          OPEN ( unit = stdout, file = TRIM(uname),status='unknown')
       ELSE
#if defined(_WIN32)
          OPEN ( unit = stdout, file='NUL:', status='unknown' )
#else
          OPEN ( unit = stdout, file='/dev/null', status='unknown' )
#endif
       END IF
       !
    END IF
    !
    CALL west_opening_message( code )
#if defined(__MPI)
    CALL report_parallel_status ( )
#else
    CALL errore(TRIM(code), 'West need MPI to run', 1 ) 
#endif
  END SUBROUTINE
  !
  !
  SUBROUTINE west_environment_end( code )
    !
    CHARACTER(LEN=*), INTENT(IN) :: code
    !
    IF ( meta_ionode ) WRITE( stdout, * )
    !
    CALL stop_clock(  TRIM(code) )
    CALL print_clock( TRIM(code) )
    !
    CALL west_closing_message( )
    !
    IF( meta_ionode ) THEN
       WRITE( stdout,'(A)')      '   JOB DONE.'
       WRITE( stdout,3335)
    END IF
3335 FORMAT('=',78('-'),'=')
    FLUSH(stdout)
    !
    RETURN
  END SUBROUTINE
  !
  !
  !
  SUBROUTINE west_opening_message( code )
    !
    USE io_global,       ONLY : stdout
    USE global_version,  ONLY : version_number, svn_revision
    USE west_version,    ONLY : west_version_number, west_svn_revision
    !
    ! I/O
    !
    CHARACTER(LEN=*), INTENT(IN) :: code
    !
    ! Workspace
    !
    CHARACTER(LEN=9)  :: cdate, ctime
    !
    CALL date_and_tim( cdate, ctime )
    !
    IF ( TRIM (west_svn_revision) /= "unknown" ) THEN
       WRITE( stdout, '(/5X,"Program ",A," v. ",A," svn rev. ",A," starts on ",A9," at ",A9)' ) &
       & TRIM(code), TRIM(west_version_number), TRIM (west_svn_revision), cdate, ctime
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
    IF ( TRIM (svn_revision) /= "unknown" ) THEN 
       WRITE( stdout, '(/5X,"Based on the Quantum ESPRESSO v. ",A," svn rev. ",A)') TRIM (version_number), TRIM (svn_revision)
    ELSE
       WRITE( stdout, '(/5X,"Based on the Quantum ESPRESSO v. ",A)') TRIM (version_number)
    ENDIF
    !
    RETURN
  END SUBROUTINE 
  !
  ! 
  !
  SUBROUTINE west_closing_message( )
    !
    CHARACTER(LEN=9)  :: cdate, ctime
    CHARACTER(LEN=80) :: time_str
    !
    CALL date_and_tim( cdate, ctime )
    !
    time_str = 'This run was terminated on:  ' // ctime // ' ' // cdate
    !
    IF( meta_ionode ) THEN
       WRITE( stdout,*)
       WRITE( stdout,3334) time_str
       WRITE( stdout,3335)
    END IF
    !
3334 FORMAT(3X,A60,/)
3335 FORMAT('=',78('-'),'=')
    !
    RETURN
  END SUBROUTINE
  !
  !
  !
  SUBROUTINE report_parallel_status()
     !-----------------------------------------------------------------------
     !
     ! ... Report the mpi/openmp status
     !
     USE io_global,        ONLY : stdout
     USE mp_global,        ONLY : nimage,npool,nbgrp,nproc_image,nproc_pool,nproc_bgrp 
     USE mp_world,         ONLY : nproc 
     USE io_push,          ONLY : io_push_title,io_push_bar
     !
     IMPLICIT NONE
     !
#if defined(__OPENMP)
     INTEGER, EXTERNAL :: omp_get_max_threads
#endif
     !
     INTEGER :: nth, ncores 
     !
#if defined(__OPENMP)
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
#if defined(__OPENMP)
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
  END SUBROUTINE
  !
END MODULE 
