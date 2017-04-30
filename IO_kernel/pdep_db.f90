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
MODULE pdep_db
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir
  !
  IMPLICIT NONE
  !
  !
  CONTAINS
    !
    !
    ! *****************************
    ! PDEP WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE pdep_db_write( )
      !------------------------------------------------------------------------
      !
      USE xml_io_base,          ONLY : create_directory
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp_global,            ONLY : my_image_id
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : wstat_calculation,n_pdep_times,n_pdep_eigen,n_pdep_maxiter,n_dfpt_maxiter, &
                                     & n_steps_write_restart,n_pdep_restart_from_itr,n_pdep_read_from_file,trev_pdep, &
                                     & tr2_dfpt,l_deflate,l_kinetic_only,ev,dvg,west_prefix,wstat_dirname,trev_pdep_rel, &
                                     & l_minimize_exx_if_active,l_use_ecutrho 
      !USE eig_distribute,       ONLY : local_npert,pdep_distr_l2g
      USE pdep_io,              ONLY : pdep_merge_and_write_G 
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : pert 
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=256)    :: fname
      CHARACTER(LEN=6)      :: my_label
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: iunout,global_j,local_j
      INTEGER :: ierr
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('pdep_db')
      time_spent(1)=get_clock('pdep_db')
      !
      ! 1)  CREATE THE INPUT FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( wstat_dirname ) // '/' // TRIM("input-file.xml") , BINARY=.FALSE.,IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      CALL errore( 'pdep_db', 'cannot open input-file.xml for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         CALL iotk_write_begin( iunout, "WSTAT_CONTROL" )
         !
         CALL iotk_write_dat( iunout, "wstat_calculation"        , wstat_calculation)
         CALL iotk_write_dat( iunout, "n_pdep_eigen"             , n_pdep_eigen)
         CALL iotk_write_dat( iunout, "n_pdep_times"             , n_pdep_times)
         CALL iotk_write_dat( iunout, "n_pdep_maxiter"           , n_pdep_maxiter)
         CALL iotk_write_dat( iunout, "n_dfpt_maxiter"           , n_dfpt_maxiter)
         CALL iotk_write_dat( iunout, "n_pdep_read_from_file"    , n_pdep_read_from_file)
         CALL iotk_write_dat( iunout, "trev_pdep"                , trev_pdep)
         CALL iotk_write_dat( iunout, "trev_pdep_rel"            , trev_pdep_rel)
         CALL iotk_write_dat( iunout, "tr2_dfpt"                 , tr2_dfpt)
         CALL iotk_write_dat( iunout, "l_kinetic_only"           , l_kinetic_only)
         CALL iotk_write_dat( iunout, "l_minimize_exx_if_active" , l_minimize_exx_if_active)
         CALL iotk_write_dat( iunout, "l_use_ecutrho"            , l_use_ecutrho)
         !
         CALL iotk_write_end( iunout, "WSTAT_CONTROL"  )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! 2) CREATE THE EIGENVALUE FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( wstat_dirname ) // '/' // TRIM("dbs_eigenvalues.xml"),BINARY=.FALSE.,IERR=ierr)
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      CALL errore( 'pdep_db', 'cannot open dbs_eigenvalues.xml file for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         CALL iotk_write_begin( iunout, "EIGENVALUES" )
         CALL iotk_write_dat( iunout, "ndim", n_pdep_eigen )
         CALL iotk_write_dat( iunout, "ev", ev(1:n_pdep_eigen))
         CALL iotk_write_end( iunout, "EIGENVALUES" )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! 3) CREATE THE EIGENVECTOR FILES
      !
      DO local_j=1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j>n_pdep_eigen) CYCLE
         ! 
         fname = TRIM( wstat_dirname ) // "/E"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_merge_and_write_G(fname,dvg(:,local_j))
         !
      ENDDO
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2)=get_clock('pdep_db')
      CALL stop_clock('pdep_db')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'Database written in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wstat_dirname )  
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
    ! *****************************
    ! PDEP READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE pdep_db_read( nglob_to_be_read )
      !------------------------------------------------------------------------
      !
      USE westcom,             ONLY : n_pdep_eigen,ev,dvg,west_prefix,npwq0x
      USE io_global,           ONLY : stdout 
      USE mp,                  ONLY : mp_bcast,mp_barrier
      USE mp_world,            ONLY : world_comm,mpime,root
      USE mp_global,           ONLY : my_image_id
      USE pdep_io,             ONLY : pdep_read_G_and_distribute
      USE io_push,             ONLY : io_push_bar
      USE distribution_center, ONLY : pert
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nglob_to_be_read  
      !
      CHARACTER(LEN=256) :: dirname,fname
      CHARACTER(LEN=6)      :: my_label
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: ierr, n_eigen_to_get
      INTEGER :: tmp_n_pdep_eigen
      INTEGER :: dime, iun, global_j, local_j
      REAL(DP),ALLOCATABLE :: tmp_ev(:)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('pdep_db')
      !
      ! TIMING
      !
      time_spent(1)=get_clock('pdep_db')
      !
      ! ... the main db directory
      !
      dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
      !
      ! 1)  READ THE INPUT FILE
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'input-file.xml' ), IERR = ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      IF ( ierr /=0 ) CALL errore( 'pdep_db', 'cannot open input-file.xml file for reading', ierr )
      !
      IF ( mpime==root ) THEN
         !
         CALL iotk_scan_begin( iun, "WSTAT_CONTROL" )
         CALL iotk_scan_dat( iun, "n_pdep_eigen"         , tmp_n_pdep_eigen)
         CALL iotk_scan_end( iun, "WSTAT_CONTROL"  )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_read( iun )
         !
      ENDIF
      !
      CALL mp_bcast( tmp_n_pdep_eigen, root, world_comm )
      !
      ! In case nglob_to_be_read is 0, overwrite it with the read value 
      !
      IF (nglob_to_be_read==0) THEN 
         n_eigen_to_get = tmp_n_pdep_eigen
         n_pdep_eigen=tmp_n_pdep_eigen
      ELSE
         n_eigen_to_get = MIN(tmp_n_pdep_eigen,nglob_to_be_read)
      ENDIF
      !
      ! 2)  READ THE EIGENVALUES FILE
      !
      IF(.NOT.ALLOCATED(ev)) ALLOCATE(ev(n_eigen_to_get))
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'dbs_eigenvalues.xml' ), IERR = ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'pdep_db', 'cannot open dbs_eigenvalues.xml file for reading', ierr )
      !
      IF ( mpime==root ) THEN
         !
         CALL iotk_scan_begin( iun, "EIGENVALUES" )
         CALL iotk_scan_dat( iun, "ndim"     , dime)
         ALLOCATE(tmp_ev(dime))
         CALL iotk_scan_dat( iun, "ev"     , tmp_ev)
         CALL iotk_scan_end( iun, "EIGENVALUES"  )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_read( iun )
         ev(1:nglob_to_be_read) = tmp_ev(1:nglob_to_be_read)
         DEALLOCATE(tmp_ev)
         !
      ENDIF
      !
      CALL mp_bcast( ev, root, world_comm )
      !
      ! 3)  READ THE EIGENVECTOR FILES
      !
      IF(.NOT.ALLOCATED(dvg)) THEN
         ALLOCATE(dvg(npwq0x,pert%nlocx))
         dvg = 0._DP
      ENDIF
      !
      DO local_j=1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j>n_eigen_to_get) CYCLE
         ! 
         fname = TRIM( dirname ) // "/E"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_read_G_and_distribute(fname,dvg(:,local_j))
         !
      ENDDO
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2)=get_clock('pdep_db')
      CALL stop_clock('pdep_db')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'Database read in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( dirname )  
      WRITE(stdout, "(5x, 'Eigen. found : ',i12)") n_eigen_to_get
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
END MODULE
