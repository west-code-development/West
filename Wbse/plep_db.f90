!
! Copyright (C) 2015-2021 M. Govoni
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
MODULE plep_db
  !----------------------------------------------------------------------------
  !   target folder wbse_save_dir
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
    SUBROUTINE plep_db_write( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp_global,            ONLY : my_image_id
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wstat_calculation,n_pdep_times,n_pdep_eigen,n_pdep_maxiter,n_dfpt_maxiter, &
                                     & n_steps_write_restart,n_pdep_restart_from_itr,n_pdep_read_from_file,trev_pdep, &
                                     & tr2_dfpt,l_deflate,l_kinetic_only,ev,west_prefix,wbse_save_dir,trev_pdep_rel, &
                                     & l_minimize_exx_if_active,l_use_ecutrho, dvg_exc
      USE plep_io,              ONLY : plep_merge_and_write_G
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
      CALL start_clock('plep_db')
      time_spent(1)=get_clock('plep_db')
      !
      ! 1)  CREATE THE INPUT FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( wbse_save_dir ) // '/' // TRIM("input-file.xml") , BINARY=.FALSE.,IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      CALL errore( 'plep_db', 'cannot open input-file.xml for writing', ierr )
      !
      IF ( mpime == root ) THEN
         !
         CALL iotk_write_begin( iunout, "WBSE_CONTROL" )
         !
         CALL iotk_write_dat( iunout, "wbse_calculation"        , wstat_calculation)
         CALL iotk_write_dat( iunout, "n_plep_eigen"             , n_pdep_eigen)
         CALL iotk_write_dat( iunout, "n_plep_times"             , n_pdep_times)
         CALL iotk_write_dat( iunout, "n_plep_maxiter"           , n_pdep_maxiter)
         CALL iotk_write_dat( iunout, "n_ldep_read_from_file"    , n_pdep_read_from_file)
         CALL iotk_write_dat( iunout, "trev_plep"                , trev_pdep)
         CALL iotk_write_dat( iunout, "trev_plep_rel"            , trev_pdep_rel)
         !
         CALL iotk_write_end( iunout, "WBSE_CONTROL"  )
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
         CALL iotk_open_write( iunout, FILE = TRIM( wbse_save_dir ) // '/' // TRIM("dbs_eigenvalues.xml"),BINARY=.FALSE.,IERR=ierr)
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      CALL errore( 'plep_db', 'cannot open dbs_eigenvalues.xml file for writing', ierr )
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
         fname = TRIM( wbse_save_dir ) // "/E"//TRIM(ADJUSTL(my_label))//".dat"
         CALL plep_merge_and_write_G(fname,dvg_exc(:,:,:,local_j))
         !
      ENDDO
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2)=get_clock('plep_db')
      CALL stop_clock('plep_db')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'Database written in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wbse_save_dir )
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
    ! *****************************
    ! PLEP READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE plep_db_read( nglob_to_be_read )
      !------------------------------------------------------------------------
      !
      USE pwcom,               ONLY : nks,npwx
      USE westcom,             ONLY : n_pdep_eigen,ev,west_prefix,wbse_save_dir,&
                                       dvg_exc,nbndval0x
      USE io_global,           ONLY : stdout
      USE mp,                  ONLY : mp_bcast,mp_barrier
      USE mp_world,            ONLY : world_comm,mpime,root
      USE mp_global,           ONLY : my_image_id
      USE plep_io,             ONLY : plep_read_G_and_distribute
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
      CALL start_clock('plep_db')
      !
      ! TIMING
      !
      time_spent(1)=get_clock('plep_db')
      !
      ! ... the main db directory
      !
      !dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wbse.save'
      dirname = wbse_save_dir
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
      IF ( ierr /=0 ) CALL errore( 'plep_db', 'cannot open input-file.xml file for reading', ierr )
      !
      IF ( mpime==root ) THEN
         !
         CALL iotk_scan_begin( iun, "WBSE_CONTROL" )
         CALL iotk_scan_dat( iun, "n_plep_eigen"         , tmp_n_pdep_eigen)
         CALL iotk_scan_end( iun, "WBSE_CONTROL"  )
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
         CALL iotk_open_read( iun, FILE = TRIM( wbse_save_dir ) // '/' // TRIM( 'dbs_eigenvalues.xml' ), IERR = ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'plep_db', 'cannot open dbs_eigenvalues.xml file for reading', ierr )
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
      IF(.NOT.ALLOCATED(dvg_exc)) THEN
         ALLOCATE(dvg_exc(npwx,nbndval0x,nks,pert%nlocx))
         dvg_exc = 0._DP
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
         fname = TRIM( wbse_save_dir ) // "/E"//TRIM(ADJUSTL(my_label))//".dat"
         CALL plep_read_G_and_distribute(fname,dvg_exc(:,:,:,local_j))
         !
      ENDDO
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2)=get_clock('plep_db')
      CALL stop_clock('plep_db')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'Database read in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wbse_save_dir )
      WRITE(stdout, "(5x, 'Eigen. found : ',i12)") n_eigen_to_get
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
END MODULE
