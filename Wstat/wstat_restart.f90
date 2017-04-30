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
MODULE wstat_restart
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER, PRIVATE :: iunout
  !
  INTERFACE wstat_restart_write
     MODULE PROCEDURE wstat_restart_write_real, wstat_restart_write_complex
  END INTERFACE
  !
  INTERFACE wstat_restart_read
     MODULE PROCEDURE wstat_restart_read_real, wstat_restart_read_complex
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_write_real( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,me_bgrp,inter_image_comm,nimage
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : n_pdep_basis,n_pdep_eigen,ev,conv,west_prefix,dvg,dng
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
      USE pdep_io,              ONLY : pdep_merge_and_write_G 
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: dav_iter, notcnv, nbase
      REAL(DP),INTENT(IN) :: ew(n_pdep_basis)
      REAL(DP),INTENT(IN) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(IN) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: ierr
      CHARACTER(LEN=256) :: dirname,fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(6) :: my_label
      INTEGER :: local_j,global_j
      INTEGER :: im
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! MKDIR 
      !
      dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.restart'
      CALL my_mkdir( dirname )
      !
      CALL start_clock('wstat_restart')
      time_spent(1)=get_clock('wstat_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("summary.xml") , BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         CALL iotk_write_begin( iunout, "R-SUMMARY" )
         CALL iotk_write_dat( iunout, "dav_iter", dav_iter)
         CALL iotk_write_dat( iunout, "notcnv", notcnv)
         CALL iotk_write_dat( iunout, "nbase", nbase)
         CALL iotk_write_dat( iunout, "conv", conv(:))
         CALL iotk_write_end( iunout, "R-SUMMARY"  )
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! CREATE THE EIG FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("eig.xml") , BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         CALL iotk_write_begin( iunout, "RESTART_EIG" )
         CALL iotk_write_dat( iunout, "ev", ev(:))
         CALL iotk_write_dat( iunout, "ew", ew(:))
         CALL iotk_write_end( iunout, "RESTART_EIG"  )
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! CREATE THE HR FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("hr.dat") , BINARY = .TRUE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         IF(me_bgrp==0) CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_write_begin( iunout, "RESTART_HR_"//TRIM(my_label) )
            CALL iotk_write_dat( iunout, "hr", tmp_distr(:,:))
            CALL iotk_write_end( iunout, "RESTART_HR"//TRIM(my_label) )
            !
         END IF
         !
      ENDDO
      !
      IF ( mpime == root ) CALL iotk_close_write( iunout )
      !
      DEALLOCATE( tmp_distr )
      !
      ! CREATE THE VR FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("vr.dat") , BINARY = .TRUE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         IF(me_bgrp==0) CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_write_begin( iunout, "RESTART_VR_"//TRIM(my_label) )
            CALL iotk_write_dat( iunout, "vr", tmp_distr(:,:))
            CALL iotk_write_end( iunout, "RESTART_VR"//TRIM(my_label) )
            !
         END IF
         !
      ENDDO
      !
      IF ( mpime == root ) CALL iotk_close_write( iunout )
      !
      DEALLOCATE( tmp_distr )
      !
      ! CREATE THE EIGENVECTOR FILES
      !
      DO local_j=1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j>nbase) CYCLE
         ! 
         fname = TRIM( dirname ) // "/V"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_merge_and_write_G(fname,dvg(:,local_j))
         fname = TRIM( dirname ) // "/N"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_merge_and_write_G(fname,dng(:,local_j))
         !
      ENDDO
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      time_spent(2)=get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout, "(5x, '[I/O] RESTART written in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, '[I/O] In location   : ',a)") TRIM( dirname )  
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_write_complex( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,me_bgrp,inter_image_comm,nimage
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : n_pdep_basis,n_pdep_eigen,ev,conv,west_prefix,dvg,dng
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
      USE pdep_io,              ONLY : pdep_merge_and_write_G 
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: dav_iter, notcnv, nbase
      REAL(DP),INTENT(IN) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(IN) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: ierr
      CHARACTER(LEN=256) :: dirname,fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(6) :: my_label
      INTEGER :: local_j,global_j
      INTEGER :: im
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! MKDIR 
      !
      dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.restart'
      CALL my_mkdir( dirname )
      !
      CALL start_clock('wstat_restart')
      time_spent(1)=get_clock('wstat_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("summary.xml") , BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         CALL iotk_write_begin( iunout, "R-SUMMARY" )
         CALL iotk_write_dat( iunout, "dav_iter", dav_iter)
         CALL iotk_write_dat( iunout, "notcnv", notcnv)
         CALL iotk_write_dat( iunout, "nbase", nbase)
         CALL iotk_write_dat( iunout, "conv", conv(:))
         CALL iotk_write_end( iunout, "R-SUMMARY"  )
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! CREATE THE EIG FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("eig.xml") , BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         CALL iotk_write_begin( iunout, "RESTART_EIG" )
         CALL iotk_write_dat( iunout, "ev", ev(:))
         CALL iotk_write_dat( iunout, "ew", ew(:))
         CALL iotk_write_end( iunout, "RESTART_EIG"  )
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! CREATE THE HR FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("hr.dat") , BINARY = .TRUE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         IF(me_bgrp==0) CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_write_begin( iunout, "RESTART_HR_"//TRIM(my_label) )
            CALL iotk_write_dat( iunout, "hr", tmp_distr(:,:))
            CALL iotk_write_end( iunout, "RESTART_HR"//TRIM(my_label) )
            !
         END IF
         !
      ENDDO
      !
      IF ( mpime == root ) CALL iotk_close_write( iunout )
      !
      DEALLOCATE( tmp_distr )
      !
      ! CREATE THE VR FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // TRIM("vr.dat") , BINARY = .TRUE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wstat_restart', 'cannot open restart file for writing', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         IF(me_bgrp==0) CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_write_begin( iunout, "RESTART_VR_"//TRIM(my_label) )
            CALL iotk_write_dat( iunout, "vr", tmp_distr(:,:))
            CALL iotk_write_end( iunout, "RESTART_VR"//TRIM(my_label) )
            !
         END IF
         !
      ENDDO
      !
      IF ( mpime == root ) CALL iotk_close_write( iunout )
      !
      DEALLOCATE( tmp_distr )
      !
      ! CREATE THE EIGENVECTOR FILES
      !
      DO local_j=1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j>nbase) CYCLE
         ! 
         fname = TRIM( dirname ) // "/V"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_merge_and_write_G(fname,dvg(:,local_j))
         fname = TRIM( dirname ) // "/N"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_merge_and_write_G(fname,dng(:,local_j))
         !
      ENDDO
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      time_spent(2)=get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout, "(5x, '[I/O] RESTART written in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, '[I/O] In location   : ',a)") TRIM( dirname )  
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_clear( )
      !------------------------------------------------------------------------
      !
      USE mp_world,             ONLY : root,mpime,world_comm
      USE mp,                   ONLY : mp_barrier,mp_bcast
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : n_pdep_basis,n_pdep_eigen,west_prefix
      USE wrappers,             ONLY : f_rmdir
      USE io_files,             ONLY : delete_if_present
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      CHARACTER(LEN=256) :: dirname,fname
      INTEGER :: ierr,ip
      CHARACTER(6) :: my_label
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! ... clear the main restart directory
      !
      dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.restart'
      !
      IF(mpime==root) THEN
         CALL delete_if_present( TRIM( dirname ) // '/' // TRIM( 'summary.xml' ) )
         CALL delete_if_present( TRIM( dirname ) // '/' // TRIM( 'eig.xml' ) )
         CALL delete_if_present( TRIM( dirname ) // '/' // TRIM( 'hr.dat' ) )
         CALL delete_if_present( TRIM( dirname ) // '/' // TRIM( 'vr.dat' ) )
         DO ip=1,n_pdep_basis
            WRITE(my_label,'(i6.6)') ip
            fname="V"//TRIM(ADJUSTL(my_label))//".dat"
            CALL delete_if_present( TRIM( dirname ) // '/' // TRIM( fname ) )
            fname="N"//TRIM(ADJUSTL(my_label))//".dat"
            CALL delete_if_present( TRIM( dirname ) // '/' // TRIM( fname ) )
         ENDDO
         ierr =  f_rmdir( TRIM( dirname ) )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm)
      !
      CALL errore( 'wstat_restart', 'cannot clear restart', ierr )
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_read_real( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,           ONLY : world_comm
      USE mp,                  ONLY : mp_barrier
      USE westcom,             ONLY : n_pdep_eigen,west_prefix,n_pdep_basis
      USE io_global,           ONLY : stdout 
      USE distribution_center, ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter, notcnv, nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      REAL(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      CHARACTER(LEN=256) :: dirname
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: ierr
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock('wstat_restart')
      time_spent(1)=get_clock('wstat_restart')
      !
      dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.restart'
      !
      CALL read_restart1_( dirname, dav_iter, notcnv, nbase )
      !
      CALL read_restart2_( dirname, ew )
      !
      CALL read_restart3d_( dirname, hr_distr, vr_distr )
      !
      CALL read_restart4_( dirname, nbase )
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      time_spent(2)=get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(1/, 5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout, "(5x, '[I/O] RESTART read in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, '[I/O] In location : ',a)") TRIM( dirname )  
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_read_complex( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,           ONLY : world_comm
      USE mp,                  ONLY : mp_barrier
      USE westcom,             ONLY : n_pdep_eigen,west_prefix,n_pdep_basis
      USE io_global,           ONLY : stdout 
      USE distribution_center, ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter, notcnv, nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      CHARACTER(LEN=256) :: dirname
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: ierr
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock('wstat_restart')
      time_spent(1)=get_clock('wstat_restart')
      !
      dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.restart'
      !
      CALL read_restart1_( dirname, dav_iter, notcnv, nbase )
      !
      CALL read_restart2_( dirname, ew )
      !
      CALL read_restart3z_( dirname, hr_distr, vr_distr )
      !
      CALL read_restart4_( dirname, nbase )
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      time_spent(2)=get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(1/, 5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout, "(5x, '[I/O] RESTART read in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, '[I/O] In location : ',a)") TRIM( dirname )  
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart1_( dirname, dav_iter, notcnv, nbase )
      !------------------------------------------------------------------------
      !
      USE westcom,       ONLY : conv,n_pdep_eigen,n_pdep_basis
      USE mp_world,      ONLY : world_comm,mpime,root
      USE mp,            ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !  
      INTEGER,INTENT(OUT) :: dav_iter, notcnv, nbase
      !
      ! Workspace
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER :: ierr,iun
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'summary.xml' ), IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'wstat_restart', 'cannot open restart file for reading', ierr )
      !
      IF ( mpime==root ) THEN
         !
         CALL iotk_scan_begin( iun, "R-SUMMARY" )
         CALL iotk_scan_dat( iun, "dav_iter", dav_iter)
         CALL iotk_scan_dat( iun, "notcnv", notcnv)
         CALL iotk_scan_dat( iun, "nbase", nbase)
         CALL iotk_scan_dat( iun, "conv", conv(1:n_pdep_eigen))
         CALL iotk_scan_end( iun, "R-SUMMARY"  )
         !
         CALL iotk_close_read( iun )
         !
      ENDIF
      !
      CALL mp_bcast( dav_iter, root, world_comm )
      CALL mp_bcast( notcnv, root, world_comm )
      CALL mp_bcast( nbase, root, world_comm )
      CALL mp_bcast( conv, root, world_comm )
      !
    END SUBROUTINE
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart2_( dirname, ew )
      !------------------------------------------------------------------------
      !
      USE westcom,       ONLY : ev,n_pdep_eigen,n_pdep_basis
      USE mp_world,      ONLY : world_comm,mpime,root
      USE mp,            ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !  
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      !
      ! Workspace
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER :: ierr,iun
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'eig.xml' ), IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'wstat_restart', 'cannot open restart file for reading', ierr )
      !
      IF ( mpime==root ) THEN
         !
         CALL iotk_scan_begin( iun, "RESTART_EIG" )
         CALL iotk_scan_dat( iun, "ev", ev(:))
         CALL iotk_scan_dat( iun, "ew", ew(:))
         CALL iotk_scan_end( iun, "RESTART_EIG"  )
         !
         CALL iotk_close_read( iun )
         !
      ENDIF
      !
      CALL mp_bcast( ev, root, world_comm )
      CALL mp_bcast( ew, root, world_comm )
      !
    END SUBROUTINE
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart3d_( dirname, hr_distr, vr_distr )
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : n_pdep_eigen,n_pdep_basis
      USE mp_world,             ONLY : world_comm,mpime,root
      USE mp,                   ONLY : mp_bcast,mp_get
      USE distribution_center,  ONLY : pert
      USE mp_global,            ONLY : nimage,me_bgrp,inter_image_comm,intra_image_comm,my_image_id
      !
      IMPLICIT NONE
      !
      ! I/O
      !  
      REAL(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      !
      ! Workspace
      !
      INTEGER :: ierr,iun
      INTEGER :: im
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      CHARACTER(6) :: my_label
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'hr.dat' ), BINARY = .TRUE., IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'wstat_restart', 'cannot open restart file for reading', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_scan_begin( iun, "RESTART_HR_"//TRIM(my_label) )
            CALL iotk_scan_dat( iun, "hr", tmp_distr(:,:))
            CALL iotk_scan_end( iun, "RESTART_HR"//TRIM(my_label) )
            !
         END IF
         !
         IF(me_bgrp==0) CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !CALL mp_bcast( tmp_distr, root, world_comm )
         !IF( my_image_id == im ) hr_distr = tmp_distr
         !
      ENDDO
      DEALLOCATE(tmp_distr)
      !
      IF ( mpime==root ) CALL iotk_close_read( iun )
      !
      CALL mp_bcast( hr_distr, 0, intra_image_comm )
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'vr.dat' ), BINARY = .TRUE., IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'wstat_restart', 'cannot open restart file for reading', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_scan_begin( iun, "RESTART_VR_"//TRIM(my_label) )
            CALL iotk_scan_dat( iun, "vr", tmp_distr(:,:))
            CALL iotk_scan_end( iun, "RESTART_VR"//TRIM(my_label) )
            !
         END IF
         !
         IF(me_bgrp==0) CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !CALL mp_bcast( tmp_distr, root, world_comm )
         !IF( my_image_id == im ) vr_distr = tmp_distr
         !
      ENDDO
      DEALLOCATE(tmp_distr)
      !
      IF ( mpime==root ) CALL iotk_close_read( iun )
      !
      CALL mp_bcast( vr_distr, 0, intra_image_comm )
      !
    END SUBROUTINE
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart3z_( dirname, hr_distr, vr_distr )
      !------------------------------------------------------------------------
      !
      USE westcom,             ONLY : n_pdep_eigen,n_pdep_basis
      USE mp_world,            ONLY : world_comm,mpime,root
      USE mp,                  ONLY : mp_bcast,mp_get
      USE distribution_center, ONLY : pert
      USE mp_global,           ONLY : nimage,me_bgrp,inter_image_comm,intra_image_comm,my_image_id
      !
      IMPLICIT NONE
      !
      ! I/O
      !  
      COMPLEX(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      !
      ! Workspace
      !
      INTEGER :: ierr,iun
      INTEGER :: im
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      CHARACTER(6) :: my_label
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'hr.dat' ), BINARY = .TRUE., IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'wstat_restart', 'cannot open restart file for reading', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_scan_begin( iun, "RESTART_HR_"//TRIM(my_label) )
            CALL iotk_scan_dat( iun, "hr", tmp_distr(:,:))
            CALL iotk_scan_end( iun, "RESTART_HR"//TRIM(my_label) )
            !
         END IF
         !
         IF(me_bgrp==0) CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !CALL mp_bcast( tmp_distr, root, world_comm )
         !IF( my_image_id == im ) hr_distr = tmp_distr
         !
      ENDDO
      DEALLOCATE(tmp_distr)
      !
      IF ( mpime==root ) CALL iotk_close_read( iun )
      !
      CALL mp_bcast( hr_distr, 0, intra_image_comm )
      !
      ierr = 0
      !
      IF ( mpime==root ) THEN
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( dirname ) // '/' // TRIM( 'vr.dat' ), BINARY = .TRUE., IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'wstat_restart', 'cannot open restart file for reading', ierr )
      !
      ALLOCATE( tmp_distr(n_pdep_basis,pert%nlocx) )
      DO im = 0, nimage-1
         !
         WRITE(my_label,'(i6.6)') im
         !
         IF ( mpime == root ) THEN  
            !
            CALL iotk_scan_begin( iun, "RESTART_VR_"//TRIM(my_label) )
            CALL iotk_scan_dat( iun, "vr", tmp_distr(:,:))
            CALL iotk_scan_end( iun, "RESTART_VR"//TRIM(my_label) )
            !
         END IF
         !
         IF(me_bgrp==0) CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !CALL mp_bcast( tmp_distr, root, world_comm )
         !IF( my_image_id == im ) vr_distr = tmp_distr
         !
      ENDDO
      DEALLOCATE(tmp_distr)
      !
      IF ( mpime==root ) CALL iotk_close_read( iun )
      !
      CALL mp_bcast( vr_distr, 0, intra_image_comm )
      !
    END SUBROUTINE
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart4_( dirname, nbase )
      !------------------------------------------------------------------------
      !
      USE westcom,             ONLY : dvg,dng,npwq0x
      USE mp_global,           ONLY : my_image_id
      USE pdep_io,             ONLY : pdep_read_G_and_distribute
      USE distribution_center, ONLY : pert
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER, INTENT(IN) :: nbase
      !
      INTEGER :: global_j, local_j, group_j
      INTEGER :: npw_g,ierr
      CHARACTER(6) :: my_label
      CHARACTER(LEN=256) :: fname
      INTEGER :: iun
      !
      IF(.NOT.ALLOCATED(dvg)) ALLOCATE(dvg(npwq0x,pert%nlocx))
      IF(.NOT.ALLOCATED(dng)) ALLOCATE(dng(npwq0x,pert%nlocx))
      dvg = 0._DP
      dng = 0._DP
      !
      DO local_j=1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j>nbase) CYCLE
         ! 
         fname = TRIM( dirname ) // "/V"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_read_G_and_distribute(fname,dvg(:,local_j))
         fname = TRIM( dirname ) // "/N"//TRIM(ADJUSTL(my_label))//".dat"
         CALL pdep_read_G_and_distribute(fname,dng(:,local_j))
         !
      ENDDO
      !
    END SUBROUTINE
    !
END MODULE
