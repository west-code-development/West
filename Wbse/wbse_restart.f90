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
MODULE wbse_restart
  !----------------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE json_module, ONLY : json_file
  !
  IMPLICIT NONE
  !
  INTERFACE wbse_restart_write
     MODULE PROCEDURE wbse_restart_write_real, wbse_restart_write_complex
  END INTERFACE
  !
  INTERFACE wbse_restart_read
     MODULE PROCEDURE wbse_restart_read_real, wbse_restart_read_complex
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE wbse_restart_write_real(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,me_bgrp,inter_image_comm,nimage
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : n_pdep_basis,ev,conv,dvg_exc,dng_exc,wbse_restart_dir
      USE mp,                   ONLY : mp_barrier,mp_get
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(IN) :: ew(n_pdep_basis)
      REAL(DP),INTENT(IN) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(IN) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(6) :: my_label
      INTEGER :: local_j,global_j
      INTEGER :: im
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      TYPE(json_file) :: json
      INTEGER :: iun
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! MKDIR
      !
      CALL my_mkdir(TRIM(wbse_restart_dir))
      !
      CALL start_clock('wbse_restart')
      time_spent(1) = get_clock('wbse_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add('dav_iter',dav_iter)
         CALL json%add('notcnv',notcnv)
         CALL json%add('nbase',nbase)
         CALL json%add('conv',conv(:))
         CALL json%add('ev',ev(:))
         CALL json%add('ew',ew(:))
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wbse_restart_dir)//'/'//TRIM('summary.json'))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! CREATE THE HR, VR FILE
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wbse_restart_dir)//'/hr_vr.dat',FORM='unformatted')
         !
      ENDIF
      !
      DO im = 0,nimage-1
         !
         IF(me_bgrp == 0) CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
         IF(me_bgrp == 0) CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      ! CREATE THE EIGENVECTOR FILES
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j > nbase) CYCLE
         !
         fname = TRIM(wbse_restart_dir)//'/V'//my_label//'.dat'
         CALL plep_merge_and_write_G(fname,dvg_exc(:,:,:,local_j))
         fname = TRIM(wbse_restart_dir)//'/N'//my_label//'.dat'
         CALL plep_merge_and_write_G(fname,dng_exc(:,:,:,local_j))
         !
      ENDDO
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      time_spent(2) = get_clock('wbse_restart')
      CALL stop_clock('wbse_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location   : ",a)') TRIM(wbse_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wbse_restart_write_complex(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,me_bgrp,inter_image_comm,nimage
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : n_pdep_basis,ev,conv,dvg_exc,dng_exc,wbse_restart_dir
      USE mp,                   ONLY : mp_barrier,mp_get
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN)  :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(IN) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(IN) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(6) :: my_label
      INTEGER :: local_j,global_j
      INTEGER :: im
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      TYPE(json_file) :: json
      INTEGER :: iun
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! MKDIR
      !
      CALL my_mkdir(TRIM(wbse_restart_dir))
      !
      CALL start_clock('wbse_restart')
      time_spent(1) = get_clock('wbse_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add('dav_iter',dav_iter)
         CALL json%add('notcnv',notcnv)
         CALL json%add('nbase',nbase)
         CALL json%add('conv',conv(:))
         CALL json%add('ev',ev(:))
         CALL json%add('ew',ew(:))
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wbse_restart_dir)//'/'//TRIM('summary.json'))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! CREATE THE HR, VR FILE
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wbse_restart_dir)//'/'//TRIM('hr_vr.dat'),FORM='unformatted')
         !
      ENDIF
      !
      DO im = 0,nimage-1
         !
         IF(me_bgrp == 0) CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
         IF(me_bgrp == 0) CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      ! CREATE THE EIGENVECTOR FILES
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j > nbase) CYCLE
         !
         fname = TRIM(wbse_restart_dir)//'/V'//my_label//'.dat'
         CALL plep_merge_and_write_G(fname,dvg_exc(:,:,:,local_j))
         fname = TRIM(wbse_restart_dir)//'/N'//my_label//'.dat'
         CALL plep_merge_and_write_G(fname,dng_exc(:,:,:,local_j))
         !
      ENDDO
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      time_spent(2) = get_clock('wbse_restart')
      CALL stop_clock('wbse_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location   : ",a)') TRIM(wbse_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wbse_restart_clear()
      !------------------------------------------------------------------------
      !
      USE mp_world,             ONLY : root,mpime,world_comm
      USE mp,                   ONLY : mp_barrier,mp_bcast
      USE westcom,              ONLY : n_pdep_basis,wbse_restart_dir
      USE clib_wrappers,        ONLY : f_rmdir
      USE west_io,              ONLY : remove_if_present
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      INTEGER :: ierr,ip
      CHARACTER(6) :: my_label
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! ... clear the main restart directory
      !
      IF(mpime == root) THEN
         CALL remove_if_present(TRIM(wbse_restart_dir)//'/summary.json')
         CALL remove_if_present(TRIM(wbse_restart_dir)//'/hr_vr.dat')
         DO ip = 1,n_pdep_basis
            WRITE(my_label,'(i6.6)') ip
            fname = 'V'//my_label//'.dat'
            CALL remove_if_present(TRIM(wbse_restart_dir)//'/'//TRIM(fname))
            fname = 'N'//my_label//'.dat'
            CALL remove_if_present(TRIM(wbse_restart_dir)//'/'//TRIM(fname))
         ENDDO
         ierr = f_rmdir(TRIM(wbse_restart_dir))
      ENDIF
      !
      CALL mp_bcast(ierr,root,world_comm)
      !
      CALL errore('wbse_restart','cannot clear restart',ierr)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wbse_restart_read_real(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : world_comm
      USE mp,                   ONLY : mp_barrier
      USE westcom,              ONLY : n_pdep_basis,wbse_restart_dir
      USE io_global,            ONLY : stdout
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      REAL(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('wbse_restart')
      time_spent(1) = get_clock('wbse_restart')
      !
      CALL read_restart12_(dav_iter,notcnv,nbase,ew)
      !
      CALL read_restart3d_(hr_distr,vr_distr)
      !
      CALL read_restart4_(nbase)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      time_spent(2) = get_clock('wbse_restart')
      CALL stop_clock('wbse_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wbse_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wbse_restart_read_complex(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : world_comm
      USE mp,                   ONLY : mp_barrier
      USE westcom,              ONLY : n_pdep_basis,wbse_restart_dir
      USE io_global,            ONLY : stdout
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('wbse_restart')
      time_spent(1) = get_clock('wbse_restart')
      !
      CALL read_restart12_(dav_iter,notcnv,nbase,ew)
      !
      CALL read_restart3z_(hr_distr,vr_distr)
      !
      CALL read_restart4_(nbase)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      time_spent(2) = get_clock('wbse_restart')
      CALL stop_clock('wbse_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wbse_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart12_(dav_iter,notcnv,nbase,ew)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : conv,n_pdep_eigen,n_pdep_basis,wbse_restart_dir,ev
      USE mp_world,             ONLY : world_comm,mpime,root
      USE mp,                   ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      !
      ! Workspace
      !
      LOGICAL :: found
      TYPE(json_file) :: json
      REAL(DP),ALLOCATABLE :: rvals(:)
      LOGICAL,ALLOCATABLE :: lvals(:)
      INTEGER :: ival
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(wbse_restart_dir)//'/'//TRIM('summary.json'))
         !
         CALL json%get('dav_iter',ival,found)
         IF(found) dav_iter = ival
         CALL json%get('notcnv',ival,found)
         IF(found) notcnv = ival
         CALL json%get('nbase',ival,found)
         IF(found) nbase = ival
         CALL json%get('conv',lvals,found)
         IF(found) conv(1:n_pdep_eigen) = lvals(1:n_pdep_eigen)
         CALL json%get('ev',rvals,found)
         IF(found) ev(:) = rvals(:)
         CALL json%get('ew',rvals,found)
         IF(found) ew(1:n_pdep_basis) = rvals(1:n_pdep_basis)
         !
         CALL json%destroy()
         !
      ENDIF
      !
      CALL mp_bcast(dav_iter,root,world_comm)
      CALL mp_bcast(notcnv,root,world_comm)
      CALL mp_bcast(nbase,root,world_comm)
      CALL mp_bcast(conv,root,world_comm)
      !
      CALL mp_bcast(ev,root,world_comm)
      CALL mp_bcast(ew,root,world_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart3d_(hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : n_pdep_basis,wbse_restart_dir
      USE mp_world,             ONLY : mpime,root
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
      !
      ! Workspace
      !
      INTEGER :: iun
      INTEGER :: im
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) OPEN(NEWUNIT=iun,FILE=TRIM(wbse_restart_dir)//'/'//TRIM('hr_vr.dat'),FORM='unformatted')
      !
      DO im = 0,nimage-1
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0) CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0) CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      CALL mp_bcast(hr_distr,0,intra_image_comm)
      CALL mp_bcast(vr_distr,0,intra_image_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart3z_(hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : n_pdep_basis,wbse_restart_dir
      USE mp_world,             ONLY : mpime,root
      USE mp,                   ONLY : mp_bcast,mp_get
      USE distribution_center,  ONLY : pert
      USE mp_global,            ONLY : nimage,me_bgrp,inter_image_comm,intra_image_comm,my_image_id
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: iun
      INTEGER :: im
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) OPEN(NEWUNIT=iun,FILE=TRIM(wbse_restart_dir)//'/'//TRIM('hr_vr.dat'),FORM='unformatted')
      !
      DO im = 0,nimage-1
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0) CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0) CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      CALL mp_bcast(hr_distr,0,intra_image_comm)
      CALL mp_bcast(vr_distr,0,intra_image_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart4_(nbase)
      !------------------------------------------------------------------------
      !
      USE pwcom,                ONLY : nks
      USE westcom,              ONLY : npwqx,dvg_exc,dng_exc,nbndval0x,wbse_restart_dir
      USE plep_io,              ONLY : plep_read_G_and_distribute
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: nbase
      !
      INTEGER :: global_j,local_j
      CHARACTER(6) :: my_label
      CHARACTER(LEN=512) :: fname
      !
      IF(.NOT. ALLOCATED(dvg_exc)) ALLOCATE(dvg_exc(npwqx,nbndval0x,nks,pert%nlocx))
      IF(.NOT. ALLOCATED(dng_exc)) ALLOCATE(dng_exc(npwqx,nbndval0x,nks,pert%nlocx))
      dvg_exc = 0._DP
      dng_exc = 0._DP
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j > nbase) CYCLE
         !
         fname = TRIM(wbse_restart_dir)//'/V'//my_label//'.dat'
         CALL plep_read_G_and_distribute(fname,dvg_exc(:,:,:,local_j))
         fname = TRIM(wbse_restart_dir)//'/N'//my_label//'.dat'
         CALL plep_read_G_and_distribute(fname,dng_exc(:,:,:,local_j))
         !
      ENDDO
      !
    END SUBROUTINE
    !
END MODULE
