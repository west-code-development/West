!
! Copyright (C) 2015-2024 M. Govoni
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
MODULE davidson_restart
  !----------------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE json_module, ONLY : json_file
  !
  IMPLICIT NONE
  !
  INTERFACE davidson_restart_write
     MODULE PROCEDURE davidson_restart_write_real, davidson_restart_write_complex
  END INTERFACE
  !
  INTERFACE davidson_restart_read
     MODULE PROCEDURE davidson_restart_read_real, davidson_restart_read_complex
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE davidson_restart_write_real(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_get
      USE mp_world,             ONLY : mpime,root
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_pool_comm,&
                                     & my_pool_id,inter_bgrp_comm,my_bgrp_id,me_bgrp
      USE io_global,            ONLY : stdout
      USE pwcom,                ONLY : npwx
      USE westcom,              ONLY : n_pdep_basis,ev,conv,dvg,dng,dvg_exc,dng_exc,nbndval0x,&
                                     & n_trunc_bands,wstat_restart_dir,wbse_restart_dir
      USE pdep_io,              ONLY : pdep_merge_and_write_G
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      USE west_mp,              ONLY : west_mp_root_sum
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
      CHARACTER(LEN=20) :: which
      CHARACTER(LEN=512) :: dirname
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(6) :: my_label
      INTEGER :: lbnd,ibnd,iks,iks_g
      INTEGER :: local_j,global_j
      INTEGER :: im
      LOGICAL :: l_bse
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      COMPLEX(DP),ALLOCATABLE :: tmp_exc(:,:,:)
      !
      TYPE(json_file) :: json
      INTEGER :: iun
      !
      IF(ALLOCATED(dvg_exc)) THEN
         !
         l_bse = .TRUE.
         which = 'wbse_restart'
         dirname = wbse_restart_dir
         !
         ALLOCATE(tmp_exc(npwx,nbndval0x-n_trunc_bands,kpt_pool%nglob))
         !
      ELSE
         !
         l_bse = .FALSE.
         which = 'wstat_restart'
         dirname = wstat_restart_dir
         !
      ENDIF
      !
      ! MKDIR
      !
      CALL my_mkdir(TRIM(dirname))
      !
      CALL start_clock(TRIM(which))
      time_spent(1) = get_clock(TRIM(which))
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
         OPEN(NEWUNIT=iun,FILE=TRIM(dirname)//'/summary.json')
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
         OPEN(NEWUNIT=iun,FILE=TRIM(dirname)//'/hr_vr.dat',FORM='unformatted')
         !
      ENDIF
      !
      DO im = 0,nimage-1
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
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
         IF(l_bse) THEN
            !
            tmp_exc(:,:,:) = (0._DP,0._DP)
            !
            DO iks = 1,kpt_pool%nloc
               iks_g = kpt_pool%l2g(iks)
               DO lbnd = 1,band_group%nloc
                  ibnd = band_group%l2g(lbnd)
                  tmp_exc(:,ibnd,iks_g) = dvg_exc(:,lbnd,iks,local_j)
               ENDDO
            ENDDO
            !
            CALL west_mp_root_sum(tmp_exc,0,inter_pool_comm)
            CALL west_mp_root_sum(tmp_exc,0,inter_bgrp_comm)
            !
            IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
               fname = TRIM(wbse_restart_dir)//'/V'//my_label//'.dat'
               CALL plep_merge_and_write_G(fname,tmp_exc)
            ENDIF
            !
            tmp_exc(:,:,:) = (0._DP,0._DP)
            !
            DO iks = 1,kpt_pool%nloc
               iks_g = kpt_pool%l2g(iks)
               DO lbnd = 1,band_group%nloc
                  ibnd = band_group%l2g(lbnd)
                  tmp_exc(:,ibnd,iks_g) = dng_exc(:,lbnd,iks,local_j)
               ENDDO
            ENDDO
            !
            CALL west_mp_root_sum(tmp_exc,0,inter_pool_comm)
            CALL west_mp_root_sum(tmp_exc,0,inter_bgrp_comm)
            !
            IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
               fname = TRIM(wbse_restart_dir)//'/N'//my_label//'.dat'
               CALL plep_merge_and_write_G(fname,tmp_exc)
            ENDIF
            !
         ELSE
            !
            IF(my_bgrp_id == 0) THEN
               fname = TRIM(dirname)//'/V'//my_label//'.dat'
               CALL pdep_merge_and_write_G(fname,dvg(:,local_j))
            ENDIF
            !
            IF(my_bgrp_id == 0) THEN
               fname = TRIM(dirname)//'/N'//my_label//'.dat'
               CALL pdep_merge_and_write_G(fname,dng(:,local_j))
            ENDIF
            !
         ENDIF
         !
      ENDDO
      !
      IF(l_bse) THEN
         DEALLOCATE(tmp_exc)
      ENDIF
      !
      time_spent(2) = get_clock(TRIM(which))
      CALL stop_clock(TRIM(which))
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART written in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location   : ",a)') TRIM(dirname)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE davidson_restart_write_complex(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr,lastdone_iq)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_get
      USE mp_world,             ONLY : mpime,root
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,my_pool_id,my_bgrp_id,&
                                     & me_bgrp
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : n_pdep_basis,ev,conv,dvg,dng,wstat_restart_dir
      USE pdep_io,              ONLY : pdep_merge_and_write_G
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(IN) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(IN) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: vr_distr(n_pdep_basis,pert%nlocx)
      INTEGER,INTENT(IN),OPTIONAL :: lastdone_iq
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
      ! MKDIR
      !
      CALL my_mkdir(wstat_restart_dir)
      !
      CALL start_clock('wstat_restart')
      time_spent(1) = get_clock('wstat_restart')
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
         IF(PRESENT(lastdone_iq)) THEN
            CALL json%add('lastdone_iq',lastdone_iq)
         ENDIF
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/summary.json')
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
         OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/hr_vr.dat',FORM='unformatted')
         !
      ENDIF
      !
      DO im = 0,nimage-1
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
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
         IF(my_bgrp_id == 0) THEN
            fname = TRIM(wstat_restart_dir)//'/V'//my_label//'.dat'
            !
            IF(PRESENT(lastdone_iq)) THEN
               CALL pdep_merge_and_write_G(fname,dvg(:,local_j),lastdone_iq)
            ELSE
               CALL pdep_merge_and_write_G(fname,dvg(:,local_j))
            ENDIF
         ENDIF
         !
         IF(my_bgrp_id == 0) THEN
            fname = TRIM(wstat_restart_dir)//'/N'//my_label//'.dat'
            !
            IF(PRESENT(lastdone_iq)) THEN
               CALL pdep_merge_and_write_G(fname,dng(:,local_j),lastdone_iq)
            ELSE
               CALL pdep_merge_and_write_G(fname,dng(:,local_j))
            ENDIF
         ENDIF
         !
      ENDDO
      !
      time_spent(2) = get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART written in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location   : ",a)') TRIM(wstat_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE davidson_restart_clear()
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : root,mpime,world_comm
      USE westcom,              ONLY : dvg_exc,n_pdep_basis,wstat_restart_dir,wbse_restart_dir
      USE clib_wrappers,        ONLY : f_rmdir
      USE west_io,              ONLY : remove_if_present
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      CHARACTER(LEN=20) :: which
      CHARACTER(LEN=512) :: dirname
      CHARACTER(LEN=512) :: fname
      INTEGER :: ierr,ip
      CHARACTER(6) :: my_label
      !
      IF(ALLOCATED(dvg_exc)) THEN
         which = 'wbse_restart'
         dirname = wbse_restart_dir
      ELSE
         which = 'wstat_restart'
         dirname = wstat_restart_dir
      ENDIF
      !
      ! ... clear the main restart directory
      !
      IF(mpime == root) THEN
         !
         CALL remove_if_present(TRIM(dirname)//'/summary.json')
         CALL remove_if_present(TRIM(dirname)//'/hr_vr.dat')
         !
         DO ip = 1,n_pdep_basis
            WRITE(my_label,'(i6.6)') ip
            fname = 'V'//my_label//'.dat'
            CALL remove_if_present(TRIM(dirname)//'/'//TRIM(fname))
            fname = 'N'//my_label//'.dat'
            CALL remove_if_present(TRIM(dirname)//'/'//TRIM(fname))
         ENDDO
         !
         ierr = f_rmdir(TRIM(dirname))
         !
      ENDIF
      !
      CALL mp_bcast(ierr,root,world_comm)
      !
      CALL errore(TRIM(which),'cannot clear restart',ierr)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE davidson_restart_read_real(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : dvg_exc,n_pdep_basis,wstat_restart_dir,wbse_restart_dir
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
      CHARACTER(LEN=20) :: which
      CHARACTER(LEN=512) :: dirname
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      IF(ALLOCATED(dvg_exc)) THEN
         which = 'wbse_restart'
         dirname = wbse_restart_dir
      ELSE
         which = 'wstat_restart'
         dirname = wstat_restart_dir
      ENDIF
      !
      CALL start_clock(TRIM(which))
      time_spent(1) = get_clock(TRIM(which))
      !
      CALL read_restart12_(dav_iter,notcnv,nbase,ew)
      !
      CALL read_restart3d_(hr_distr,vr_distr)
      !
      CALL read_restart4_(nbase)
      !
      time_spent(2) = get_clock(TRIM(which))
      CALL stop_clock(TRIM(which))
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(dirname)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE davidson_restart_read_complex(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr,lastdone_iq,iq)
      !------------------------------------------------------------------------
      !
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : n_pdep_basis,wstat_restart_dir
      USE distribution_center,  ONLY : pert
      USE types_bz_grid,        ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      INTEGER,INTENT(OUT) :: lastdone_iq
      INTEGER,INTENT(IN) :: iq
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: ipol
      !
      CALL start_clock('wstat_restart')
      time_spent(1) = get_clock('wstat_restart')
      !
      CALL read_restart12_(dav_iter,notcnv,nbase,ew,lastdone_iq)
      !
      CALL read_restart3z_(hr_distr,vr_distr)
      !
      CALL read_restart4_(nbase,lastdone_iq)
      !
      time_spent(2) = get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      IF(iq == lastdone_iq) THEN
         WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------------------")')
         WRITE(stdout,'(5x,"[I/O] Restarting from q(",i5,") = (",3f12.7,")")') &
              lastdone_iq,(q_grid%p_cryst(ipol,lastdone_iq),ipol=1,3)
         WRITE(stdout,'(5x,"[I/O] RESTART read in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
         WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wstat_restart_dir)
         WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------------------")')
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart12_(dav_iter,notcnv,nbase,ew,iq)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm,mpime,root
      USE westcom,              ONLY : dvg_exc,conv,n_pdep_eigen,n_pdep_basis,ev,wstat_restart_dir,&
                                     & wbse_restart_dir
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      INTEGER,INTENT(OUT),OPTIONAL :: iq
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: dirname
      LOGICAL :: found
      TYPE(json_file) :: json
      REAL(DP),ALLOCATABLE :: rvals(:)
      LOGICAL,ALLOCATABLE :: lvals(:)
      INTEGER :: ival
      !
      IF(mpime == root) THEN
         !
         IF(ALLOCATED(dvg_exc)) THEN
            dirname = wbse_restart_dir
         ELSE
            dirname = wstat_restart_dir
         ENDIF
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(dirname)//'/summary.json')
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
         IF(PRESENT(iq)) THEN
            CALL json%get('lastdone_iq',ival,found)
            IF(found) iq = ival
         ENDIF
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
      IF(PRESENT(iq)) THEN
         CALL mp_bcast(iq,root,world_comm)
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart3d_(hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_get
      USE mp_world,             ONLY : mpime,root
      USE mp_global,            ONLY : inter_image_comm,intra_image_comm,nimage,my_image_id,&
                                     & my_pool_id,my_bgrp_id,me_bgrp
      USE westcom,              ONLY : dvg_exc,n_pdep_basis,wstat_restart_dir,wbse_restart_dir
      USE distribution_center,  ONLY : pert
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
      CHARACTER(LEN=512) :: dirname
      INTEGER :: iun
      INTEGER :: im
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) THEN
         !
         IF(ALLOCATED(dvg_exc)) THEN
            dirname = wbse_restart_dir
         ELSE
            dirname = wstat_restart_dir
         ENDIF
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(dirname)//'/hr_vr.dat',FORM='unformatted')
         !
      ENDIF
      !
      DO im = 0,nimage-1
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
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
      USE mp,                   ONLY : mp_bcast,mp_get
      USE mp_world,             ONLY : mpime,root
      USE mp_global,            ONLY : inter_image_comm,intra_image_comm,nimage,my_image_id,&
                                     & my_pool_id,my_bgrp_id,me_bgrp
      USE westcom,              ONLY : n_pdep_basis,wstat_restart_dir
      USE distribution_center,  ONLY : pert
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
      IF(mpime == root) OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/hr_vr.dat',FORM='unformatted')
      !
      DO im = 0,nimage-1
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
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
    SUBROUTINE read_restart4_(nbase,iq)
      !------------------------------------------------------------------------
      !
      USE pwcom,                ONLY : npwx
      USE westcom,              ONLY : dvg,dng,dvg_exc,dng_exc,nbndval0x,n_trunc_bands,&
                                     & wstat_restart_dir,wbse_restart_dir
      USE pdep_io,              ONLY : pdep_read_G_and_distribute
      USE plep_io,              ONLY : plep_read_G_and_distribute
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nbase
      INTEGER,INTENT(IN),OPTIONAL :: iq
      !
      ! Workspace
      !
      LOGICAL :: l_bse
      INTEGER :: lbnd,ibnd,iks,iks_g
      INTEGER :: global_j,local_j
      CHARACTER(6) :: my_label
      CHARACTER(LEN=512) :: fname
      COMPLEX(DP),ALLOCATABLE :: tmp_exc(:,:,:)
      !
      IF(ALLOCATED(dvg_exc)) THEN
         !
         l_bse = .TRUE.
         dvg_exc(:,:,:,:) = (0._DP,0._DP)
         dng_exc(:,:,:,:) = (0._DP,0._DP)
         !
         ALLOCATE(tmp_exc(npwx,nbndval0x-n_trunc_bands,kpt_pool%nglob))
         !
      ELSE
         !
         l_bse = .FALSE.
         dvg(:,:) = (0._DP,0._DP)
         dng(:,:) = (0._DP,0._DP)
         !
      ENDIF
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j > nbase) CYCLE
         !
         IF(l_bse) THEN
            !
            fname = TRIM(wbse_restart_dir)//'/V'//my_label//'.dat'
            CALL plep_read_G_and_distribute(fname,tmp_exc)
            !
            DO iks = 1,kpt_pool%nloc
               iks_g = kpt_pool%l2g(iks)
               DO lbnd = 1,band_group%nloc
                  ibnd = band_group%l2g(lbnd)
                  dvg_exc(:,lbnd,iks,local_j) = tmp_exc(:,ibnd,iks_g)
               ENDDO
            ENDDO
            !
            fname = TRIM(wbse_restart_dir)//'/N'//my_label//'.dat'
            CALL plep_read_G_and_distribute(fname,tmp_exc)
            !
            DO iks = 1,kpt_pool%nloc
               iks_g = kpt_pool%l2g(iks)
               DO lbnd = 1,band_group%nloc
                  ibnd = band_group%l2g(lbnd)
                  dng_exc(:,lbnd,iks,local_j) = tmp_exc(:,ibnd,iks_g)
               ENDDO
            ENDDO
            !
         ELSE
            !
            fname = TRIM(wstat_restart_dir)//'/V'//my_label//'.dat'
            IF(PRESENT(iq)) THEN
               CALL pdep_read_G_and_distribute(fname,dvg(:,local_j),iq)
            ELSE
               CALL pdep_read_G_and_distribute(fname,dvg(:,local_j))
            ENDIF
            !
            fname = TRIM(wstat_restart_dir)//'/N'//my_label//'.dat'
            IF(PRESENT(iq)) THEN
               CALL pdep_read_G_and_distribute(fname,dng(:,local_j),iq)
            ELSE
               CALL pdep_read_G_and_distribute(fname,dng(:,local_j))
            ENDIF
            !
         ENDIF
         !
      ENDDO
      !
      IF(l_bse) THEN
         DEALLOCATE(tmp_exc)
      ENDIF
      !
    END SUBROUTINE
    !
END MODULE
