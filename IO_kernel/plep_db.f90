!
! Copyright (C) 2015-2023 M. Govoni
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
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    ! *****************************
    ! PDEP WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE plep_db_write()
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp_global,            ONLY : inter_pool_comm,my_pool_id,inter_bgrp_comm,my_bgrp_id
      USE io_global,            ONLY : stdout
      USE pwcom,                ONLY : npwx
      USE westcom,              ONLY : n_pdep_eigen,ev,wbse_save_dir,dvg_exc,nbndval0x,n_trunc_bands
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      USE json_module,          ONLY : json_file
      USE west_mp,              ONLY : west_mp_root_sum
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: lbnd,ibnd,iks,iks_g
      INTEGER :: iun,global_j,local_j
      CHARACTER(LEN=6) :: label_j
      CHARACTER(LEN=256) :: fname
      TYPE(json_file) :: json
      LOGICAL :: lexists
      COMPLEX(DP),ALLOCATABLE :: dvg_tmp(:,:,:)
      !
      ! MPI barrier
      !
      CALL mp_barrier(world_comm)
      !
      ! Timing
      !
      CALL start_clock('plep_db')
      time_spent(1) = get_clock('plep_db')
      !
      fname = TRIM(wbse_save_dir)//'/summary.json'
      !
      ! Create summary file if it does not exist
      !
      IF(mpime == root) THEN
         !
         INQUIRE(FILE=TRIM(fname),EXIST=lexists)
         IF(.NOT. lexists) THEN
           CALL json%initialize()
           CALL json%add('plep.n_plep_eigen',n_pdep_eigen)
           !
           OPEN(NEWUNIT=iun,FILE=TRIM(fname))
           CALL json%print(iun)
           CLOSE(iun)
           !
           CALL json%destroy()
         ENDIF
         !
      ENDIF
      !
      ! Update summary file with current structure
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(fname))
         !
         CALL json%add('plep.eigenval',ev(1:n_pdep_eigen))
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! Dump eigenvectors
      !
      ALLOCATE(dvg_tmp(npwx,nbndval0x-n_trunc_bands,kpt_pool%nglob))
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         IF(global_j > n_pdep_eigen) CYCLE
         !
         WRITE(label_j,'(i6.6)') global_j
         !
         dvg_tmp(:,:,:) = (0._DP,0._DP)
         !
         DO iks = 1,kpt_pool%nloc
            iks_g = kpt_pool%l2g(iks)
            DO lbnd = 1,band_group%nloc
               ibnd = band_group%l2g(lbnd)
               dvg_tmp(:,ibnd,iks_g) = dvg_exc(:,lbnd,iks,local_j)
            ENDDO
         ENDDO
         !
         CALL west_mp_root_sum(dvg_tmp,0,inter_pool_comm)
         CALL west_mp_root_sum(dvg_tmp,0,inter_bgrp_comm)
         !
         IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
            fname = TRIM(wbse_save_dir)//'/E'//label_j//'.dat'
            CALL plep_merge_and_write_G(TRIM(fname),dvg_tmp)
         ENDIF
         !
      ENDDO
      !
      DEALLOCATE(dvg_tmp)
      !
      ! MPI barrier
      !
      CALL mp_barrier(world_comm)
      !
      ! Timing
      !
      time_spent(2) = get_clock('plep_db')
      CALL stop_clock('plep_db')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"SAVE written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wbse_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    ! *****************************
    ! PLEP READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE plep_db_read(nglob_to_be_read)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE pwcom,                ONLY : npwx
      USE westcom,              ONLY : n_pdep_eigen,ev,wbse_save_dir,dvg_exc,nbndval0x,n_trunc_bands
      USE plep_io,              ONLY : plep_read_G_and_distribute
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      USE json_module,          ONLY : json_file
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nglob_to_be_read
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: lbnd,ibnd,iks,iks_g
      INTEGER :: n_eigen_to_get
      INTEGER :: tmp_n_pdep_eigen
      INTEGER :: global_j,local_j
      CHARACTER(LEN=6) :: label_j
      REAL(DP),ALLOCATABLE :: tmp_ev(:)
      TYPE(json_file) :: json
      CHARACTER(LEN=256) :: fname
      COMPLEX(DP),ALLOCATABLE :: dvg_tmp(:,:,:)
      !
      ! MPI barrier
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('plep_db')
      !
      ! Timing
      !
      time_spent(1) = get_clock('plep_db')
      !
      ! 1) READ THE INPUT FILE
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(wbse_save_dir)//'/summary.json')
         IF(json%failed()) THEN
            CALL errore('plep_db_read','Cannot open file: '//TRIM(wbse_save_dir)//'/summary.json',1)
         ENDIF
         !
         CALL json%get('plep.eigenval',tmp_ev)
         tmp_n_pdep_eigen = SIZE(tmp_ev,1)
         !
         CALL json%destroy()
         !
      ENDIF
      !
      CALL mp_bcast(tmp_n_pdep_eigen,root,world_comm)
      !
      ! In case nglob_to_be_read is 0, overwrite it with the read value
      !
      IF(nglob_to_be_read == 0) THEN
         n_eigen_to_get = tmp_n_pdep_eigen
         n_pdep_eigen = tmp_n_pdep_eigen
      ELSE
         n_eigen_to_get = MIN(tmp_n_pdep_eigen,nglob_to_be_read)
      ENDIF
      !
      ! 2) READ THE EIGENVALUES FILE
      !
      IF(.NOT. ALLOCATED(ev)) ALLOCATE(ev(n_eigen_to_get))
      IF(mpime == root) ev(1:nglob_to_be_read) = tmp_ev(1:nglob_to_be_read)
      CALL mp_bcast(ev,root,world_comm)
      !
      ! 3) READ THE EIGENVECTOR FILES
      !
      IF(.NOT. ALLOCATED(dvg_exc)) THEN
         ALLOCATE(dvg_exc(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx))
      ENDIF
      !
      ALLOCATE(dvg_tmp(npwx,nbndval0x-n_trunc_bands,kpt_pool%nglob))
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         IF(global_j > n_eigen_to_get) CYCLE
         !
         WRITE(label_j,'(i6.6)') global_j
         fname = TRIM(wbse_save_dir)//'/E'//label_j//'.dat'
         !
         CALL plep_read_G_and_distribute(TRIM(fname),dvg_tmp)
         !
         DO iks = 1,kpt_pool%nloc
            iks_g = kpt_pool%l2g(iks)
            DO lbnd = 1,band_group%nloc
               ibnd = band_group%l2g(lbnd)
               dvg_exc(:,lbnd,iks,local_j) = dvg_tmp(:,ibnd,iks_g)
            ENDDO
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE(dvg_tmp)
      !
      ! MPI barrier
      !
      CALL mp_barrier(world_comm)
      !
      ! Timing
      !
      time_spent(2) = get_clock('plep_db')
      CALL stop_clock('plep_db')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"SAVE read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wbse_save_dir)
      WRITE(stdout,'(5x,"Eigen. found : ",i12)') n_eigen_to_get
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
END MODULE
