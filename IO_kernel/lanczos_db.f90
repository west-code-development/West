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
MODULE lanczos_db
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    ! *****************************
    ! D0PSI WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanczos_d0psi_write()
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : world_comm
      USE mp_global,            ONLY : inter_image_comm,my_image_id
      USE io_global,            ONLY : stdout
      USE pwcom,                ONLY : npwx,nks
      USE westcom,              ONLY : wbse_save_dir,d0psi,nbndval0x,n_trunc_bands
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : aband
      USE west_mp,              ONLY : mp_root_sum_c16_3d
      !
      IMPLICIT NONE
      !
      INTEGER :: ipol
      INTEGER :: lbnd,ibnd
      CHARACTER(LEN=256) :: fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20), EXTERNAL :: human_readable_time
      CHARACTER(LEN=6) :: my_label
      INTEGER, PARAMETER :: n_ipol = 3
      COMPLEX(DP), ALLOCATABLE :: d0psi_tmp(:,:,:)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('lan_d0psi_write')
      time_spent(1) = get_clock('lan_d0psi_write')
      !
      ! WRITE D0PSI TO DISK
      !
      ALLOCATE(d0psi_tmp(npwx,nbndval0x-n_trunc_bands,nks))
      !
      DO ipol = 1,n_ipol
         !
         WRITE(my_label,'(i6.6)') ipol
         !
         d0psi_tmp(:,:,:) = (0._DP,0._DP)
         !
         DO lbnd = 1,aband%nloc
            ibnd = aband%l2g(lbnd)
            d0psi_tmp(:,ibnd,:) = d0psi(:,lbnd,:,ipol)
         ENDDO
         !
         CALL mp_root_sum_c16_3d(d0psi_tmp,0,inter_image_comm)
         !
         IF(my_image_id == 0) THEN
            fname = TRIM(wbse_save_dir)//'/d0psi_'//my_label//'.dat'
            CALL plep_merge_and_write_G(fname,d0psi_tmp)
         ENDIF
         !
      ENDDO
      !
      DEALLOCATE(d0psi_tmp)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL stop_clock('lan_d0psi_write')
      time_spent(2) = get_clock('lan_d0psi_write')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"d0psi written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wbse_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    ! *****************************
    ! D0PSI READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanczos_d0psi_read()
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : world_comm
      USE io_global,            ONLY : stdout
      USE pwcom,                ONLY : npwx,nks
      USE westcom,              ONLY : wbse_save_dir,d0psi,nbndval0x,n_trunc_bands
      USE plep_io,              ONLY : plep_read_G_and_distribute
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : aband
      !
      IMPLICIT NONE
      !
      INTEGER :: ipol
      INTEGER :: lbnd,ibnd
      CHARACTER(LEN=256) :: fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20), EXTERNAL :: human_readable_time
      CHARACTER(LEN=6) :: my_label
      INTEGER, PARAMETER :: n_ipol = 3
      COMPLEX(DP), ALLOCATABLE :: d0psi_tmp(:,:,:)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('lan_d0psi_read')
      time_spent(1) = get_clock('lan_d0psi_read')
      !
      ALLOCATE(d0psi_tmp(npwx,nbndval0x-n_trunc_bands,nks))
      !
      DO ipol = 1,n_ipol
         !
         WRITE(my_label,'(i6.6)') ipol
         fname = TRIM(wbse_save_dir)//'/d0psi_'//my_label//'.dat'
         !
         CALL plep_read_G_and_distribute(fname,d0psi_tmp)
         !
         DO lbnd = 1,aband%nloc
            ibnd = aband%l2g(lbnd)
            d0psi(:,lbnd,:,ipol) = d0psi_tmp(:,ibnd,:)
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE(d0psi_tmp)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL stop_clock('lan_d0psi_read')
      time_spent(2) = get_clock('lan_d0psi_read')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"d0psi read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wbse_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanczos_evcs_write(evc1, evc1_old)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : world_comm
      USE mp_global,            ONLY : inter_image_comm,my_image_id
      USE io_global,            ONLY : stdout
      USE pwcom,                ONLY : npwx,nks
      USE westcom,              ONLY : wbse_save_dir,nbndval0x,n_trunc_bands
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : aband
      USE west_mp,              ONLY : mp_root_sum_c16_3d
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(IN) :: evc1(npwx,aband%nlocx,nks)
      COMPLEX(DP), INTENT(IN) :: evc1_old(npwx,aband%nlocx,nks)
      !
      INTEGER :: lbnd,ibnd
      CHARACTER(LEN=256) :: fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20), EXTERNAL :: human_readable_time
      COMPLEX(DP), ALLOCATABLE :: evc1_tmp(:,:,:)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('lan_evc_write')
      time_spent(1) = get_clock('lan_evc_write')
      !
      ! WRITE EVC1 & EVC1_OLD TO DISK
      !
      ALLOCATE(evc1_tmp(npwx,nbndval0x-n_trunc_bands,nks))
      !
      evc1_tmp(:,:,:) = (0._DP,0._DP)
      !
      DO lbnd = 1,aband%nloc
         ibnd = aband%l2g(lbnd)
         evc1_tmp(:,ibnd,:) = evc1(:,lbnd,:)
      ENDDO
      !
      CALL mp_root_sum_c16_3d(evc1_tmp,0,inter_image_comm)
      !
      IF(my_image_id == 0) THEN
         fname = TRIM(wbse_save_dir)//'/evc1.dat'
         CALL plep_merge_and_write_G(fname,evc1_tmp)
      ENDIF
      !
      evc1_tmp(:,:,:) = (0._DP,0._DP)
      !
      DO lbnd = 1,aband%nloc
         ibnd = aband%l2g(lbnd)
         evc1_tmp(:,ibnd,:) = evc1_old(:,lbnd,:)
      ENDDO
      !
      CALL mp_root_sum_c16_3d(evc1_tmp,0,inter_image_comm)
      !
      IF(my_image_id == 0) THEN
         fname = TRIM(wbse_save_dir)//'/evc1_old.dat'
         CALL plep_merge_and_write_G(fname,evc1_tmp)
      ENDIF
      !
      DEALLOCATE(evc1_tmp)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL stop_clock('lan_evc_write')
      time_spent(2) = get_clock('lan_evc_write')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"evc1 evc1_old written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wbse_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    ! *****************************
    ! D0PSI READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanczos_evcs_read(evc1, evc1_old)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : world_comm
      USE io_global,            ONLY : stdout
      USE pwcom,                ONLY : npwx,nks
      USE westcom,              ONLY : wbse_save_dir,nbndval0x,n_trunc_bands
      USE plep_io,              ONLY : plep_read_G_and_distribute
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : aband
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(OUT) :: evc1(npwx,aband%nlocx,nks)
      COMPLEX(DP), INTENT(OUT) :: evc1_old(npwx,aband%nlocx,nks)
      !
      INTEGER :: lbnd,ibnd
      CHARACTER(LEN=256) :: fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20), EXTERNAL :: human_readable_time
      COMPLEX(DP), ALLOCATABLE :: evc1_tmp(:,:,:)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('lan_evc_read')
      time_spent(1) = get_clock('lan_evc_read')
      !
      ALLOCATE(evc1_tmp(npwx,nbndval0x-n_trunc_bands,nks))
      !
      fname = TRIM(wbse_save_dir)//'/evc1.dat'
      CALL plep_read_G_and_distribute(fname,evc1_tmp)
      !
      DO lbnd = 1,aband%nloc
         ibnd = aband%l2g(lbnd)
         evc1(:,lbnd,:) = evc1_tmp(:,ibnd,:)
      ENDDO
      !
      fname = TRIM(wbse_save_dir)//'/evc1_old.dat'
      CALL plep_read_G_and_distribute(fname,evc1_tmp)
      !
      DO lbnd = 1,aband%nloc
         ibnd = aband%l2g(lbnd)
         evc1_old(:,lbnd,:) = evc1_tmp(:,ibnd,:)
      ENDDO
      !
      DEALLOCATE(evc1_tmp)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL stop_clock('lan_evc_read')
      time_spent(2) = get_clock('lan_evc_read')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"evc1 evc1_old read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wbse_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
END MODULE
