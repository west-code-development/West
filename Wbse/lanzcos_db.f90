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
MODULE lanzcos_db
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  CONTAINS
    !
    !
    ! *****************************
    ! D0PSI WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanzcos_d0psi_write ()
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wbse_save_dir,  d0psi
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE io_push,              ONLY : io_push_bar
      !
      IMPLICIT NONE
      !
      INTEGER                :: ipol
      CHARACTER(LEN=256)     :: fname
      REAL(DP), EXTERNAL     :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(LEN=6)       :: my_label
      INTEGER, PARAMETER     :: n_ipol = 3
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('lanzcos_d0psi_write')
      time_spent(1)=get_clock('lanzcos_d0psi_write')
      !
      ! 1) WRITE TO DISK THE D0PSI
      !
      DO ipol = 1, n_ipol
         !
         WRITE (my_label,'(i6.6)')  ipol
         fname = TRIM( wbse_save_dir ) // "/D0PSI_"//TRIM( my_label )//".dat"
         CALL plep_merge_and_write_G(fname,d0psi(:,:,:,ipol))
         !
      ENDDO
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      CALL stop_clock('lanzcos_d0psi_write')
      time_spent(2)=get_clock('lanzcos_d0psi_write')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'D0PSI written in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wbse_save_dir )
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
    ! *****************************
    ! D0PSI READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanzcos_d0psi_read ()
      !------------------------------------------------------------------------
      !
      USE westcom,             ONLY : wbse_save_dir, d0psi
      USE io_global,           ONLY : stdout
      USE mp,                  ONLY : mp_bcast,mp_barrier
      USE mp_world,            ONLY : world_comm,mpime,root
      USE plep_io,             ONLY : plep_read_G_and_distribute
      USE io_push,             ONLY : io_push_bar
      !
      IMPLICIT NONE
      !
      INTEGER            :: ipol
      CHARACTER(LEN=256) :: dirname,fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(LEN=6)   :: my_label
      INTEGER, PARAMETER :: n_ipol = 3
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('lanzcos_d0psi_read')
      time_spent(1)=get_clock('lanzcos_d0psi_read')
      !
      DO ipol = 1, n_ipol
         !
         WRITE (my_label,'(i6.6)')  ipol
         fname = TRIM( wbse_save_dir ) // "/D0PSI_"//TRIM( my_label )//".dat"
         CALL plep_read_G_and_distribute(fname,d0psi(:,:,:,ipol))
         !
      ENDDO
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      CALL stop_clock('lanzcos_d0psi_read')
      time_spent(2)=get_clock('lanzcos_d0psi_read')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'D0PSI read in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wbse_save_dir )
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanzcos_evcs_write( evc1, evc1_old )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wbse_save_dir,d0psi
      USE plep_io,              ONLY : plep_merge_and_write_G
      USE io_push,              ONLY : io_push_bar
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(IN) :: evc1(:,:,:)
      COMPLEX(DP), INTENT(IN) :: evc1_old(:,:,:)
      !
      CHARACTER(LEN=256)    :: fname
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('lanzcos_evcs_write')
      time_spent(1)=get_clock('lanzcos_evcs_write')
      !
      ! 1) WRITE TO DISK THE D0PSI
      !
      fname = TRIM( wbse_save_dir ) // "/EVC1.dat"
      CALL plep_merge_and_write_G(fname,evc1)
      !
      fname = TRIM( wbse_save_dir ) // "/EVC1_OLD.dat"
      CALL plep_merge_and_write_G(fname,evc1_old)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      CALL stop_clock('lanzcos_evcs_write')
      time_spent(2)=get_clock('lanzcos_evcs_write')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'EVC1 EVC1_OLD written in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wbse_save_dir )
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
    ! *****************************
    ! D0PSI READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanzcos_evcs_read(evc1, evc1_old)
      !------------------------------------------------------------------------
      !
      USE westcom,             ONLY : wbse_save_dir, d0psi
      USE io_global,           ONLY : stdout
      USE mp,                  ONLY : mp_bcast,mp_barrier
      USE mp_world,            ONLY : world_comm,mpime,root
      USE plep_io,             ONLY : plep_read_G_and_distribute
      USE io_push,             ONLY : io_push_bar
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(INOUT) :: evc1(:,:,:)
      COMPLEX(DP), INTENT(INOUT) :: evc1_old(:,:,:)
      !
      CHARACTER(LEN=256)  :: fname
      REAL(DP), EXTERNAL  :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('lanzcos_evcs_read')
      time_spent(1)=get_clock('lanzcos_evcs_read')
      !
      fname = TRIM( wbse_save_dir ) // "/EVC1.dat"
      CALL plep_read_G_and_distribute(fname,evc1)
      fname = TRIM( wbse_save_dir ) // "/EVC1_OLD.dat"
      CALL plep_read_G_and_distribute(fname,evc1_old)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      CALL stop_clock('lanzcos_evcs_read')
      time_spent(2)=get_clock('lanzcos_evcs_read')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'EVC1 EVC1_OLD read in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wbse_save_dir )
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
END MODULE
