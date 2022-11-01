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
MODULE lanczos_restart
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir
  USE lsda_mod,  ONLY : nspin
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER, PRIVATE :: iunout
  !
  PUBLIC :: lanczos_restart_write
  PUBLIC :: lanczos_restart_read
  PUBLIC :: lanczos_postpro_write
  !
  CONTAINS
    !
    !------------------------------------------------------------------
    SUBROUTINE lanczos_restart_write (nipol_input, pliter_stop, lriter_stop)
      !----------------------------------------------------------------
      !
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wbse_save_dir,ipol_input,n_lanczos,alpha_store,&
                                     & beta_store,gamma_store,zeta_store
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN)  :: nipol_input, pliter_stop, lriter_stop
      !
      ! Workspace
      !
      INTEGER :: ierr, ipol, ipol2, is
      CHARACTER(LEN=256) :: dirname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('wbse_lanczos_restart')
      time_spent(1)=get_clock('wbse_lanczos_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      dirname = TRIM( wbse_save_dir ) // '/' // TRIM("summary.xml")
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( wbse_save_dir ) // '/' // TRIM("summary.xml") , BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wbse_restart', 'cannot open restart file for writing', ierr )
      !
      IF ( mpime == root ) THEN
         !
         CALL iotk_write_begin( iunout, "SUMMARY" )
         CALL iotk_write_dat( iunout, "nipol_input", nipol_input)
         CALL iotk_write_dat( iunout, "ipol_input", ipol_input)
         CALL iotk_write_dat( iunout, "n_lanczos", n_lanczos)
         CALL iotk_write_dat( iunout, "pliter_stop", pliter_stop)
         CALL iotk_write_dat( iunout, "lriter_stop", lriter_stop)
         CALL iotk_write_end( iunout, "SUMMARY"  )
         !
         CALL iotk_write_begin( iunout, "ALPHA_STORE" )
         DO ipol = 1, nipol_input
            !
            DO is = 1, nspin
               CALL iotk_write_dat( iunout, "alpha_store_k", alpha_store(ipol,:,is))
            ENDDO
            !
         ENDDO
         CALL iotk_write_end( iunout, "ALPHA_STORE" )
         !
         CALL iotk_write_begin( iunout, "BETA_STORE" )
         DO ipol = 1, nipol_input
            !
            DO is = 1, nspin
               CALL iotk_write_dat( iunout, "beta_store_k", beta_store(ipol,:,is))
            ENDDO
            !
         ENDDO
         CALL iotk_write_end( iunout, "BETA_STORE" )
         !
         CALL iotk_write_begin( iunout, "GAMMA_STORE" )
         DO ipol = 1, nipol_input
            !
            DO is = 1, nspin
               CALL iotk_write_dat( iunout, "gamma_store_k", gamma_store(ipol,:,is))
            ENDDO
            !
         ENDDO
         CALL iotk_write_end( iunout, "GAMMA_STORE" )
         !
         CALL iotk_write_begin( iunout, "ZETA_STORE" )
         DO ipol = 1, nipol_input
            CALL iotk_write_dat( iunout, "zeta_store_ipol_i", ipol )
            DO ipol2 = 1, 3
               CALL iotk_write_dat( iunout, "zeta_store_ipol_j", ipol2 )
               DO is = 1, nspin
                  CALL iotk_write_dat( iunout, "zeta_store_k", zeta_store(ipol,ipol2,:,is))
               ENDDO
            ENDDO
         ENDDO
         CALL iotk_write_end( iunout, "ZETA_STORE" )
         !
         CALL iotk_close_write( iunout )
         !
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL stop_clock('wbse_lanczos_restart')
      time_spent(2)=get_clock('wbse_lanczos_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout, "(5x, '[I/O] RESTART written in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, '[I/O] In location   : ',a)") TRIM( dirname )
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    SUBROUTINE lanczos_restart_read (nipol_input, pliter_stop, lriter_stop)
      !
      USE mp_global,           ONLY : world_comm
      USE mp_world,            ONLY : mpime,root,world_comm
      USE mp,                  ONLY : mp_barrier,mp_bcast
      USE io_global,           ONLY : stdout
      USE westcom,             ONLY : wbse_save_dir,ipol_input,n_lanczos,alpha_store,&
                                    & beta_store,gamma_store,zeta_store
      USE mp_global,           ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN)  :: nipol_input
      INTEGER, INTENT(OUT) :: pliter_stop, lriter_stop

      !
      ! Workspace
      !
      CHARACTER(LEN=256) :: dirname
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(LEN=3)   :: ipol_input_tmp
      INTEGER :: iun, ierr, ipol, ipol2, is
      INTEGER :: nipol_input_tmp, n_lanczos_tmp
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock('wbse_lanczos_restart')
      time_spent(1)=get_clock('wbse_lanczos_restart')
      !
      dirname = TRIM( wbse_save_dir ) // '/' // TRIM( 'summary.xml' )
      IF ( mpime==root ) THEN
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_read( iun, FILE = TRIM( wbse_save_dir ) // '/' // TRIM( 'summary.xml' ), IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      IF ( ierr /=0 ) CALL errore( 'wbse_lanczos_restart', 'cannot open restart file for reading', ierr )
      !
      IF ( mpime==root ) THEN
         !
         CALL iotk_scan_begin( iun, "SUMMARY")
         CALL iotk_scan_dat( iun, "nipol_input", nipol_input_tmp )
         CALL iotk_scan_dat( iun, "ipol_input",  ipol_input_tmp )
         CALL iotk_scan_dat( iun, "n_lanczos",    n_lanczos_tmp )
         CALL iotk_scan_dat( iun, "pliter_stop", pliter_stop )
         CALL iotk_scan_dat( iun, "lriter_stop", lriter_stop )
         CALL iotk_scan_end( iun, "SUMMARY"  )
         !
      ENDIF
      !
      CALL mp_bcast(nipol_input_tmp, root, world_comm )
      CALL mp_bcast(ipol_input_tmp, root, world_comm )
      CALL mp_bcast(n_lanczos_tmp, root, world_comm )
      CALL mp_bcast(pliter_stop, root, world_comm )
      CALL mp_bcast(lriter_stop, root, world_comm )
      !
      IF (ipol_input_tmp /= ipol_input) THEN
         CALL errore( 'wbse_lanczos_restart', 'there is inconsistent between previous ipol and ipol in input para', 1)
      ENDIF
      !
      IF (n_lanczos_tmp > n_lanczos) THEN
         CALL errore( 'wbse_lanczos_restart', 'last n_lanczos > n_lanczos', 1)
      ENDIF
      !
      IF (pliter_stop > nipol_input) THEN
         CALL errore( 'wbse_lanczos_restart', 'ipol stopped > nipol_input', 1)
      ENDIF
      !
      IF (lriter_stop > n_lanczos) THEN
         CALL errore( 'wbse_lanczos_restart', ' lriter_stop > n_lanczos', 1)
      ENDIF
      !
      IF ((lriter_stop == n_lanczos).AND.(pliter_stop+1 <= nipol_input)) THEN
         !
         pliter_stop = pliter_stop+1
         lriter_stop = 0
         !
      ENDIF
      !
      IF ( mpime==root ) THEN
         !
         alpha_store(:,:,:) = 0.0_DP
         CALL iotk_scan_begin( iun, "ALPHA_STORE")
         DO ipol = 1, nipol_input
            DO is = 1, nspin
               CALL iotk_scan_dat( iun, "alpha_store_k", alpha_store(ipol,1:n_lanczos_tmp,is) )
            ENDDO
         ENDDO
         CALL iotk_scan_end( iun, "ALPHA_STORE"  )
         !
         beta_store(:,:,:) = 0.0_DP
         CALL iotk_scan_begin( iun, "BETA_STORE")
         DO ipol = 1, nipol_input
            DO is = 1, nspin
               CALL iotk_scan_dat( iun, "beta_store_k", beta_store(ipol,1:n_lanczos_tmp,is) )
            ENDDO
         ENDDO
         CALL iotk_scan_end( iun, "BETA_STORE"  )
         !
         gamma_store(:,:,:) = 0.0_DP
         CALL iotk_scan_begin( iun, "GAMMA_STORE")
         DO ipol = 1, nipol_input
            DO is = 1, nspin
               CALL iotk_scan_dat( iun, "gamma_store_k", gamma_store(ipol,1:n_lanczos_tmp,is) )
            ENDDO
         ENDDO
         CALL iotk_scan_end( iun, "GAMMA_STORE"  )
         !
         zeta_store(:,:,:,:) = 0.0_DP
         CALL iotk_scan_begin( iun, "ZETA_STORE")
         DO ipol = 1, nipol_input
            CALL iotk_scan_dat( iun, "zeta_store_ipol_i", ipol )
            DO ipol2 = 1, 3
               CALL iotk_scan_dat( iun, "zeta_store_ipol_j", ipol2 )
               DO is = 1, nspin
                  CALL iotk_scan_dat( iun, "zeta_store_k", zeta_store(ipol,ipol2,1:n_lanczos_tmp,is) )
               ENDDO
            ENDDO
         ENDDO
         CALL iotk_scan_end( iun, "ZETA_STORE"  )
         !
         CALL iotk_close_read( iun )
         !
      ENDIF
      !
      CALL mp_bcast( alpha_store, 0, intra_image_comm )
      CALL mp_bcast( beta_store, 0, intra_image_comm )
      CALL mp_bcast( gamma_store, 0, intra_image_comm )
      CALL mp_bcast( zeta_store, 0, intra_image_comm )
      !
      ! BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      CALL stop_clock('wbse_lanczos_restart')
      time_spent(2)=get_clock('wbse_lanczos_restart')
      !
      WRITE(stdout,'(1/, 5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout, "(5x, '[I/O] RESTART read in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout, "(5x, '[I/O] In location : ',a)") TRIM( dirname )
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------
    SUBROUTINE lanczos_postpro_write (nipol_input, ipol_iter, ipol_label)
      !----------------------------------------------------------------
      !
      USE mp_world,             ONLY : mpime,root,world_comm
      USE westcom,              ONLY : n_lanczos,alpha_store,beta_store,gamma_store,zeta_store
      USE mp,                   ONLY : mp_barrier,mp_bcast,mp_get
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN)  :: nipol_input, ipol_iter
      CHARACTER(LEN=3),INTENT(IN)   :: ipol_label
      !
      ! Workspace
      !
      INTEGER :: ierr, ipol2,is
      CHARACTER(LEN=256) :: fname
      CHARACTER(LEN=4)   :: my_ip
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! CREATE THE SUMMARY FILE
      !
      WRITE(my_ip,'(i1)') ipol_iter
      fname = "summary."//TRIM(my_ip)//".xml"
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         CALL iotk_open_write( iunout, FILE = TRIM( fname ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      !
      CALL errore( 'wbse_restart', 'cannot open restart file for writing', ierr )
      !
      IF ( mpime == root ) THEN
         !
         CALL iotk_write_begin( iunout, "SUMMARY" )
         CALL iotk_write_dat( iunout, "nspin", nspin)
         CALL iotk_write_dat( iunout, "nipol_input", nipol_input)
         CALL iotk_write_dat( iunout, "ipol_label", ipol_label)
         CALL iotk_write_dat( iunout, "n_lanczos", n_lanczos)
         CALL iotk_write_end( iunout, "SUMMARY"  )
         !
         CALL iotk_write_begin( iunout, "ALPHA_STORE" )
         DO is = 1, nspin
            CALL iotk_write_dat  ( iunout, "alpha_store_k", alpha_store(ipol_iter,:,is))
         ENDDO
         CALL iotk_write_end  ( iunout, "ALPHA_STORE" )
         !
         CALL iotk_write_begin( iunout, "BETA_STORE" )
         DO is = 1, nspin
            CALL iotk_write_dat  ( iunout, "beta_store_k", beta_store(ipol_iter,:,is))
         ENDDO
         CALL iotk_write_end  ( iunout, "BETA_STORE" )
         !
         CALL iotk_write_begin( iunout, "GAMMA_STORE" )
         DO is = 1, nspin
            CALL iotk_write_dat  ( iunout, "gamma_store_k", gamma_store(ipol_iter,:,is))
         ENDDO
         CALL iotk_write_end  ( iunout, "GAMMA_STORE" )
         !
         CALL iotk_write_begin( iunout, "ZETA_STORE" )
         CALL iotk_write_dat  ( iunout, "zeta_store_ipol_i", ipol_iter )
         DO ipol2 = 1, 3
            CALL iotk_write_dat( iunout, "zeta_store_ipol_j", ipol2 )
            DO is = 1, nspin
               CALL iotk_write_dat( iunout, "zeta_store_k", zeta_store(ipol_iter,ipol2,:,is))
            ENDDO
         ENDDO
         CALL iotk_write_end( iunout, "ZETA_STORE" )
         !
         CALL iotk_close_write( iunout )
         !
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
END MODULE
