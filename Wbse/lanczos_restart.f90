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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
MODULE lanczos_restart
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !------------------------------------------------------------------
    SUBROUTINE lanczos_restart_write(nipol_input,ipol_stopped,ilan_stopped)
      !----------------------------------------------------------------
      !
      USE kinds,               ONLY : DP,i8b
      USE mp_world,            ONLY : mpime,root,world_comm
      USE mp,                  ONLY : mp_barrier
      USE io_global,           ONLY : stdout
      USE westcom,             ONLY : wbse_restart_dir,n_lanczos,beta_store,zeta_store
      USE lsda_mod,            ONLY : nspin
      USE json_module,         ONLY : json_file
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN
      USE base64_module,       ONLY : islittleendian
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: nipol_input,ipol_stopped,ilan_stopped
      !
      ! Workspace
      !
      INTEGER :: iun
      CHARACTER(LEN=256) :: fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20), EXTERNAL :: human_readable_time
      TYPE(json_file) :: json
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! MKDIR
      !
      CALL my_mkdir(TRIM(wbse_restart_dir))
      !
      CALL start_clock('lan_restart')
      time_spent(1) = get_clock('lan_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add('nipol_input',nipol_input)
         CALL json%add('n_lanczos',n_lanczos)
         CALL json%add('nspin',nspin)
         CALL json%add('ipol_stopped',ipol_stopped)
         CALL json%add('ilan_stopped',ilan_stopped)
         !
         fname = TRIM(wbse_restart_dir)//'/summary.json'
         OPEN(NEWUNIT=iun,FILE=TRIM(fname))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
         header = 0
         header(HD_ID_VERSION) = HD_VERSION
         IF(islittleendian()) THEN
            header(HD_ID_LITTLE_ENDIAN) = 1
         ENDIF
         !
         fname = TRIM(wbse_restart_dir)//'/lanczos.dat'
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
         offset = 1
         WRITE(iun,POS=offset) header
         offset = offset+SIZEOF(header)
         WRITE(iun,POS=offset) beta_store(1:n_lanczos,1:nipol_input,1:nspin)
         offset = offset+SIZEOF(beta_store)
         WRITE(iun,POS=offset) zeta_store(1:n_lanczos,1:3,1:nipol_input,1:nspin)
         CLOSE(iun)
         !
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL stop_clock('lan_restart')
      time_spent(2) = get_clock('lan_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location   : ",a)') TRIM(wbse_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    SUBROUTINE lanczos_restart_read(nipol_input,ipol_stopped,ilan_stopped)
      !
      USE kinds,               ONLY : DP,i8b
      USE mp_world,            ONLY : mpime,root,world_comm
      USE mp,                  ONLY : mp_barrier,mp_bcast
      USE io_global,           ONLY : stdout
      USE westcom,             ONLY : wbse_restart_dir,n_lanczos,beta_store,zeta_store
      USE lsda_mod,            ONLY : nspin
      USE json_module,         ONLY : json_file
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN
      USE base64_module,       ONLY : islittleendian
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: nipol_input
      INTEGER, INTENT(OUT) :: ipol_stopped,ilan_stopped
      !
      ! Workspace
      !
      INTEGER :: iun,ierr
      INTEGER :: nipol_input_tmp,n_lanczos_tmp,nspin_tmp
      LOGICAL :: found
      CHARACTER(LEN=256) :: fname
      REAL(DP), EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20), EXTERNAL :: human_readable_time
      TYPE(json_file) :: json
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('lan_restart')
      time_spent(1) = get_clock('lan_restart')
      !
      IF(mpime == root) THEN
         !
         fname = TRIM(wbse_restart_dir)//'/summary.json'
         !
         CALL json%initialize()
         CALL json%load(filename=fname)
         !
         CALL json%get('nipol_input',nipol_input_tmp,found)
         IF(.NOT. found) CALL errore('lanczos_restart_read','nipol_input not found',1)
         CALL json%get('n_lanczos',n_lanczos_tmp,found)
         IF(.NOT. found) CALL errore('lanczos_restart_read','n_lanczos not found',1)
         CALL json%get('nspin',nspin_tmp,found)
         IF(.NOT. found) CALL errore('lanczos_restart_read','nspin not found',1)
         CALL json%get('ipol_stopped',ipol_stopped,found)
         IF(.NOT. found) CALL errore('lanczos_restart_read','ipol_stopped not found',1)
         CALL json%get('ilan_stopped',ilan_stopped,found)
         IF(.NOT. found) CALL errore('lanczos_restart_read','ilan_stopped not found',1)
         !
         CALL json%destroy()
         !
         IF(nipol_input_tmp /= nipol_input) CALL errore('lanczos_restart_read','inconsistent nipol_input',1)
         IF(n_lanczos_tmp > n_lanczos) CALL errore('lanczos_restart_read','last n_lanczos > n_lanczos',1)
         IF(nspin_tmp /= nspin) CALL errore('lanczos_restart_read','inconsistent nspin',1)
         IF(ipol_stopped > nipol_input) CALL errore('lanczos_restart_read','ipol_stopped > nipol_input',1)
         IF(ilan_stopped > n_lanczos) CALL errore('lanczos_restart_read','ilan_stopped > n_lanczos',1)
         !
         IF(ilan_stopped == n_lanczos .AND. ipol_stopped+1 <= nipol_input) THEN
            ipol_stopped = ipol_stopped+1
            ilan_stopped = 0
         ENDIF
         !
         fname = TRIM(wbse_restart_dir)//'/lanczos.dat'
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
         IF(ierr /= 0) THEN
            CALL errore('lanczos_restart_read','Cannot read file: '//TRIM(fname),1)
         ENDIF
         !
         offset = 1
         READ(iun,POS=offset) header
         IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
            CALL errore('lanczos_restart_read','Unknown file format: '//TRIM(fname),1)
         ENDIF
         IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
            .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
            CALL errore('lanczos_restart_read','Endianness mismatch: '//TRIM(fname),1)
         ENDIF
         !
         offset = offset+SIZEOF(header)
         READ(iun,POS=offset) beta_store(1:n_lanczos,1:nipol_input,1:nspin)
         offset = offset+SIZEOF(beta_store)
         READ(iun,POS=offset) zeta_store(1:n_lanczos,1:3,1:nipol_input,1:nspin)
         CLOSE(iun)
         !
      ENDIF
      !
      CALL mp_bcast(ipol_stopped,root,world_comm)
      CALL mp_bcast(ilan_stopped,root,world_comm)
      CALL mp_bcast(beta_store,root,world_comm)
      CALL mp_bcast(zeta_store,root,world_comm)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL stop_clock('lan_restart')
      time_spent(2) = get_clock('lan_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wbse_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------
    SUBROUTINE lanczos_log(ipol_iter,ipol_label)
      !----------------------------------------------------------------
      !
      USE kinds,               ONLY : DP
      USE mp_world,            ONLY : mpime,root,world_comm
      USE mp,                  ONLY : mp_barrier
      USE westcom,             ONLY : logfile,n_lanczos,beta_store,zeta_store
      USE lsda_mod,            ONLY : nspin
      USE json_module,         ONLY : json_file
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: ipol_iter
      CHARACTER(LEN=3), INTENT(IN) :: ipol_label
      !
      ! Workspace
      !
      CHARACTER(LEN=6) :: labels
      INTEGER :: iun,is
      TYPE(json_file) :: json
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(logfile))
         !
         DO is = 1,nspin
            WRITE(labels,'(I6.6)') is
            CALL json%add('output.lanczos.K'//labels//'.'//TRIM(ipol_label)//'.beta',&
            & beta_store(:,ipol_iter,is))
            CALL json%add('output.lanczos.K'//labels//'.'//TRIM(ipol_label)//'.zeta',&
            & RESHAPE(zeta_store(:,:,ipol_iter,is),(/n_lanczos*3/)))
         ENDDO
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(logfile))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE lanczos_restart_clear()
      !------------------------------------------------------------------------
      !
      USE mp_world,             ONLY : root,mpime,world_comm
      USE mp,                   ONLY : mp_barrier,mp_bcast
      USE westcom,              ONLY : wbse_restart_dir
      USE clib_wrappers,        ONLY : f_rmdir
      USE west_io,              ONLY : remove_if_present
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      INTEGER :: ierr
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! ... clear the main restart directory
      !
      IF(mpime == root) THEN
         CALL remove_if_present(TRIM(wbse_restart_dir)//'/summary.json')
         ierr = f_rmdir(TRIM(wbse_restart_dir))
      ENDIF
      !
      CALL mp_bcast(ierr,root,world_comm)
      !
      CALL errore('lanczos_restart','cannot clear restart',ierr)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
END MODULE
