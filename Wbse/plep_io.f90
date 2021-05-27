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
MODULE plep_io
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,       ONLY : DP
  USE mp_global,   ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm,my_pool_id,my_bgrp_id,inter_bgrp_comm,inter_pool_comm
  USE mp,          ONLY : mp_max
  USE westcom,     ONLY : npwq,npwq_g,nbndval0x
  !wbsecom combined into westcom
  !USE wbsecom,     ONLY : nbndval0x
  USE pwcom,       ONLY : nks,npwx
  USE wvfct,       ONLY : npwx
  USE gvect,       ONLY : ig_l2g
  !
  IMPLICIT NONE
  !
  PUBLIC plep_merge_and_write_G_wfc
  PUBLIC plep_read_G_and_distribute_wfc
  !
  INTERFACE plep_merge_and_write_G
     MODULE PROCEDURE plep_merge_and_write_G_2d, &
          plep_merge_and_write_G_3d
  END INTERFACE
  !
  INTERFACE plep_read_G_and_distribute
     MODULE PROCEDURE plep_read_G_and_distribute_2d, &
          plep_read_G_and_distribute_3d
  END INTERFACE
  !
  CONTAINS
    !
    ! ******************************************
    ! WRITE IN G SPACE
    !       wfc is passed distributed in G space
    !       then merged and written in R space
    ! ******************************************
    !
    SUBROUTINE plep_merge_and_write_G_2d(fname,plepg,nbndval)
      !
      USE mp_wave,      ONLY : mergewf
      USE mp,           ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,      INTENT(IN) :: nbndval
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP),  INTENT(IN) :: plepg(npwx,nbndval)
      !
      ! Scratch
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ibnd
      !
      !
      IF(my_pool_id.NE.0) RETURN
      IF(my_bgrp_id.NE.0) RETURN
      !
      ! Resume all components
      !
      IF(me_bgrp==root_bgrp) THEN
        !
        ! ... open XML descriptor
        !
        CALL iotk_free_unit( iun, ierr )
        CALL iotk_open_write( iun, FILE = TRIM(fname), BINARY = .TRUE.)
        CALL iotk_write_begin( iun, 'PLEP_GSPACE' )
        !
      ENDIF
      !
      ALLOCATE( tmp_vec(npwq_g) )
      !
      DO ibnd = 1, nbndval
         !
         tmp_vec=0._DP
         !
         CALL mergewf( plepg(:,ibnd), tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
         !
         ! ONLY ROOT W/IN BGRP WRITES
         !
         IF(me_bgrp==root_bgrp) THEN
           !
           ! ... open XML descriptor
           !
           CALL iotk_write_dat( iun, "ndim" , npwq_g )
           CALL iotk_write_dat( iun, "plep" // iotk_index( ibnd ), tmp_vec(1:npwq_g) )
           !
         ENDIF
         !
      ENDDO
      !
      IF(me_bgrp==root_bgrp) THEN
        !
        CALL iotk_write_end( iun, 'PLEP_GSPACE' )
        CALL iotk_close_write( iun )
        !
      ENDIF
      !
      DEALLOCATE( tmp_vec )
      !
    END SUBROUTINE
    !
    ! ******************************************
    ! READ IN G SPACE
    !       wfc is read merged in G space
    !       then split in G space
    ! ******************************************
    !
    SUBROUTINE plep_read_G_and_distribute_2d(fname,plepg,nbndval)
      !
      USE mp_wave,      ONLY : splitwf
      USE mp,           ONLY : mp_bcast
      USE mp_global,    ONLY : intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,      INTENT(IN) :: nbndval
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP), INTENT(OUT) :: plepg(npwx,nbndval)
      !
      ! Scratch
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ig,ibnd
      !
      ! Resume all components
      !
      ALLOCATE( tmp_vec(npwq_g) )
      tmp_vec=0._DP
      plepg=0._DP
      IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            ! ... open XML descriptor
            !
            CALL iotk_free_unit( iun, ierr )
            CALL iotk_open_read( iun, FILE = TRIM(fname), BINARY = .TRUE., IERR = ierr)
            !
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(ierr,0,inter_bgrp_comm)
      !CALL errore( 'read_plep ', &
      !             'cannot open restart file for reading', ierr )
      !
      IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            CALL iotk_scan_begin( iun, 'PLEP_GSPACE' )
            !
         ENDIF
         !
      ENDIF
      !
      !
      DO ibnd = 1, nbndval
         !
         !
         IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            IF(me_bgrp==root_bgrp) THEN
               !
               CALL iotk_scan_dat( iun, &
                                   "plep" // iotk_index( ibnd ), tmp_vec(1:npwq_g) )
               !
            ENDIF
            !
            CALL splitwf( plepg(:,ibnd), tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
            !
         ENDIF
         !
      ENDDO
      !
      IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            CALL iotk_scan_end( iun, 'PLEP_GSPACE' )
            CALL iotk_close_read( iun )
            !
         ENDIF
         !
      ENDIF
      !
      DEALLOCATE( tmp_vec )
      !
      CALL mp_bcast(plepg,0,inter_bgrp_comm)
      CALL mp_bcast(plepg,0,inter_pool_comm)
      !
    END SUBROUTINE
    !
    !
    ! ******************************************
    ! WRITE IN G SPACE
    !       wfc is passed distributed in G space
    !       then merged and written in R space
    ! ******************************************
    !
    SUBROUTINE plep_merge_and_write_G_3d(fname,plepg)
      !
      USE mp_wave,      ONLY : mergewf
      USE mp,           ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP), INTENT(IN)  :: plepg(npwx,nbndval0x,nks)
      !
      ! Scratch
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ibnd,ik,npwx_g
      !
      npwx_g=MAXVAL(ig_l2g(1:npwx))
      CALL mp_max(npwx_g,intra_bgrp_comm)
      !
      IF(my_pool_id.NE.0) RETURN
      IF(my_bgrp_id.NE.0) RETURN
      !
      ! Resume all components
      !
      IF(me_bgrp==root_bgrp) THEN
        !
        ! ... open XML descriptor
        !
        CALL iotk_free_unit( iun, ierr )
        CALL iotk_open_write( iun, FILE = TRIM(fname), BINARY = .TRUE.)
        CALL iotk_write_begin( iun, 'PLEP_GSPACE' )
        !
      ENDIF
      !
      ALLOCATE( tmp_vec(npwx_g) )
      !
      DO ik = 1, nks
         !
         DO ibnd = 1, nbndval0x
            !
            tmp_vec=0._DP
            !
            CALL mergewf( plepg(:,ibnd,ik), tmp_vec, npwx, ig_l2g(1:npwx), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
            !
            ! ONLY ROOT W/IN BGRP WRITES
            !
            IF(me_bgrp==root_bgrp) THEN
              !
              ! ... open XML descriptor
              !
              CALL iotk_write_dat( iun, "ik" , ik )
              CALL iotk_write_dat( iun, "ibnd" , ibnd )
              CALL iotk_write_dat( iun, "ndim" , npwx_g )
              CALL iotk_write_dat( iun, "plep", tmp_vec(1:npwx_g) )
              !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      !
      IF(me_bgrp==root_bgrp) THEN
        !
        CALL iotk_write_end( iun, 'PLEP_GSPACE' )
        CALL iotk_close_write( iun )
        !
      ENDIF
      !
      DEALLOCATE( tmp_vec )
      !
    END SUBROUTINE
    !
    ! ******************************************
    ! READ IN G SPACE
    !       wfc is read merged in G space
    !       then split in G space
    ! ******************************************
    !
    SUBROUTINE plep_read_G_and_distribute_3d(fname,plepg)
      !
      USE mp_wave,      ONLY : splitwf
      USE mp,           ONLY : mp_bcast
      USE mp_global,    ONLY : intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP), INTENT(OUT) :: plepg(npwx,nbndval0x,nks)
      !
      ! Scratch
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ig,ibnd,ik,npwx_g
      !
      npwx_g=MAXVAL(ig_l2g(1:npwx))
      CALL mp_max(npwx_g,intra_bgrp_comm)
      !
      ! Resume all components
      !
      ALLOCATE( tmp_vec(npwx_g) )
      tmp_vec=0._DP
      plepg=0._DP
      IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            ! ... open XML descriptor
            !
            CALL iotk_free_unit( iun, ierr )
            CALL iotk_open_read( iun, FILE = TRIM(fname), BINARY = .TRUE., IERR = ierr)
            !
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(ierr,0,inter_bgrp_comm)
      !CALL errore( 'read_plep ', &
      !             'cannot open restart file for reading', ierr )
      !
      IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            CALL iotk_scan_begin( iun, 'PLEP_GSPACE' )
            !
         ENDIF
         !
      ENDIF
      !
      !
      DO ik = 1, nks
         !
         DO ibnd = 1, nbndval0x
            !
            !
            IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
              !
              ! ONLY ROOT W/IN BGRP READS
              !
              IF(me_bgrp==root_bgrp) THEN
                !
                CALL iotk_scan_dat( iun, &
                                   "plep", tmp_vec(1:npwx_g) )
                !
              ENDIF
              !
              CALL splitwf( plepg(:,ibnd,ik), tmp_vec, npwx, ig_l2g(1:npwx), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
              !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      !
      IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            CALL iotk_scan_end( iun, 'PLEP_GSPACE' )
            CALL iotk_close_read( iun )
            !
         ENDIF
         !
      ENDIF
      !
      DEALLOCATE( tmp_vec )
      !
      CALL mp_bcast(plepg,0,inter_bgrp_comm)
      CALL mp_bcast(plepg,0,inter_pool_comm)
      !
    END SUBROUTINE
    !
    !
    ! ******************************************
    ! WRITE IN G SPACE
    !       wfc is passed distributed in G space
    !       then merged and written in R space
    ! ******************************************
    !
    SUBROUTINE plep_merge_and_write_G_wfc(fname,plepg,nbnd)
      !
      USE mp_wave,      ONLY : mergewf
      USE mp,           ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*), INTENT(IN) :: fname
      INTEGER,      INTENT(IN) :: nbnd
      COMPLEX(DP),  INTENT(IN) :: plepg(npwx,nbnd)
      !
      ! Scratch
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ibnd
      !
      !
      IF(my_pool_id.NE.0) RETURN
      IF(my_bgrp_id.NE.0) RETURN
      !
      ! Resume all components
      !
      IF(me_bgrp==root_bgrp) THEN
        !
        ! ... open XML descriptor
        !
        CALL iotk_free_unit( iun, ierr )
        CALL iotk_open_write( iun, FILE = TRIM(fname), BINARY = .TRUE.)
        CALL iotk_write_begin( iun, 'PLEP_GSPACE' )
        !
      ENDIF
      !
      ALLOCATE( tmp_vec(npwq_g) )
      !
      DO ibnd = 1, nbnd
         !
         tmp_vec=0._DP
         !
         CALL mergewf( plepg(:,ibnd), tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
         !
         ! ONLY ROOT W/IN BGRP WRITES
         !
         IF(me_bgrp==root_bgrp) THEN
           !
           ! ... open XML descriptor
           !
           CALL iotk_write_dat( iun, "ndim" , npwq_g )
           CALL iotk_write_dat( iun, "plep" // iotk_index( ibnd ), tmp_vec(1:npwq_g) )
           !
         ENDIF
         !
      ENDDO
      !
      IF(me_bgrp==root_bgrp) THEN
        !
        CALL iotk_write_end( iun, 'PLEP_GSPACE' )
        CALL iotk_close_write( iun )
        !
      ENDIF
      !
      DEALLOCATE( tmp_vec )
      !
    END SUBROUTINE
    !
    ! ******************************************
    ! READ IN G SPACE
    !       wfc is read merged in G space
    !       then split in G space
    ! ******************************************
    !
    SUBROUTINE plep_read_G_and_distribute_wfc(fname,plepg,nbnd)
      !
      USE io_global,            ONLY : stdout
      USE mp_wave,      ONLY : splitwf
      USE mp,           ONLY : mp_bcast
      USE mp_global,    ONLY : intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*), INTENT(IN) :: fname
      INTEGER,      INTENT(IN) :: nbnd
      COMPLEX(DP), INTENT(OUT) :: plepg(npwx,nbnd)
      !
      ! Scratch
      !
      COMPLEX(DP), ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun, ierr, ig, ibnd
      INTEGER :: ngw_, nbnd_, ispin_, nspin_, igwx_, ik_, nk_
      REAL(DP):: scalef_
      CHARACTER(iotk_attlenx)  :: attr
      !
      ! Resume all components
      !
      ALLOCATE( tmp_vec(npwq_g) )
      !
      tmp_vec=0._DP
      plepg=0._DP
      !
      IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF (me_bgrp==root_bgrp) THEN
            !
            ! ... open XML descriptor
            !
            CALL iotk_free_unit( iun, ierr )
            CALL iotk_open_read( iun, FILE = TRIM(fname), BINARY = .TRUE., IERR = ierr)
            !
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(ierr,0,inter_bgrp_comm)
      !
      IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF (me_bgrp==root_bgrp) THEN
            !
            ! ... open XML descriptor
            !
            CALL iotk_scan_empty( iun, "INFO", attr )
            !
            CALL iotk_scan_attr( attr, "ngw",          ngw_  )
            CALL iotk_scan_attr( attr, "nbnd",         nbnd_ )
            CALL iotk_scan_attr( attr, "ik",           ik_   )
            CALL iotk_scan_attr( attr, "nk",           nk_   )
            CALL iotk_scan_attr( attr, "ispin",        ispin_ )
            CALL iotk_scan_attr( attr, "nspin",        nspin_ )
            CALL iotk_scan_attr( attr, "igwx",         igwx_ )
            CALL iotk_scan_attr( attr, "scale_factor", scalef_ )
            !
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast( ngw_,    root_bgrp, intra_bgrp_comm)
      CALL mp_bcast( nbnd_,   root_bgrp, intra_bgrp_comm)
      CALL mp_bcast( ik_,     root_bgrp, intra_bgrp_comm)
      CALL mp_bcast( nk_,     root_bgrp, intra_bgrp_comm)
      CALL mp_bcast( ispin_,  root_bgrp, intra_bgrp_comm)
      CALL mp_bcast( nspin_,  root_bgrp, intra_bgrp_comm)
      CALL mp_bcast( igwx_,   root_bgrp, intra_bgrp_comm)
      CALL mp_bcast( scalef_, root_bgrp, intra_bgrp_comm)
      !
      DO ibnd = 1, nbnd_
         !
         IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            IF (me_bgrp==root_bgrp) THEN
               !
               CALL iotk_scan_dat( iun, &
                                   "evc" // iotk_index( ibnd ), tmp_vec(1:igwx_) )
               IF ( npwq_g > igwx_ ) tmp_vec((igwx_+1):npwq_g) = 0.0_DP
               !
            ENDIF
            !
            CALL splitwf( plepg(:,ibnd), tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
            !
         ENDIF
         !
      ENDDO
      !
      IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
         !
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            CALL iotk_close_read( iun )
            !
         ENDIF
         !
      ENDIF
      !
      DEALLOCATE( tmp_vec )
      !
      CALL mp_bcast(plepg,0,inter_bgrp_comm)
      CALL mp_bcast(plepg,0,inter_pool_comm)
      !
    END SUBROUTINE
    !
END MODULE
