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
MODULE wfreq_restart
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! BKS: band, k
  !
  TYPE :: bks_type
     INTEGER :: lastdone_ks
     INTEGER :: lastdone_band
     INTEGER :: old_ks
     INTEGER :: old_band
     INTEGER :: max_ks
     INTEGER :: max_band
     INTEGER :: min_ks
     INTEGER :: min_band
  END TYPE bks_type
  !
  ! BKSQ: band, k, q (solve_wfreq_q)
  !
  TYPE :: bksq_type
     INTEGER :: lastdone_q
     INTEGER :: lastdone_ks
     INTEGER :: lastdone_band
     INTEGER :: old_q
     INTEGER :: old_ks
     INTEGER :: old_band
     INTEGER :: max_q
     INTEGER :: max_ks
     INTEGER :: max_band
     INTEGER :: min_q
     INTEGER :: min_ks
     INTEGER :: min_band
  END TYPE bksq_type
  !
  ! BKSKS: band, k, k1 (solve_gfreq_q)
  !
  TYPE :: bksks_type
     INTEGER :: lastdone_ks
     INTEGER :: lastdone_kks
     INTEGER :: lastdone_band
     INTEGER :: old_ks
     INTEGER :: old_kks
     INTEGER :: old_band
     INTEGER :: max_ks
     INTEGER :: max_kks
     INTEGER :: max_band
     INTEGER :: min_ks
     INTEGER :: min_kks
     INTEGER :: min_band
  END TYPE bksks_type
  !
  INTERFACE solvewfreq_restart_write
     MODULE PROCEDURE solvewfreq_restart_write_real, solvewfreq_restart_write_complex, solvewfreq_restart_write_complex_q
  END INTERFACE
  !
  INTERFACE solvewfreq_restart_read
     MODULE PROCEDURE solvewfreq_restart_read_real, solvewfreq_restart_read_complex, solvewfreq_restart_read_complex_q
  END INTERFACE
  !
  CONTAINS
    !
    ! BKS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_bks(bks,dirname,fname)
      !------------------------------------------------------------------------
      !
      USE json_module,          ONLY : json_file
      USE mp_world,             ONLY : mpime,root
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(IN) :: bks
      CHARACTER(LEN=*),INTENT(IN) :: dirname,fname
      !
      ! Workspace
      !
      INTEGER :: iunit
      TYPE(json_file) :: json
      CHARACTER(LEN=*),PARAMETER :: what = 'BKS-SUMMARY.'
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add(what//'lastdone_ks',bks%lastdone_ks)
         CALL json%add(what//'lastdone_band',bks%lastdone_band)
         CALL json%add(what//'old_ks',bks%old_ks)
         CALL json%add(what//'old_band',bks%old_band)
         CALL json%add(what//'max_ks',bks%max_ks)
         CALL json%add(what//'max_band',bks%max_band)
         CALL json%add(what//'min_ks',bks%min_ks)
         CALL json%add(what//'min_band',bks%min_band)
         !
         OPEN(NEWUNIT=iunit,FILE=TRIM(dirname)//'/'//TRIM(fname))
         CALL json%print(iunit)
         CLOSE(iunit)
         CALL json%destroy()
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_bksq(bksq,dirname,fname)
      !------------------------------------------------------------------------
      !
      USE json_module,          ONLY : json_file
      USE mp_world,             ONLY : mpime,root
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksq_type),INTENT(IN) :: bksq
      CHARACTER(LEN=*),INTENT(IN) :: dirname,fname
      !
      ! Workspace
      !
      INTEGER :: iunit
      TYPE(json_file) :: json
      CHARACTER(LEN=*),PARAMETER :: what = 'BKSQ-SUMMARY.'
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add(what//'lastdone_q',bksq%lastdone_q)
         CALL json%add(what//'lastdone_ks',bksq%lastdone_ks)
         CALL json%add(what//'lastdone_band',bksq%lastdone_band)
         CALL json%add(what//'old_q',bksq%old_q)
         CALL json%add(what//'old_ks',bksq%old_ks)
         CALL json%add(what//'old_band',bksq%old_band)
         CALL json%add(what//'max_q',bksq%max_q)
         CALL json%add(what//'max_ks',bksq%max_ks)
         CALL json%add(what//'max_band',bksq%max_band)
         CALL json%add(what//'min_q',bksq%min_q)
         CALL json%add(what//'min_ks',bksq%min_ks)
         CALL json%add(what//'min_band',bksq%min_band)
         !
         OPEN(NEWUNIT=iunit,FILE=TRIM(dirname)//'/'//TRIM(fname))
         CALL json%print(iunit)
         CLOSE(iunit)
         CALL json%destroy()
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_bksks(bksks,dirname,fname)
      !------------------------------------------------------------------------
      !
      USE json_module,          ONLY : json_file
      USE mp_world,             ONLY : mpime,root
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksks_type),INTENT(IN) :: bksks
      CHARACTER(LEN=*),INTENT(IN) :: dirname,fname
      !
      ! Workspace
      !
      INTEGER :: iunit
      TYPE(json_file) :: json
      CHARACTER(LEN=*),PARAMETER :: what = 'BKSKS-SUMMARY.'
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add(what//'lastdone_ks',bksks%lastdone_ks)
         CALL json%add(what//'lastdone_kks',bksks%lastdone_kks)
         CALL json%add(what//'lastdone_band',bksks%lastdone_band)
         CALL json%add(what//'old_ks',bksks%old_ks)
         CALL json%add(what//'old_kks',bksks%old_kks)
         CALL json%add(what//'old_band',bksks%old_band)
         CALL json%add(what//'max_ks',bksks%max_ks)
         CALL json%add(what//'max_kks',bksks%max_kks)
         CALL json%add(what//'max_band',bksks%max_band)
         CALL json%add(what//'min_ks',bksks%min_ks)
         CALL json%add(what//'min_kks',bksks%min_kks)
         CALL json%add(what//'min_band',bksks%min_band)
         !
         OPEN(NEWUNIT=iunit,FILE=TRIM(dirname)//'/'//TRIM(fname))
         CALL json%print(iunit)
         CLOSE(iunit)
         CALL json%destroy()
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_bks(bks,dirname,fname)
      !------------------------------------------------------------------------
      !
      USE json_module,          ONLY : json_file
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp,                   ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(OUT) :: bks
      CHARACTER(LEN=*),INTENT(IN) :: dirname,fname
      !
      ! Workspace
      !
      LOGICAL :: found
      INTEGER :: ival
      INTEGER :: tmp(8)
      TYPE(json_file) :: json
      CHARACTER(LEN=*),PARAMETER :: what = 'BKS-SUMMARY.'
      !
      tmp = 0
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(dirname)//'/'//TRIM(fname))
         !
         CALL json%get(what//'lastdone_ks',ival,found)
         IF(found) tmp(1) = ival
         CALL json%get(what//'lastdone_band',ival,found)
         IF(found) tmp(2) = ival
         CALL json%get(what//'old_ks',ival,found)
         IF(found) tmp(3) = ival
         CALL json%get(what//'old_band',ival,found)
         IF(found) tmp(4) = ival
         CALL json%get(what//'max_ks',ival,found)
         IF(found) tmp(5) = ival
         CALL json%get(what//'max_band',ival,found)
         IF(found) tmp(6) = ival
         CALL json%get(what//'min_ks',ival,found)
         IF(found) tmp(7) = ival
         CALL json%get(what//'min_band',ival,found)
         IF(found) tmp(8) = ival
         !
         CALL json%destroy()
         !
      ENDIF
      !
      CALL mp_bcast(tmp,root,world_comm)
      !
      bks%lastdone_ks = tmp(1)
      bks%lastdone_band = tmp(2)
      bks%old_ks = tmp(3)
      bks%old_band = tmp(4)
      bks%max_ks = tmp(5)
      bks%max_band = tmp(6)
      bks%min_ks = tmp(7)
      bks%min_band = tmp(8)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_bksq(bksq,dirname,fname)
      !------------------------------------------------------------------------
      !
      USE json_module,          ONLY : json_file
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp,                   ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksq_type),INTENT(OUT) :: bksq
      CHARACTER(LEN=*),INTENT(IN) :: dirname,fname
      !
      ! Workspace
      !
      LOGICAL :: found
      INTEGER :: ival
      INTEGER :: tmp(12)
      TYPE(json_file) :: json
      CHARACTER(LEN=*),PARAMETER :: what = 'BKSQ-SUMMARY.'
      !
      tmp = 0
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(dirname)//'/'//TRIM(fname))
         !
         CALL json%get(what//'lastdone_q',ival,found)
         IF(found) tmp(1) = ival
         CALL json%get(what//'lastdone_ks',ival,found)
         IF(found) tmp(2) = ival
         CALL json%get(what//'lastdone_band',ival,found)
         IF(found) tmp(3) = ival
         CALL json%get(what//'old_q',ival,found)
         IF(found) tmp(4) = ival
         CALL json%get(what//'old_ks',ival,found)
         IF(found) tmp(5) = ival
         CALL json%get(what//'old_band',ival,found)
         IF(found) tmp(6) = ival
         CALL json%get(what//'max_q',ival,found)
         IF(found) tmp(7) = ival
         CALL json%get(what//'max_ks',ival,found)
         IF(found) tmp(8) = ival
         CALL json%get(what//'max_band',ival,found)
         IF(found) tmp(9) = ival
         CALL json%get(what//'min_q',ival,found)
         IF(found) tmp(10) = ival
         CALL json%get(what//'min_ks',ival,found)
         IF(found) tmp(11) = ival
         CALL json%get(what//'min_band',ival,found)
         IF(found) tmp(12) = ival
         !
         CALL json%destroy()
         !
      ENDIF
      !
      CALL mp_bcast(tmp,root,world_comm)
      !
      bksq%lastdone_q = tmp(1)
      bksq%lastdone_ks = tmp(2)
      bksq%lastdone_band = tmp(3)
      bksq%old_q = tmp(4)
      bksq%old_ks = tmp(5)
      bksq%old_band = tmp(6)
      bksq%max_q = tmp(7)
      bksq%max_ks = tmp(8)
      bksq%max_band = tmp(9)
      bksq%min_q = tmp(10)
      bksq%min_ks = tmp(11)
      bksq%min_band = tmp(12)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_bksks(bksks,dirname,fname)
      !------------------------------------------------------------------------
      !
      USE json_module,          ONLY : json_file
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp,                   ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksks_type),INTENT(OUT) :: bksks
      CHARACTER(LEN=*),INTENT(IN) :: dirname,fname
      !
      ! Workspace
      !
      LOGICAL :: found
      INTEGER :: ival
      INTEGER :: tmp(12)
      TYPE(json_file) :: json
      CHARACTER(LEN=*),PARAMETER :: what = 'BKSKS-SUMMARY.'
      !
      tmp = 0
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(dirname)//'/'//TRIM(fname))
         !
         CALL json%get(what//'lastdone_ks',ival,found)
         IF(found) tmp(1) = ival
         CALL json%get(what//'lastdone_kks',ival,found)
         IF(found) tmp(2) = ival
         CALL json%get(what//'lastdone_band',ival,found)
         IF(found) tmp(3) = ival
         CALL json%get(what//'old_ks',ival,found)
         IF(found) tmp(4) = ival
         CALL json%get(what//'old_kks',ival,found)
         IF(found) tmp(5) = ival
         CALL json%get(what//'old_band',ival,found)
         IF(found) tmp(6) = ival
         CALL json%get(what//'max_ks',ival,found)
         IF(found) tmp(7) = ival
         CALL json%get(what//'max_kks',ival,found)
         IF(found) tmp(8) = ival
         CALL json%get(what//'max_band',ival,found)
         IF(found) tmp(9) = ival
         CALL json%get(what//'min_ks',ival,found)
         IF(found) tmp(10) = ival
         CALL json%get(what//'min_kks',ival,found)
         IF(found) tmp(11) = ival
         CALL json%get(what//'min_band',ival,found)
         IF(found) tmp(12) = ival
         !
         CALL json%destroy()
         !
      ENDIF
      !
      CALL mp_bcast(tmp,root,world_comm)
      !
      bksks%lastdone_ks = tmp(1)
      bksks%lastdone_kks = tmp(2)
      bksks%lastdone_band = tmp(3)
      bksks%old_ks = tmp(4)
      bksks%old_kks = tmp(5)
      bksks%old_band = tmp(6)
      bksks%max_ks = tmp(7)
      bksks%max_kks = tmp(8)
      bksks%max_band = tmp(9)
      bksks%min_ks = tmp(10)
      bksks%min_kks = tmp(11)
      bksks%min_band = tmp(12)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE clear_bks(dirname,fname)
      !------------------------------------------------------------------------
      !
      USE mp_world,             ONLY : mpime,root
      USE west_io,              ONLY : remove_if_present
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(LEN=*),INTENT(IN) :: dirname,fname
      !
      IF(mpime == root) THEN
         !
         CALL remove_if_present(TRIM(dirname)//'/'//TRIM(fname))
         !
      ENDIF
      !
    END SUBROUTINE
    !
    ! SOLVEWFREQ
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvewfreq_restart_write_real(bks,dmat,zmat,npg,npl)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : my_image_id,me_bgrp,intra_bgrp_comm
      USE westcom,              ONLY : wfreq_restart_dir
      USE mp,                   ONLY : mp_sum
      USE distribution_center,  ONLY : ifr,rfr
      USE west_io,              ONLY : serial_data_write,remove_if_present
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(IN) :: bks
      INTEGER,INTENT(IN) :: npg,npl
      REAL(DP),INTENT(IN) :: dmat(npg,npl,ifr%nloc)
      COMPLEX(DP),INTENT(IN) :: zmat(npg,npl,rfr%nloc)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      CHARACTER(LEN=31) :: my_label
      CHARACTER(LEN=35) :: my_label2
      INTEGER :: ip,ip_glob
      REAL(DP),ALLOCATABLE :: tmp_dmat(:,:,:)
      COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:)
      LOGICAL :: lproc
      !
      CALL start_clock('sw_restart')
      !
      ! MKDIR
      !
      CALL my_mkdir(wfreq_restart_dir)
      !
      ! DMAT
      !
      ALLOCATE(tmp_dmat(npg,npl,ifr%nglob))
      tmp_dmat = 0._DP
      DO ip = 1,ifr%nloc
         ip_glob = ifr%l2g(ip)
         tmp_dmat(:,:,ip_glob) = dmat(:,:,ip)
      ENDDO
      CALL mp_sum(tmp_dmat,intra_bgrp_comm)
      !
      WRITE(my_label,'("dmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp == 0)
      CALL serial_data_write(lproc,fname,tmp_dmat,npg,npl,ifr%nglob)
      !
      DEALLOCATE(tmp_dmat)
      !
      ! ZMAT
      !
      ALLOCATE(tmp_zmat(npg,npl,rfr%nglob))
      tmp_zmat = 0._DP
      DO ip = 1,rfr%nloc
         ip_glob = rfr%l2g(ip)
         tmp_zmat(:,:,ip_glob) = zmat(:,:,ip)
      ENDDO
      CALL mp_sum(tmp_zmat,intra_bgrp_comm)
      !
      WRITE(my_label,'("zmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp == 0)
      CALL serial_data_write(lproc,fname,tmp_zmat,npg,npl,rfr%nglob)
      !
      DEALLOCATE(tmp_zmat)
      !
      ! CLEAR
      !
      IF(bks%old_ks /= 0 .AND. bks%old_band /= 0) THEN
         IF(me_bgrp == 0) THEN
            WRITE(my_label2,'("dmat_iks",i6.6,"_iv",i6.6,"_I",i6.6,".dat")') &
            & bks%old_ks,bks%old_band,my_image_id
            fname = TRIM(my_label2)
            CALL remove_if_present(TRIM(wfreq_restart_dir)//'/'//TRIM(fname))
            WRITE(my_label2,'("zmat_iks",i6.6,"_iv",i6.6,"_I",i6.6,".dat")') &
            & bks%old_ks,bks%old_band,my_image_id
            fname = TRIM(my_label2)
            CALL remove_if_present(TRIM(wfreq_restart_dir)//'/'//TRIM(fname))
         ENDIF
      ENDIF
      !
      ! CREATE THE SUMMARY FILE
      !
      fname = 'summary_w.json'
      CALL write_bks(bks,wfreq_restart_dir,fname)
      !
      CALL stop_clock('sw_restart')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvewfreq_restart_write_complex(bks,dmat,zmat,npg,npl)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : my_image_id,me_bgrp,intra_bgrp_comm
      USE westcom,              ONLY : wfreq_restart_dir
      USE mp,                   ONLY : mp_sum
      USE distribution_center,  ONLY : ifr,rfr
      USE west_io,              ONLY : serial_data_write,remove_if_present
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(IN) :: bks
      INTEGER,INTENT(IN) :: npg,npl
      COMPLEX(DP),INTENT(IN) :: dmat(npg,npl,ifr%nloc)
      COMPLEX(DP),INTENT(IN) :: zmat(npg,npl,rfr%nloc)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      CHARACTER(LEN=31) :: my_label
      CHARACTER(LEN=35) :: my_label2
      INTEGER :: ip,ip_glob
      COMPLEX(DP),ALLOCATABLE :: tmp_dmat(:,:,:)
      COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:)
      LOGICAL :: lproc
      !
      CALL start_clock('sw_restart')
      !
      ! MKDIR
      !
      CALL my_mkdir(wfreq_restart_dir)
      !
      ! DMAT
      !
      ALLOCATE(tmp_dmat(npg,npl,ifr%nglob))
      tmp_dmat = 0._DP
      DO ip = 1,ifr%nloc
         ip_glob = ifr%l2g(ip)
         tmp_dmat(:,:,ip_glob) = dmat(:,:,ip)
      ENDDO
      CALL mp_sum(tmp_dmat,intra_bgrp_comm)
      !
      WRITE(my_label,'("dmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp == 0)
      CALL serial_data_write(lproc,fname,tmp_dmat,npg,npl,ifr%nglob)
      !
      DEALLOCATE(tmp_dmat)
      !
      ! ZMAT
      !
      ALLOCATE(tmp_zmat(npg,npl,rfr%nglob))
      tmp_zmat = 0._DP
      DO ip = 1,rfr%nloc
         ip_glob = rfr%l2g(ip)
         tmp_zmat(:,:,ip_glob) = zmat(:,:,ip)
      ENDDO
      CALL mp_sum(tmp_zmat,intra_bgrp_comm)
      !
      WRITE(my_label,'("zmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp == 0)
      CALL serial_data_write(lproc,fname,tmp_zmat,npg,npl,rfr%nglob)
      !
      DEALLOCATE(tmp_zmat)
      !
      ! CLEAR
      !
      IF(bks%old_ks /= 0 .AND. bks%old_band /= 0) THEN
         IF(me_bgrp == 0) THEN
            WRITE(my_label2,'("dmat_iks",i6.6,"_iv",i6.6,"_I",i6.6,".dat")') &
            & bks%old_ks,bks%old_band,my_image_id
            fname = TRIM(my_label2)
            CALL remove_if_present(TRIM(wfreq_restart_dir)//'/'//TRIM(fname))
            WRITE(my_label2,'("zmat_iks",i6.6,"_iv",i6.6,"_I",i6.6,".dat")') &
            & bks%old_ks,bks%old_band,my_image_id
            fname = TRIM(my_label2)
            CALL remove_if_present(TRIM(wfreq_restart_dir)//'/'//TRIM(fname))
         ENDIF
      ENDIF
      !
      ! CREATE THE SUMMARY FILE
      !
      fname = 'summary_w.json'
      CALL write_bks(bks,wfreq_restart_dir,fname)
      !
      CALL stop_clock('sw_restart')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvewfreq_restart_write_complex_q(bksq,dmat,zmat,npg,npl)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : my_image_id,me_bgrp,intra_bgrp_comm
      USE westcom,              ONLY : wfreq_restart_dir
      USE mp,                   ONLY : mp_sum
      USE distribution_center,  ONLY : ifr,rfr
      USE west_io,              ONLY : serial_data_write,remove_if_present
      USE types_bz_grid,        ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksq_type),INTENT(IN) :: bksq
      INTEGER,INTENT(IN) :: npg,npl
      COMPLEX(DP),INTENT(IN) :: dmat(npg,npl,ifr%nloc,q_grid%np)
      COMPLEX(DP),INTENT(IN) :: zmat(npg,npl,rfr%nloc,q_grid%np)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      CHARACTER(LEN=40) :: my_label
      CHARACTER(LEN=44) :: my_label2
      INTEGER :: ip,ip_glob,iq
      COMPLEX(DP),ALLOCATABLE :: tmp_dmat(:,:,:,:)
      COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:,:)
      LOGICAL :: lproc
      !
      CALL start_clock('sw_restart')
      !
      ! MKDIR
      !
      CALL my_mkdir(wfreq_restart_dir)
      !
      ! DMAT
      !
      ALLOCATE(tmp_dmat(npg,npl,ifr%nglob,q_grid%np))
      tmp_dmat = 0._DP
      DO iq = 1,q_grid%np
         DO ip = 1,ifr%nloc
            ip_glob = ifr%l2g(ip)
            tmp_dmat(:,:,ip_glob,iq) = dmat(:,:,ip,iq)
         ENDDO
      ENDDO
      CALL mp_sum(tmp_dmat,intra_bgrp_comm)
      !
      WRITE(my_label,'("dmat_iq",i6.6,"_iks",i6.6,"_iv",i6.6,"_I",i6.6)') &
      & bksq%lastdone_q,bksq%lastdone_ks,bksq%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp == 0)
      CALL serial_data_write(lproc,fname,tmp_dmat,npg,npl,ifr%nglob,q_grid%np)
      !
      DEALLOCATE(tmp_dmat)
      !
      ! ZMAT
      !
      ALLOCATE(tmp_zmat(npg,npl,rfr%nglob,q_grid%np))
      tmp_zmat = 0._DP
      DO iq = 1,q_grid%np
         DO ip = 1,rfr%nloc
            ip_glob = rfr%l2g(ip)
            tmp_zmat(:,:,ip_glob,iq) = zmat(:,:,ip,iq)
         ENDDO
      ENDDO
      CALL mp_sum(tmp_zmat,intra_bgrp_comm)
      !
      WRITE(my_label,'("zmat_iq",i6.6,"_iks",i6.6,"_iv",i6.6,"_I",i6.6)') &
      & bksq%lastdone_q,bksq%lastdone_ks,bksq%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp == 0)
      CALL serial_data_write(lproc,fname,tmp_zmat,npg,npl,rfr%nglob,q_grid%np)
      !
      DEALLOCATE(tmp_zmat)
      !
      ! CLEAR
      !
      IF(bksq%old_q /= 0 .AND. bksq%old_ks /= 0 .AND. bksq%old_band /= 0) THEN
         IF(me_bgrp == 0) THEN
            WRITE(my_label2,'("dmat_iq",i6.6,"_iks",i6.6,"_iv",i6.6,"_I",i6.6,".dat")') &
            & bksq%old_q,bksq%old_ks,bksq%old_band,my_image_id
            fname = TRIM(my_label2)
            CALL remove_if_present(TRIM(wfreq_restart_dir)//'/'//TRIM(fname))
            WRITE(my_label2,'("zmat_iq",i6.6,"_iks",i6.6,"_iv",i6.6,"_I",i6.6,".dat")') &
            & bksq%old_q,bksq%old_ks,bksq%old_band,my_image_id
            fname = TRIM(my_label2)
            CALL remove_if_present(TRIM(wfreq_restart_dir)//'/'//TRIM(fname))
         ENDIF
      ENDIF
      !
      ! CREATE THE SUMMARY FILE
      !
      fname = 'summary_w.json'
      CALL write_bksq(bksq,wfreq_restart_dir,fname)
      !
      CALL stop_clock('sw_restart')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvewfreq_restart_read_real(bks,dmat,zmat,npg,npl)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : my_image_id,me_bgrp,intra_bgrp_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_restart_dir
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : ifr,rfr
      USE west_io,              ONLY : serial_data_read
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(OUT) :: bks
      INTEGER,INTENT(IN) :: npg,npl
      REAL(DP),INTENT(OUT) :: dmat(npg,npl,ifr%nloc)
      COMPLEX(DP),INTENT(OUT) :: zmat(npg,npl,rfr%nloc)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: get_clock
      REAL(DP) :: time_spent(2)
      CHARACTER(LEN=20),EXTERNAL :: human_readable_time
      CHARACTER(LEN=31) :: my_label
      INTEGER :: ip,ip_glob
      REAL(DP),ALLOCATABLE :: tmp_dmat(:,:,:)
      COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:)
      LOGICAL :: lproc
      !
      CALL start_clock('sw_restart')
      time_spent(1) = get_clock('sw_restart')
      !
      ! READ THE SUMMARY FILE
      !
      fname = 'summary_w.json'
      CALL read_bks(bks,wfreq_restart_dir,fname)
      !
      ! READ
      !
      ALLOCATE(tmp_dmat(npg,npl,ifr%nglob))
      !
      WRITE(my_label,'("dmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp==0)
      CALL serial_data_read(lproc,fname,tmp_dmat,npg,npl,ifr%nglob)
      !
      CALL mp_bcast(tmp_dmat,0,intra_bgrp_comm)
      DO ip = 1,ifr%nloc
         ip_glob = ifr%l2g(ip)
         dmat(:,:,ip) = tmp_dmat(:,:,ip_glob)
      ENDDO
      DEALLOCATE(tmp_dmat)
      !
      ALLOCATE(tmp_zmat(npg,npl,rfr%nglob))
      !
      WRITE(my_label,'("zmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp==0)
      CALL serial_data_read(lproc,fname,tmp_zmat,npg,npl,rfr%nglob)
      !
      CALL mp_bcast(tmp_zmat,0,intra_bgrp_comm)
      DO ip = 1,rfr%nloc
         ip_glob = rfr%l2g(ip)
         zmat(:,:,ip) = tmp_zmat(:,:,ip_glob)
      ENDDO
      DEALLOCATE(tmp_zmat)
      !
      time_spent(2) = get_clock('sw_restart')
      CALL stop_clock('sw_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wfreq_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvewfreq_restart_read_complex(bks,dmat,zmat,npg,npl)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : my_image_id,me_bgrp,intra_bgrp_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_restart_dir
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : ifr,rfr
      USE west_io,              ONLY : serial_data_read
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(OUT) :: bks
      INTEGER,INTENT(IN) :: npg,npl
      COMPLEX(DP),INTENT(OUT) :: dmat(npg,npl,ifr%nloc)
      COMPLEX(DP),INTENT(OUT) :: zmat(npg,npl,rfr%nloc)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: get_clock
      REAL(DP) :: time_spent(2)
      CHARACTER(LEN=20),EXTERNAL :: human_readable_time
      CHARACTER(LEN=31) :: my_label
      INTEGER :: ip,ip_glob
      COMPLEX(DP),ALLOCATABLE :: tmp_dmat(:,:,:)
      COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:)
      LOGICAL :: lproc
      !
      CALL start_clock('sw_restart')
      time_spent(1) = get_clock('sw_restart')
      !
      ! READ THE SUMMARY FILE
      !
      fname = 'summary_w.json'
      CALL read_bks(bks,wfreq_restart_dir,fname)
      !
      ! READ
      !
      ALLOCATE(tmp_dmat(npg,npl,ifr%nglob))
      !
      WRITE(my_label,'("dmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp==0)
      CALL serial_data_read(lproc,fname,tmp_dmat,npg,npl,ifr%nglob)
      !
      CALL mp_bcast(tmp_dmat,0,intra_bgrp_comm)
      DO ip = 1,ifr%nloc
         ip_glob = ifr%l2g(ip)
         dmat(:,:,ip) = tmp_dmat(:,:,ip_glob)
      ENDDO
      DEALLOCATE(tmp_dmat)
      !
      ALLOCATE(tmp_zmat(npg,npl,rfr%nglob))
      !
      WRITE(my_label,'("zmat_iks",i6.6,"_iv",i6.6,"_I",i6.6)') bks%lastdone_ks,bks%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp==0)
      CALL serial_data_read(lproc,fname,tmp_zmat,npg,npl,rfr%nglob)
      !
      CALL mp_bcast(tmp_zmat,0,intra_bgrp_comm)
      DO ip = 1,rfr%nloc
         ip_glob = rfr%l2g(ip)
         zmat(:,:,ip) = tmp_zmat(:,:,ip_glob)
      ENDDO
      DEALLOCATE(tmp_zmat)
      !
      time_spent(2) = get_clock('sw_restart')
      CALL stop_clock('sw_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wfreq_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvewfreq_restart_read_complex_q(bksq,dmat,zmat,npg,npl)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : my_image_id,me_bgrp,intra_bgrp_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_restart_dir
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : ifr,rfr
      USE west_io,              ONLY : serial_data_read
      USE types_bz_grid,        ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksq_type),INTENT(OUT) :: bksq
      INTEGER,INTENT(IN) :: npg,npl
      COMPLEX(DP),INTENT(OUT) :: dmat(npg,npl,ifr%nloc,q_grid%np)
      COMPLEX(DP),INTENT(OUT) :: zmat(npg,npl,rfr%nloc,q_grid%np)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: get_clock
      REAL(DP) :: time_spent(2)
      CHARACTER(LEN=20),EXTERNAL :: human_readable_time
      CHARACTER(LEN=40) :: my_label
      INTEGER :: ip,ip_glob,iq
      COMPLEX(DP),ALLOCATABLE :: tmp_dmat(:,:,:,:)
      COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:,:)
      LOGICAL :: lproc
      !
      CALL start_clock('sw_restart')
      time_spent(1) = get_clock('sw_restart')
      !
      ! READ THE SUMMARY FILE
      !
      fname = 'summary_w.json'
      CALL read_bksq(bksq,wfreq_restart_dir,fname)
      !
      ! READ
      !
      ALLOCATE(tmp_dmat(npg,npl,ifr%nglob,q_grid%np))
      !
      WRITE(my_label,'("dmat_iq",i6.6,"_iks",i6.6,"_iv",i6.6,"_I",i6.6)') &
      & bksq%lastdone_q,bksq%lastdone_ks,bksq%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp==0)
      CALL serial_data_read(lproc,fname,tmp_dmat,npg,npl,ifr%nglob,q_grid%np)
      !
      CALL mp_bcast(tmp_dmat,0,intra_bgrp_comm)
      DO iq = 1,q_grid%np
         DO ip = 1,ifr%nloc
            ip_glob = ifr%l2g(ip)
            dmat(:,:,ip,iq) = tmp_dmat(:,:,ip_glob,iq)
         ENDDO
      ENDDO
      DEALLOCATE(tmp_dmat)
      !
      ALLOCATE(tmp_zmat(npg,npl,rfr%nglob,q_grid%np))
      !
      WRITE(my_label,'("zmat_iq",i6.6,"_iks",i6.6,"_iv",i6.6,"_I",i6.6)') &
      & bksq%lastdone_q,bksq%lastdone_ks,bksq%lastdone_band,my_image_id
      fname = TRIM(wfreq_restart_dir)//'/'//TRIM(my_label)
      lproc = (me_bgrp==0)
      CALL serial_data_read(lproc,fname,tmp_zmat,npg,npl,rfr%nglob,q_grid%np)
      !
      CALL mp_bcast(tmp_zmat,0,intra_bgrp_comm)
      DO iq = 1,q_grid%np
         DO ip = 1,rfr%nloc
            ip_glob = rfr%l2g(ip)
            zmat(:,:,ip,iq) = tmp_zmat(:,:,ip_glob,iq)
         ENDDO
      ENDDO
      DEALLOCATE(tmp_zmat)
      !
      time_spent(2) = get_clock('sw_restart')
      CALL stop_clock('sw_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wfreq_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    ! GFREQ
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvegfreq_restart_write(bks)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : wfreq_restart_dir
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(IN) :: bks
      !
      ! Workspace
      !
      CHARACTER(LEN=*),PARAMETER :: fname = 'summary_g.json'
      !
      ! MKDIR
      !
      CALL my_mkdir(wfreq_restart_dir)
      !
      CALL start_clock('sg_restart')
      !
      CALL write_bks(bks,wfreq_restart_dir,fname)
      !
      CALL stop_clock('sg_restart')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvegfreq_restart_write_q(bksks)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : wfreq_restart_dir
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksks_type),INTENT(IN) :: bksks
      !
      ! Workspace
      !
      CHARACTER(LEN=*),PARAMETER :: fname = 'summary_g.json'
      !
      ! MKDIR
      !
      CALL my_mkdir(wfreq_restart_dir)
      !
      CALL start_clock('sg_restart')
      !
      CALL write_bksks(bksks,wfreq_restart_dir,fname)
      !
      CALL stop_clock('sg_restart')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvegfreq_restart_read(bks)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_restart_dir
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bks_type),INTENT(OUT) :: bks
      !
      ! Workspace
      !
      CHARACTER(LEN=*),PARAMETER :: fname = 'summary_g.json'
      REAL(DP),EXTERNAL :: get_clock
      REAL(DP) :: time_spent(2)
      CHARACTER(LEN=20),EXTERNAL :: human_readable_time
      !
      CALL start_clock('sg_restart')
      time_spent(1) = get_clock('sg_restart')
      !
      CALL read_bks(bks,wfreq_restart_dir,fname)
      !
      time_spent(2) = get_clock('sg_restart')
      CALL stop_clock('sg_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wfreq_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE solvegfreq_restart_read_q(bksks)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_restart_dir
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bksks_type),INTENT(OUT) :: bksks
      !
      ! Workspace
      !
      CHARACTER(LEN=*),PARAMETER :: fname = 'summary_g.json'
      REAL(DP),EXTERNAL :: get_clock
      REAL(DP) :: time_spent(2)
      CHARACTER(LEN=20),EXTERNAL :: human_readable_time
      !
      CALL start_clock('sg_restart')
      time_spent(1) = get_clock('sg_restart')
      !
      CALL read_bksks(bksks,wfreq_restart_dir,fname)
      !
      time_spent(2) = get_clock('sg_restart')
      CALL stop_clock('sg_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wfreq_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
END MODULE
