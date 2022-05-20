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
MODULE wfreq_io
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTERFACE writeout_overlap
     MODULE PROCEDURE writeout_overlap_real, writeout_overlap_complex, writeout_overlap_complex_q
  END INTERFACE
  !
  INTERFACE readin_overlap
     MODULE PROCEDURE readin_overlap_real, readin_overlap_complex, readin_overlap_complex_q
  END INTERFACE
  !
  INTERFACE writeout_solvegfreq
     MODULE PROCEDURE writeout_solvegfreq_real, writeout_solvegfreq_complex, writeout_solvegfreq_complex_q
  END INTERFACE
  !
  INTERFACE readin_solvegfreq
     MODULE PROCEDURE readin_solvegfreq_real, readin_solvegfreq_complex, readin_solvegfreq_complex_q
  END INTERFACE
  !
  CONTAINS
  !
  SUBROUTINE writeout_solvegfreq_real(glob_iks,glob_ib,diago,braket,nloc,nglob,myoffset)
    !
    ! WHO WRITES? ... root_bgrp
    ! WHAT?       ... diago(        n_lanczos )
    ! WHAT?       ... braket( ndim, n_lanczos )
    !
    USE kinds,               ONLY : DP,i8b
    USE westcom,             ONLY : n_lanczos,wfreq_save_dir,l_enable_off_diagonal
    USE mod_mpiio,           ONLY : mp_write_dmsg_at
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(IN) :: diago(n_lanczos,nloc)
    REAL(DP),INTENT(IN) :: braket(nglob,n_lanczos,nloc)
    INTEGER,INTENT(IN) :: nloc,nglob,myoffset
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    INTEGER(i8b) :: myoffset_i8b
    !
    CALL start_clock('write_gfreq')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    IF (l_enable_off_diagonal) THEN
      myoffset_i8b = 1_i8b*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'B'//c_glob_ib//'full.dat'
      CALL mp_write_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
      !
      myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'B'//c_glob_ib//'full.dat'
      CALL mp_write_dmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    ELSE
      myoffset_i8b = 1_i8b*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
      CALL mp_write_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
      !
      myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
      CALL mp_write_dmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    ENDIF
    !
    CALL stop_clock('write_gfreq')
    !
  END SUBROUTINE
  !
  SUBROUTINE writeout_solvegfreq_complex(glob_iks,glob_ib,diago,braket,nloc,nglob,myoffset)
    !
    ! WHO WRITES? ... root_bgrp
    ! WHAT?       ... diago(        n_lanczos )
    ! WHAT?       ... braket( ndim, n_lanczos )
    !
    USE kinds,               ONLY : DP,i8b
    USE westcom,             ONLY : n_lanczos,wfreq_save_dir
    USE mod_mpiio,           ONLY : mp_write_dmsg_at,mp_write_zmsg_at
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(IN) :: diago(n_lanczos,nloc)
    COMPLEX(DP),INTENT(IN) :: braket(nglob,n_lanczos,nloc)
    INTEGER,INTENT(IN) :: nloc,nglob,myoffset
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    INTEGER(i8b) :: myoffset_i8b
    !
    CALL start_clock('write_gfreq')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    myoffset_i8b = 1_i8b*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
    CALL mp_write_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
    !
    myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
    CALL mp_write_zmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    !
    CALL stop_clock('write_gfreq')
    !
  END SUBROUTINE
  !
  SUBROUTINE writeout_solvegfreq_complex_q(glob_iks,glob_ikks,glob_ib,diago,braket,nloc,nglob,myoffset)
    !
    ! WHO WRITES? ... root_bgrp
    ! WHAT?       ... diago(        n_lanczos )
    ! WHAT?       ... braket( ndim, n_lanczos )
    !
    USE kinds,               ONLY : DP,i8b
    USE westcom,             ONLY : n_lanczos,wfreq_save_dir
    USE mod_mpiio,           ONLY : mp_write_dmsg_at,mp_write_zmsg_at
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ikks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(IN) :: diago(n_lanczos,nloc)
    COMPLEX(DP),INTENT(IN) :: braket(nglob,n_lanczos,nloc)
    INTEGER,INTENT(IN) :: nloc,nglob,myoffset
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ikks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    INTEGER(i8b) :: myoffset_i8b
    !
    CALL start_clock('write_gfreq')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ikks,'(i5.5)') glob_ikks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    myoffset_i8b = 1_i8b*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib//'.dat'
    CALL mp_write_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
    !
    myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib//'.dat'
    CALL mp_write_zmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    !
    CALL stop_clock('write_gfreq')
    !
  END SUBROUTINE
  !
  SUBROUTINE readin_solvegfreq_real(glob_iks,glob_ib,diago,braket,nloc,nglob,myoffset)
    !
    ! WHO READS?  ... root_bgrp
    ! WHAT?       ... diago(        n_lanczos )
    ! WHAT?       ... braket( ndim, n_lanczos )
    !
    USE kinds,          ONLY : DP,i8b
    USE westcom,        ONLY : n_lanczos,wfreq_save_dir,l_enable_off_diagonal
    USE mod_mpiio,      ONLY : mp_read_dmsg_at
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(OUT) :: diago(n_lanczos,nloc)
    REAL(DP),INTENT(OUT) :: braket(nglob,n_lanczos,nloc)
    INTEGER,INTENT(IN) :: nloc,nglob,myoffset
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    INTEGER(i8b) :: myoffset_i8b
    !
    CALL start_clock('read_gfreq')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    IF (l_enable_off_diagonal) THEN
      myoffset_i8b = 1_i8b*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'B'//c_glob_ib//'_full.dat'
      CALL mp_read_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
      myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'B'//c_glob_ib//'_full.dat'
      CALL mp_read_dmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    ELSE
      myoffset_i8b = 1_i8b*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
      CALL mp_read_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
      myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
      fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
      CALL mp_read_dmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    ENDIF
    !
    CALL stop_clock('read_gfreq')
    !
  END SUBROUTINE
  !
  SUBROUTINE readin_solvegfreq_complex(glob_iks,glob_ib,diago,braket,nloc,nglob,myoffset)
    !
    ! WHO READS?  ... root_bgrp
    ! WHAT?       ... diago(        n_lanczos )
    ! WHAT?       ... braket( ndim, n_lanczos )
    !
    USE kinds,          ONLY : DP,i8b
    USE westcom,        ONLY : n_lanczos,wfreq_save_dir
    USE mod_mpiio,      ONLY : mp_read_dmsg_at,mp_read_zmsg_at
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(OUT) :: diago(n_lanczos,nloc)
    COMPLEX(DP),INTENT(OUT) :: braket(nglob,n_lanczos,nloc)
    INTEGER,INTENT(IN) :: nloc,nglob,myoffset
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    INTEGER(i8b) :: myoffset_i8b
    !
    CALL start_clock('read_gfreq')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    myoffset_i8b = 1_i8b*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
    CALL mp_read_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
    !
    myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'B'//c_glob_ib//'.dat'
    CALL mp_read_zmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    !
    CALL stop_clock('read_gfreq')
    !
  END SUBROUTINE
  !
  SUBROUTINE readin_solvegfreq_complex_q(glob_iks,glob_ikks,glob_ib,diago,braket,nloc,nglob,myoffset)
    !
    ! WHO READS?  ... root_bgrp
    ! WHAT?       ... diago(        n_lanczos )
    ! WHAT?       ... braket( ndim, n_lanczos )
    !
    USE kinds,          ONLY : DP,i8b
    USE westcom,        ONLY : n_lanczos,wfreq_save_dir
    USE mod_mpiio,      ONLY : mp_read_dmsg_at,mp_read_zmsg_at
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ikks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(OUT) :: diago(n_lanczos,nloc)
    COMPLEX(DP),INTENT(OUT) :: braket(nglob,n_lanczos,nloc)
    INTEGER,INTENT(IN) :: nloc,nglob,myoffset
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ikks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    INTEGER(i8b) :: myoffset_i8b
    !
    CALL start_clock('read_gfreq')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ikks,'(i5.5)') glob_ikks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    myoffset_i8b = 1_i8b*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_diag_K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib//'.dat'
    CALL mp_read_dmsg_at(fname,diago,n_lanczos*nloc,myoffset_i8b)
    !
    myoffset_i8b = 1_i8b*nglob*n_lanczos*myoffset
    fname = TRIM(wfreq_save_dir)//'/g_brak_K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib//'.dat'
    CALL mp_read_zmsg_at(fname,braket,nglob*n_lanczos*nloc,myoffset_i8b)
    !
    CALL stop_clock('read_gfreq')
    !
  END SUBROUTINE
  !
  SUBROUTINE writeout_overlap_real(labellina,glob_iks,glob_ib,overlap,no1,no2)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... overlap( no1, no2 )
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_global,      ONLY : my_image_id,root_image,me_bgrp,root_bgrp
    USE west_io,        ONLY : serial_data_write
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: labellina
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(IN) :: overlap(no1,no2)
    INTEGER,INTENT(IN) :: no1,no2
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('write_over')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//labellina//'_K'//c_glob_iks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_write(lproc,fname,overlap,no1,no2)
    !
    CALL stop_clock('write_over')
    !
  END SUBROUTINE
  !
  SUBROUTINE writeout_overlap_complex(labellina,glob_iks,glob_ib,overlap,no1,no2)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... overlap( no1, no2 )
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_global,      ONLY : my_image_id,root_image,me_bgrp,root_bgrp
    USE west_io,        ONLY : serial_data_write
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: labellina
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    COMPLEX(DP),INTENT(IN) :: overlap(no1,no2)
    INTEGER,INTENT(IN) :: no1,no2
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('write_over')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//labellina//'_K'//c_glob_iks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_write(lproc,fname,overlap,no1,no2)
    !
    CALL stop_clock('write_over')
    !
  END SUBROUTINE
  !
  SUBROUTINE writeout_overlap_complex_q(labellina,glob_iks,glob_ikks,glob_ib,overlap,no1,no2)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... overlap( no1, no2 )
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_global,      ONLY : my_image_id,root_image,me_bgrp,root_bgrp
    USE west_io,        ONLY : serial_data_write
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: labellina
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ikks
    INTEGER,INTENT(IN) :: glob_ib
    COMPLEX(DP),INTENT(IN) :: overlap(no1,no2)
    INTEGER,INTENT(IN) :: no1,no2
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ikks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('write_over')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ikks,'(i5.5)') glob_ikks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//labellina//'K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_write(lproc,fname,overlap,no1,no2)
    !
    CALL stop_clock('write_over')
    !
  END SUBROUTINE
  !
  SUBROUTINE readin_overlap_real(labellina,glob_iks,glob_ib,overlap,no1,no2)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... overlap( no1, no2 )
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_global,      ONLY : my_image_id,root_image,inter_image_comm,me_bgrp,&
                             & root_bgrp,intra_bgrp_comm
    USE mp,             ONLY : mp_bcast
    USE west_io,        ONLY : serial_data_read
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: labellina
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    REAL(DP),INTENT(OUT) :: overlap(no1,no2)
    INTEGER,INTENT(IN) :: no1,no2
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('read_over')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//labellina//'_K'//c_glob_iks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_read(lproc,fname,overlap,no1,no2)
    !
    CALL mp_bcast(overlap,root_bgrp,intra_bgrp_comm)
    CALL mp_bcast(overlap,root_image,inter_image_comm)
    !
    CALL stop_clock('read_over')
    !
  END SUBROUTINE
  !
  SUBROUTINE readin_overlap_complex(labellina,glob_iks,glob_ib,overlap,no1,no2)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... overlap( no1, no2 )
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_global,      ONLY : my_image_id,root_image,inter_image_comm,me_bgrp,&
                             & root_bgrp,intra_bgrp_comm
    USE mp,             ONLY : mp_bcast
    USE west_io,        ONLY : serial_data_read
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: labellina
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ib
    COMPLEX(DP),INTENT(OUT) :: overlap(no1,no2)
    INTEGER,INTENT(IN) :: no1,no2
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('read_over')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//labellina//'_K'//c_glob_iks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_read(lproc,fname,overlap,no1,no2)
    !
    CALL mp_bcast(overlap,root_bgrp,intra_bgrp_comm)
    CALL mp_bcast(overlap,root_image,inter_image_comm)
    !
    CALL stop_clock('read_over')
    !
  END SUBROUTINE
  !
  SUBROUTINE readin_overlap_complex_q(labellina,glob_iks,glob_ikks,glob_ib,overlap,no1,no2)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... overlap( no1, no2 )
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_global,      ONLY : my_image_id,root_image,inter_image_comm,me_bgrp,&
                             & root_bgrp,intra_bgrp_comm
    USE mp,             ONLY : mp_bcast
    USE west_io,        ONLY : serial_data_read
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: labellina
    INTEGER,INTENT(IN) :: glob_iks
    INTEGER,INTENT(IN) :: glob_ikks
    INTEGER,INTENT(IN) :: glob_ib
    COMPLEX(DP),INTENT(OUT) :: overlap(no1,no2)
    INTEGER,INTENT(IN) :: no1,no2
    !
    ! Workspace
    !
    CHARACTER(LEN=5) :: c_glob_iks
    CHARACTER(LEN=5) :: c_glob_ikks
    CHARACTER(LEN=5) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('read_over')
    !
    ! Generate the filename
    !
    WRITE(c_glob_iks,'(i5.5)') glob_iks
    WRITE(c_glob_ikks,'(i5.5)') glob_ikks
    WRITE(c_glob_ib,'(i5.5)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//labellina//'K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_read(lproc,fname,overlap,no1,no2)
    !
    CALL mp_bcast(overlap,root_bgrp,intra_bgrp_comm)
    CALL mp_bcast(overlap,root_image,inter_image_comm)
    !
    CALL stop_clock('read_over')
    !
  END SUBROUTINE
  !
  SUBROUTINE writeout_solvehf(sigma_hf,nb,nk)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... sigma_hf(nb, nk)
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_world,       ONLY : mpime,root
    USE west_io,        ONLY : serial_data_write
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    REAL(DP),INTENT(IN) :: sigma_hf(nb,nk)
    INTEGER,INTENT(IN) :: nb,nk
    !
    ! Workspace
    !
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('write_hf')
    !
    ! Generate the filename
    !
    fname = TRIM(wfreq_save_dir)//'/hf'
    lproc = (mpime == root)
    !
    CALL serial_data_write(lproc,fname,sigma_hf,nb,nk)
    !
    CALL stop_clock('write_hf')
    !
  END SUBROUTINE
  !
  SUBROUTINE readin_solvehf(sigma_hf,nb,nk)
    !
    ! WHO WRITES? ... root
    ! WHAT?       ... sigma_hf(nb, nk)
    !
    USE kinds,          ONLY : DP
    USE westcom,        ONLY : wfreq_save_dir
    USE mp_world,       ONLY : mpime,root,world_comm
    USE mp,             ONLY : mp_bcast
    USE west_io,        ONLY : serial_data_read
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    REAL(DP),INTENT(OUT) :: sigma_hf(nb,nk)
    INTEGER,INTENT(IN) :: nb,nk
    !
    ! Workspace
    !
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('read_hf')
    !
    ! Generate the filename
    !
    fname = TRIM(wfreq_save_dir)//'/hf'
    lproc = (mpime == root)
    !
    CALL serial_data_read(lproc,fname,sigma_hf,nb,nk)
    !
    CALL mp_bcast(sigma_hf,root,world_comm)
    !
    CALL stop_clock('read_hf')
    !
  END SUBROUTINE
END MODULE
