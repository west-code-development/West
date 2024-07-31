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
MODULE wfreq_io
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTERFACE write_overlap
     MODULE PROCEDURE write_overlap_real, write_overlap_complex_q
  END INTERFACE
  !
  INTERFACE read_overlap
     MODULE PROCEDURE read_overlap_real, read_overlap_complex_q
  END INTERFACE
  !
  INTERFACE write_wfreq
     MODULE PROCEDURE write_wfreq_real, write_wfreq_complex_q
  END INTERFACE
  !
  INTERFACE read_wfreq
     MODULE PROCEDURE read_wfreq_real, read_wfreq_complex_q
  END INTERFACE
  !
  INTERFACE write_gfreq
     MODULE PROCEDURE write_gfreq_real, write_gfreq_complex_q
  END INTERFACE
  !
  INTERFACE read_gfreq
     MODULE PROCEDURE read_gfreq_real, read_gfreq_complex_q
  END INTERFACE
  !
  CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_overlap_real(label,glob_iks,glob_ib,overlap,no1,no2)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,root_image,me_bgrp,root_bgrp
    USE west_io,              ONLY : serial_data_write
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: label
    INTEGER,INTENT(IN) :: glob_iks,glob_ib
    INTEGER,INTENT(IN) :: no1,no2
    REAL(DP),INTENT(IN) :: overlap(no1,no2)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_glob_iks
    CHARACTER(LEN=9) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('write_over')
    !
    WRITE(c_glob_iks,'(i6.6)') glob_iks
    WRITE(c_glob_ib,'(i9.9)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//label//'_K'//c_glob_iks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_write(lproc,fname,overlap,no1,no2)
    !
    CALL stop_clock('write_over')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_overlap_complex_q(label,glob_iks,glob_ikks,glob_ib,overlap,no1,no2)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,root_image,me_bgrp,root_bgrp
    USE west_io,              ONLY : serial_data_write
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: label
    INTEGER,INTENT(IN) :: glob_iks,glob_ikks,glob_ib
    INTEGER,INTENT(IN) :: no1,no2
    COMPLEX(DP),INTENT(IN) :: overlap(no1,no2)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_glob_iks
    CHARACTER(LEN=6) :: c_glob_ikks
    CHARACTER(LEN=9) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('write_over')
    !
    WRITE(c_glob_iks,'(i6.6)') glob_iks
    WRITE(c_glob_ikks,'(i6.6)') glob_ikks
    WRITE(c_glob_ib,'(i9.9)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//label//'K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib
    lproc = (my_image_id == root_image .AND. me_bgrp == root_bgrp)
    !
    CALL serial_data_write(lproc,fname,overlap,no1,no2)
    !
    CALL stop_clock('write_over')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_overlap_real(label,glob_iks,glob_ib,overlap,no1,no2)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,root_image,inter_image_comm,me_bgrp,&
                                   & root_bgrp,intra_bgrp_comm
    USE mp,                   ONLY : mp_bcast
    USE west_io,              ONLY : serial_data_read
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: label
    INTEGER,INTENT(IN) :: glob_iks,glob_ib
    INTEGER,INTENT(IN) :: no1,no2
    REAL(DP),INTENT(OUT) :: overlap(no1,no2)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_glob_iks
    CHARACTER(LEN=9) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('read_over')
    !
    WRITE(c_glob_iks,'(i6.6)') glob_iks
    WRITE(c_glob_ib,'(i9.9)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//label//'_K'//c_glob_iks//'B'//c_glob_ib
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
  !------------------------------------------------------------------------
  SUBROUTINE read_overlap_complex_q(label,glob_iks,glob_ikks,glob_ib,overlap,no1,no2)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,root_image,inter_image_comm,me_bgrp,&
                                   & root_bgrp,intra_bgrp_comm
    USE mp,                   ONLY : mp_bcast
    USE west_io,              ONLY : serial_data_read
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    CHARACTER(LEN=1),INTENT(IN) :: label
    INTEGER,INTENT(IN) :: glob_iks,glob_ikks,glob_ib
    INTEGER,INTENT(IN) :: no1,no2
    COMPLEX(DP),INTENT(OUT) :: overlap(no1,no2)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_glob_iks
    CHARACTER(LEN=6) :: c_glob_ikks
    CHARACTER(LEN=9) :: c_glob_ib
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('read_over')
    !
    WRITE(c_glob_iks,'(i6.6)') glob_iks
    WRITE(c_glob_ikks,'(i6.6)') glob_ikks
    WRITE(c_glob_ib,'(i9.9)') glob_ib
    !
    fname = TRIM(wfreq_save_dir)//'/over_'//label//'K'//c_glob_iks//'KK'//c_glob_ikks//'B'//c_glob_ib
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
  !------------------------------------------------------------------------
  SUBROUTINE write_hf(sigma_hf,nb)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_world,             ONLY : mpime,root
    USE west_io,              ONLY : serial_data_write
    USE types_bz_grid,        ONLY : k_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: nb
    REAL(DP),INTENT(IN) :: sigma_hf(nb,k_grid%nps)
    !
    ! Workspace
    !
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('write_hf')
    !
    fname = TRIM(wfreq_save_dir)//'/hf'
    lproc = (mpime == root)
    !
    CALL serial_data_write(lproc,fname,sigma_hf,nb,k_grid%nps)
    !
    CALL stop_clock('write_hf')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_hf(sigma_hf,nb)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_world,             ONLY : mpime,root,world_comm
    USE mp,                   ONLY : mp_bcast
    USE west_io,              ONLY : serial_data_read
    USE types_bz_grid,        ONLY : k_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: nb
    REAL(DP),INTENT(OUT) :: sigma_hf(nb,k_grid%nps)
    !
    ! Workspace
    !
    CHARACTER(LEN=512) :: fname
    LOGICAL :: lproc
    !
    CALL start_clock('read_hf')
    !
    fname = TRIM(wfreq_save_dir)//'/hf'
    lproc = (mpime == root)
    !
    CALL serial_data_read(lproc,fname,sigma_hf,nb,k_grid%nps)
    !
    CALL mp_bcast(sigma_hf,root,world_comm)
    !
    CALL stop_clock('read_hf')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_wfreq_real(dmat,zmat,npg,npl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,my_pool_id,my_bgrp_id,me_bgrp,intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE distribution_center,  ONLY : ifr,rfr
    USE west_io,              ONLY : serial_data_write
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npg,npl
    REAL(DP),INTENT(IN) :: dmat(npg,npl,ifr%nloc)
    COMPLEX(DP),INTENT(IN) :: zmat(npg,npl,rfr%nloc)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_image
    CHARACTER(LEN=512) :: fname
    INTEGER :: ip,ip_glob
    REAL(DP),ALLOCATABLE :: tmp_dmat(:,:,:)
    COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:)
    LOGICAL :: lproc
    !
    CALL start_clock('write_w')
    !
    IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
       !
       WRITE(c_image,'(i6.6)') my_image_id
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
       fname = TRIM(wfreq_save_dir)//'/dmat_I'//c_image
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
       fname = TRIM(wfreq_save_dir)//'/zmat_I'//c_image
       lproc = (me_bgrp == 0)
       CALL serial_data_write(lproc,fname,tmp_zmat,npg,npl,rfr%nglob)
       !
       DEALLOCATE(tmp_zmat)
       !
    ENDIF
    !
    CALL stop_clock('write_w')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_wfreq_complex_q(dmat,zmat,npg,npl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,my_pool_id,my_bgrp_id,me_bgrp,intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE distribution_center,  ONLY : ifr,rfr
    USE west_io,              ONLY : serial_data_write
    USE types_bz_grid,        ONLY : q_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npg,npl
    COMPLEX(DP),INTENT(IN) :: dmat(npg,npl,ifr%nloc,q_grid%np)
    COMPLEX(DP),INTENT(IN) :: zmat(npg,npl,rfr%nloc,q_grid%np)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_image
    CHARACTER(LEN=512) :: fname
    INTEGER :: ip,ip_glob,iq
    COMPLEX(DP),ALLOCATABLE :: tmp_dmat(:,:,:,:)
    COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:,:)
    LOGICAL :: lproc
    !
    CALL start_clock('write_w')
    !
    IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
       !
       WRITE(c_image,'(i6.6)') my_image_id
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
       fname = TRIM(wfreq_save_dir)//'/dmat_I'//c_image
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
       fname = TRIM(wfreq_save_dir)//'/zmat_I'//c_image
       lproc = (me_bgrp == 0)
       CALL serial_data_write(lproc,fname,tmp_zmat,npg,npl,rfr%nglob,q_grid%np)
       !
       DEALLOCATE(tmp_zmat)
       !
    ENDIF
    !
    CALL stop_clock('write_w')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_wfreq_real(dmat,zmat,npg,npl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,my_pool_id,inter_pool_comm,my_bgrp_id,me_bgrp,&
                                   & inter_bgrp_comm,intra_bgrp_comm
    USE mp,                   ONLY : mp_bcast
    USE distribution_center,  ONLY : ifr,rfr
    USE west_io,              ONLY : serial_data_read
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npg,npl
    REAL(DP),INTENT(OUT) :: dmat(npg,npl,ifr%nloc)
    COMPLEX(DP),INTENT(OUT) :: zmat(npg,npl,rfr%nloc)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_image
    CHARACTER(LEN=512) :: fname
    INTEGER :: ip,ip_glob
    REAL(DP),ALLOCATABLE :: tmp_dmat(:,:,:)
    COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:)
    LOGICAL :: lproc
    !
    CALL start_clock('read_w')
    !
    IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
       !
       WRITE(c_image,'(i6.6)') my_image_id
       !
       ! DMAT
       !
       ALLOCATE(tmp_dmat(npg,npl,ifr%nglob))
       !
       fname = TRIM(wfreq_save_dir)//'/dmat_I'//c_image
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
       ! ZMAT
       !
       ALLOCATE(tmp_zmat(npg,npl,rfr%nglob))
       !
       fname = TRIM(wfreq_save_dir)//'/zmat_I'//c_image
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
    ENDIF
    !
    CALL mp_bcast(dmat,0,inter_bgrp_comm)
    CALL mp_bcast(dmat,0,inter_pool_comm)
    CALL mp_bcast(zmat,0,inter_bgrp_comm)
    CALL mp_bcast(zmat,0,inter_pool_comm)
    !
    CALL stop_clock('read_w')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_wfreq_complex_q(dmat,zmat,npg,npl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : wfreq_save_dir
    USE mp_global,            ONLY : my_image_id,my_pool_id,inter_pool_comm,my_bgrp_id,me_bgrp,&
                                   & inter_bgrp_comm,intra_bgrp_comm
    USE mp,                   ONLY : mp_bcast
    USE distribution_center,  ONLY : ifr,rfr
    USE west_io,              ONLY : serial_data_read
    USE types_bz_grid,        ONLY : q_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npg,npl
    COMPLEX(DP),INTENT(OUT) :: dmat(npg,npl,ifr%nloc,q_grid%np)
    COMPLEX(DP),INTENT(OUT) :: zmat(npg,npl,rfr%nloc,q_grid%np)
    !
    ! Workspace
    !
    CHARACTER(LEN=6) :: c_image
    CHARACTER(LEN=512) :: fname
    INTEGER :: ip,ip_glob,iq
    COMPLEX(DP),ALLOCATABLE :: tmp_dmat(:,:,:,:)
    COMPLEX(DP),ALLOCATABLE :: tmp_zmat(:,:,:,:)
    LOGICAL :: lproc
    !
    CALL start_clock('read_w')
    !
    IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
       !
       WRITE(c_image,'(i6.6)') my_image_id
       !
       ! DMAT
       !
       ALLOCATE(tmp_dmat(npg,npl,ifr%nglob,q_grid%np))
       !
       fname = TRIM(wfreq_save_dir)//'/dmat_I'//c_image
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
       ! ZMAT
       !
       ALLOCATE(tmp_zmat(npg,npl,rfr%nglob,q_grid%np))
       !
       fname = TRIM(wfreq_save_dir)//'/zmat_I'//c_image
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
    ENDIF
    !
    CALL mp_bcast(dmat,0,inter_bgrp_comm)
    CALL mp_bcast(dmat,0,inter_pool_comm)
    CALL mp_bcast(zmat,0,inter_bgrp_comm)
    CALL mp_bcast(zmat,0,inter_pool_comm)
    !
    CALL stop_clock('read_w')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_gfreq_real(d_diago,d_body2,npl,nbl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : n_lanczos
    USE distribution_center,  ONLY : ifr
    USE types_bz_grid,        ONLY : k_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npl,nbl
    REAL(DP),INTENT(IN) :: d_diago(n_lanczos,npl,nbl,k_grid%nps)
    REAL(DP),INTENT(IN) :: d_body2(n_lanczos,npl,ifr%nloc,nbl,k_grid%nps)
    !
    CALL start_clock('write_g')
    !
    CALL stop_clock('write_g')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_gfreq_complex_q(d_diago,z_body2,npl,nbl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : n_lanczos
    USE distribution_center,  ONLY : ifr
    USE types_bz_grid,        ONLY : k_grid,q_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npl,nbl
    REAL(DP),INTENT(IN) :: d_diago(n_lanczos,npl,nbl,k_grid%nps,q_grid%nps)
    COMPLEX(DP),INTENT(IN) :: z_body2(n_lanczos,npl,ifr%nloc,nbl,k_grid%nps,q_grid%nps)
    !
    CALL start_clock('write_g')
    !
    CALL stop_clock('write_g')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_gfreq_real(d_diago,d_body2,npl,nbl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : n_lanczos
    USE distribution_center,  ONLY : ifr
    USE types_bz_grid,        ONLY : k_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npl,nbl
    REAL(DP),INTENT(OUT) :: d_diago(n_lanczos,npl,nbl,k_grid%nps)
    REAL(DP),INTENT(OUT) :: d_body2(n_lanczos,npl,ifr%nloc,nbl,k_grid%nps)
    !
    CALL start_clock('read_g')
    !
    CALL stop_clock('read_g')
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_gfreq_complex_q(d_diago,z_body2,npl,nbl)
    !------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE westcom,              ONLY : n_lanczos
    USE distribution_center,  ONLY : ifr
    USE types_bz_grid,        ONLY : k_grid,q_grid
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER,INTENT(IN) :: npl,nbl
    REAL(DP),INTENT(OUT) :: d_diago(n_lanczos,npl,nbl,k_grid%nps,q_grid%nps)
    COMPLEX(DP),INTENT(OUT) :: z_body2(n_lanczos,npl,ifr%nloc,nbl,k_grid%nps,q_grid%nps)
    !
    CALL start_clock('read_g')
    !
    CALL stop_clock('read_g')
    !
  END SUBROUTINE
  !
END MODULE
