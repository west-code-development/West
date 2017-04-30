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
MODULE pdep_io
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,       ONLY : DP
  USE mp_global,   ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm,my_pool_id,my_bgrp_id,inter_bgrp_comm,inter_pool_comm
  USE westcom,     ONLY : npwq0, npwq0_g, npwq0x
  USE gvect,       ONLY : ig_l2g
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    ! ******************************************
    ! WRITE IN G SPACE 
    !       wfc is passed distributed in G space
    !       then merged and written in R space
    ! ******************************************
    !
    SUBROUTINE pdep_merge_and_write_G(fname,pdepg)
      !
      USE mp_wave,      ONLY : mergewf
      USE mp,           ONLY : mp_bcast
      !
      ! I/O
      !    
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP), INTENT(IN) :: pdepg(npwq0x)
      !
      ! Scratch
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr
      !
      !
      IF(my_pool_id.NE.0) RETURN
      IF(my_bgrp_id.NE.0) RETURN
      !
      ! Resume all components 
      !
      ALLOCATE( tmp_vec(npwq0_g) )
      tmp_vec=0._DP
      !
      CALL mergewf( pdepg(:), tmp_vec, npwq0, ig_l2g(1:npwq0), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm) 
      !
      ! ONLY ROOT W/IN BGRP WRITES
      !
      IF(me_bgrp==root_bgrp) THEN 
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iun, ierr )
         CALL iotk_open_write( iun, FILE = TRIM(fname), BINARY = .TRUE.)
         CALL iotk_write_begin( iun, 'PDEP_GSPACE' )
         CALL iotk_write_dat( iun, "ndim" , npwq0_g )
         CALL iotk_write_dat( iun, "pdep" , tmp_vec(1:npwq0_g) )
         CALL iotk_write_end( iun, 'PDEP_GSPACE' )
         !
         CALL iotk_close_write( iun )
         !
      END IF
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
    SUBROUTINE pdep_read_G_and_distribute(fname,pdepg)
      !
      USE mp_wave,      ONLY : splitwf
      USE mp,           ONLY : mp_bcast
      USE mp_global,    ONLY : intra_bgrp_comm
      !
      ! I/O
      !    
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP), INTENT(OUT) :: pdepg(npwq0x)
      !
      ! Scratch
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ig
      !
      ! Resume all components 
      !
      ALLOCATE( tmp_vec(npwq0_g) )
      tmp_vec=0._DP
      pdepg=0._DP
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
            CALL iotk_scan_begin( iun, 'PDEP_GSPACE' )
            CALL iotk_scan_dat( iun, "pdep" , tmp_vec(1:npwq0_g) )
            CALL iotk_scan_end( iun, 'PDEP_GSPACE' )
            CALL iotk_close_read( iun )
            !
         END IF
         !
         CALL splitwf( pdepg, tmp_vec, npwq0, ig_l2g(1:npwq0), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm) 
         !
      ENDIF
      !
      DEALLOCATE( tmp_vec )
      !
      CALL mp_bcast(pdepg,0,inter_bgrp_comm)
      CALL mp_bcast(pdepg,0,inter_pool_comm)
      !
    END SUBROUTINE
    !
END MODULE
