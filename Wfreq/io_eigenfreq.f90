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
MODULE io_eigenfreq 
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir
  !
  IMPLICIT NONE
  !
  !
  CONTAINS
    !
    !
    ! *****************************
    ! EIGENFREQ WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE eigenfreq_write(e, ne, ind)
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast
      USE mp_global,            ONLY : my_image_id
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : wfreq_dirname 
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN) :: e(ne)
      INTEGER,INTENT(IN) :: ind
      INTEGER,INTENT(IN) :: ne
      !
      CHARACTER(LEN=256)    :: fname
      CHARACTER(LEN=6)      :: my_label
      INTEGER :: iunout
      INTEGER :: ierr
      !
      WRITE(my_label,'(i6.6)') ind
      fname = TRIM( wfreq_dirname ) // "/EVF"//TRIM(ADJUSTL(my_label))//".dat"
      !
      IF ( my_image_id == 0 ) THEN
         !
         ! ... open XML descriptor
         !
         !CALL iotk_free_unit( iunout, ierr )
         iunout = 2000 + ind
         CALL iotk_open_write( iunout, FILE = TRIM( fname ), BINARY=.FALSE., IERR=ierr)
         !
         ! ... dump
         !
         CALL iotk_write_begin( iunout, "EIGENVALUES" )
         CALL iotk_write_dat( iunout, "ndim", ne )
         CALL iotk_write_dat( iunout, "ev", e(1:ne ))
         CALL iotk_write_end( iunout, "EIGENVALUES" )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
    END SUBROUTINE
    !
    !
    ! *****************************
    ! EIGENFREQ READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE eigenfreq_read( e, ne, ind )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast
      USE mp_global,            ONLY : my_image_id,inter_image_comm
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : wfreq_dirname 
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(OUT) :: e(ne)
      INTEGER,INTENT(IN) :: ind
      INTEGER,INTENT(IN) :: ne
      !
      CHARACTER(LEN=256)    :: fname
      CHARACTER(LEN=6)      :: my_label
      INTEGER :: iunout
      INTEGER :: ierr
      !
      WRITE(my_label,'(i6.6)') ind
      fname = TRIM( wfreq_dirname ) // "/EVF"//TRIM(ADJUSTL(my_label))//".dat"
      !
      IF ( my_image_id == 0 ) THEN
         !
         ! ... open XML descriptor
         !
         iunout = 2000 + ind
         CALL iotk_open_read( iunout, FILE = TRIM( fname ), BINARY=.FALSE., IERR=ierr)
         !
         ! ... pick 
         !
         CALL iotk_scan_begin( iunout, "EIGENVALUES" )
         CALL iotk_scan_dat( iunout, "ev", e(1:ne ))
         CALL iotk_scan_end( iunout, "EIGENVALUES" )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_read( iunout )
         !
      END IF
      !
      CALL mp_bcast( e, 0, inter_image_comm ) 
      !
    END SUBROUTINE
    !
END MODULE
