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
MODULE pdep_io
  !----------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP,i8b
  USE mp_global,     ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm,my_pool_id,&
                          & my_bgrp_id,inter_bgrp_comm,inter_pool_comm
  USE westcom,       ONLY : npwq,npwq_g,npwqx,ngq,ngq_g,igq_q
  USE gvect,         ONLY : ig_l2g
  USE control_flags, ONLY : gamma_only
  USE base64_module, ONLY : islittleendian
  !
  IMPLICIT NONE
  !
  ! Base64 was changed to binary in order to improve I/O performance.
  !
  ! A simple format is used here:
  ! a header consisting of 32 (HD_LENGTH) integers, followed by raw data.
  ! Currently only 3 integers are used in the header, storing the version
  ! of the binary format, which could change in the future, the endianness,
  ! and the length of the raw data which is useful when reading the file.
  !
  INTEGER, PARAMETER :: HD_LENGTH = 32
  INTEGER, PARAMETER :: HD_VERSION = 210405
  INTEGER, PARAMETER :: HD_ID_VERSION = 1
  INTEGER, PARAMETER :: HD_ID_LITTLE_ENDIAN = 2
  INTEGER, PARAMETER :: HD_ID_DIMENSION = 3
  !
  CONTAINS
    !
    ! ******************************************
    ! WRITE IN G SPACE
    !       wfc is passed distributed in G space
    !       then merged and written in R space
    ! ******************************************
    !
    SUBROUTINE pdep_merge_and_write_G(fname,pdepg,iq)
      !
      USE mp_wave,      ONLY : mergewf
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP), INTENT(IN) :: pdepg(npwqx)
      INTEGER, INTENT(IN), OPTIONAL :: iq
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: ig
      INTEGER :: ndim, iunit
      INTEGER :: npwqx_g
      INTEGER, ALLOCATABLE :: igq_l2g_kdip(:), igq_l2g(:)
      INTEGER, PARAMETER :: default_iq = 1
      INTEGER :: iq_
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock( 'pdep_write' )
      !
      IF(PRESENT(iq)) THEN
         iq_ = iq
      ELSE
         iq_ = default_iq
      ENDIF
      !
      IF(.NOT. gamma_only) THEN
         !
         ! Resume all components
         !
         ndim = ngq_g(iq_)
         !
         npwqx_g = MAXVAL( ngq_g(:) )
         ALLOCATE( igq_l2g_kdip(npwqx_g) )
         igq_l2g_kdip(:) = 0
         !
         ALLOCATE( igq_l2g(ngq(iq_)) )
         DO ig = 1, ngq(iq_)
            igq_l2g(ig) = ig_l2g( igq_q(ig,iq_) )
         ENDDO
         CALL gq_l2gmap_kdip( npwq_g, ngq_g(iq_), ngq(iq_), igq_l2g, igq_l2g_kdip )
         DEALLOCATE( igq_l2g )
         !
         ALLOCATE( tmp_vec(npwq_g) )
         tmp_vec = 0._DP
         !
         CALL mergewf( pdepg(:), tmp_vec, npwq, igq_l2g_kdip, me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm )
         DEALLOCATE( igq_l2g_kdip )
         !
         ! ONLY ROOT W/IN BGRP WRITES
         !
         IF(me_bgrp == root_bgrp) THEN
            !
            header = 0
            header(HD_ID_VERSION) = HD_VERSION
            header(HD_ID_DIMENSION) = ndim
            IF(islittleendian()) THEN
               header(HD_ID_LITTLE_ENDIAN) = 1
            ENDIF
            !
            OPEN( NEWUNIT=iunit, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED' )
            offset = 1
            WRITE( iunit, POS=offset ) header
            offset = 1+HD_LENGTH*4
            WRITE( iunit, POS=offset ) tmp_vec(1:ndim)
            CLOSE( iunit )
            !
         END IF
         !
         DEALLOCATE( tmp_vec )
      !
      ELSE
         !
         ! Resume all components
         !
         ALLOCATE( tmp_vec(npwq_g) )
         tmp_vec = 0._DP
         !
         CALL mergewf( pdepg(:), tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm )
         !
         ! ONLY ROOT W/IN BGRP WRITES
         !
         IF(me_bgrp == root_bgrp) THEN
            !
            ndim = npwq_g
            header = 0
            header(HD_ID_VERSION) = HD_VERSION
            header(HD_ID_DIMENSION) = ndim
            IF(islittleendian()) THEN
               header(HD_ID_LITTLE_ENDIAN) = 1
            ENDIF
            !
            OPEN( NEWUNIT=iunit, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED' )
            offset = 1
            WRITE( iunit, POS=offset ) header
            offset = 1+HD_LENGTH*4
            WRITE( iunit, POS=offset ) tmp_vec(1:ndim)
            CLOSE( iunit )
            !
         ENDIF
         !
         DEALLOCATE( tmp_vec )
         !
      ENDIF
      !
      CALL stop_clock( 'pdep_write' )
      !
    END SUBROUTINE
    !
    ! ******************************************
    ! READ IN G SPACE
    !       wfc is read merged in G space
    !       then split in G space
    ! ******************************************
    !
    SUBROUTINE pdep_read_G_and_distribute(fname,pdepg,iq)
      !
      USE mp_wave,      ONLY : splitwf
      USE mp,           ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*), INTENT(IN) :: fname
      COMPLEX(DP), INTENT(OUT) :: pdepg(npwqx)
      INTEGER, INTENT(IN), OPTIONAL :: iq
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: ig
      INTEGER :: ndim, iunit
      INTEGER :: npwqx_g
      INTEGER, ALLOCATABLE :: igq_l2g_kdip(:), igq_l2g(:)
      INTEGER, PARAMETER :: default_iq = 1
      INTEGER :: iq_
      INTEGER :: ierr
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock( 'pdep_read' )
      !
      IF(PRESENT(iq)) THEN
         iq_ = iq
      ELSE
         iq_ = default_iq
      ENDIF
      !
      IF(.NOT. gamma_only) THEN
         !
         ! Resume all components
         !
         ndim = ngq_g(iq_)
         !
         ALLOCATE( tmp_vec(npwq_g) )
         tmp_vec = 0._DP
         pdepg = 0._DP
         !
         IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            IF(me_bgrp == root_bgrp) THEN
               !
               OPEN( NEWUNIT=iunit, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ierr )
               IF(ierr /= 0) THEN
                  CALL errore('pdep_read', 'Cannot RD F:'//TRIM(ADJUSTL(fname)),1)
               ENDIF
               !
               offset = 1
               READ( iunit, POS=offset ) header
               IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
                  CALL errore('pdep_read', 'Unknown file format:'//TRIM(ADJUSTL(fname)),1)
               ENDIF
               IF(ndim /= header(HD_ID_DIMENSION)) THEN
                  CALL errore('pdep_read', 'Dimension mismatch:'//TRIM(ADJUSTL(fname)),1)
               ENDIF
               IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
                  .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
                  CALL errore('pdep_read', 'Endianness mismatch:'//TRIM(ADJUSTL(fname)),1)
               ENDIF
               !
               offset = 1+HD_LENGTH*4
               READ( iunit, POS=offset ) tmp_vec(1:ndim)
               CLOSE( iunit )
               !
            END IF
            !
            npwqx_g = MAXVAL( ngq_g(:) )
            ALLOCATE( igq_l2g_kdip(npwqx_g) )
            igq_l2g_kdip(:) = 0
            !
            ALLOCATE( igq_l2g(ngq(iq_)) )
            DO ig = 1, ngq(iq_)
               igq_l2g(ig) = ig_l2g( igq_q(ig,iq_) )
            ENDDO
            CALL gq_l2gmap_kdip( npwq_g, ngq_g(iq_), ngq(iq_), igq_l2g, igq_l2g_kdip )
            DEALLOCATE( igq_l2g )
            !
            CALL splitwf( pdepg, tmp_vec, npwq, igq_l2g_kdip, me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm )
            DEALLOCATE( igq_l2g_kdip )
            !
         ENDIF
         !
         DEALLOCATE( tmp_vec )
         !
         CALL mp_bcast( pdepg, 0, inter_bgrp_comm )
         CALL mp_bcast( pdepg, 0, inter_pool_comm )
         !
      ELSE
         !
         ! Resume all components
         !
         ALLOCATE( tmp_vec(npwq_g) )
         tmp_vec = 0._DP
         pdepg = 0._DP
         !
         IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            ndim = npwq_g
            !
            IF(me_bgrp == root_bgrp) THEN
               !
               OPEN( NEWUNIT=iunit, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ierr )
               IF(ierr /= 0) THEN
                  CALL errore('pdep_read', 'Cannot RD F:'//TRIM(ADJUSTL(fname)), 1)
               ENDIF
               !
               offset = 1
               READ( iunit, POS=offset ) header
               IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
                  CALL errore('pdep_read', 'Unknown file format:'//TRIM(ADJUSTL(fname)),1)
               ENDIF
               IF(ndim /= header(HD_ID_DIMENSION)) THEN
                  CALL errore('pdep_read', 'Dimension mismatch:'//TRIM(ADJUSTL(fname)),1)
               ENDIF
               IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
                  .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
                  CALL errore('pdep_read', 'Endianness mismatch:'//TRIM(ADJUSTL(fname)),1)
               ENDIF
               !
               offset = 1+HD_LENGTH*4
               READ( iunit, POS=offset ) tmp_vec(1:ndim)
               CLOSE( iunit )
               !
            END IF
            !
            CALL splitwf( pdepg, tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm )
            !
         ENDIF
         !
         DEALLOCATE( tmp_vec )
         !
         CALL mp_bcast( pdepg, 0, inter_bgrp_comm )
         CALL mp_bcast( pdepg, 0, inter_pool_comm )
         !
      ENDIF
      !
      CALL stop_clock( 'pdep_read' )
      !
    END SUBROUTINE
    !
END MODULE
