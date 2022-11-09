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
MODULE pdep_io
  !----------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP,i8b
  USE mp_global,     ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm
  USE westcom,       ONLY : npwq,npwq_g,npwqx,ngq,ngq_g,igq_q
  USE gvect,         ONLY : ig_l2g
  USE control_flags, ONLY : gamma_only
  USE base64_module, ONLY : islittleendian
  USE west_io,       ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,HD_ID_DIMENSION
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
    SUBROUTINE pdep_merge_and_write_G(fname,pdepg,iq)
      !
      USE mp_wave,      ONLY : mergewf
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(LEN=*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(IN) :: pdepg(npwqx)
      INTEGER,INTENT(IN),OPTIONAL :: iq
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: ig
      INTEGER :: ndim
      INTEGER :: iun
      INTEGER :: npwqx_g
      INTEGER,ALLOCATABLE :: igq_l2g_kdip(:)
      INTEGER,ALLOCATABLE :: igq_l2g(:)
      INTEGER,PARAMETER :: default_iq = 1
      INTEGER :: iq_
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock('pdep_write')
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
         npwqx_g = MAXVAL(ngq_g)
         ALLOCATE(igq_l2g_kdip(npwqx_g))
         igq_l2g_kdip = 0
         !
         ALLOCATE(igq_l2g(ngq(iq_)))
         DO ig = 1,ngq(iq_)
            igq_l2g(ig) = ig_l2g(igq_q(ig,iq_))
         ENDDO
         CALL gq_l2gmap_kdip(npwq_g,ngq_g(iq_),ngq(iq_),igq_l2g,igq_l2g_kdip)
         DEALLOCATE(igq_l2g)
         !
         ALLOCATE(tmp_vec(npwq_g))
         tmp_vec = 0._DP
         !
         CALL mergewf(pdepg,tmp_vec,npwq,igq_l2g_kdip,me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
         DEALLOCATE(igq_l2g_kdip)
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
            OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
            offset = 1
            WRITE(iun,POS=offset) header
            offset = 1+HD_LENGTH*SIZEOF(header(1))
            WRITE(iun,POS=offset) tmp_vec(1:ndim)
            CLOSE(iun)
            !
         ENDIF
         !
         DEALLOCATE(tmp_vec)
         !
      ELSE
         !
         ! Resume all components
         !
         ALLOCATE(tmp_vec(npwq_g))
         tmp_vec = 0._DP
         !
         CALL mergewf(pdepg,tmp_vec,npwq,ig_l2g(1:npwq),me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
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
            OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
            offset = 1
            WRITE(iun,POS=offset) header
            offset = 1+HD_LENGTH*SIZEOF(header(1))
            WRITE(iun,POS=offset) tmp_vec(1:ndim)
            CLOSE(iun)
            !
         ENDIF
         !
         DEALLOCATE(tmp_vec)
         !
      ENDIF
      !
      CALL stop_clock('pdep_write')
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
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(LEN=*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(OUT) :: pdepg(npwqx)
      INTEGER,INTENT(IN),OPTIONAL :: iq
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: ig
      INTEGER :: ndim
      INTEGER :: iun
      INTEGER :: npwqx_g
      INTEGER,ALLOCATABLE :: igq_l2g_kdip(:)
      INTEGER,ALLOCATABLE :: igq_l2g(:)
      INTEGER,PARAMETER :: default_iq = 1
      INTEGER :: iq_
      INTEGER :: ierr
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock('pdep_read')
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
         ALLOCATE(tmp_vec(npwq_g))
         tmp_vec = 0._DP
         pdepg = 0._DP
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         IF(me_bgrp == root_bgrp) THEN
            !
            OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
            IF(ierr /= 0) THEN
               CALL errore('pdep_read','Cannot read file:'//TRIM(fname),1)
            ENDIF
            !
            offset = 1
            READ(iun,POS=offset) header
            IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
               CALL errore('pdep_read','Unknown file format:'//TRIM(fname),1)
            ENDIF
            IF(ndim /= header(HD_ID_DIMENSION)) THEN
               CALL errore('pdep_read','Dimension mismatch:'//TRIM(fname),1)
            ENDIF
            IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
               .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
               CALL errore('pdep_read','Endianness mismatch:'//TRIM(fname),1)
            ENDIF
            !
            offset = 1+HD_LENGTH*SIZEOF(header(1))
            READ(iun,POS=offset) tmp_vec(1:ndim)
            CLOSE(iun)
            !
         ENDIF
         !
         npwqx_g = MAXVAL(ngq_g)
         ALLOCATE(igq_l2g_kdip(npwqx_g))
         igq_l2g_kdip = 0
         !
         ALLOCATE(igq_l2g(ngq(iq_)))
         DO ig = 1,ngq(iq_)
            igq_l2g(ig) = ig_l2g(igq_q(ig,iq_))
         ENDDO
         CALL gq_l2gmap_kdip(npwq_g,ngq_g(iq_),ngq(iq_),igq_l2g,igq_l2g_kdip)
         DEALLOCATE(igq_l2g)
         !
         CALL splitwf(pdepg,tmp_vec,npwq,igq_l2g_kdip,me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
         DEALLOCATE(igq_l2g_kdip)
         !
         DEALLOCATE(tmp_vec)
         !
      ELSE
         !
         ! Resume all components
         !
         ALLOCATE(tmp_vec(npwq_g))
         tmp_vec = 0._DP
         pdepg = 0._DP
         !
         ! ONLY ROOT W/IN BGRP READS
         !
         ndim = npwq_g
         !
         IF(me_bgrp == root_bgrp) THEN
            !
            OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
            IF(ierr /= 0) THEN
               CALL errore('pdep_read','Cannot read file:'//TRIM(fname),1)
            ENDIF
            !
            offset = 1
            READ(iun,POS=offset) header
            IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
               CALL errore('pdep_read','Unknown file format:'//TRIM(fname),1)
            ENDIF
            IF(ndim /= header(HD_ID_DIMENSION)) THEN
               CALL errore('pdep_read','Dimension mismatch:'//TRIM(fname),1)
            ENDIF
            IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
               .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
               CALL errore('pdep_read','Endianness mismatch:'//TRIM(fname),1)
            ENDIF
            !
            offset = 1+HD_LENGTH*SIZEOF(header(1))
            READ(iun,POS=offset) tmp_vec(1:ndim)
            CLOSE(iun)
            !
         ENDIF
         !
         CALL splitwf(pdepg,tmp_vec,npwq,ig_l2g(1:npwq),me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
         !
         DEALLOCATE(tmp_vec)
         !
      ENDIF
      !
      CALL stop_clock('pdep_read')
      !
    END SUBROUTINE
    !
END MODULE
