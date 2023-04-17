!
! Copyright (C) 2015-2023 M. Govoni
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
    SUBROUTINE plep_merge_and_write_G(fname,plepg)
      !
      USE kinds,         ONLY : DP,i8b
      USE mp_global,     ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm
      USE westcom,       ONLY : nbndval0x,n_trunc_bands
      USE gvect,         ONLY : ig_l2g
      USE pwcom,         ONLY : nks,npwx
      USE base64_module, ONLY : islittleendian
      USE west_io,       ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,HD_ID_DIMENSION
      USE mp_wave,       ONLY : mergewf
      USE mp,            ONLY : mp_max
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(LEN=*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(IN) :: plepg(npwx,nbndval0x-n_trunc_bands,nks)
      !
      ! Workspae
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ibnd,ik,npwx_g
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock('plep_write')
      !
      npwx_g = MAXVAL(ig_l2g(1:npwx))
      CALL mp_max(npwx_g,intra_bgrp_comm)
      !
      IF(me_bgrp == root_bgrp) THEN
         !
         header = 0
         header(HD_ID_VERSION) = HD_VERSION
         header(HD_ID_DIMENSION) = npwx_g
         IF(islittleendian()) THEN
            header(HD_ID_LITTLE_ENDIAN) = 1
         ENDIF
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
         offset = 1
         WRITE(iun,POS=offset) header
         offset = offset+SIZEOF(header)
         !
      ENDIF
      !
      ! Resume all components
      !
      ALLOCATE(tmp_vec(npwx_g))
      !
      DO ik = 1,nks
         DO ibnd = 1,nbndval0x-n_trunc_bands
            !
            tmp_vec = 0._DP
            !
            CALL mergewf(plepg(:,ibnd,ik),tmp_vec,npwx,ig_l2g(1:npwx),me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
            !
            ! ONLY ROOT W/IN BGRP WRITES
            !
            IF(me_bgrp == root_bgrp) THEN
               WRITE(iun,POS=offset) tmp_vec(1:npwx_g)
               offset = offset+SIZEOF(tmp_vec)
            ENDIF
            !
         ENDDO
      ENDDO
      !
      IF(me_bgrp == root_bgrp) THEN
         CLOSE(iun)
      ENDIF
      !
      DEALLOCATE(tmp_vec)
      !
      CALL stop_clock('plep_write')
      !
    END SUBROUTINE
    !
    ! ******************************************
    ! READ IN G SPACE
    !       wfc is read merged in G space
    !       then split in G space
    ! ******************************************
    !
    SUBROUTINE plep_read_G_and_distribute(fname,plepg)
      !
      USE kinds,         ONLY : DP,i8b
      USE mp_global,     ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm
      USE westcom,       ONLY : nbndval0x,n_trunc_bands
      USE gvect,         ONLY : ig_l2g
      USE pwcom,         ONLY : nks,npwx
      USE base64_module, ONLY : islittleendian
      USE west_io,       ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,HD_ID_DIMENSION
      USE mp_wave,       ONLY : splitwf
      USE mp,            ONLY : mp_max
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(LEN=*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(OUT) :: plepg(npwx,nbndval0x-n_trunc_bands,nks)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ibnd,ik,npwx_g
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock('plep_read')
      !
      npwx_g = MAXVAL(ig_l2g(1:npwx))
      CALL mp_max(npwx_g,intra_bgrp_comm)
      !
      ! Resume all components
      !
      ALLOCATE(tmp_vec(npwx_g))
      tmp_vec = 0._DP
      plepg = 0._DP
      !
      IF(me_bgrp == root_bgrp) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
         IF(ierr /= 0) THEN
            CALL errore('plep_read','Cannot read file: '//TRIM(fname),1)
         ENDIF
         !
         offset = 1
         READ(iun,POS=offset) header
         IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
            CALL errore('plep_read','Unknown file format: '//TRIM(fname),1)
         ENDIF
         IF(npwx_g /= header(HD_ID_DIMENSION)) THEN
            CALL errore('plep_read','Dimension mismatch: '//TRIM(fname),1)
         ENDIF
         IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
            .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
            CALL errore('plep_read','Endianness mismatch: '//TRIM(fname),1)
         ENDIF
         !
         offset = offset+SIZEOF(header)
         !
      ENDIF
      !
      DO ik = 1,nks
         DO ibnd = 1,nbndval0x-n_trunc_bands
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            IF(me_bgrp == root_bgrp) THEN
               READ(iun,POS=offset) tmp_vec(1:npwx_g)
               offset = offset+SIZEOF(tmp_vec)
            ENDIF
            !
            CALL splitwf(plepg(:,ibnd,ik),tmp_vec,npwx,ig_l2g(1:npwx),me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
            !
         ENDDO
      ENDDO
      !
      IF(me_bgrp == root_bgrp) THEN
         CLOSE(iun)
      ENDIF
      !
      DEALLOCATE(tmp_vec)
      !
      CALL stop_clock('plep_read')
      !
    END SUBROUTINE
    !
END MODULE
