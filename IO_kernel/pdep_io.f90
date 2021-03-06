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
  USE kinds,         ONLY : DP
  USE mp_global,     ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm,my_pool_id,&
                          & my_bgrp_id,inter_bgrp_comm,inter_pool_comm,intra_pool_comm
  USE westcom,       ONLY : npwq,npwq_g,npwqx,ngq,ngq_g,igq_q
  USE gvect,         ONLY : ig_l2g
  USE json_module,   ONLY : json_file
  USE control_flags, ONLY : gamma_only
  USE base64_module
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
      USE mp,           ONLY : mp_bcast,mp_max
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
      CHARACTER(LEN=:),ALLOCATABLE :: charbase64
      INTEGER :: nbytes, ndim, iunit, nlen
      CHARACTER(LEN=30) :: endian
      INTEGER :: npwqx_g
      INTEGER, ALLOCATABLE :: igq_l2g_kdip(:), igq_l2g(:)
      INTEGER, PARAMETER :: default_iq = 1
      INTEGER :: iq_
      !
      CALL start_clock('pdep_write')
      !
      IF( PRESENT(iq) ) THEN
         iq_ = iq
      ELSE
         iq_ = default_iq
      ENDIF
      !
      IF ( .NOT. gamma_only) THEN
         !
         ! Resume all components
         !
         ndim = ngq_g(iq_)
         !
         ! <NEW>
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
         ! </NEW>
         !
         ! npwq_g = MAXVAL(igq_l2g_kdip(1:ndim,iq))
         ! CALL mp_max(npwq_g,intra_pool_comm)
         ! CALL mp_max(npwq_g,intra_bgrp_comm)
         !
         ALLOCATE( tmp_vec(npwq_g) )
         tmp_vec=0._DP
         !
         CALL mergewf( pdepg(:), tmp_vec, npwq, igq_l2g_kdip, me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
         DEALLOCATE( igq_l2g_kdip )
         !
         ! ONLY ROOT W/IN BGRP WRITES
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            nbytes = SIZEOF(tmp_vec(1)) * ndim
            nlen = lenbase64(nbytes)
            ALLOCATE(CHARACTER(LEN=nlen) :: charbase64)
            CALL base64_encode_complex(tmp_vec(1:ndim), ndim, charbase64)
            !
            IF( islittleendian() ) THEN
               endian = '"islittleendian" : true'
            ELSE
               endian = '"islittleendian" : false'
            ENDIF
            !
            OPEN( NEWUNIT=iunit, FILE = TRIM(fname) )
            WRITE( iunit, '(a)' ) '{'
            WRITE( iunit, '(a,i0,a)' ) '"meta" : { "name" : "eigenpotential", "type" : "complex double", "space" : "G",&
                         "ndim" : ', ndim, ', "encoding" : "base64", '//TRIM(endian)//' }'
            WRITE( iunit, '(a)') ', "data" : '
            WRITE( iunit, '(a)' ) '"'//charbase64//'"'
            WRITE( iunit, '(a)' ) '}'
            CLOSE( iunit )
            !
            DEALLOCATE( charbase64 )
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
         tmp_vec=0._DP
         !
         CALL mergewf( pdepg(:), tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
         !
         ! ONLY ROOT W/IN BGRP WRITES
         !
         IF(me_bgrp==root_bgrp) THEN
            !
            ndim = npwq_g
            nbytes = SIZEOF(tmp_vec(1)) * ndim
            nlen = lenbase64(nbytes)
            ALLOCATE(CHARACTER(LEN=nlen) :: charbase64)
            CALL base64_encode_complex(tmp_vec(1:ndim), ndim, charbase64)
            !
            IF( islittleendian() ) THEN
               endian = '"islittleendian" : true'
            ELSE
               endian = '"islittleendian" : false'
            ENDIF
            !
            OPEN( NEWUNIT=iunit, FILE = TRIM(fname) )
            WRITE( iunit, '(a)' ) '{'
            WRITE( iunit, '(a,i0,a)' ) '"meta" : { "name" : "eigenpotential", "type" : "complex double", "space" : "G",&
                         "ndim" : ', ndim, ', "encoding" : "base64", '//TRIM(endian)//' }'
            WRITE( iunit, '(a)') ', "data" : '
            WRITE( iunit, '(a)' ) '"'//charbase64//'"'
            WRITE( iunit, '(a)' ) '}'
            CLOSE( iunit )
            !
            DEALLOCATE( charbase64 )
            !
         ENDIF
         !
         DEALLOCATE( tmp_vec )
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
      USE mp,           ONLY : mp_bcast, mp_max
      USE mp_global,    ONLY : intra_bgrp_comm
      USE base64_module
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
      TYPE(json_file) :: json
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: ig
      CHARACTER(LEN=1000) :: line
      CHARACTER(LEN=:),ALLOCATABLE :: charbase64
      CHARACTER(LEN=:),ALLOCATABLE :: endian
      INTEGER :: nbytes, ndim, iunit, nlen
      LOGICAL :: found, isle
      INTEGER :: npwqx_g
      INTEGER, ALLOCATABLE :: igq_l2g_kdip(:), igq_l2g(:)
      INTEGER, PARAMETER :: default_iq = 1
      INTEGER :: iq_
      !
      CALL start_clock('pdep_read')
      !
      IF( PRESENT(iq) ) THEN
         iq_ = iq
      ELSE
         iq_ = default_iq
      ENDIF
      !
      IF ( .NOT. gamma_only ) THEN
         !
         ! Resume all components
         !
         ndim = ngq_g(iq_)
         ! npwq_g = MAXVAL(igq_l2g_kdip(1:ndim,iq))
         ! CALL mp_max(npwq_g,intra_pool_comm)
         ! CALL mp_max(npwq_g,intra_bgrp_comm)
         !
         ALLOCATE( tmp_vec(npwq_g) )
         tmp_vec=0._DP
         pdepg=0._DP
         !
         IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            nbytes = SIZEOF(tmp_vec(1)) * ndim
            nlen = lenbase64(nbytes)
            !
            IF(me_bgrp==root_bgrp) THEN
               !
               ALLOCATE(CHARACTER(LEN=(nlen+2)) :: charbase64)
               !
               OPEN( NEWUNIT=iunit, FILE = TRIM(fname) )
               READ( iunit, * )
               READ( iunit, '(a)' ) line
               READ( iunit, * )
               READ( iunit, '(a)' ) charbase64
               CLOSE( iunit )
               CALL base64_decode_complex(charbase64(2:(nlen+1)), ndim, tmp_vec(1:ndim))
               DEALLOCATE( charbase64 )
               !
               CALL json%load_from_string("{"//TRIM(line)//"}")
               CALL json%get('meta.islittleendian', isle, found)
               CALL json%destroy()
               !
               IF (islittleendian() .NEQV. isle) CALL base64_byteswap_complex(nbytes,tmp_vec(1:ndim))
               !
            END IF
            !
            !
            ! <NEW>
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
            ! </NEW>
            !
            CALL splitwf( pdepg, tmp_vec, npwq, igq_l2g_kdip, me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
            DEALLOCATE( igq_l2g_kdip )
            !
         ENDIF
         !
         DEALLOCATE( tmp_vec )
         !
         CALL mp_bcast(pdepg,0,inter_bgrp_comm)
         CALL mp_bcast(pdepg,0,inter_pool_comm)
      !
      ELSE
         !
         ! Resume all components
         !
         ALLOCATE( tmp_vec(npwq_g) )
         tmp_vec=0._DP
         pdepg=0._DP
         !
         IF(my_pool_id==0.AND.my_bgrp_id==0) THEN
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            ndim = npwq_g
            nbytes = SIZEOF(tmp_vec(1)) * ndim
            nlen = lenbase64(nbytes)
            !
            IF(me_bgrp==root_bgrp) THEN
               !
               ALLOCATE(CHARACTER(LEN=(nlen+2)) :: charbase64)
               !
               OPEN( NEWUNIT=iunit, FILE = TRIM(fname) )
               READ( iunit, * )
               READ( iunit, '(a)' ) line
               READ( iunit, * )
               READ( iunit, '(a)' ) charbase64
               CLOSE( iunit )
               CALL base64_decode_complex(charbase64(2:(nlen+1)), ndim, tmp_vec(1:ndim))
               DEALLOCATE( charbase64 )
               !
               CALL json%load_from_string("{"//TRIM(line)//"}")
               CALL json%get('meta.islittleendian', isle, found)
               CALL json%destroy()
               !
               IF (islittleendian() .NEQV. isle) CALL base64_byteswap_complex(nbytes,tmp_vec(1:ndim))
               !
            END IF
            !
            CALL splitwf( pdepg, tmp_vec, npwq, ig_l2g(1:npwq), me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
            !
         ENDIF
         !
         DEALLOCATE( tmp_vec )
         !
         CALL mp_bcast(pdepg,0,inter_bgrp_comm)
         CALL mp_bcast(pdepg,0,inter_pool_comm)
      !
      ENDIF
      !
      CALL stop_clock('pdep_read')
      !
    END SUBROUTINE
    !
END MODULE
