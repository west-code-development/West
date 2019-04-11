!
! Copyright (C) 2015-2017 M. Govoni 
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
! -------------------------------------------------------------------
MODULE function3d 
 ! -----------------------------------------------------------------
 !
 IMPLICIT NONE
 !
 CONTAINS
 ! 
 !-----------------------------------------------------------------
   SUBROUTINE write_function3d ( dfft, fname, descriptor, ng, ngx, funct3d_g, nmaps, nl )
   ! -----------------------------------------------------------------
   !
   USE kinds,                       ONLY : DP
   USE cell_base,                   ONLY : celldm, at
   USE control_flags,               ONLY : gamma_only
   USE mp_bands,                    ONLY : me_bgrp
   USE westcom,                     ONLY : fftdriver
   USE scatter_mod,                 ONLY : gather_grid
   USE fft_types,                   ONLY : fft_type_descriptor
   USE base64_module
   USE fourier_interpolation
   !
   IMPLICIT NONE
   !
   ! I/O 
   !
   TYPE(fft_type_descriptor), INTENT(IN) :: dfft
   CHARACTER(LEN=*),INTENT(IN) :: fname
   CHARACTER(LEN=*),INTENT(IN) :: descriptor
   INTEGER, INTENT(IN) :: ng, ngx, nmaps
   COMPLEX(DP),INTENT(IN) :: funct3d_g(ngx)
   INTEGER,INTENT(IN) :: nl(nmaps,ngx)
   ! 
   ! Workspace
   !
   INTEGER              :: iu, ndim, nbytes, nlen, i
   COMPLEX(DP),ALLOCATABLE :: funct3d_r_complex(:)
   COMPLEX(DP),ALLOCATABLE :: funct3d_r_complex_gathered(:)
   REAL(DP),ALLOCATABLE :: funct3d_r_double(:)
   CHARACTER(LEN=14)    :: lab(3)
   CHARACTER(LEN=:),ALLOCATABLE :: charbase64
   CHARACTER(LEN=:),ALLOCATABLE :: ctype
   !
   ALLOCATE( funct3d_r_complex(dfft%nnr) )
   IF ( gamma_only ) THEN 
      SELECT CASE(fftdriver) 
      CASE('Wave') 
         CALL single_interp_invfft_gamma(dfft,ng,ngx,funct3d_g,funct3d_r_complex,'Wave',nl)
      CASE('Dense')
         CALL single_interp_invfft_k(dfft,ng,ngx,funct3d_g,funct3d_r_complex,'Dense',nl)
      END SELECT
   ELSE
      CALL single_interp_invfft_k(dfft,ng,ngx,funct3d_g,funct3d_r_complex,'Wave',nl)
   ENDIF
   !
   IF ( mpime == 1 ) WRITE(*,*) "funct3d_r_complex when writing: ", funct3d_r_complex(1:10)
   !
   ALLOCATE( funct3d_r_complex_gathered(dfft%nr1x*dfft%nr2x*dfft%nr3x) )
   funct3d_r_complex_gathered = 0.0_DP
   CALL gather_grid(dfft,funct3d_r_complex,funct3d_r_complex_gathered)
   !IF ( mpime == 0 ) WRITE(*,*) "funct3d_r_complex_gathered when writing: ", funct3d_r_complex_gathered(:)
   !
   IF( dfft%mype == dfft%root ) THEN
      !
      ! 2) Encode 
      !
      IF( gamma_only ) THEN 
         ALLOCATE( funct3d_r_double( dfft%nr1x*dfft%nr2x*dfft%nr3x ) ) 
         funct3d_r_double = REAL(funct3d_r_complex_gathered(:),KIND=DP)
         DEALLOCATE(funct3d_r_complex_gathered)
         ndim = dfft%nr1x*dfft%nr2x*dfft%nr3x
         nbytes = SIZEOF(funct3d_r_double(1)) * ndim
         nlen = lenbase64(nbytes)
         ALLOCATE(CHARACTER(LEN=nlen) :: charbase64)
         IF (.NOT. islittleendian()) CALL base64_byteswap_double(nbytes,funct3d_r_double(1:ndim))
         CALL base64_encode_double(funct3d_r_double(1:ndim), ndim, charbase64)
         DEALLOCATE(funct3d_r_double)
         ctype = "double"
      ELSE
         ndim = dfft%nr1x*dfft%nr2x*dfft%nr3x
         nbytes = SIZEOF(funct3d_r_complex_gathered(1)) * ndim
         nlen = lenbase64(nbytes)
         ALLOCATE(CHARACTER(LEN=nlen) :: charbase64)
         IF (.NOT. islittleendian()) CALL base64_byteswap_complex(nbytes,funct3d_r_complex_gathered(1:ndim))
         CALL base64_encode_complex(funct3d_r_complex_gathered(1:ndim), ndim, charbase64)
         DEALLOCATE(funct3d_r_complex_gathered)
         ctype = "complex"
      ENDIF
      !
      ! 3) Write 
      !
      OPEN(NEWUNIT=iu,FILE=TRIM(ADJUSTL(fname)))
      !
      WRITE(iu,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
      WRITE(iu,'(a)') '<fpmd:function3d xmlns:fpmd="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0"'
      WRITE(iu,'(a)') 'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
      WRITE(iu,'(a)') 'xsi:schemaLocation="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0 function3d.xsd"'
      WRITE(iu,'(a)') 'name="',TRIM(ADJUSTL(descriptor)),'">'
      DO i = 1, 3
         WRITE(lab(i),'(f14.6)') celldm(1) * at(i,1) 
      ENDDO     
      WRITE(iu,'(a)') '<domain a="'//TRIM(ADJUSTL(lab(1)))//" "//TRIM(ADJUSTL(lab(2)))//" "//TRIM(ADJUSTL(lab(3)))//'"'
      DO i = 1, 3
         WRITE(lab(i),'(f14.6)') celldm(1) * at(i,2) 
      ENDDO     
      WRITE(iu,'(a)') 'b="'//TRIM(ADJUSTL(lab(1)))//" "//TRIM(ADJUSTL(lab(2)))//" "//TRIM(ADJUSTL(lab(3)))//'"'
      DO i = 1, 3
         WRITE(lab(i),'(f14.6)') celldm(1) * at(i,3) 
      ENDDO     
      WRITE(iu,'(a)') 'c="'//TRIM(ADJUSTL(lab(1)))//" "//TRIM(ADJUSTL(lab(2)))//" "//TRIM(ADJUSTL(lab(3)))//'"/>'
      WRITE(lab(1),'(i14)') dfft%nr1x
      WRITE(lab(2),'(i14)') dfft%nr2x
      WRITE(lab(3),'(i14)') dfft%nr3x
      WRITE(iu,'(a)') '<grid nx="'//TRIM(ADJUSTL(lab(1)))//'" ny="'//TRIM(ADJUSTL(lab(2)))//'" nz="'//TRIM(ADJUSTL(lab(3)))//'"/>'
      WRITE(iu,'(a)') '<grid_function type="'//ctype//'" nx="'//TRIM(ADJUSTL(lab(1)))//'" ny="'//TRIM(ADJUSTL(lab(2)))// &
            &'" nz="'//TRIM(ADJUSTL(lab(3)))//'" encoding="base64">'
      CALL write_long_string(iu,charbase64) 
      WRITE(iu,'(a)') '</grid_function>'
      WRITE(iu,'(a)') '</fpmd:function3d>'
      !
      CLOSE(iu)
      !
      DEALLOCATE(charbase64)
      DEALLOCATE(ctype)
      !
   ELSE 
      DEALLOCATE( funct3d_r_complex )
   ENDIF 
   !
   END SUBROUTINE 
 ! 
 !-----------------------------------------------------------------
   SUBROUTINE read_function3d ( dfft, fname, ng, ngx, funct3d_g, nmaps, nl)
   ! -----------------------------------------------------------------
   !
   USE kinds,                       ONLY : DP
   USE control_flags,               ONLY : gamma_only
   USE mp,                          ONLY : mp_bcast
   USE mp_bands,                    ONLY : me_bgrp, intra_bgrp_comm
   USE mp_world,                    ONLY : mpime
   USE fft_types,                   ONLY : fft_type_descriptor
   USE scatter_mod,                 ONLY : scatter_grid
   USE westcom,                     ONLY : fftdriver
   USE base64_module
   USE fourier_interpolation
   !
   IMPLICIT NONE
   !
   ! I/O 
   !
   TYPE(fft_type_descriptor), INTENT(IN) :: dfft 
   CHARACTER(LEN=*),INTENT(IN) :: fname
   INTEGER,INTENT(IN) :: ng, ngx, nmaps
   COMPLEX(DP),INTENT(OUT) :: funct3d_g(ngx)
   INTEGER,INTENT(IN) :: nl(nmaps,ngx)
   ! 
   ! Workspace
   !
   INTEGER :: iu, i, nlen, ndim, nbytes, nx, ny, nz 
   COMPLEX(DP),ALLOCATABLE :: funct3d_r_complex_gathered(:)
   COMPLEX(DP),ALLOCATABLE :: funct3d_r_complex(:)
   REAL(DP),ALLOCATABLE :: funct3d_r_double(:)
   INTEGER :: nlines, j, is, ie, ios, nline_start, nline_end
   CHARACTER(LEN=256) :: buffline
   CHARACTER(LEN=:),ALLOCATABLE :: bufftag
   CHARACTER(LEN=:),ALLOCATABLE :: buff2
   CHARACTER(LEN=:),ALLOCATABLE :: charbase64
   LOGICAL :: lread 
   CHARACTER(LEN=:),ALLOCATABLE :: ctype 
   !
   IF( dfft%mype == dfft%root ) THEN
      !
      OPEN(NEWUNIT=iu,FILE=TRIM(ADJUSTL(fname)))
      !
      DO
         READ(iu,'(a)',IOSTAT=ios) buffline
         IF( ios /=0 ) EXIT
         IF( INDEX(buffline,"<grid_function ") /= 0 ) THEN
            bufftag = TRIM(buffline)
            EXIT
         ENDIF
      ENDDO
      !
      DO WHILE ( INDEX(bufftag,">") == 0 )
         READ(iu,'(a)',IOSTAT=ios) buffline
         IF( ios /=0 ) EXIT
         bufftag = bufftag // TRIM(buffline)
      ENDDO
      !
      IF( INDEX(bufftag,">") == 0 ) THEN
        lread = .false.
      ELSE
        is = INDEX(bufftag,'nx="') + 4
        buff2 = bufftag(is:)
        ie = INDEX(buff2,'"') + is - 2
        READ(bufftag(is:ie),*) nx
        is = INDEX(bufftag,'ny="') + 4
        buff2 = bufftag(is:)
        ie = INDEX(buff2,'"') + is - 2
        READ(bufftag(is:ie),*) ny
        is = INDEX(bufftag,'nz="') + 4
        buff2 = bufftag(is:)
        ie = INDEX(buff2,'"') + is - 2
        READ(bufftag(is:ie),*) nz
        is = INDEX(bufftag,'type="') + 6
        buff2 = bufftag(is:)
        ie = INDEX(buff2,'"') + is - 2
        ctype = bufftag(is:ie)
        lread = .true.
      ENDIF
      !
      IF( lread ) THEN
         ndim = nx*ny*nz
         SELECT CASE (ctype)
         CASE("double") 
            nbytes = SIZEOF(1._DP) * ndim
         CASE("complex")
            nbytes = SIZEOF(CMPLX(1._DP,0._DP,KIND=DP)) * ndim
         CASE DEFAULT
         END SELECT 
         nlen = lenbase64(nbytes)
         ALLOCATE( CHARACTER(LEN=nlen) :: charbase64 )
         !
         CALL read_long_string(iu,charbase64)
         !
      ELSE
         CALL errore("","Could not start tag",1)
      ENDIF
      !
      !
      CLOSE(iu)
      !
   ENDIF
   !
   CALL mp_bcast(ndim,0,intra_bgrp_comm)
   CALL mp_bcast(nx,0,intra_bgrp_comm)
   CALL mp_bcast(ny,0,intra_bgrp_comm)
   CALL mp_bcast(nz,0,intra_bgrp_comm)
   !
   IF( nx /= dfft%nr1x ) CALL errore('read','Wrong nx',1) 
   IF( ny /= dfft%nr2x ) CALL errore('read','Wrong ny',1) 
   IF( nz /= dfft%nr3x ) CALL errore('read','Wrong nz',1) 
   !
   ALLOCATE( funct3d_r_complex_gathered(1:ndim) )
   !
   IF ( dfft%mype == dfft%root ) THEN
      SELECT CASE(ctype) 
      CASE("double")
         ALLOCATE( funct3d_r_double(1:ndim) )
         CALL base64_decode_double(charbase64(1:nlen), ndim, funct3d_r_double(1:ndim))
         DEALLOCATE( charbase64 )
         IF (.NOT. islittleendian()) CALL base64_byteswap_double(nbytes,funct3d_r_double(1:ndim)) 
         funct3d_r_complex_gathered(:) = CMPLX( funct3d_r_double(:), 0._DP, KIND=DP)
         DEALLOCATE( funct3d_r_double )
      CASE("complex")
         CALL base64_decode_complex(charbase64(1:nlen), ndim, funct3d_r_complex_gathered(1:ndim))
         DEALLOCATE( charbase64 )
         IF (.NOT. islittleendian()) CALL base64_byteswap_complex(nbytes,funct3d_r_complex_gathered(1:ndim)) 
      CASE DEFAULT
      END SELECT
   ENDIF
   !
   ALLOCATE( funct3d_r_complex(dfft%nnr) )
   !
   !IF ( me_bgrp == 0 ) WRITE(*,*) fname,"size of gathered : ", SIZE(funct3d_r_complex_gathered)
   !IF ( me_bgrp == 0 ) WRITE(*,*) fname,"size of distributed : ", SIZE(funct3d_r_complex)
   !
   !IF ( mpime == 0 ) WRITE(*,*) "funct3d_r_complex_gathered when reading: ", funct3d_r_complex_gathered(:)
   CALL scatter_grid(dfft,funct3d_r_complex_gathered,funct3d_r_complex)
   !
   !IF ( mpime == 1 ) WRITE(*,*) "funct3d_r_complex: ",  funct3d_r_complex(:)
   !
   IF ( gamma_only ) THEN
      SELECT CASE(fftdriver)
      CASE('Wave')
         CALL single_interp_fwfft_gamma(dfft,ng,ngx,funct3d_r_complex,funct3d_g,fftdriver,nl)
      CASE('Dense')
         CALL single_interp_fwfft_k(dfft,ng,ngx,funct3d_r_complex,funct3d_g,fftdriver,nl)
      END SELECT
   ELSE
      CALL single_interp_fwfft_k(dfft,ng,ngx,funct3d_r_complex,funct3d_g,fftdriver,nl)
   ENDIF
   !
   !IF ( mpime == 1 ) WRITE(*,*) "funct3d_g: ",  funct3d_g(1:10)
   !
   DEALLOCATE( funct3d_r_complex_gathered )
   DEALLOCATE( funct3d_r_complex )
   !
 END SUBROUTINE
 !
 !
 SUBROUTINE write_long_string(iu,longstring) 
   !
   ! Write a long string on multiple lines (each line has a max of 72 charachter)
   ! The unit "iu" is NOT opened and closed here 
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER,INTENT(IN) :: iu
   CHARACTER(LEN=*),INTENT(IN) :: longstring
   !
   ! Workspace
   !
   INTEGER :: j, nlines, thislen
   INTEGER, PARAMETER :: maxlen = 72
   !
   thislen = LEN(longstring)
   nlines = thislen / maxlen
   IF( MOD( thislen, maxlen ) > 0 ) nlines = nlines + 1
   DO j = 1, nlines
      WRITE(iu,'(a)') longstring((j-1)*maxlen+1:MIN(j*maxlen,thislen))
   ENDDO
   !
 END SUBROUTINE
 !
 !
 SUBROUTINE read_long_string(iu,longstring) 
   !
   ! Read a long string on multiple lines (each line has a max of 72 charachter)
   ! The unit "iu" is NOT opened and closed here 
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER,INTENT(IN) :: iu
   CHARACTER(LEN=*),INTENT(INOUT) :: longstring
   !
   ! Workspace
   !
   INTEGER :: j, nlines, thislen
   INTEGER, PARAMETER :: maxlen = 72
   !
   thislen = LEN(longstring)
   nlines = thislen / maxlen
   IF( MOD( thislen, maxlen ) > 0 ) nlines = nlines + 1
   !
   DO j = 1, nlines 
      READ(iu,'(a)') longstring((j-1)*maxlen+1:MIN(j*maxlen,thislen))
   ENDDO
   !
 END SUBROUTINE
 !
END MODULE
