!
! Copyright (C) 2015 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file:
! He Ma
!
!#define DEBUG_BASE64
#define C_BINDING
!
!-----------------------------------------------------------------------
MODULE qbox_interface
  !----------------------------------------------------------------------------
  !
  USE fft_base,       ONLY : dffts
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout
  USE mp,             ONLY : mp_barrier, mp_bcast, mp_sum
  USE mp_world,       ONLY : mpime, root
  USE mp_global,      ONLY : my_image_id, me_image, intra_image_comm
  USE io_files,       ONLY : delete_if_present
  USE fft_types,      ONLY : fft_type_descriptor
  USE scatter_mod,    ONLY : gather_grid, scatter_grid
  USE cubefile,       ONLY : write_wfc_cube_r, read_wfc_cube_r
  USE io_push,        ONLY : io_push_title, io_push_value
  USE iso_c_binding,  ONLY : c_double, c_int, c_char, c_bool, c_null_char
  USE function3d,     ONLY : write_function3d, read_function3d
  !USE base64_module,  ONLY : base64_type, c_base64__new, &
  !                           & c_base64__encode_double, c_base64__decode_double
  USE parallel_include
  USE wrappers,       ONLY : f_copy
  USE cpp_wrappers,   ONLY : c_sleep, c_wait_for_file
  USE pwcom,          ONLY : lsda
  USE westcom,        ONLY : nrowmax, xml_file, xc, alpha_pbe0, amplitude, wf_dyn, &
                             nitscf, nite, blHF, btHF, qbox_bisec_wfc_filename
  !
  IMPLICIT NONE
  !
  !----------------------------------------------------------------------------
  !
  ! currently only consider nks = nbgrb = 1 case
  !
  PRIVATE
  SAVE
  !
  !PUBLIC :: qbox_ks_wfc_filename, qbox_bisec_wfc_filename, load_qbox_wfc
  PUBLIC :: load_qbox_wfc, init_qbox_interface, finalize_qbox_interface, apply_kernel_by_qbox, add_debug_log
  !PUBLIC :: nrowmax, xml_file, xc, alpha_pbe0, amplitude, wf_dyn, nitscf, nite, io, blHF, btHF
  !
  ! parameters for using qbox calculations
  !LOGICAL            :: l_load_qbox_wfc
  !CHARACTER(LEN=256) :: qbox_ks_wfc_filename
  !CHARACTER(LEN=256) :: qbox_bisec_wfc_filename
  !
  ! parameters for doing qbox calculations
  !
  CHARACTER(LEN=256) :: ready_for_kernel = ''
  !
  !CHARACTER(LEN=256) :: xc      ! xc functional
  !CHARACTER(LEN=256) :: wf_dyn  ! wavefunction update algorithm
  !CHARACTER(LEN=256) :: blHF    ! bisection levels for HF exchange computation
  !REAL(DP) :: btHF              ! bisection threshold for HF exchange computation
  !REAL(DP) :: alpha_pbe0        ! alpha for PBE0 calculation
  !REAL(DP) :: amplitude         ! amplitude for vext
  !INTEGER :: nrowmax
  !INTEGER :: nitscf
  !INTEGER :: nite
  !
  CHARACTER(LEN=8)    :: path                    ! path to qbox working dir, image dependent
  !CHARACTER(LEN=256)  :: xml_file                ! xml file from qbox ground state calculation
  CHARACTER(LEN=256)  :: lock_file
  CHARACTER(LEN=256)  :: server_input_file
  CHARACTER(LEN=256)  :: server_output_file
  CHARACTER(LEN=256)  :: vext_file
  CHARACTER(LEN=256)  :: resp_file
  !
  CHARACTER(LEN=256)  :: resp_command
  !
  LOGICAL             :: onroot                  ! dffts%mype == dffts%root
  INTEGER             :: np                      ! global size of vextr and drhor array
  INTEGER             :: nploc                   ! size of vextr and drhor array on current processor
  INTEGER,ALLOCATABLE :: nplocs(:)               ! and other processors, determined by dffts
  INTEGER             :: lastproc                ! index of last processor to have FFT grid
  !
  INTEGER             :: nsec_max      = 7200    ! maximum waiting time for one qbox response calculation
  INTEGER             :: iu = 100                ! iu for all serial I/O
  !
  !TYPE(Base64_type)   :: xcdr                    ! Base64 transcoder
  !
  CHARACTER(len=10)   :: mype
  !
  !
  INTERFACE init_qbox_interface
  MODULE PROCEDURE init_qbox, init_qbox_for_kernel
  END INTERFACE
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
    SUBROUTINE load_qbox_wfc(evc, qbox_wfc_filename, current_spin, nbndval)
      !----------------------------------------------------------------------------
      !
      USE pwcom,                  ONLY : nelec,npw,npwx,omega
      USE westcom,                ONLY : fftdriver
      USE iotk_module
      USE fft_base,               ONLY : dffts
      USE fft_at_gamma,           ONLY : single_fwfft_gamma,single_invfft_gamma
      USE fourier_interpolation,  ONLY : ft_interpolate
      USE mp_global,              ONLY : my_image_id,inter_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)             :: current_spin,nbndval
      CHARACTER(LEN=256), INTENT(IN)  :: qbox_wfc_filename
      COMPLEX(DP),        INTENT(OUT) :: evc(npwx, nbndval)
      !
      INTEGER                         :: iwfc, nwfcs, nspin, ispin, is, ir
      INTEGER                         :: nx, ny, nz
      INTEGER                         :: iunit, ierr
      CHARACTER(iotk_attlenx)         :: attr
      CHARACTER( LEN = 20 )           :: namespin
      !
      REAL(DP), ALLOCATABLE           :: psir(:)
      LOGICAL                         :: exst
      !
      CHARACTER(iotk_attlenx)         :: slen
      !
      ! read qbox xml file and overwrite current evc array
      !
      CALL io_push_title('Loading Qbox wavefunction...')
      !
      ! root of first image read qbox xml file, fourier interpolate to evc of whole image,
      ! then first image bcast evc to other images
      !
      evc(:,:) = (0._DP,0.0_DP)
      !
      IF ( me_image == 0 ) THEN
         !
         WRITE(stdout,*)  "    QboxInterface: reading from file ", TRIM(qbox_wfc_filename)
         !
         CALL iotk_free_unit( iunit )
         CALL iotk_open_read( iunit, FILE = TRIM(qbox_wfc_filename) )
         !
         CALL iotk_scan_begin( iunit, "wavefunction",  attr )
         CALL iotk_scan_attr( attr, "nspin", nspin )
         !
         ! read grid size
         !
         CALL iotk_scan_empty( iunit, "grid", attr )
         CALL iotk_scan_attr( attr, "nx", nx )
         CALL iotk_scan_attr( attr, "ny", ny )
         CALL iotk_scan_attr( attr, "nz", nz )
         !
         WRITE(stdout,*)  "    QboxInterface: grid size: ", nx, ny, nz
         !
         ALLOCATE( psir(nx*ny*nz) )
         !
      ENDIF
      !
      CALL mp_bcast( nx, 0, intra_image_comm )
      CALL mp_bcast( ny, 0, intra_image_comm )
      CALL mp_bcast( nz, 0, intra_image_comm )
      CALL mp_bcast( nspin, 0, intra_image_comm )
      !
      ! read KS orbitals
      !
      WRITE(stdout,*) "    QboxInterface: loading KS orbitals"
      !
      DO is = 1, nspin
         !
         IF (nspin > 1) THEN
            !
            IF ( me_image == 0 ) THEN
               !
               CALL iotk_scan_begin( iunit, "slater_determinant", attr )
               CALL iotk_scan_attr( attr, "spin", namespin )
               CALL iotk_scan_attr( attr, "size", nwfcs )
               !
            ENDIF
            !
            CALL mp_bcast( namespin, 0, intra_image_comm )
            CALL mp_bcast( nwfcs, 0, intra_image_comm )
            !
            IF (namespin == 'up' )   ispin = 1
            IF (namespin == 'down' ) ispin = 2
            !
         ELSE
            !
            IF ( me_image == 0 ) THEN
               !
               CALL iotk_scan_begin( iunit, "slater_determinant", attr )
               CALL iotk_scan_attr( attr, "size", nwfcs )
               !
            ENDIF
            !
            CALL mp_bcast( nwfcs, 0, intra_image_comm )
            !
            namespin = 'none'
            ispin = 1
            !
         ENDIF
         !
         IF (ispin == current_spin) THEN
            !
            WRITE(stdout,*)  "          ",  namespin, nwfcs, nspin, ispin
            !
            DO iwfc=1, nwfcs
               !
               IF ( me_image == 0 ) THEN
                  !
                  CALL iotk_scan_begin(iunit, "grid_function", slen)
                  READ(iunit,*) ( psir(ir), ir = 1, nx*ny*nz )
                  CALL iotk_scan_end(iunit, "grid_function")
                  !
               ENDIF
               !
               CALL ft_interpolate(psir, nx, ny, nz, nwfcs, evc(:,iwfc), npw, npwx)
               !
            ENDDO
            !
         ENDIF
         !
         IF ( me_image == 0 ) CALL iotk_scan_end( iunit, "slater_determinant" )
         !
      ENDDO
      !
      IF ( me_image == 0 ) THEN
         !
         CALL iotk_scan_end( iunit, "wavefunction" )
         !
         CALL iotk_close_read( iunit )
         !
      ENDIF
      !
      IF ( me_image == 0 ) DEALLOCATE(psir)
      !
    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
    SUBROUTINE init_qbox()
      !----------------------------------------------------------------------------
      !
      INTEGER                         ::  proc, ierr, numsp
      !
      ! initialize Qbox to load wavefunction, set parameters, etc.
      !
      WRITE(mype, '(I0.6)') dffts%mype
      !
      CALL mp_barrier(intra_image_comm)
      !
      ! first time to run init_qbox_interface
      !
      CALL io_push_title("Qbox interface setup")
      !
#ifndef C_BINDING
      IF ( DP /= C_DOUBLE ) call errore( 'qbox_interface', 'DP /= C_DOUBLE but C_BINDING is not used', 1 )
#endif
      !
      !IF ( TRIM(io) == "base64_serial" .OR. TRIM(io) == "base64_parallel")  xcdr%object = c_base64__new()
      !
      onroot = ( dffts%mype == dffts%root )
      !
      ! compute size of vextr and drhor arrays
      !
      !ALLOCATE(nplocs(0:dffts%nproc - 1))
      !lastproc = 0
      !DO proc = 0, ( dffts%nproc - 1 )
      !   nplocs(proc) = dffts%nnp * dffts%npp(proc+1)
      !   IF ( nplocs(proc) > 0 ) lastproc = proc
      !ENDDO
      !nploc = nplocs(dffts%mype)
      !
      !CALL MPI_ALLREDUCE(nploc, np, 1, MPI_INT, MPI_SUM, intra_image_comm, ierr)
      !IF ( np /= dffts%nr1 * dffts%nr2 * dffts%nr3 ) CALL errore('qbox_interface', 'np /= dffts%nr1 * dffts%nr2 * dffts%nr3',1)
      !
      ! define qbox working directory and file names
      !
      WRITE(path,'(I6.6)') my_image_id
      path = 'I' // TRIM(path) // '/'
      server_input_file       = 'qb.in'
      lock_file = TRIM(server_input_file) // '.lock'
      server_output_file      = 'qb.out'
      vext_file               = 'vext.dat'
      resp_file = TRIM(vext_file) // '.response'
      !
      numsp = 30
      CALL io_push_value("Qbox input file", TRIM(server_input_file), numsp)
      CALL io_push_value("Qbox output file", TRIM(server_output_file), numsp)
      CALL io_push_value("potential file", TRIM(vext_file), numsp)
      CALL io_push_value("response file", TRIM(resp_file), numsp)
      !
      ! initialize qbox
      !
      IF(me_image == 0) THEN
         OPEN(UNIT=iu, FILE=path//TRIM(server_input_file))
         IF ( nrowmax /= 0 ) WRITE(iu,*) 'set nrowmax ', nrowmax
         WRITE(iu,'(a)') 'load ../' // TRIM(xml_file)
         WRITE(iu,'(a)') ' '
         WRITE(iu,'(a)') 'set xc ' // TRIM( xc )
         IF ( TRIM(xc) == 'PBE0' ) WRITE(iu, *) 'set alpha_PBE0 ', alpha_pbe0
         WRITE(iu,'(a)') 'set wf_dyn ' // TRIM( wf_dyn )
         WRITE(iu,'(a)') 'set blHF ' // TRIM( blHF )
         WRITE(iu, *) 'set btHF ', btHF
         CLOSE(iu)
      ENDIF
      !
      CALL wait_for_lock_file()  ! wait for qbox to start
      !
      CALL delete_lock_file()    ! now qbox should start initialization
      !
      CALL io_push_title("Initializing Qbox")
      !
      CALL wait_for_lock_file()
      !
      ! archive input and output files for initialization
      !
      IF(me_image == 0) THEN
         ierr = f_copy( path // TRIM( server_input_file), &
                      & path // (TRIM( server_input_file ) ) // '.init')
         IF ( ierr /= 0 ) CALL errore("qbox_interface", "fail to archive qbox input")
      ENDIF
      CALL archive_qbox_output()
      !
      CALL mp_barrier(intra_image_comm)
      !
    END SUBROUTINE
    !
    !
    !----------------------------------------------------------------------------
    SUBROUTINE init_qbox_for_kernel( kernel )
      !----------------------------------------------------------------------------
      !
      CHARACTER(LEN=256), INTENT(IN)    ::  kernel
      INTEGER                           ::  proc, ierr, numsp
      CHARACTER(LEN=256)                ::  tmp
      !
      ! generate response command for given kernel
      !
      IF ( TRIM(ready_for_kernel) == TRIM(kernel) ) RETURN
      !
      resp_command = 'response -vext ' // TRIM( vext_file )
      !
      ! response mode
      !
      SELECT CASE( TRIM(kernel) )
      CASE('CHI')
      CASE('CHI_RPA')
         resp_command = TRIM( resp_command ) // ' -RPA'
      CASE('CHI0')
         resp_command = TRIM( resp_command ) // ' -IPA'
      CASE DEFAULT
         CALL errore("init_qbox_for_kernel", "qbox_interface: wrong kernel: " // TRIM(kernel), 1)
      END SELECT
      !
      ! amplitude
      !
      WRITE(tmp,'(E12.5)') amplitude
      resp_command = TRIM( resp_command ) // ' -amplitude ' // TRIM(tmp)
      !
      ! io scheme
      !
      !resp_command = TRIM( resp_command ) // ' -io ' // TRIM(io)
      !IF ( TRIM(io) == 'base64_serial' .OR. TRIM(io) == 'base64_parallel') THEN
      !   WRITE(tmp,'(I5)') dffts%nr1
      !   resp_command = TRIM( resp_command ) // ' -nx ' // TRIM(tmp)
      !   WRITE(tmp,'(I5)') dffts%nr2
      !   resp_command = TRIM( resp_command ) // ' -ny ' // TRIM(tmp)
      !   WRITE(tmp,'(I5)') dffts%nr3
      !   resp_command = TRIM( resp_command ) // ' -nz ' // TRIM(tmp)
      !ENDIF
      !
      ! nitscf and nite
      !
      WRITE(tmp,'(I5)') nitscf
      resp_command = TRIM( resp_command ) // ' ' // TRIM(tmp)
      WRITE(tmp,'(I5)') nite
      resp_command = TRIM( resp_command ) // ' ' // TRIM(tmp)
      !
      ! write response command to file
      !
      IF(me_image == 0) THEN
         !
         OPEN(UNIT=iu, FILE=path//TRIM(server_input_file))
         WRITE(iu,'(A)') TRIM(resp_command)
         CLOSE(iu)
         !
      ENDIF
      !
      ready_for_kernel = kernel
      !
      CALL mp_barrier(intra_image_comm)
      !
    END SUBROUTINE
    !
    !
    !----------------------------------------------------------------------------
    SUBROUTINE apply_kernel_by_qbox(kernel, fr)
      !----------------------------------------------------------------------------
      !
      ! Applying kernel to real space function fr using response command of Qbox
      !
      ! kernel: "CHI", "CHI_RPA" or "CHI0"
      ! fr: real space function distributed according to dffts
      !
      CHARACTER(LEN=256), INTENT(IN)  ::  kernel
      REAL(DP), INTENT(INOUT)         ::  fr(dffts%nnr)
      REAL(DP), ALLOCATABLE           ::  frspin(:,:)
      !
      ! Generate response command to qbox
      !
      CALL init_qbox_interface(kernel)
      !
      fr = 0.5_DP * fr  ! change from rydberg to hartree
      !
      ! Write fr --> v.dat
      !
      CALL write_vext_file(fr)
      !
      ! Remove lock file, let qbox resume running
      !
      CALL delete_lock_file()
      !
      ! Wait for lock file to appear again, archive qbox output
      !
      CALL start_clock( 'wait_qbox' )
      CALL wait_for_lock_file()
      CALL stop_clock( 'wait_qbox' )
      CALL archive_qbox_output()
      !
      ! Read v.dat.response --> fr
      !
      IF (lsda) THEN
         !
         ALLOCATE(frspin(dffts%nnr,2))
         !
         CALL read_resp_file(frspin(:,1), 0)
         CALL read_resp_file(frspin(:,2), 1)
         !
         fr(:) = frspin(:,1) + frspin(:,2)
         !
         DEALLOCATE(frspin)
         !
      ELSE
         !
         CALL read_resp_file(fr)
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !
    !----------------------------------------------------------------------------
    SUBROUTINE write_vext_file(vextr)
      !----------------------------------------------------------------------------
      !
      REAL(DP), INTENT(OUT)           :: vextr(:)
      !
      !REAL(DP),ALLOCATABLE            :: vextr_gathered(:)              ! gathered vextr on processor 0
      !REAL(DP),ALLOCATABLE            :: sbuf(:), rbuf(:)               ! send buf, receive buf
      !REAL(DP),ALLOCATABLE            :: vextr_adj(:)                   ! vextr with adjusted size for base64 encoding
      !
      !CHARACTER,ALLOCATABLE           :: wbuf(:)                        ! write buf (of size nchars)
      !
      !INTEGER              :: sbufsize, rbufsize, sbufsize0, rbufsize0
      !INTEGER              :: nploc_adj                                 ! local size of vextr after adjusting
      !INTEGER              :: nchars                                    ! size of base64 encoded vextr
      !INTEGER              :: offset                                    ! offset for MPI_File_write_at_all
      !
      !INTEGER              :: stat(MPI_STATUS_SIZE)
      !INTEGER              :: proc, fh, info, ierr, i
      !
      CALL mp_barrier(intra_image_comm)
      !
      CALL start_clock( 'write_vext' )


      CALL write_function3d(path//TRIM(vext_file),vextr,dffts,onroot)
      !
      !SELECT CASE (TRIM(io))
      !CASE ("base64_parallel")
         !
      !   IF ( dffts%nproc == 1 )  GOTO 100
!#ifdef DEBUG_BASE64
         !
      !   IF( onroot ) THEN
      !      ALLOCATE( vextr_gathered(dffts%nr1x*dffts%nr2x*dffts%nr3x) )
      !   ENDIF
         !
      !  CALL gather_grid(dffts,vextr,vextr_gathered)
      !  IF( onroot ) THEN
      !     !
      !     OPEN(UNIT=iu, FILE=path//TRIM(vext_file)//".unencoded")
      !     DO i=1, size(vextr_gathered)
      !        WRITE(iu, *) vextr_gathered(i), ' '
      !     END DO
      !     CLOSE(iu)
      !     !
      !  ENDIF
      !  IF( onroot ) DEALLOCATE( vextr_gathered )
!#endif
         !
      !   IF ( size(vextr) < nploc ) CALL errore("write_vext_file: size(vextr) < nploc", 1)
         !
         ! adjust size of vextr array on each processor such that it is multiple of 3
         !
!         sbufsize = 0
!         rbufsize = 0
!         DO proc = 0, lastproc
!            IF ( proc == 0 ) THEN
!               sbufsize0 = MOD(nplocs(proc), 3)
!               rbufsize0 = 0
!            ELSE IF ( proc == lastproc ) THEN
!               rbufsize0 = sbufsize0
!               sbufsize0 = 0
!            ELSE
!               rbufsize0 = sbufsize0
!               sbufsize0 = rbufsize0 + MOD(nplocs(proc), 3)
!            ENDIF
!            IF ( dffts%mype == proc ) THEN
!               sbufsize = sbufsize0
!               rbufsize = rbufsize0
!            ENDIF
!         ENDDO
!         !
!         IF ( sbufsize >= nploc ) CALL errore("write_vext_file: sbufsize >= nploc", 1)
!         IF ( sbufsize > 0 ) THEN
!           ALLOCATE(sbuf(sbufsize))
!           sbuf = vextr(nploc-sbufsize+1:nploc)
!         ENDIF
!         !
!         IF ( rbufsize > 0 ) THEN
!           ALLOCATE(rbuf(rbufsize))
!           rbuf = 0.
!         ENDIF
!         !
!         IF ( MOD(dffts%mype, 2) == 0 .AND. sbufsize > 0 ) THEN
!            CALL MPI_SEND(sbuf, sbufsize, MPI_DOUBLE_PRECISION, dffts%mype + 1, 0, intra_image_comm, ierr)
!            IF ( ierr /= 0 ) CALL errore("write_vext_file: error in MPI_SEND", 1)
!         ELSE IF ( MOD(dffts%mype, 2) == 1 .AND. rbufsize > 0 ) THEN
!            CALL MPI_RECV(rbuf, rbufsize, MPI_DOUBLE_PRECISION, dffts%mype - 1, 0, intra_image_comm, stat, ierr)
!            IF ( ierr /= 0 ) CALL errore("write_vext_file: error in MPI_RECV", 1)
!         ENDIF
!         !
!         IF ( MOD(dffts%mype, 2) == 1 .AND. sbufsize > 0 ) THEN
!            CALL MPI_SEND(sbuf, sbufsize, MPI_DOUBLE_PRECISION, dffts%mype + 1, 0, intra_image_comm, ierr)
!            IF ( ierr /= 0 ) CALL errore("write_vext_file: error in MPI_SEND", 1)
!         ELSE IF ( MOD(dffts%mype, 2) == 0 .AND. rbufsize > 0 ) THEN
!            CALL MPI_RECV(rbuf, rbufsize, MPI_DOUBLE_PRECISION, dffts%mype - 1, 0, intra_image_comm, stat, ierr)
!            IF ( ierr /= 0 ) CALL errore("write_vext_file: error in MPI_RECV", 1)
!         ENDIF
!         !
!         ! assemble vextr_adj with vextr and rbuf
!         !
!         nploc_adj = nploc + rbufsize - sbufsize
!         ALLOCATE(vextr_adj(nploc_adj))
!         IF ( dffts%mype < lastproc .AND. MOD(nploc_adj, 3) /= 0 )  CALL errore("write_vext_file: MOD(nploc_adj, 3) /= 0", 1)
!         IF ( rbufsize > 0 ) vextr_adj(1:rbufsize) = rbuf(:)
!         IF ( dffts%mype <= lastproc) vextr_adj(rbufsize+1:) = vextr(1:nploc-sbufsize)
!         !
!         ! encode vextr_adj into wbuf
!         !
!         nchars = 4 * ( ( nploc_adj * DP + 2 ) / 3 )
!         ALLOCATE(wbuf(nchars))
!         !
!         CALL c_base64__encode_double(xcdr%object, vextr_adj, nploc_adj, wbuf)
!         !
!         ! compute offset and write to file
!         !
!         CALL MPI_SCAN(nchars, offset, 1, MPI_INT, MPI_SUM, intra_image_comm, ierr)
!         offset = offset - nchars
!         !
!         CALL MPI_FILE_OPEN(intra_image_comm, path//TRIM(vext_file), &
!                            & IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)
!         IF ( ierr /=0 ) CALL errore('write_vext_file', 'MPI_File_open error', 1)
!         !
!         CALL MPI_FILE_SET_SIZE(fh, INT(0, KIND=MPI_OFFSET_KIND), ierr)
!         IF ( ierr /=0 ) CALL errore('write_vext_file', 'MPI_FILE_SET_SIZE error', 1)
!         !
!         CALL MPI_FILE_WRITE_AT_ALL(fh, INT(offset, KIND=MPI_OFFSET_KIND), &
!                                    & wbuf, nchars, MPI_CHARACTER, stat, ierr)
!         IF ( ierr /=0 ) CALL errore('write_vext_file', 'MPI_FILE_WRITE_AT_ALL error', 1)
!         !
!         CALL MPI_FILE_CLOSE(fh, ierr)
!         IF ( ierr /=0 ) CALL errore('write_vext_file', 'MPI_FILE_CLOSE error', 1)
!         !
!         IF (sbufsize > 0) DEALLOCATE(sbuf)
!         IF (rbufsize > 0) DEALLOCATE(rbuf)
!         DEALLOCATE(vextr_adj)
!         DEALLOCATE(wbuf)
!         !
!      CASE ("base64_serial")
!         !
!100      ALLOCATE( vextr_gathered(dffts%nr1x*dffts%nr2x*dffts%nr3x) )
!         !IF( onroot )  ALLOCATE( vextr_gathered(dffts%nr1x*dffts%nr2x*dffts%nr3x) )
!         !
!         CALL gather_grid(dffts,vextr,vextr_gathered)
!         !
!         IF( onroot ) THEN
!            !
!#ifdef DEBUG_BASE64
!            OPEN(UNIT=iu, FILE=path//TRIM(vext_file)//".unencoded")
!            DO i=1, size(vextr_gathered)
!               WRITE(iu, *) vextr_gathered(i), ' '
!            END DO
!            CLOSE(iu)
!#endif
!            !
!            nchars = 4 * ( ( np * DP + 2 ) / 3 )
!            ALLOCATE(wbuf(nchars))
!            !
!            CALL c_base64__encode_double(xcdr%object, vextr_gathered, np, wbuf)
!            !
!            OPEN(UNIT=iu, FILE=path//TRIM(vext_file), ACCESS="STREAM")
!            WRITE(iu) wbuf
!            CLOSE(iu)
!            !
!            DEALLOCATE(wbuf)
!            !
!         ENDIF
!         !
!         !IF( onroot ) DEALLOCATE(vextr_gathered)
!         DEALLOCATE(vextr_gathered)
!         !
!      CASE ("cube")
!         !
!         CALL write_wfc_cube_r(dffts,iu,path//TRIM(vext_file),vextr)
!         !
!      CASE DEFAULT
!         !
!         CALL errore('write_vext_file', 'wrong io mode', 1)
!         !
!      END SELECT
      !
      CALL mp_barrier(intra_image_comm)
      !
      CALL stop_clock( 'write_vext' )
      !
    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
    SUBROUTINE read_resp_file(drhor, ispin)
      !----------------------------------------------------------------------------
      !
      INTEGER, OPTIONAL, INTENT(IN)  :: ispin
      REAL(DP), INTENT(OUT)  :: drhor(:)
      !
      !REAL(DP),ALLOCATABLE   :: drhor_gathered(:)
      !REAL(DP),ALLOCATABLE   :: tmpr(:)
!#ifdef C_BINDING
!      REAL(C_DOUBLE),ALLOCATABLE   :: tmpr_C(:)
!#endif
      !
      !CHARACTER,ALLOCATABLE  :: rbuf(:)
      !
      !INTEGER                :: fh
      !INTEGER                :: nstart, nend, gstart, gend, istart, iend
      !INTEGER                :: nfloats, offset, nbytes, nchars, rbufsize
      !INTEGER(KIND=MPI_OFFSET_KIND)   :: filesize
      !INTEGER                :: stat(MPI_STATUS_SIZE)
      !INTEGER                :: ierr, i
      !
      CHARACTER(len=10)      :: mype
      CHARACTER(LEN=256)     :: resp_file_to_read
      !
      IF (PRESENT(ispin)) THEN
         !
         IF (.NOT. lsda) CALL errore("present ispin, but not lsda calculation ", 1)
         !
         SELECT CASE (ispin)
         !
         CASE(0)
             !
             resp_file_to_read = path//TRIM(resp_file)//".spin0"
             !
         CASE(1)
             !
             resp_file_to_read = path//TRIM(resp_file)//".spin1"
             !
         CASE DEFAULT
             !
             CALL errore("present wrong ispin, ", 1)
             !
         END SELECT
         !
      ELSE
         !
         resp_file_to_read = path//TRIM(resp_file)
         !
      ENDIF
      !
      WRITE(mype, '(I0.6)') dffts%mype
      !
      CALL mp_barrier(intra_image_comm)
      !
      CALL start_clock( 'read_resp' )
      !
      CALL read_function3d(resp_file_to_read,drhor,dffts,onroot)
!      SELECT CASE (TRIM(io))
!      CASE ('base64_parallel')
!         !
!         IF ( dffts%nproc == 1 )  GOTO 200   ! go to serial case
!         !
!         ! compute offset for MPI_File_read_at_all
!         ! for simplicity, let the size of file (in bytes) read by each processor
!         ! to be N*sizeof(double)*4, which correspond to N*3 float numbers
!         ! denote each 3 floats as a group
!         !
!         CALL MPI_SCAN(nploc, nend, 1, MPI_INT, MPI_SUM, intra_image_comm, ierr)
!         IF ( ierr /=0 ) CALL errore('read_resp_file', 'MPI_SCAN error',1)
!         !
!         nstart = nend - nploc + 1
!         gstart = (nstart - 1) / 3 + 1
!         gend = (nend - 1) / 3 + 1
!         offset =  4 * C_DOUBLE * (gstart - 1)  ! assuming the input file is encoded by C program
!         !
!         ! allocate read buffer
!         !
!         nfloats = 0
!         IF ( dffts%mype < lastproc ) THEN
!            nfloats = 3 * ( gend - gstart + 1 )
!         ELSE IF ( dffts%mype == lastproc ) THEN
!            nfloats = 3 * ( gend - gstart ) + ( MOD(np - 1, 3) + 1 )
!         ENDIF
!         IF ( nfloats < nploc ) CALL errore('read_resp_file', 'nfloats < nploc',1)
!         !
!         nbytes = nfloats * C_DOUBLE
!         nchars = 4 * ( ( nbytes + 2 ) / 3 )
!         IF ( dffts%mype < lastproc .AND. MOD(nchars, 4*C_DOUBLE) /= 0 ) THEN
!            CALL errore('read_resp_file', 'MOD(nchars, 4*=C_DOUBLE) /= 0',1)
!         ENDIF
!         !
!         ALLOCATE(rbuf(nchars))
!         !
!         ! read file
!         !
!         CALL MPI_File_open(intra_image_comm, resp_file_to_read, MPI_MODE_RDONLY , MPI_INFO_NULL, fh, ierr)
!         IF ( ierr /=0 ) CALL errore('read_resp_file', 'MPI_File_open error',1)
!         !
!         !CALL MPI_FILE_GET_SIZE(fh, filesize, ierr)
!         !
!         CALL MPI_File_read_at_all(fh, INT(offset, KIND=MPI_OFFSET_KIND), rbuf, nchars, MPI_CHARACTER, stat, ierr)
!         IF ( ierr /=0 ) CALL errore('read_resp_file', 'MPI_File_read_at_all error',1)
!         !
!         CALL MPI_File_close(fh, ierr)
!         IF ( ierr /=0 ) CALL errore('read_resp_file', 'MPI_File_close error',1)
!         !
!         ! decode rbuf into tmpr, then dump to drhor
!         !
!#ifdef C_BINDING
!         ALLOCATE(tmpr_C(nfloats))
!         ALLOCATE(tmpr(nfloats))
!         !
!         CALL c_base64__decode_double(xcdr%object, rbuf, nfloats, tmpr_C)
!         !
!         tmpr = REAL(tmpr_C, KIND=DP)
!         DEALLOCATE(tmpr_C)
!#else
!         ALLOCATE(tmpr(nfloats))
!         !
!         CALL c_base64__decode_double(xcdr%object, rbuf, nfloats, tmpr)
!#endif
!         !
!         istart = MOD((nstart - 1), 3) + 1
!         iend = istart + nploc - 1
!         drhor(1:nploc) = tmpr(istart:iend)
!         !
!#ifdef DEBUG_BASE64
!         OPEN(UNIT=iu, FILE=resp_file_to_read//'.decoded_by_west_pe'//mype)
!         DO i=1, nploc
!            WRITE(iu, *) drhor(i), ' '
!         END DO
!         CLOSE(iu)
!#endif
!         DEALLOCATE(rbuf)
!         DEALLOCATE(tmpr)
!         !
!      CASE ('base64_serial')
!         !
!200      ALLOCATE( drhor_gathered(dffts%nr1x*dffts%nr2x*dffts%nr3x) )
!         IF( onroot ) THEN
!            !
!            !ALLOCATE( drhor_gathered(dffts%nr1x*dffts%nr2x*dffts%nr3x) )
!            drhor_gathered = 0._DP
!            nfloats = np
!            nbytes = nfloats * C_DOUBLE
!            nchars = 4 * ( ( nbytes + 2 ) / 3 )
!            ALLOCATE(rbuf(nchars))
!            !
!            OPEN(UNIT=iu, FILE=resp_file_to_read, ACCESS='STREAM')
!            !OPEN(UNIT=iu, FILE=resp_file_to_read, ACCESS='DIRECT')
!            READ(iu) rbuf
!            CLOSE(iu)
!            !
!#ifdef C_BINDING
!            ALLOCATE(tmpr_C(dffts%nr1x*dffts%nr2x*dffts%nr3x))
!            CALL c_base64__decode_double(xcdr%object, rbuf, nfloats, tmpr_C)
!            drhor_gathered = REAL(tmpr_C, KIND=DP)
!            DEALLOCATE(tmpr_C)
!#else
!            CALL c_base64__decode_double(xcdr%object, rbuf, nfloats, drhor_gathered)
!#endif
!            !
!#ifdef DEBUG_BASE64
!            OPEN(UNIT=iu, FILE=resp_file_to_read//'.decoded_by_west')
!            DO i=1, dffts%nr1x*dffts%nr2x*dffts%nr3x
!               WRITE(iu, *) drhor_gathered(i)
!            END DO
!            CLOSE(iu)
!#endif
!            !
!         ENDIF
!         !
!         CALL scatter_grid(dffts,drhor_gathered,drhor)
!         !
!         IF( onroot ) THEN
!            DEALLOCATE(rbuf)
!            !DEALLOCATE(drhor_gathered)
!         ENDIF
!         !
!         DEALLOCATE(drhor_gathered)
!         !
!      CASE ('cube')
!         !
!         CALL read_wfc_cube_r(dffts,iu,resp_file_to_read,drhor)
!         !
!      CASE DEFAULT
!         CALL errore('read_resp_file', 'Wrong io mode',1)
!      END SELECT
      !
      CALL mp_barrier(intra_image_comm)
      !
      CALL stop_clock( 'read_resp' )
      !
    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
    SUBROUTINE delete_lock_file()
      !----------------------------------------------------------------------------
      !
      IF(me_image == 0) THEN
         CALL delete_if_present( path // TRIM( lock_file ) )
      ENDIF
      !
    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
    SUBROUTINE wait_for_lock_file()
      !----------------------------------------------------------------------------
      !
      ! wait for qbox to finish (indicated by existance of a lock file)
      !
      LOGICAL              :: file_exists = .FALSE.
      LOGICAL(KIND=C_BOOL) :: file_existss = .FALSE.
      INTEGER :: isec
      !
      IF( me_image == 0 ) THEN
        !
        CALL c_wait_for_file(INT(nsec_max, KIND=C_INT), file_existss, &
                             & path // TRIM( lock_file ) // C_NULL_CHAR)
        !
        ! DO isec = 1,nsec_max
        !    !
        !    INQUIRE( FILE= path // TRIM( lock_file ) , EXIST = file_exists )
        !    IF( file_exists ) THEN
        !       EXIT
        !    ELSE
        !       CALL c_sleep(1000000_C_INT)
        !    ENDIF
        !    !
        ! ENDDO
        !
        file_exists = file_existss
        !
      ENDIF
      !
      CALL mp_bcast( file_exists, 0, intra_image_comm )
      IF ( .NOT. file_exists ) CALL errore( 'qbox_interface', 'Lock file timeout', 1 )
      !
    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
    SUBROUTINE archive_qbox_output()
      !----------------------------------------------------------------------------
      !
      ! rename qbox output file to prevent it from being overwritten by next qbox run
      !
      LOGICAL             :: file_exists
      CHARACTER(LEN=256)  :: new_server_output_file
      CHARACTER (len=6)   :: file_index
      INTEGER             :: i, ierr
      !
      IF( me_image == 0 ) THEN
         !
         file_exists = .FALSE.
         INQUIRE( FILE= path // TRIM( server_output_file ) , EXIST = file_exists )
         IF( .NOT. file_exists ) CALL errore( 'qbox_interface', 'file to achive not exist', 1 )
         !
         i = 0
         DO
            WRITE(file_index,'(i6.6)') i
            new_server_output_file = TRIM(server_output_file) // "." // file_index
            INQUIRE( FILE= path // TRIM( new_server_output_file ) , EXIST = file_exists )
            IF( file_exists ) THEN
               i = i + 1
            ELSE
               ierr = f_copy( path // TRIM( server_output_file ), path // TRIM( new_server_output_file ) )
               IF ( ierr /= 0 ) CALL errore("qbox_interface", "fail to archive qbox output")
               EXIT
            ENDIF
         ENDDO
         !
      ENDIF
      !
    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
    SUBROUTINE finalize_qbox_interface()
      !----------------------------------------------------------------------------
      !
      ! send quit command to qbox
      !
      IF(me_image == 0) THEN
         OPEN(UNIT=iu, FILE=path//TRIM(server_input_file))
         WRITE(iu,'(A)') 'quit'
         CLOSE(iu)
      ENDIF
      !
      CALL delete_lock_file()
      !
      !DEALLOCATE(nplocs)
      !
    END SUBROUTINE
    !
  ! !----------------------------------------------------------------------------
  ! SUBROUTINE movefile(from, to)
  !   !----------------------------------------------------------------------------
  !   !

  !   !
  !   CHARACTER(len=*), INTENT(IN)         :: from, to
  !   CHARACTER,ALLOCATABLE                :: tmp(:)
  !   INTEGER                              :: filesize
  !   !
  !   INQUIRE(FILE=from, SIZE=filesize)
  !   ALLOCATE(tmp(filesize))
  !   !
  !   OPEN(iu, FILE=from, STATUS="OLD", ACCESS="STREAM")
  !   READ(iu) tmp
  !   CLOSE(iu)
  !   !
  !   OPEN(iu, FILE=to, STATUS="NEW", ACCESS="STREAM")
  !   WRITE(iu) tmp
  !   CLOSE(iu)
  !   !
  !   DEALLOCATE(tmp)
  !   !
  !   CALL delete_if_present(from)
  !   !
  ! END SUBROUTINE
  !
  !----------------------------------------------------------------------------
  SUBROUTINE add_debug_log(message)
    !----------------------------------------------------------------------------
    !
    CHARACTER(LEN=*), INTENT(IN) ::  message
    CHARACTER(LEN=256)           ::  logfile
    LOGICAL                      ::  file_exists = .FALSE.
    INTEGER                      ::  iu_dbg = 101
    !
    logfile = path//'dbg_west_pe'//mype
    INQUIRE( FILE=logfile, EXIST=file_exists )
    IF ( file_exists ) THEN
       OPEN(UNIT=iu_dbg, FILE=logfile, STATUS="OLD", POSITION="APPEND")
    ELSE
       OPEN(UNIT=iu_dbg, FILE=logfile)
    ENDIF
    !
    WRITE(iu_dbg, *) TRIM(message)
    CLOSE(iu_dbg)
    CALL c_sleep(100000_C_INT)
    !
  END SUBROUTINE
END MODULE
