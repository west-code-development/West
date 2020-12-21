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
!#define C_BINDING
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
  !USE iso_c_binding,  ONLY : c_double, c_int, c_char, c_bool, c_null_char
  USE function3d,     ONLY : write_function3d, read_function3d
  USE parallel_include
  USE wrappers,       ONLY : f_copy
  !USE cpp_wrappers,   ONLY : c_sleep, c_wait_for_file
  USE pwcom,          ONLY : lsda
  !USE westcom,        ONLY : nrowmax, xml_file, xc, alpha_pbe0, amplitude, wf_dyn, &
  !                           nitscf, nite, blHF, btHF, qbox_bisec_wfc_filename
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
  PUBLIC :: load_qbox_wfc
  !, init_qbox_interface, finalize_qbox_interface, apply_kernel_by_qbox !, add_debug_log
  !
  ! parameters for doing qbox calculations
  !
  CHARACTER(LEN=256) :: ready_for_kernel = ''
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
  INTEGER             :: np                      ! global size of vextr and drhor array
  INTEGER             :: nploc                   ! size of vextr and drhor array on current processor
  INTEGER,ALLOCATABLE :: nplocs(:)               ! and other processors, determined by dffts
  INTEGER             :: lastproc                ! index of last processor to have FFT grid
  !
  INTEGER             :: nsec_max      = 7200    ! maximum waiting time for one qbox response calculation
  INTEGER             :: iu = 100                ! iu for all serial I/O
  !
  CHARACTER(len=10)   :: mype
  !
  !
  !INTERFACE init_qbox_interface
  !MODULE PROCEDURE init_qbox, init_qbox_for_kernel
  !END INTERFACE
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
!    SUBROUTINE init_qbox()
!      !----------------------------------------------------------------------------
!      !
!      INTEGER                         ::  proc, ierr, numsp
!      !
!      ! initialize Qbox to load wavefunction, set parameters, etc.
!      !
!      WRITE(mype, '(I0.6)') dffts%mype
!      !
!      CALL mp_barrier(intra_image_comm)
!      !
!      ! first time to run init_qbox_interface
!      !
!      CALL io_push_title("Qbox interface setup")
!      !
!!#ifndef C_BINDING
!!      IF ( DP /= C_DOUBLE ) call errore( 'qbox_interface', 'DP /= C_DOUBLE but C_BINDING is not used', 1 )
!!#endif
!      !
!      WRITE(path,'(I6.6)') my_image_id
!      path = 'I' // TRIM(path) // '/'
!      server_input_file       = 'qb.in'
!      lock_file = TRIM(server_input_file) // '.lock'
!      server_output_file      = 'qb.out'
!      vext_file               = 'vext.dat'
!      resp_file = TRIM(vext_file) // '.response'
!      !
!      numsp = 30
!      CALL io_push_value("Qbox input file", TRIM(server_input_file), numsp)
!      CALL io_push_value("Qbox output file", TRIM(server_output_file), numsp)
!      CALL io_push_value("potential file", TRIM(vext_file), numsp)
!      CALL io_push_value("response file", TRIM(resp_file), numsp)
!      !
!      ! initialize qbox
!      !
!      IF(me_image == 0) THEN
!         OPEN(UNIT=iu, FILE=path//TRIM(server_input_file))
!         IF ( nrowmax /= 0 ) WRITE(iu,*) 'set nrowmax ', nrowmax
!         WRITE(iu,'(a)') 'load ../' // TRIM(xml_file)
!         WRITE(iu,'(a)') ' '
!         WRITE(iu,'(a)') 'set xc ' // TRIM( xc )
!         IF ( TRIM(xc) == 'PBE0' ) WRITE(iu, *) 'set alpha_PBE0 ', alpha_pbe0
!         WRITE(iu,'(a)') 'set wf_dyn ' // TRIM( wf_dyn )
!         WRITE(iu,'(a)') 'set blHF ' // TRIM( blHF )
!         WRITE(iu, *) 'set btHF ', btHF
!         CLOSE(iu)
!      ENDIF
!      !
!      CALL wait_for_lock_file()  ! wait for qbox to start
!      !
!      CALL delete_lock_file()    ! now qbox should start initialization
!      !
!      CALL io_push_title("Initializing Qbox")
!      !
!      CALL wait_for_lock_file()
!      !
!      ! archive input and output files for initialization
!      !
!      IF(me_image == 0) THEN
!         ierr = f_copy( path // TRIM( server_input_file), &
!                      & path // (TRIM( server_input_file ) ) // '.init')
!         IF ( ierr /= 0 ) CALL errore("qbox_interface", "fail to archive qbox input")
!      ENDIF
!      CALL archive_qbox_output()
!      !
!      CALL mp_barrier(intra_image_comm)
!      !
!    END SUBROUTINE
!    !
!    !
!    !----------------------------------------------------------------------------
!    SUBROUTINE init_qbox_for_kernel( kernel )
!      !----------------------------------------------------------------------------
!      !
!      CHARACTER(LEN=256), INTENT(IN)    ::  kernel
!      INTEGER                           ::  proc, ierr, numsp
!      CHARACTER(LEN=256)                ::  tmp
!      !
!      ! generate response command for given kernel
!      !
!      IF ( TRIM(ready_for_kernel) == TRIM(kernel) ) RETURN
!      !
!      resp_command = 'response -vext ' // TRIM( vext_file )
!      !
!      ! response mode
!      !
!      SELECT CASE( TRIM(kernel) )
!      CASE('CHI')
!      CASE('CHI_RPA')
!         resp_command = TRIM( resp_command ) // ' -RPA'
!      CASE('CHI0')
!         resp_command = TRIM( resp_command ) // ' -IPA'
!      CASE DEFAULT
!         CALL errore("init_qbox_for_kernel", "qbox_interface: wrong kernel: " // TRIM(kernel), 1)
!      END SELECT
!      !
!      ! amplitude
!      !
!      WRITE(tmp,'(E12.5)') amplitude
!      resp_command = TRIM( resp_command ) // ' -amplitude ' // TRIM(tmp)
!      !
!      ! nitscf and nite
!      !
!      WRITE(tmp,'(I5)') nitscf
!      resp_command = TRIM( resp_command ) // ' ' // TRIM(tmp)
!      WRITE(tmp,'(I5)') nite
!      resp_command = TRIM( resp_command ) // ' ' // TRIM(tmp)
!      !
!      ! write response command to file
!      !
!      IF(me_image == 0) THEN
!         !
!         OPEN(UNIT=iu, FILE=path//TRIM(server_input_file))
!         WRITE(iu,'(A)') TRIM(resp_command)
!         CLOSE(iu)
!         !
!      ENDIF
!      !
!      ready_for_kernel = kernel
!      !
!      CALL mp_barrier(intra_image_comm)
!      !
!    END SUBROUTINE
!    !
!    !
!    !----------------------------------------------------------------------------
!    SUBROUTINE apply_kernel_by_qbox(kernel, fr)
!      !----------------------------------------------------------------------------
!      !
!      ! Applying kernel to real space function fr using response command of Qbox
!      !
!      ! kernel: "CHI", "CHI_RPA" or "CHI0"
!      ! fr: real space function distributed according to dffts
!      !
!      CHARACTER(LEN=256), INTENT(IN)  ::  kernel
!      REAL(DP), INTENT(INOUT)         ::  fr(dffts%nnr)
!      REAL(DP), ALLOCATABLE           ::  frspin(:,:)
!      !
!      ! Generate response command to qbox
!      !
!      CALL init_qbox_interface(kernel)
!      !
!      fr = 0.5_DP * fr  ! change from rydberg to hartree
!      !
!      ! Write fr --> v.dat
!      !
!      CALL write_vext_file(fr)
!      !
!      ! Remove lock file, let qbox resume running
!      !
!      CALL delete_lock_file()
!      !
!      ! Wait for lock file to appear again, archive qbox output
!      !
!      CALL start_clock( 'wait_qbox' )
!      CALL wait_for_lock_file()
!      CALL stop_clock( 'wait_qbox' )
!      CALL archive_qbox_output()
!      !
!      ! Read v.dat.response --> fr
!      !
!      IF (lsda) THEN
!         !
!         ALLOCATE(frspin(dffts%nnr,2))
!         !
!         CALL read_resp_file(frspin(:,1), 0)
!         CALL read_resp_file(frspin(:,2), 1)
!         !
!         fr(:) = frspin(:,1) + frspin(:,2)
!         !
!         DEALLOCATE(frspin)
!         !
!      ELSE
!         !
!         CALL read_resp_file(fr)
!         !
!      ENDIF
!      !
!    END SUBROUTINE
!    !
!    !
!    !----------------------------------------------------------------------------
!    SUBROUTINE write_vext_file(vextr)
!      !----------------------------------------------------------------------------
!      !
!      REAL(DP), INTENT(OUT)           :: vextr(:)
!      !
!      CALL mp_barrier(intra_image_comm)
!      !
!      CALL start_clock( 'write_vext' )
!      !
!      CALL write_function3d(path//TRIM(vext_file),vextr,dffts)
!      !
!      CALL mp_barrier(intra_image_comm)
!      !
!      CALL stop_clock( 'write_vext' )
!      !
!    END SUBROUTINE
!    !
!    !----------------------------------------------------------------------------
!    SUBROUTINE read_resp_file(drhor, ispin)
!      !----------------------------------------------------------------------------
!      !
!      INTEGER, OPTIONAL, INTENT(IN)  :: ispin
!      REAL(DP), INTENT(OUT)  :: drhor(:)
!      !
!      CHARACTER(len=10)      :: mype
!      CHARACTER(LEN=256)     :: resp_file_to_read
!      !
!      IF (PRESENT(ispin)) THEN
!         !
!         IF (.NOT. lsda) CALL errore("present ispin, but not lsda calculation ", 1)
!         !
!         SELECT CASE (ispin)
!         !
!         CASE(0)
!             !
!             resp_file_to_read = path//TRIM(resp_file)//".spin0"
!             !
!         CASE(1)
!             !
!             resp_file_to_read = path//TRIM(resp_file)//".spin1"
!             !
!         CASE DEFAULT
!             !
!             CALL errore("present wrong ispin, ", 1)
!             !
!         END SELECT
!         !
!      ELSE
!         !
!         resp_file_to_read = path//TRIM(resp_file)
!         !
!      ENDIF
!      !
!      WRITE(mype, '(I0.6)') dffts%mype
!      !
!      CALL mp_barrier(intra_image_comm)
!      !
!      CALL start_clock( 'read_resp' )
!      !
!      CALL read_function3d(resp_file_to_read,drhor,dffts)
!      !
!      CALL mp_barrier(intra_image_comm)
!      !
!      CALL stop_clock( 'read_resp' )
!      !
!    END SUBROUTINE
!    !
!    !----------------------------------------------------------------------------
!    SUBROUTINE delete_lock_file()
!      !----------------------------------------------------------------------------
!      !
!      IF(me_image == 0) THEN
!         CALL delete_if_present( path // TRIM( lock_file ) )
!      ENDIF
!      !
!    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
!    SUBROUTINE wait_for_lock_file()
!      !----------------------------------------------------------------------------
!      !
!      ! wait for qbox to finish (indicated by existance of a lock file)
!      !
!      USE forpy_mod,  ONLY: call_py, call_py_noret, import_py, module_py
!      USE forpy_mod,  ONLY: tuple, tuple_create
!      USE forpy_mod,  ONLY: dict, dict_create
!      USE forpy_mod,  ONLY: list, list_create
!      USE forpy_mod,  ONLY: object, cast
!      USE forpy_mod,  ONLY: exception_matches, KeyError, err_clear, err_print
!      !
!      LOGICAL              :: file_exists = .FALSE.
!      !LOGICAL(KIND=C_BOOL) :: file_existss = .FALSE.
!      !
!      INTEGER :: IERR
!      TYPE(tuple) :: args
!      TYPE(dict) :: kwargs
!      TYPE(module_py) :: pymod
!      TYPE(object) :: return_obj
!      INTEGER :: return_int
!      !
!      IF( me_image == 0 ) THEN
!        !
!        IERR = import_py(pymod, "west_clientserver")
!        !
!        IERR = tuple_create(args, 1)
!        IERR = args%setitem(0, TRIM(path // TRIM( lock_file )) )
!        !IERR = dict_create(kwargs)
!        !
!        IERR = call_py(return_obj, pymod, "wait_for_lockfile", args)
!        !
!        IERR = cast(return_int, return_obj)
!        !
!        IF( return_int /= 0 ) CALL errore("sleep","wait for lockfile error",return_int)
!        !
!        !CALL kwargs%destroy
!        CALL args%destroy
!        CALL return_obj%destroy
!        CALL pymod%destroy
!        !CALL c_wait_for_file(INT(nsec_max, KIND=C_INT), file_existss, &
!        !                     & path // TRIM( lock_file ) // C_NULL_CHAR)
!        !
!        ! DO isec = 1,nsec_max
!        !    !
!        !    INQUIRE( FILE= path // TRIM( lock_file ) , EXIST = file_exists )
!        !    IF( file_exists ) THEN
!        !       EXIT
!        !    ELSE
!        !       CALL c_sleep(1000000_C_INT)
!        !    ENDIF
!        !    !
!        ! ENDDO
!        !
!        file_exists = (return_int == 0)
!        !file_exists = file_existss
!        !
!      ENDIF
!      !
!      CALL mp_bcast( file_exists, 0, intra_image_comm )
!      IF ( .NOT. file_exists ) CALL errore( 'qbox_interface', 'Lock file timeout', 1 )
!      !
!    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
!    SUBROUTINE archive_qbox_output()
!      !----------------------------------------------------------------------------
!      !
!      ! rename qbox output file to prevent it from being overwritten by next qbox run
!      !
!      LOGICAL             :: file_exists
!      CHARACTER(LEN=256)  :: new_server_output_file
!      CHARACTER (len=6)   :: file_index
!      INTEGER             :: i, ierr
!      !
!      IF( me_image == 0 ) THEN
!         !
!         file_exists = .FALSE.
!         INQUIRE( FILE= path // TRIM( server_output_file ) , EXIST = file_exists )
!         IF( .NOT. file_exists ) CALL errore( 'qbox_interface', 'file to achive not exist', 1 )
!         !
!         i = 0
!         DO
!            WRITE(file_index,'(i6.6)') i
!            new_server_output_file = TRIM(server_output_file) // "." // file_index
!            INQUIRE( FILE= path // TRIM( new_server_output_file ) , EXIST = file_exists )
!            IF( file_exists ) THEN
!               i = i + 1
!            ELSE
!               ierr = f_copy( path // TRIM( server_output_file ), path // TRIM( new_server_output_file ) )
!               IF ( ierr /= 0 ) CALL errore("qbox_interface", "fail to archive qbox output")
!               EXIT
!            ENDIF
!         ENDDO
!         !
!      ENDIF
!      !
!    END SUBROUTINE
    !
    !----------------------------------------------------------------------------
!    SUBROUTINE finalize_qbox_interface()
!      !----------------------------------------------------------------------------
!      !
!      ! send quit command to qbox
!      !
!      IF(me_image == 0) THEN
!         OPEN(UNIT=iu, FILE=path//TRIM(server_input_file))
!         WRITE(iu,'(A)') 'quit'
!         CLOSE(iu)
!      ENDIF
!      !
!      CALL delete_lock_file()
!      !
!      !DEALLOCATE(nplocs)
!      !
!    END SUBROUTINE
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
  !SUBROUTINE add_debug_log(message)
  !  !----------------------------------------------------------------------------
  !  !
  !  CHARACTER(LEN=*), INTENT(IN) ::  message
  !  CHARACTER(LEN=256)           ::  logfile
  !  LOGICAL                      ::  file_exists = .FALSE.
  !  INTEGER                      ::  iu_dbg = 101
  !  !
  !  logfile = path//'dbg_west_pe'//mype
  !  INQUIRE( FILE=logfile, EXIST=file_exists )
  !  IF ( file_exists ) THEN
  !     OPEN(UNIT=iu_dbg, FILE=logfile, STATUS="OLD", POSITION="APPEND")
  !  ELSE
  !     OPEN(UNIT=iu_dbg, FILE=logfile)
  !  ENDIF
  !  !
  !  WRITE(iu_dbg, *) TRIM(message)
  !  CLOSE(iu_dbg)
  !  CALL c_sleep(100000_C_INT)
  !  !
  !END SUBROUTINE
END MODULE
