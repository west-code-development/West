!
! Copyright (C) 2015-2022 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file:
! He Ma
!
!-----------------------------------------------------------------------
MODULE qbox_interface
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! currently only consider nks = nbgrp = 1 case
  !
  PRIVATE
  !
  PUBLIC :: load_qbox_wfc
  PUBLIC :: init_qbox
  PUBLIC :: finalize_qbox
  PUBLIC :: sleep_and_wait_for_lock_to_be_removed
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE load_qbox_wfc(evc, qbox_wfc_filename, current_spin, nbndval)
    !----------------------------------------------------------------------------
    !
    USE iotk_module
    USE kinds,                  ONLY : DP
    USE mp,                     ONLY : mp_bcast
    USE mp_global,              ONLY : intra_image_comm
    USE io_global,              ONLY : stdout,ionode
    USE pwcom,                  ONLY : npw,npwx
    USE fourier_interpolation,  ONLY : ft_interpolate
    USE io_push,                ONLY : io_push_title
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: current_spin,nbndval
    CHARACTER(LEN=*),INTENT(IN) :: qbox_wfc_filename
    COMPLEX(DP),INTENT(OUT) :: evc(npwx,nbndval)
    !
    INTEGER :: iwfc, nwfcs, nspin, ispin, is, ir
    INTEGER :: nx, ny, nz
    INTEGER :: iun
    CHARACTER(iotk_attlenx) :: attr
    CHARACTER(LEN=20) :: namespin
    REAL(DP),ALLOCATABLE :: psir(:)
    CHARACTER(iotk_attlenx) :: slen
    !
    ! read Qbox XML file and overwrite current evc array
    !
    CALL io_push_title('Loading Qbox wavefunction...')
    !
    ! root of first image read Qbox XML file, Fourier interpolate to evc of whole image,
    ! then first image bcast evc to other images
    !
    evc(:,:) = (0._DP,0._DP)
    !
    IF(ionode) THEN
       !
       WRITE(stdout,'(5X,2A)') 'Qbox interface: reading from file ', TRIM(qbox_wfc_filename)
       !
       CALL iotk_free_unit(iun)
       CALL iotk_open_read(iun, FILE=TRIM(qbox_wfc_filename))
       !
       CALL iotk_scan_begin(iun, 'wavefunction', attr)
       CALL iotk_scan_attr(attr, 'nspin', nspin)
       !
       ! read grid size
       !
       CALL iotk_scan_empty(iun, 'grid', attr)
       CALL iotk_scan_attr(attr, 'nx', nx)
       CALL iotk_scan_attr(attr, 'ny', ny)
       CALL iotk_scan_attr(attr, 'nz', nz)
       !
       WRITE(stdout,'(5X,A,3I5)') 'Qbox interface: grid size: ', nx, ny, nz
       !
       ALLOCATE(psir(nx*ny*nz))
       !
    ENDIF
    !
    CALL mp_bcast(nx, 0, intra_image_comm)
    CALL mp_bcast(ny, 0, intra_image_comm)
    CALL mp_bcast(nz, 0, intra_image_comm)
    CALL mp_bcast(nspin, 0, intra_image_comm)
    !
    ! read KS orbitals
    !
    WRITE(stdout,'(5X,A)') 'Qbox interface: loading KS orbitals'
    !
    DO is = 1,nspin
       !
       IF(nspin > 1) THEN
          IF(ionode) THEN
             CALL iotk_scan_begin(iun, 'slater_determinant', attr)
             CALL iotk_scan_attr(attr, 'spin', namespin)
             CALL iotk_scan_attr(attr, 'size', nwfcs)
          ENDIF
          !
          CALL mp_bcast(namespin, 0, intra_image_comm)
          CALL mp_bcast(nwfcs, 0, intra_image_comm)
          !
          IF(namespin == 'up') ispin = 1
          IF(namespin == 'down') ispin = 2
       ELSE
          IF(ionode) THEN
             CALL iotk_scan_begin(iun, 'slater_determinant', attr)
             CALL iotk_scan_attr(attr, 'size', nwfcs)
          ENDIF
          !
          CALL mp_bcast(nwfcs, 0, intra_image_comm)
          !
          namespin = 'none'
          ispin = 1
       ENDIF
       !
       IF(ispin == current_spin) THEN
          !
          WRITE(stdout,'(10X,A,I5,2I3)') namespin, nwfcs, nspin, ispin
          !
          DO iwfc = 1,nwfcs
             !
             IF(ionode) THEN
                CALL iotk_scan_begin(iun, 'grid_function', slen)
                READ(iun,*) (psir(ir), ir = 1, nx*ny*nz)
                CALL iotk_scan_end(iun, 'grid_function')
             ENDIF
             !
             CALL ft_interpolate(psir, nx, ny, nz, evc(:,iwfc), npw, npwx)
             !
          ENDDO
          !
       ENDIF
       !
       IF(ionode) CALL iotk_scan_end(iun, 'slater_determinant')
       !
    ENDDO
    !
    IF(ionode) THEN
       CALL iotk_scan_end(iun, 'wavefunction')
       CALL iotk_close_read(iun)
       !
       DEALLOCATE(psir)
    ENDIF
    !
  END SUBROUTINE
  !
  !----------------------------------------------------------------------------
  SUBROUTINE init_qbox()
    !----------------------------------------------------------------------------
    !
    USE mp,              ONLY : mp_barrier
    USE mp_global,       ONLY : intra_image_comm,my_image_id,me_bgrp
    USE conversions,     ONLY : itoa
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=:),ALLOCATABLE :: lockfile
    INTEGER :: iun
    !
    ! DUMP A LOCK FILE
    !
    IF(me_bgrp == 0) THEN
       lockfile = 'I.'//itoa(my_image_id)//'.lock'
       OPEN(NEWUNIT=iun, FILE=lockfile)
       CLOSE(iun)
       !
       ! SLEEP AND WAIT FOR LOCKFILE TO BE REMOVED
       !
       CALL sleep_and_wait_for_lock_to_be_removed(lockfile, '["script"]')
    ENDIF
    !
    CALL mp_barrier(intra_image_comm)
    !
  END SUBROUTINE
  !
  !----------------------------------------------------------------------------
  SUBROUTINE finalize_qbox()
    !----------------------------------------------------------------------------
    !
    ! send quit command to qbox
    !
    USE mp_global,              ONLY : my_image_id,me_bgrp
    USE west_io,                ONLY : remove_if_present
    USE conversions,            ONLY : itoa
    !
    IMPLICIT NONE
    !
    INTEGER :: iun
    !
    IF(me_bgrp == 0) THEN
       OPEN(NEWUNIT=iun, FILE='qb.'//itoa(my_image_id)//'.in')
       WRITE(iun, '(A)') 'quit'
       CLOSE(iun)
       CALL remove_if_present('qb.'//itoa(my_image_id)//'.in.lock')
    ENDIF
    !
  END SUBROUTINE
  !
  !----------------------------------------------------------------------------
  SUBROUTINE sleep_and_wait_for_lock_to_be_removed(lockfile, consider_only)
    !----------------------------------------------------------------------------
    !
    USE westcom,                ONLY : document
    USE forpy_mod,              ONLY : call_py,import_py,module_py,tuple,tuple_create,dict,&
                                     & dict_create,object,cast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),INTENT(IN) :: lockfile
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: consider_only
    !
    INTEGER :: ierr
    TYPE(tuple) :: args
    TYPE(dict) :: kwargs
    TYPE(module_py) :: pymod
    TYPE(object) :: return_obj
    INTEGER :: return_int
    !
    ierr = import_py(pymod, 'west_clientserver')
    !
    ierr = tuple_create(args, 1)
    ierr = args%setitem(0, TRIM(ADJUSTL(lockfile)))
    ierr = dict_create(kwargs)
    ierr = kwargs%setitem('document', document)
    !
    IF(PRESENT(consider_only)) THEN
       ierr = kwargs%setitem('consider_only', consider_only)
    ENDIF
    !
    ierr = call_py(return_obj, pymod, 'sleep_and_wait', args, kwargs)
    !
    ierr = cast(return_int, return_obj)
    !
    IF(return_int /= 0) CALL errore('sleep', 'Did not wake well', return_int)
    !
    CALL kwargs%destroy
    CALL args%destroy
    CALL return_obj%destroy
    CALL pymod%destroy
    !
  END SUBROUTINE
  !
END MODULE
