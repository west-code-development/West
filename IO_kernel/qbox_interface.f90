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
  ! currently only consider nks = nbgrp = 1 case
  !
  IMPLICIT NONE
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE load_qbox_wfc(ispin,nbndval,evc)
    !----------------------------------------------------------------------------
    !
    USE kinds,                  ONLY : DP
    USE mp,                     ONLY : mp_bcast
    USE mp_global,              ONLY : intra_image_comm
    USE io_global,              ONLY : stdout,ionode
    USE pwcom,                  ONLY : npw,npwx
    USE fourier_interpolation,  ONLY : ft_interpolate
    USE io_push,                ONLY : io_push_title
    USE westcom,                ONLY : wfc_from_qbox
    USE forpy_mod,              ONLY : call_py,import_py,module_py,tuple,tuple_create,dict,&
                                     & dict_create,ndarray,object,cast
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: ispin,nbndval
    COMPLEX(DP),INTENT(OUT) :: evc(npwx,nbndval)
    !
    CHARACTER(LEN=256) :: fname
    INTEGER :: ierr
    INTEGER :: iwfc,nwfcs
    INTEGER :: nx,ny,nz
    REAL(DP),POINTER :: psir(:)
    TYPE(module_py) :: pymod
    TYPE(tuple) :: args
    TYPE(dict) :: kwargs,return_dict
    TYPE(ndarray) :: return_ndarray
    TYPE(object) :: return_obj
    INTEGER,PARAMETER :: DUMMY_DEFAULT = -1210
    !
    ! read Qbox XML file and overwrite current evc array
    !
    CALL io_push_title('Loading Qbox wavefunction...')
    !
    CALL start_clock('load_qb')
    !
    ! root of first image read Qbox XML file, Fourier interpolate to evc of whole image,
    ! then first image bcast evc to other images
    !
    evc(:,:) = (0._DP,0._DP)
    !
    IF(ionode) THEN
       !
       WRITE(fname,'(A,I1)') TRIM(wfc_from_qbox)//'.',ispin
       WRITE(stdout,'(5X,2A)') 'Qbox interface: reading from file ',TRIM(fname)
       !
       ierr = import_py(pymod,'west_read_bse')
       ierr = tuple_create(args,1)
       ierr = args%setitem(0,TRIM(fname))
       ierr = dict_create(kwargs)
       ierr = call_py(return_obj,pymod,'bse_read_qb_param',args,kwargs)
       ierr = cast(return_dict,return_obj)
       ierr = return_dict%get(nwfcs,'nwfcs',DUMMY_DEFAULT)
       ierr = return_dict%get(nx,'nx',DUMMY_DEFAULT)
       ierr = return_dict%get(ny,'ny',DUMMY_DEFAULT)
       ierr = return_dict%get(nz,'nz',DUMMY_DEFAULT)
       !
       IF(nwfcs == DUMMY_DEFAULT) CALL errore('load_qbox_wfc','cannot read nwfcs',1)
       IF(nx == DUMMY_DEFAULT) CALL errore('load_qbox_wfc','cannot read nx',1)
       IF(ny == DUMMY_DEFAULT) CALL errore('load_qbox_wfc','cannot read ny',1)
       IF(nz == DUMMY_DEFAULT) CALL errore('load_qbox_wfc','cannot read nz',1)
       !
       CALL args%destroy
       CALL kwargs%destroy
       CALL return_obj%destroy
       CALL return_dict%destroy
       !
       WRITE(stdout,'(5X,A,3I5)') 'Qbox interface: grid size ',nx,ny,nz
       !
       ALLOCATE(psir(nx*ny*nz))
       !
    ELSE
       !
       ! ft_interpolate segfaults without this allocation
       !
       ALLOCATE(psir(1))
       !
    ENDIF
    !
    CALL mp_bcast(nwfcs,0,intra_image_comm)
    CALL mp_bcast(nx,0,intra_image_comm)
    CALL mp_bcast(ny,0,intra_image_comm)
    CALL mp_bcast(nz,0,intra_image_comm)
    !
    ! read KS orbitals
    !
    WRITE(stdout,'(5X,A,I5,A)') 'Qbox interface: loading ',nwfcs,' KS orbitals'
    !
    DO iwfc = 1,nwfcs
       !
       IF(ionode) THEN
          !
          ierr = tuple_create(args,2)
          ierr = args%setitem(0,TRIM(fname))
          ierr = args%setitem(1,iwfc)
          ierr = dict_create(kwargs)
          ierr = call_py(return_obj,pymod,'bse_read_qb_wfc',args,kwargs)
          ierr = cast(return_ndarray,return_obj)
          ierr = return_ndarray%get_data(psir)
          !
          CALL args%destroy
          CALL kwargs%destroy
          CALL return_obj%destroy
          CALL return_ndarray%destroy
          !
       ENDIF
       !
       CALL ft_interpolate(psir,nx,ny,nz,evc(:,iwfc),npw,npwx)
       !
    ENDDO
    !
    IF(ionode) THEN
       !
       CALL pymod%destroy
       !
       NULLIFY(psir)
       !
    ELSE
       !
       DEALLOCATE(psir)
       !
    ENDIF
    !
    CALL stop_clock('load_qb')
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
       OPEN(NEWUNIT=iun,FILE=lockfile)
       CLOSE(iun)
       CALL sleep_and_wait_for_lock_to_be_removed(lockfile,'["script"]')
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
  SUBROUTINE sleep_and_wait_for_lock_to_be_removed(lockfile,consider_only)
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
    ierr = import_py(pymod,'west_clientserver')
    ierr = tuple_create(args,1)
    ierr = args%setitem(0,TRIM(ADJUSTL(lockfile)))
    ierr = dict_create(kwargs)
    ierr = kwargs%setitem('document',document)
    !
    IF(PRESENT(consider_only)) THEN
       ierr = kwargs%setitem('consider_only',consider_only)
    ENDIF
    !
    ierr = call_py(return_obj,pymod,'sleep_and_wait',args,kwargs)
    ierr = cast(return_int,return_obj)
    !
    IF(return_int /= 0) CALL errore('sleep','Did not wake well',return_int)
    !
    CALL kwargs%destroy
    CALL args%destroy
    CALL return_obj%destroy
    CALL pymod%destroy
    !
  END SUBROUTINE
  !
END MODULE
