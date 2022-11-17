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
  SUBROUTINE load_qbox_wfc(evc,qbox_wfc_fname,current_spin,nbndval)
    !----------------------------------------------------------------------------
    !
    USE kinds,                  ONLY : DP
    USE mp,                     ONLY : mp_bcast
    USE mp_global,              ONLY : intra_image_comm
    USE io_global,              ONLY : stdout,ionode
    USE pwcom,                  ONLY : npw,npwx
    USE fourier_interpolation,  ONLY : ft_interpolate
    USE io_push,                ONLY : io_push_title
    USE forpy_mod,              ONLY : call_py,import_py,module_py,tuple,tuple_create,dict,&
                                     & dict_create,list,object,cast
    USE base64_module,          ONLY : lenbase64,base64_byteswap_double,base64_decode_double,&
                                     & islittleendian
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: current_spin,nbndval
    CHARACTER(LEN=*),INTENT(IN) :: qbox_wfc_fname
    COMPLEX(DP),INTENT(OUT) :: evc(npwx,nbndval)
    !
    INTEGER :: iwfc,nwfcs,nspin,qbox_ispin,ispin
    INTEGER :: nx,ny,nz
    CHARACTER(LEN=20) :: namespin
    REAL(DP),ALLOCATABLE :: psir(:)
    CHARACTER(LEN=:),ALLOCATABLE :: charbase64
    INTEGER :: ndim,nbytes,nlen,ierr
    TYPE(tuple) :: args
    TYPE(dict) :: kwargs,return_dict
    TYPE(list) :: tmp_list
    TYPE(module_py) :: pymod
    TYPE(object) :: return_obj,tmp_obj
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
       WRITE(stdout,'(5X,2A)') 'Qbox interface: reading from file ',TRIM(qbox_wfc_fname)
       !
       ierr = import_py(pymod,'west_function3d')
       ierr = tuple_create(args,1)
       ierr = args%setitem(0,TRIM(qbox_wfc_fname))
       ierr = dict_create(kwargs)
       ierr = call_py(return_obj,pymod,'qb_wfc_info',args,kwargs)
       ierr = cast(return_dict,return_obj)
       ierr = return_dict%get(nspin,'nspin',DUMMY_DEFAULT)
       ierr = return_dict%get(nwfcs,'nwfcs',DUMMY_DEFAULT)
       ierr = return_dict%getitem(tmp_obj,'grid')
       ierr = cast(tmp_list,tmp_obj)
       ierr = tmp_list%getitem(nx,0)
       ierr = tmp_list%getitem(ny,1)
       ierr = tmp_list%getitem(nz,2)
       !
       IF(nspin == DUMMY_DEFAULT) CALL errore('load_qbox_wfc','cannot read nspin',1)
       IF(nwfcs == DUMMY_DEFAULT) CALL errore('load_qbox_wfc','cannot read nwfcs',1)
       !
       CALL args%destroy
       CALL kwargs%destroy
       CALL return_obj%destroy
       CALL return_dict%destroy
       CALL tmp_list%destroy
       CALL tmp_obj%destroy
       !
       WRITE(stdout,'(5X,A,3I5)') 'Qbox interface: grid size: ',nx,ny,nz
       !
       ndim = nx*ny*nz
       ALLOCATE(psir(ndim))
       nbytes = SIZEOF(psir(1)) * ndim
       nlen = lenbase64(nbytes)
       ALLOCATE(CHARACTER(LEN=nlen) :: charbase64)
       !
    ENDIF
    !
    CALL mp_bcast(nspin,0,intra_image_comm)
    CALL mp_bcast(nwfcs,0,intra_image_comm)
    CALL mp_bcast(nx,0,intra_image_comm)
    CALL mp_bcast(ny,0,intra_image_comm)
    CALL mp_bcast(nz,0,intra_image_comm)
    !
    ! read KS orbitals
    !
    WRITE(stdout,'(5X,A)') 'Qbox interface: loading KS orbitals'
    !
    DO ispin = 1,nspin
       !
       IF(nspin > 1) THEN
          CALL errore('load_qbox_wfc','only nspin = 1 supported',nspin)
       ELSE
          namespin = 'none'
          qbox_ispin = 1
       ENDIF
       !
       IF(qbox_ispin == current_spin) THEN
          !
          WRITE(stdout,'(10X,A,I5,2I3)') namespin,nwfcs,nspin,qbox_ispin
          !
          DO iwfc = 1,nwfcs
             !
             IF(ionode) THEN
                !
                ierr = tuple_create(args,2)
                ierr = args%setitem(0,TRIM(qbox_wfc_fname))
                ierr = args%setitem(1,iwfc)
                ierr = dict_create(kwargs)
                ierr = call_py(return_obj,pymod,'qb_wfc_to_base64',args,kwargs)
                ierr = cast(return_dict,return_obj)
                ierr = return_dict%getitem(charbase64,'grid_function')
                !
                CALL args%destroy
                CALL kwargs%destroy
                CALL return_obj%destroy
                CALL return_dict%destroy
                !
                CALL base64_decode_double(charbase64(1:nlen),ndim,psir)
                IF(.NOT. islittleendian()) CALL base64_byteswap_double(nbytes,psir)
                !
             ENDIF
             !
             CALL ft_interpolate(psir,nx,ny,nz,evc(:,iwfc),npw,npwx)
             !
          ENDDO
          !
       ENDIF
       !
    ENDDO
    !
    IF(ionode) THEN
       !
       CALL pymod%destroy
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
