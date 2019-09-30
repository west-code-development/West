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
!-------------------------------------------------------------------------
SUBROUTINE apply_operator (m,dvg,dng,tr2,iq)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE mp,                    ONLY : mp_barrier
  USE westcom,               ONLY : npwqx,npwq,wstat_calculation
  USE mp_world,              ONLY : world_comm
  USE types_coulomb,         ONLY : pot3D
  USE dfpt_module,           ONLY : dfpt
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m
  COMPLEX(DP), INTENT(IN) :: dvg(npwqx,m)
  COMPLEX(DP), INTENT(OUT) :: dng(npwqx,m)
  REAL(DP),INTENT(IN), OPTIONAL :: tr2
  INTEGER, INTENT(IN), OPTIONAL :: iq
  !
  ! Workspace
  !
  INTEGER :: ipert, ig, i
  INTEGER :: iq_
  !
  REAL(DP) :: tr2_
  !
  COMPLEX(DP), ALLOCATABLE ::aux_g(:,:)
  !
  LOGICAL :: l_outsource
  !
  CALL mp_barrier( world_comm )
  !
  l_outsource = .FALSE.
  DO i = 1,2 
     IF( wstat_calculation(i:i) == 'E' ) l_outsource = .TRUE.
  ENDDO
  !
  IF(PRESENT(iq)) THEN 
     iq_ = iq
  ELSE
     iq_ = 1
  ENDIF
  IF(PRESENT(tr2)) THEN 
     tr2_ = tr2
  ELSE
     tr2_ = 1.d-8
  ENDIF
  !
  dng=0._DP
  !
  ALLOCATE( aux_g(npwqx,m) ); aux_g=0._DP
  !
  DO ipert = 1, m
     DO ig = 1, npwq
        aux_g(ig,ipert) = dvg(ig,ipert) * pot3D%sqvc(ig) ! perturbation acts only on body 
     ENDDO
  ENDDO
  !
  IF( l_outsource ) THEN
     CALL calc_outsourced(m,aux_g,dng,iq_)
  ELSE 
     CALL dfpt(m,aux_g,dng,tr2_,iq_)
  ENDIF
  !
  DEALLOCATE( aux_g ) 
  !
  DO ipert = 1, m
     DO ig = 1, npwq
        dng(ig,ipert) = dng(ig,ipert) * pot3D%sqvc(ig) ! perturbation acts only on body  
     ENDDO
  ENDDO
  !
  CALL mp_barrier( world_comm )
  !
END SUBROUTINE
!
!
!-------------------------------------------------------------------------
SUBROUTINE calc_outsourced (m,dvg,dng,iq)
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE mp,              ONLY : mp_barrier
  USE westcom,         ONLY : npwq,npwqx,fftdriver,igq_q      
  USE mp_global,       ONLY : intra_image_comm,inter_pool_comm,my_image_id,me_bgrp
  USE fft_at_k,        ONLY : single_fwfft_k,single_invfft_k
  USE fft_at_gamma,    ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_base,        ONLY : dfftp,dffts
  USE control_flags,   ONLY : gamma_only
  USE function3d,      ONLY : write_function3d,read_function3d
  USE conversions,     ONLY : itoa
  USE bar,             ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,         ONLY : io_push_title
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: iq
  INTEGER, INTENT(IN) :: m
  COMPLEX(DP), INTENT(IN) :: dvg(npwqx,m)
  COMPLEX(DP), INTENT(OUT) :: dng(npwqx,m)
  COMPLEX(DP), ALLOCATABLE :: aux_r(:)
  REAL(DP), ALLOCATABLE :: aux_r_double(:)
  CHARACTER(LEN=:),ALLOCATABLE :: filename
  CHARACTER(LEN=:),ALLOCATABLE :: lockfile
  !
  INTEGER :: ipert, iu, stat 
  TYPE(bar_type) :: barra
  !
  IF(iq/=1) CALL errore("outsourced","iq /= 1 not allowed",iq)
  !
  CALL io_push_title("Calculation outsourced")
  !
  CALL start_bar_type( barra, 'outsourced', MAX(m,1) )
  !
  IF (m>0) THEN 
     !
     ALLOCATE(aux_r(dffts%nnr)); aux_r=0._DP
     ALLOCATE(aux_r_double(dffts%nnr)); aux_r=0._DP
     !
     ! WRITE PERTURBATIONS TO FILE
     !
     DO ipert = 1, m
        !
        IF (gamma_only) THEN
          CALL single_invfft_gamma(dffts,npwq,npwqx,dvg(:,ipert),aux_r,TRIM(fftdriver))
        ELSE
           CALL single_invfft_k(dffts,npwq,npwqx,dvg(:,ipert),aux_r,'Wave',igq_q(1,iq))
        ENDIF
        !
        filename = "I."//itoa(my_image_id)//"_P."//itoa(ipert)//".xml"
        aux_r_double(:) = DBLE(aux_r(:)) / 2._DP ! The output must be in Ha Atomic units 
        CALL write_function3d(filename,aux_r_double,dffts)
        !
     ENDDO
     !
     ! DUMP A LOCK FILE 
     !
     IF( me_bgrp == 0 ) THEN  
        lockfile = "I."//itoa(my_image_id)//".lock"
        OPEN(NEWUNIT=iu,FILE=lockfile) 
        DO ipert = 1, m 
           filename = "I."//itoa(my_image_id)//"_P."//itoa(ipert)//".xml"
           WRITE(iu,'(A)') filename
        ENDDO
        CLOSE(iu)
        !
        ! SLEEP AND WAIT FOR LOCKFILE TO BE REMOVED 
        !
        CALL sleep_and_wait_for_lock_to_be_removed(lockfile)
        !
     ENDIF
     !
     CALL mp_barrier(intra_image_comm)
     !
     ! READ RESPONSES
     !
     DO ipert = 1, m
        !
        filename = "I."//itoa(my_image_id)//"_P."//itoa(ipert)//".xml.response"
        CALL read_function3d(filename,aux_r_double,dffts)
        aux_r(:) = CMPLX(aux_r_double(:),0._DP) 
        !       
        IF(gamma_only) THEN
           CALL single_fwfft_gamma(dffts,npwq,npwqx,aux_r,dng(:,ipert),TRIM(fftdriver))
        ELSE
           CALL single_fwfft_k(dffts,npwq,npwqx,aux_r,dng(:,ipert),'Wave',igq_q(1,iq))
        ENDIF
        !
     ENDDO
     !
     ! CLEANUP
     !
     IF( me_bgrp == 0 ) THEN  
        DO ipert = 1, m
           filename = "I."//itoa(my_image_id)//"_P."//itoa(ipert)//".xml" 
           OPEN(NEWUNIT=iu, IOSTAT=stat, FILE=filename, STATUS='OLD')
           IF (stat == 0) CLOSE(iu, STATUS='DELETE')
           filename = "I."//itoa(my_image_id)//"_P."//itoa(ipert)//".xml.response"
           OPEN(NEWUNIT=iu, IOSTAT=stat, FILE=filename, STATUS='OLD')
           IF (stat == 0) CLOSE(iu, STATUS='DELETE')
        ENDDO
     ENDIF
     CALL update_bar_type( barra, 'outsourced', m ) 
     DEALLOCATE(aux_r)
     DEALLOCATE(aux_r_double)
  ELSE 
     CALL update_bar_type( barra, 'outsourced', 1 )
  ENDIF 
  CALL stop_bar_type( barra, 'outsourced' )
  !
  !
END SUBROUTINE


SUBROUTINE sleep_and_wait_for_lock_to_be_removed(lockfile)
    !
    USE forpy_mod,  ONLY: call_py, call_py_noret, import_py, module_py
    USE forpy_mod,  ONLY: tuple, tuple_create 
    USE forpy_mod,  ONLY: dict, dict_create 
    USE forpy_mod,  ONLY: list, list_create 
    USE forpy_mod,  ONLY: object, cast
    USE forpy_mod,  ONLY: exception_matches, KeyError, err_clear, err_print 
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),INTENT(IN) :: lockfile
    !
    INTEGER :: IERR
    TYPE(tuple) :: args
    TYPE(dict) :: kwargs
    TYPE(module_py) :: pymod
    TYPE(object) :: return_obj
    INTEGER :: return_int
    !
    IERR = import_py(pymod, "west_clientserver")
    !  
    IERR = tuple_create(args, 1)
    IERR = args%setitem(0, TRIM(ADJUSTL(lockfile)) )
    IERR = dict_create(kwargs)
    !
    IERR = call_py(return_obj, pymod, "sleep_and_wait", args, kwargs)
    !
    IERR = cast(return_int, return_obj)
    !
    IF( return_int /= 0 ) CALL errore("sleep","Did not wake well",return_int)
    !
    CALL kwargs%destroy
    CALL args%destroy
    CALL return_obj%destroy
    CALL pymod%destroy 
    !
END SUBROUTINE
