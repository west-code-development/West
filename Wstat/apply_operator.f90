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
!-------------------------------------------------------------------------
SUBROUTINE apply_operator (m,dvg,dng,tr2,iq)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE westcom,               ONLY : npwqx,npwq,wstat_calculation
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
  REAL(DP),INTENT(IN) :: tr2
  INTEGER, INTENT(IN) :: iq
  !
  ! Workspace
  !
  INTEGER :: ipert, ig, i
  !
  COMPLEX(DP), ALLOCATABLE ::aux_g(:,:)
  !
  LOGICAL :: l_outsource
  !
  l_outsource = .FALSE.
  DO i = 1,2
     IF( wstat_calculation(i:i) == 'E' ) l_outsource = .TRUE.
  ENDDO
  !
  dng=(0.0_DP,0.0_DP)
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
     CALL calc_outsourced(m,aux_g,dng,iq)
  ELSE
     CALL dfpt(m,aux_g,dng,tr2,iq)
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
  USE mp_global,       ONLY : intra_image_comm,my_image_id,me_bgrp
  USE fft_at_k,        ONLY : single_fwfft_k,single_invfft_k
  USE fft_at_gamma,    ONLY : single_fwfft_gamma,single_invfft_gamma
  USE fft_base,        ONLY : dffts
  USE control_flags,   ONLY : gamma_only
  USE function3d,      ONLY : write_function3d,read_function3d
  USE conversions,     ONLY : itoa
  USE bar,             ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,         ONLY : io_push_title
  USE qbox_interface,  ONLY : sleep_and_wait_for_lock_to_be_removed
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
#if defined(__CUDA)
  CALL errore("outsourced","GPU not implemented",1)
#endif
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
#if !defined(__CUDA)
        IF (gamma_only) THEN
           CALL single_invfft_gamma(dffts,npwq,npwqx,dvg(:,ipert),aux_r,TRIM(fftdriver))
        ELSE
           CALL single_invfft_k(dffts,npwq,npwqx,dvg(:,ipert),aux_r,'Wave',igq_q(:,iq))
        ENDIF
#endif
        !
        filename = "I."//itoa(my_image_id)//"_P."//itoa(ipert)//".xml"
        aux_r_double(:) = REAL(aux_r(:),KIND=DP) / 2._DP ! The output must be in Ha Atomic units
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
        aux_r(:) = CMPLX(aux_r_double(:),0._DP,KIND=DP)
        !
#if !defined(__CUDA)
        IF(gamma_only) THEN
           CALL single_fwfft_gamma(dffts,npwq,npwqx,aux_r,dng(:,ipert),TRIM(fftdriver))
        ELSE
           CALL single_fwfft_k(dffts,npwq,npwqx,aux_r,dng(:,ipert),'Wave',igq_q(:,iq))
        ENDIF
#endif
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
END SUBROUTINE
