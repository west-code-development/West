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
!----------------------------------------------------------------------------
SUBROUTINE do_wfc2 ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE pwcom,                 ONLY : igk_k,npw,npwx,current_k,ngk
  USE io_push,               ONLY : io_push_title
  USE westcom,               ONLY : westpp_sign,iuwfc,lrwfc,westpp_range,westpp_save_dir
  USE mp_global,             ONLY : inter_image_comm,my_image_id
  USE mp,                    ONLY : mp_bcast
  USE fft_base,              ONLY : dffts
  USE buffers,               ONLY : get_buffer
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : aband
  USE class_idistribute,     ONLY : idistribute
  USE control_flags,         ONLY : gamma_only
  USE types_bz_grid,         ONLY : k_grid
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
  USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ir,iks,local_ib,global_ib,dffts_nnr
  REAL(DP),ALLOCATABLE :: auxr(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  CHARACTER(LEN=6) :: labelb,labelk
  !
  aband = idistribute()
  CALL aband%init(westpp_range(2)-westpp_range(1)+1,'i','westpp_range',.TRUE.)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  ALLOCATE(auxr(dffts%nnr))
  !$acc enter data create(auxr)
  !
  dffts_nnr = dffts%nnr
  !
  CALL io_push_title('(W)avefunctions')
  !
  CALL start_bar_type( barra, 'westpp', k_grid%nps * MAX(aband%nloc,1) )
  !
  DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(k_grid%nps > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     DO local_ib=1,aband%nloc
        !
        ! local -> global
        !
        global_ib = aband%l2g(local_ib)+westpp_range(1)-1
        !
        IF( gamma_only ) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,global_ib),psic,'Wave')
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc_work(:,global_ib),psic,'Wave',igk_k(:,current_k))
        ENDIF
        IF( westpp_sign ) THEN
           !$acc parallel loop present(auxr)
           DO ir = 1, dffts_nnr
              auxr(ir) = REAL( psic(ir), KIND=DP) *  ABS( REAL( psic(ir), KIND=DP) )
           ENDDO
           !$acc end parallel loop
        ELSE
           !$acc parallel loop present(auxr)
           DO ir = 1, dffts_nnr
              auxr(ir) = REAL( psic(ir), KIND=DP) *  REAL( psic(ir), KIND=DP)
           ENDDO
           !$acc end parallel loop
        ENDIF
        !
        WRITE(labelb,'(i6.6)') global_ib
        WRITE(labelk,'(i6.6)') iks
        fname = TRIM( westpp_save_dir ) // '/wfcK'//labelk//'B'//labelb
        !$acc update host(auxr)
        CALL dump_r( auxr, fname )
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  !$acc exit data delete(auxr)
  DEALLOCATE( auxr )
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
