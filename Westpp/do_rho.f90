!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE do_rho ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE pwcom,                 ONLY : igk_k,npw,npwx,current_k,ngk,nbnd
  USE io_push,               ONLY : io_push_title
  USE westcom,               ONLY : iuwfc,lrwfc,westpp_save_dir,nbnd_occ
  USE mp_global,             ONLY : inter_image_comm,my_image_id
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_base,              ONLY : dffts
  USE buffers,               ONLY : get_buffer
  USE wavefunctions,         ONLY : evc,psic
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : aband
  USE class_idistribute,     ONLY : idistribute
  USE control_flags,         ONLY : gamma_only
  USE types_bz_grid,         ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ir, iks, local_ib, global_ib
  REAL(DP),ALLOCATABLE :: auxr(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  !
  CALL io_push_title('(R)ho')
  !
  aband = idistribute()
  CALL aband%init(nbnd,'i','nbnd',.TRUE.)
  !
  ALLOCATE(auxr(dffts%nnr))
  !
  auxr = 0._DP
  psic = 0._DP
  !
  CALL start_bar_type( barra, 'westpp', k_grid%nps )
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
     IF(k_grid%nps>1) THEN
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     DO local_ib=1,aband%nloc
        !
        ! local -> global
        !
        global_ib = aband%l2g(local_ib)
        IF( global_ib > nbnd_occ(iks) ) CYCLE
        !
        IF( gamma_only ) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,global_ib),psic,'Wave')
           DO ir = 1, dffts%nnr
              auxr(ir) = auxr(ir) + REAL( psic(ir), KIND=DP) *  REAL( psic(ir), KIND=DP) * k_grid%weight(iks)
           ENDDO
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(:,global_ib),psic,'Wave',igk_k(:,current_k))
           DO ir = 1, dffts%nnr
              auxr(ir) = auxr(ir) + REAL( CONJG( psic(ir) ) * psic(ir), KIND=DP) * k_grid%weight(iks)
           ENDDO
        ENDIF
        !
     ENDDO
     !
     CALL update_bar_type( barra,'westpp', 1 )
     !
  ENDDO
  !
  CALL mp_sum( auxr, inter_image_comm )
  !
  fname = TRIM( westpp_save_dir ) // '/rho'
  IF(my_image_id==0) CALL dump_r( auxr, fname)
  !
  DEALLOCATE( auxr )
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
END SUBROUTINE
