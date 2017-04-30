!
! Copyright (C) 2015-2016 M. Govoni 
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
SUBROUTINE do_eigenpot2 ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE uspp,                  ONLY : vkb,nkb
  USE io_global,             ONLY : stdout
  USE pwcom,                 ONLY : current_spin,wk,nks,nelup,neldw,isk,g,igk_k,ngm,tpiba2,xk,npw,npwx,lsda,nkstot,&
                                  & current_k,ngk
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE westcom,               ONLY : westpp_sign,iuwfc,lrwfc,westpp_calculation,westpp_range,westpp_dirname,fftdriver,&
                                  & npwq0,npwq0x,dvg 
  USE mp_global,             ONLY : inter_image_comm,my_image_id
  USE mp,                    ONLY : mp_bcast
  USE fft_base,              ONLY : dfftp,dffts
  USE wvfct,                 ONLY : nbnd
  USE buffers,               ONLY : get_buffer
  USE wavefunctions_module,  ONLY : evc,psic
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : pert
  USE control_flags,         ONLY : gamma_only 
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: i1,i2, ipol, ir, local_j, global_j, i, ig, iks, ibnd, local_ib, global_ib
  REAL(DP),ALLOCATABLE :: auxr(:)
  CHARACTER(LEN=256)    :: fname
  TYPE(bar_type) :: barra
  CHARACTER(LEN=6) :: labelq,labeli
  !
  CALL io_push_title("(E)igenpotentials")
  !
  ALLOCATE(auxr(dffts%nnr))
  !
  auxr = 0._DP
  psic = 0._DP
  !
  CALL start_bar_type( barra, 'westpp', pert%nloc ) 
  !
  DO local_j=1,pert%nloc
     !
     ! local -> global
     !
     global_j = pert%l2g(local_j)
     IF( global_j < westpp_range(1) .OR. global_j > westpp_range(2) ) CYCLE
     !
     auxr = 0._DP
     IF( gamma_only ) THEN 
        CALL single_invfft_gamma(dffts,npwq0,npwq0x,dvg(1,local_j),psic,TRIM(fftdriver))
     ELSE
        CALL single_invfft_k(dffts,npwq0,npwq0x,dvg(1,local_j),psic,TRIM(fftdriver))
     ENDIF
     IF( westpp_sign ) THEN
        DO ir = 1, dffts%nnr
           auxr(ir) = REAL( psic(ir), KIND=DP) *  ABS( REAL( psic(ir), KIND=DP) )
        ENDDO 
     ELSE
        DO ir = 1, dffts%nnr
           auxr(ir) = REAL( psic(ir), KIND=DP) *  REAL( psic(ir), KIND=DP) 
        ENDDO 
     ENDIF
     !
     WRITE(labeli,'(i6.6)') global_j 
     WRITE(labelq,'(i6.6)') 1 
     fname = TRIM( westpp_dirname ) // "/eigQ"//TRIM(labelq)//"I"//TRIM(labeli)
     CALL dump_r( auxr, fname)
     !
     CALL update_bar_type( barra,'westpp', 1 )
     !
  ENDDO
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  DEALLOCATE( auxr )
  !
END SUBROUTINE
