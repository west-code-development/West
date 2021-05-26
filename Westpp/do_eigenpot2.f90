!
! Copyright (C) 2015-2021 M. Govoni 
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
  USE pwcom,                 ONLY : current_spin,wk,nks,nelup,neldw,isk,igk_k,xk,lsda,nkstot,&
                                  & current_k,ngk
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE westcom,               ONLY : westpp_sign,iuwfc,lrwfc,westpp_calculation,westpp_range,westpp_save_dir,fftdriver,&
                                  & npwq,npwqx,dvg,ngq,igq_q,westpp_n_pdep_eigen_to_use
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
  USE pdep_db,               ONLY : pdep_db_read
  USE types_bz_grid,         ONLY : q_grid
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: i1,i2, ipol, ir, local_j, global_j, i, ig, iks, ibnd, local_ib, global_ib, iq
  REAL(DP),ALLOCATABLE :: auxr(:)
  CHARACTER(LEN=512)    :: fname
  TYPE(bar_type) :: barra
  CHARACTER(LEN=6) :: labeli
  CHARACTER(LEN=5) :: labelq
  LOGICAL :: l_print_pdep_read
  !
  CALL io_push_title("(E)igenpotentials")
  !
  ALLOCATE(auxr(dffts%nnr))
  !
  CALL start_bar_type( barra, 'westpp', pert%nloc*q_grid%np) 
  !
  DO iq = 1, q_grid%np
     !
     IF (gamma_only) THEN
        CALL pdep_db_read(westpp_n_pdep_eigen_to_use)
     ELSE
        !l_print_pdep_read = .FALSE.
        IF (iq==1) THEN 
           l_print_pdep_read = .TRUE.
        ELSE 
           l_print_pdep_read = .FALSE.
        ENDIF
        CALL pdep_db_read(westpp_n_pdep_eigen_to_use,iq,l_print_pdep_read)
        npwq = ngq(iq)
     ENDIF
     !
     auxr = 0._DP
     psic = 0._DP
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
           CALL single_invfft_gamma(dffts,npwq,npwqx,dvg(1,local_j),psic,TRIM(fftdriver))
        ELSE
           CALL single_invfft_k(dffts,npwq,npwqx,dvg(1,local_j),psic,'Wave',igq_q(1,iq))
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
        WRITE(labelq,'(i5.5)') iq
        fname = TRIM( westpp_save_dir ) // "/eigQ"//TRIM(labelq)//"I"//TRIM(labeli)
        CALL dump_r( auxr, fname)
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  DEALLOCATE( auxr )
  !
END SUBROUTINE
