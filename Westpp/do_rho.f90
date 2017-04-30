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
SUBROUTINE do_rho ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE uspp,                  ONLY : vkb,nkb
  USE io_global,             ONLY : stdout
  USE pwcom,                 ONLY : current_spin,wk,nks,nelup,neldw,isk,g,igk_k,ngm,tpiba2,xk,npw,npwx,lsda,nkstot,&
                                  & current_k,ngk
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE westcom,               ONLY : westpp_sign,iuwfc,lrwfc,westpp_calculation,westpp_range,westpp_dirname,nbnd_occ 
  USE mp_global,             ONLY : inter_image_comm,my_image_id,intra_image_comm
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_base,              ONLY : dfftp,dffts
  USE wvfct,                 ONLY : nbnd
  USE buffers,               ONLY : get_buffer
  USE wavefunctions_module,  ONLY : evc,psic
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : aband
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
  CHARACTER(LEN=6) :: label
  !
  CALL io_push_title("(R)ho")
  !
  ALLOCATE(auxr(dffts%nnr))
  !
  auxr = 0._DP
  psic = 0._DP
  !
  CALL start_bar_type( barra, 'westpp', nks ) 
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     !IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks>1) THEN
        !iuwfc = 20
        !lrwfc = nbnd * npwx * npol 
        !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        !CALL mp_bcast(evc,0,inter_image_comm)
        !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     !nbndval = nbnd_occ(iks)
     !
     DO local_ib=1,aband%nloc
        !
        ! local -> global
        !
        global_ib = aband%l2g(local_ib)
        IF( global_ib > nbnd_occ(iks) ) CYCLE
        !
        IF( gamma_only ) THEN 
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,global_ib),psic,'Wave')
           DO ir = 1, dffts%nnr
              auxr(ir) = auxr(ir) + REAL( psic(ir), KIND=DP) *  REAL( psic(ir), KIND=DP) * wk(iks)  
           ENDDO 
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,global_ib),psic,'Wave',igk_k(1,current_k))
           DO ir = 1, dffts%nnr
              auxr(ir) = auxr(ir) + REAL( CONJG( psic(ir) ) * psic(ir), KIND=DP) * wk(iks)
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
  fname = TRIM( westpp_dirname ) // "/rho"
  IF(my_image_id==0) CALL dump_r( auxr, fname)
  !
  DEALLOCATE( auxr )
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
END SUBROUTINE
