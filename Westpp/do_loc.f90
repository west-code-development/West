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
SUBROUTINE do_loc ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE uspp,                  ONLY : vkb,nkb
  USE io_global,             ONLY : stdout
  USE pwcom,                 ONLY : current_spin,wk,nelup,neldw,isk,igk_k,xk,npw,npwx,lsda,current_k,ngk
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE westcom,               ONLY : westpp_sign,iuwfc,lrwfc,westpp_calculation,westpp_range,westpp_save_dir 
  USE mp_global,             ONLY : inter_image_comm,my_image_id
  USE mp,                    ONLY : mp_bcast
  USE fft_base,              ONLY : dfftp,dffts
  USE wvfct,                 ONLY : nbnd
  USE buffers,               ONLY : get_buffer
  USE wavefunctions_module,  ONLY : evc,psic
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : aband
  USE control_flags,         ONLY : gamma_only 
  USE types_bz_grid,         ONLY : k_grid
  USE cell_base,             ONLY : alat, at, omega
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: i1,i2, ipol, ir, local_j, global_j, i, ig, iks, ibnd, local_ib, &
   & global_ib, ir1, ir2, ir3, index1, index2, n_points
  REAL(DP),ALLOCATABLE :: aux_loc(:, :)
  REAL(DP) :: r_vec(3)
  CHARACTER(LEN=512)    :: fname
  TYPE(bar_type) :: barra
  CHARACTER(LEN=6) :: labelb,labelk

  
  REAL(DP),ALLOCATABLE :: density_loc(:)
  REAL(DP) :: density_gat(dfft%nr1x*dfft%nr2x*dfft%nr3x)
  !
  CALL io_push_title("(L)ocalization Factor")
  !
  ALLOCATE(density_loc(dffts%nnr))
  ALLOCATE(aux_loc(k_grid%nps, aband%nloc))
  !
  auxr = 0._DP
  psic = 0._DP
  !
  CALL start_bar_type( barra, 'westpp', k_grid%nps * MAX(aband%nloc,1) ) 
  !
  DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
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
     IF(k_grid%nps>1) THEN
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
        IF( global_ib < westpp_range(1) .OR. global_ib > westpp_range(2) ) CYCLE
        !
        auxr = 0._DP
        IF( gamma_only ) THEN 
           CALL singe_invfft_gamma(dffts,npw,npwx,evc(1,global_ib),psic,'Wave')
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,global_ib),psic,'Wave',igk_k(1,current_k))
        ENDIF
        ! create density
        DO ir = 1, dffts%nnr
          density_loc(ir) = REAL( psic(ir), KIND=DP) *  REAL( psic(ir), KIND=DP) 
        ENDDO 
        ! gather distributed FFT data
        CALL gather_grid(dfft, density_loc, density_gat)
        ! only the root has all the data and will continue 
        IF( dfft%mype == dfft%root ) THEN
          n_points = 0
          index2 = 0
          index1 = global_ib - westpp_range(1) + 1
          
          DO ir1 = 1, dfft%nr1x
            DO ir2 = 1, dfft%nr2x
              Do ir3 = 1, dfft%nr3x
                ! update current index
                index2 = index2 + 1
                ! find r-vector for this point
                r_vec(:) = dble(ir1 - 1)/dble(dfft%nr1x) * alat*at(:,1) &
                  & + dble(ir2 - 1)/dble(dfft%nr2x) * alat*at(:,2) &
                  & + dble(ir3 - 1)/dble(dfft%nr3x) *alat*at(:,3)

                ! check whether r-vector is in the box
                IF ((r_vec(1) .GT. westpp_box(1)) .AND. &
                 & r_vec(1) .LT. westpp_box(2)) THEN
                  IF ((r_vec(2) .GT. westpp_box(3)) .AND. &
                   & r_vec(2) .LT. westpp_box(4)) THEN
                    IF ((r_vec(3) .GT. westpp_box(5)) .AND. &
                    & r_vec(3) .LT. westpp_box(6)) THEN
                      
                      n_points = n_points + 1
                      aux_loc(iks, index1) = aux_loc(iks, index1) + density_gat(index2)
                    END IF
                  END IF
                END IF
              END DO
            END DO
          ENDDO
          ! Post-processing
          aux_loc(:,:) = omega/n_points * aux_loc(:,:)
          ! print to file
          print *, iks, global_ib aux_loc(iks, index1)

        ENDIF
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
