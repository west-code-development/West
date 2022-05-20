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
  USE io_push,               ONLY : io_push_title,io_push_bar,io_push_value
  USE westcom,               ONLY : westpp_sign,iuwfc,lrwfc,westpp_calculation, &
    &                               westpp_range,westpp_save_dir, westpp_box
  USE mp_global,             ONLY : inter_image_comm,my_image_id
  USE mp_world,              ONLY : mpime, root, world_comm
  USE mp,                    ONLY : mp_bcast, mp_sum
  USE fft_base,              ONLY : dfftp,dffts
  USE scatter_mod,           ONLY : gather_grid
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
  USE json_module,           ONLY : json_core, json_value
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: i1,i2, ipol, ir, local_j, global_j, i, ig, iks, ibnd, local_ib, &
   & global_ib, ir1, ir2, ir3, index1, index2, n_points, nbnd_, iunit
  REAL(DP),ALLOCATABLE :: aux_loc(:, :), density_loc(:), density_gat(:)
  REAL(DP) :: r_vec(3), norm, volume
  CHARACTER(LEN=512)    :: fname, aname, ikstring
  TYPE(bar_type) :: barra
  CHARACTER(LEN=6) :: labelb,labelk
  !REAL(DP) :: alat
  LOGICAL :: is_it_in
  TYPE(json_core) :: json
  TYPE(json_value), POINTER :: json_root, localization

  
  !
  CALL io_push_title("(L)ocalization Factor")
  ! determine integration volume
  volume = (westpp_box(2) - westpp_box(1))*(westpp_box(4) - westpp_box(3))*(westpp_box(6) - westpp_box(5))
  nbnd_ = westpp_range(2) - westpp_range(1) + 1
  !
  ALLOCATE(density_loc(dffts%nnr))
  ALLOCATE(aux_loc(nbnd_,k_grid%nps))
  ALLOCATE(density_gat(dffts%nr1x*dffts%nr2x*dffts%nr3x))
  !
  psic = 0._DP
  !
  CALL start_bar_type( barra, 'westpp', k_grid%nps * MAX(aband%nloc,1) ) 
  !

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
       norm = 0.0
        !
        ! local -> global
        !
        !WRITE(stdout,*) 'Made it here for ', local_ib, '/', aband%nloc
        global_ib = aband%l2g(local_ib)
        IF( global_ib < westpp_range(1) .OR. global_ib > westpp_range(2) ) CYCLE
        !
        IF( gamma_only ) THEN 
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,global_ib),psic,'Wave')
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,global_ib),psic,'Wave',igk_k(1,current_k))
        ENDIF
        ! create density
        DO ir = 1, dffts%nnr
          density_loc(ir) = REAL( psic(ir), KIND=DP) *  REAL( psic(ir), KIND=DP) 
        ENDDO 
        ! gather distributed FFT data
        CALL gather_grid(dffts, density_loc, density_gat)
        
        ! only the root has all the data and will continue 
        IF( dffts%mype == dffts%root ) THEN
          n_points = 0
          index2 = 0
          index1 = global_ib - westpp_range(1) + 1

          
          DO ir3 = 1, dffts%nr3
            Do ir2 = 1, dffts%nr2
              DO ir1 = 1, dffts%nr1
         ! DO ir3 = 1, dffts%nr3
         !   DO ir2 = 1, dffts%nr2
         !     Do ir1 = 1, dffts%nr3
                ! update current index
                is_it_in = .FALSE.
                index2 = index2 + 1
                norm = norm + density_gat(index2)
                ! find r-vector for this point
                DO i = 1, 3
                  r_vec(i) = dble(ir1 - 1)/dble(dffts%nr1) * alat*at(i,1) &
                    & + dble(ir2 - 1)/dble(dffts%nr2) * alat*at(i,2) &
                    & + dble(ir3 - 1)/dble(dffts%nr3) *alat*at(i,3)
                END DO
                
                !WRITE(stdout, '(5I7, E15.8)') ir1, ir2, ir3, index2, global_ib, &
                ! & density_gat(index2)

                ! check whether r-vector is in the box
                IF ((r_vec(1) .GT. westpp_box(1)) .AND. &
                 & (r_vec(1) .LT. westpp_box(2)) .AND. &
                  & (r_vec(2) .GT. westpp_box(3)) .AND. &
                   & (r_vec(2) .LT. westpp_box(4)) .AND. &
                    & (r_vec(3) .GT. westpp_box(5)) .AND. &
                    & (r_vec(3) .LT. westpp_box(6))) THEN
                      
                      n_points = n_points + 1
                      aux_loc(index1, iks) = aux_loc(index1, iks) + density_gat(index2)
                      is_it_in = .TRUE.
                END IF
                !WRITE(stdout, '(5I7, 3E15.8, L5)') ir1, ir2, ir3, index2, global_ib, &
                 !& r_vec(:), is_it_in
              END DO
            END DO
          ENDDO
          ! Post-processing
          ! print to file
          !print *, iks, global_ib, aux_loc(iks, index1)
          !write(stdout, *) iks, global_ib, norm/(dffts%nr1x*dffts%nr2x*dffts%nr3x)
          !CALL io_push_value("index", index1, 20)
          !CALL io_push_value("localization", aux_loc(index1, iks),20)
          !CALL io_push_value("Npoints", n_points, 20)
        
        ENDIF
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO
     !
  ENDDO
  
  aux_loc(:,:) = volume/omega * aux_loc(:,:)/dble(n_points)
  ! gather localization functions for all bands
  CALL mp_sum(aux_loc,inter_image_comm)
  ! root writes JSON file
  IF (mpime == root) THEN
    CALL json%initialize
    ! initialize the structure
    CALL json%create_object(json_root, '')

    IF (k_grid%nps == 1) THEN
      CALL json%add(json_root, 'localization', aux_loc(:,1))
    ELSE
      ! add "localization" object to "root"
      CALL json%create_object(localization, 'localization')
      CALL json%add(json_root, localization)

      DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
        WRITE(ikstring, '(I4)') iks
        CALL json%add(localization, TRIM(ADJUSTL(ikstring)), aux_loc(:,iks))
      ENDDO 
      ! don't need pointer anymore
      nullify(localization)
    ENDIF
    
    OPEN( NEWUNIT=iunit, FILE=TRIM(westpp_save_dir)//"/localization.json" )
    CALL json%print(json_root,iunit)
    CLOSE( iunit )
      !
    CALL json%destroy()
  END IF
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  !
  DEALLOCATE(density_loc,density_gat,aux_loc)
END SUBROUTINE
