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
  USE pwcom,                 ONLY : current_spin,isk,igk_k,npw,npwx,lsda,current_k,ngk
  USE io_push,               ONLY : io_push_title
  USE westcom,               ONLY : iuwfc,lrwfc,westpp_range,westpp_save_dir,westpp_box
  USE mp_global,             ONLY : inter_image_comm,my_image_id
  USE mp_world,              ONLY : mpime,root
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_base,              ONLY : dffts
  USE scatter_mod,           ONLY : scatter_grid
  USE buffers,               ONLY : get_buffer
  USE wavefunctions,         ONLY : evc,psic
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : aband
  USE class_idistribute,     ONLY : idistribute
  USE control_flags,         ONLY : gamma_only
  USE types_bz_grid,         ONLY : k_grid
  USE cell_base,             ONLY : alat,at,omega
  USE json_module,           ONLY : json_core,json_value
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ir, i, iks, local_ib, global_ib, ir1, ir2, ir3, index1, index2, n_points, nbnd_, iunit
  REAL(DP), ALLOCATABLE :: aux_loc(:, :), ipr(:, :), density_loc(:), density_gat(:), ipr_loc(:)
  REAL(DP), ALLOCATABLE :: filter(:), filter_loc(:)
  REAL(DP) :: r_vec(3), volume
  CHARACTER(LEN=512) :: ikstring
  TYPE(bar_type) :: barra
  TYPE(json_core) :: json
  TYPE(json_value), POINTER :: json_root, localization_object, ipr_object
  !
  CALL io_push_title('(L)ocalization Factor')
  !
  ! determine integration volume
  volume = (westpp_box(2) - westpp_box(1))*(westpp_box(4) - westpp_box(3))*(westpp_box(6) - westpp_box(5))
  nbnd_ = westpp_range(2) - westpp_range(1) + 1
  aband = idistribute()
  CALL aband%init(nbnd_,'i','westpp_range',.TRUE.)
  !
  ALLOCATE(density_loc(dffts%nnr))
  ALLOCATE(filter_loc(dffts%nnr))
  ALLOCATE(ipr_loc(dffts%nnr))
  ALLOCATE(aux_loc(nbnd_,k_grid%nps))
  ALLOCATE(ipr(nbnd_,k_grid%nps))
  ALLOCATE(density_gat(dffts%nr1x*dffts%nr2x*dffts%nr3x))
  ALLOCATE(filter(dffts%nr1x*dffts%nr2x*dffts%nr3x))
  !
  ipr = 0._DP
  aux_loc = 0._DP
  psic = 0._DP
  !
  CALL start_bar_type( barra, 'westpp', k_grid%nps * MAX(aband%nloc,1) )
  !
  ! only the FFT root generates the filter
  IF( dffts%mype == dffts%root ) THEN
    !
    ! create filter: for each point in space, the filter is ONE when point
    ! is in box, ZERO if not. Loop over all points on FFT grid
    filter(:) = 0._DP
    n_points = 0
    index2 = 0
    DO ir3 = 1, dffts%nr3
      Do ir2 = 1, dffts%nr2
        DO ir1 = 1, dffts%nr1
          ! update current index
          index2 = index2 + 1
          ! create real-space vector
          DO i = 1, 3
            r_vec(i) = REAL(ir1-1,KIND=DP)/REAL(dffts%nr1,KIND=DP) * alat*at(i,1) &
              & + REAL(ir2-1,KIND=DP)/REAL(dffts%nr2,KIND=DP) * alat*at(i,2) &
              & + REAL(ir3-1,KIND=DP)/REAL(dffts%nr3,KIND=DP) *alat*at(i,3)
          END DO
          ! check whether r-vector is in the box
          IF ((r_vec(1) > westpp_box(1)) .AND. &
          &   (r_vec(1) < westpp_box(2)) .AND. &
          &   (r_vec(2) > westpp_box(3)) .AND. &
          &   (r_vec(2) < westpp_box(4)) .AND. &
          &   (r_vec(3) > westpp_box(5)) .AND. &
          &   (r_vec(3) < westpp_box(6))) THEN
            ! if condition fulfilled, set filter to 1
            n_points = n_points + 1
            filter(index2) = 1._DP
          END IF
        END DO
      END DO
    END DO
  END IF
  !
  ! scatter the filter function to all FFT processes
  CALL scatter_grid(dffts, filter, filter_loc)
  ! broadcast the number of points in rectangle to all FFT processes
  CALL mp_bcast(n_points, dffts%root, dffts%comm)
  !
  DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     npw = ngk(iks)
     IF(lsda) current_spin = isk(iks)
     call g2_kin(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(k_grid%nps>1) THEN
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     !
     DO local_ib=1,aband%nloc
        !
        ! local -> global
        !
        global_ib = aband%l2g(local_ib)+westpp_range(1)-1
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
        !
        index1 = global_ib - westpp_range(1) + 1
        ! generate Inverse Participation Ratio for each point
        ipr_loc(:) = density_loc(:)**2
        ! sum over all points
        ipr(index1, iks) = SUM(ipr_loc)
        ! filter density
        density_loc(:) = filter_loc(:) * density_loc(:)
        ! sum localization factor
        aux_loc(index1, iks) = SUM(density_loc)
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO
     !
  ENDDO
  ! Post processing
  aux_loc(:,:) = volume/omega * aux_loc(:,:)/REAL(n_points,KIND=DP)
  ipr(:,:) = ipr(:,:)/REAL(dffts%nr1*dffts%nr2*dffts%nr3,KIND=DP)/omega
  !
  ! sum ipr and localization factor over all FFT processes
  CALL mp_sum(aux_loc, dffts%comm)
  CALL mp_sum(ipr, dffts%comm)
  !
  ! gather localization functions for all bands
  CALL mp_sum(aux_loc,inter_image_comm)
  ! gather IPR for all bands
  CALL mp_sum(ipr,inter_image_comm)
  !
  ! root writes JSON file
  IF (mpime == root) THEN
    CALL json%initialize
    ! initialize the structure
    CALL json%create_object(json_root, '')
    !
    IF (k_grid%nps == 1) THEN
      ! write localization factor and IPR to JSON file directly
      CALL json%add(json_root, 'localization', aux_loc(:,1))
      CALL json%add(json_root, 'ipr', ipr(:,1))
    ELSE
      ! add localization object to root
      CALL json%create_object(localization_object, 'localization')
      CALL json%create_object(ipr_object, 'ipr')
      CALL json%add(json_root, localization_object)
      CALL json%add(json_root, ipr_object)
      !
      DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
        WRITE(ikstring, '(I4)') iks
        ! write localization factor and IPR for each iks-point
        CALL json%add(localization_object, TRIM(ADJUSTL(ikstring)), aux_loc(:,iks))
        CALL json%add(ipr_object, TRIM(ADJUSTL(ikstring)), ipr(:,iks))
      ENDDO
      ! don't need pointer anymore
      NULLIFY(localization_object)
      NULLIFY(ipr_object)
    ENDIF
    !
    OPEN( NEWUNIT=iunit, FILE=TRIM(westpp_save_dir)//'/localization.json' )
    CALL json%print(json_root,iunit)
    CLOSE( iunit )
    !
    CALL json%destroy()
  END IF
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  DEALLOCATE(density_loc)
  DEALLOCATE(filter_loc)
  DEALLOCATE(ipr_loc)
  DEALLOCATE(aux_loc)
  DEALLOCATE(ipr)
  DEALLOCATE(density_gat)
  DEALLOCATE(filter)
  !
END SUBROUTINE
