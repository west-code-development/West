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
SUBROUTINE do_loc ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE constants,             ONLY : fpi
  USE pwcom,                 ONLY : igk_k,npw,npwx,current_k,ngk
  USE gvect,                 ONLY : ngm
  USE io_push,               ONLY : io_push_title
  USE westcom,               ONLY : iuwfc,lrwfc,logfile,westpp_format,westpp_range,westpp_box,&
                                  & westpp_r0,westpp_nr,westpp_rmax
  USE mp_global,             ONLY : intra_bgrp_comm,me_bgrp,root_bgrp,inter_image_comm,my_image_id
  USE mp_world,              ONLY : mpime,root
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_base,              ONLY : dffts
  USE scatter_mod,           ONLY : scatter_grid
  USE buffers,               ONLY : get_buffer
  USE wavefunctions,         ONLY : evc,psic
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : aband
  USE class_idistribute,     ONLY : idistribute
  USE control_flags,         ONLY : gamma_only
  USE types_bz_grid,         ONLY : k_grid
  USE cell_base,             ONLY : alat,at,omega
  USE json_module,           ONLY : json_file
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  LOGICAL :: l_box
  INTEGER :: ir, i, iks, local_ib, global_ib, global_ib2, ir1, ir2, ir3, n_points, nbnd_, iunit
  REAL(DP), ALLOCATABLE :: local_fac(:,:), ipr(:,:)
  REAL(DP), ALLOCATABLE :: filter(:), filter_loc(:)
  REAL(DP), ALLOCATABLE :: spav(:)
  COMPLEX(DP), ALLOCATABLE :: auxc(:), auxg(:)
  REAL(DP) :: rho, r_vec(3)
  REAL(DP) :: r, dr
  CHARACTER(LEN=5) :: label_k
  TYPE(bar_type) :: barra
  TYPE(json_file) :: json
  !
  nbnd_ = westpp_range(2) - westpp_range(1) + 1
  aband = idistribute()
  CALL aband%init(nbnd_,'i','westpp_range',.TRUE.)
  !
  l_box = .TRUE.
  DO i = 1, 7
     IF( westpp_format(i:i) == 's' .OR. westpp_format(i:i) == 'S' ) l_box = .FALSE.
  ENDDO
  !
  IF(l_box) THEN
     ALLOCATE(filter(dffts%nr1x*dffts%nr2x*dffts%nr3x))
     ALLOCATE(filter_loc(dffts%nnr))
  ELSE
     ALLOCATE(auxc(dffts%nnr))
     ALLOCATE(auxg(ngm))
     ALLOCATE(spav(westpp_nr+1))
  ENDIF
  ALLOCATE(local_fac(nbnd_,k_grid%nps))
  ALLOCATE(ipr(nbnd_,k_grid%nps))
  !
  CALL io_push_title('(L)ocalization Factor')
  !
  CALL start_bar_type( barra, 'westpp', k_grid%nps * MAX(aband%nloc,1) )
  !
  IF(l_box) THEN
     !
     ! only the FFT root generates the filter
     !
     IF(me_bgrp == root_bgrp) THEN
        !
        ! for each point of FFT grid, the filter is ONE when point is in box, ZERO if not
        !
        filter(:) = 0._DP
        n_points = 0
        ir = 0
        DO ir3 = 1, dffts%nr3
           DO ir2 = 1, dffts%nr2
              DO ir1 = 1, dffts%nr1
                 ir = ir + 1
                 !
                 ! create real-space vector
                 !
                 DO i = 1, 3
                    r_vec(i) = alat * ( REAL(ir1-1,KIND=DP)/REAL(dffts%nr1,KIND=DP)*at(i,1) &
                                    & + REAL(ir2-1,KIND=DP)/REAL(dffts%nr2,KIND=DP)*at(i,2) &
                                    & + REAL(ir3-1,KIND=DP)/REAL(dffts%nr3,KIND=DP)*at(i,3) )
                 ENDDO
                 !
                 ! check point is in box
                 !
                 IF( (r_vec(1) > westpp_box(1)) .AND. (r_vec(1) < westpp_box(2)) .AND. &
                   & (r_vec(2) > westpp_box(3)) .AND. (r_vec(2) < westpp_box(4)) .AND. &
                   & (r_vec(3) > westpp_box(5)) .AND. (r_vec(3) < westpp_box(6)) ) THEN
                    n_points = n_points + 1
                    filter(ir) = 1._DP
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     !
     ! scatter filter to all FFT processes
     !
     CALL scatter_grid(dffts, filter, filter_loc)
     !
     ! broadcast the number of points in box to all FFT processes
     !
     CALL mp_bcast(n_points, root_bgrp, intra_bgrp_comm)
     !
     IF( n_points == 0 ) CALL errore('do_loc','no point found in integration volume',1)
     !
  ELSE
     !
     dr = westpp_rmax/REAL(westpp_nr,KIND=DP)
     !
  ENDIF
  !
  ipr = 0._DP
  local_fac = 0._DP
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
        global_ib = aband%l2g(local_ib)+westpp_range(1)-1
        global_ib2 = aband%l2g(local_ib)
        !
        IF( gamma_only ) THEN
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,global_ib),psic,'Wave')
           !
           IF(l_box) THEN
              !
              DO ir = 1, dffts%nnr
                 rho = REAL(psic(ir),KIND=DP)**2
                 local_fac(global_ib2,iks) = local_fac(global_ib2,iks)+filter_loc(ir)*rho
                 ipr(global_ib2,iks) = ipr(global_ib2,iks)+rho**2
              ENDDO
              !
           ELSE
              !
              DO ir = 1, dffts%nnr
                 rho = REAL(psic(ir),KIND=DP)**2
                 auxc(ir) = CMPLX(rho,KIND=DP)
                 ipr(global_ib2,iks) = ipr(global_ib2,iks)+rho**2
              ENDDO
              !
              CALL single_fwfft_gamma(dffts,ngm,ngm,auxc,auxg,'Rho')
              !
              CALL wfc_spav(auxg,westpp_r0,westpp_nr,westpp_rmax,spav)
              !
              DO ir = 1, westpp_nr+1
                 r = REAL(ir-1,KIND=DP) * dr
                 local_fac(global_ib2,iks) = local_fac(global_ib2,iks)+(r**2)*spav(ir)
              ENDDO
              !
           ENDIF
           !
        ELSE
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(:,global_ib),psic,'Wave',igk_k(:,current_k))
           !
           DO ir = 1, dffts%nnr
              rho = REAL(CONJG(psic(ir))*psic(ir),KIND=DP)
              local_fac(global_ib2,iks) = local_fac(global_ib2,iks)+filter_loc(ir)*rho
              ipr(global_ib2,iks) = ipr(global_ib2,iks)+rho**2
           ENDDO
           !
        ENDIF
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO
     !
  ENDDO
  !
  ! Post processing
  !
  IF(l_box) THEN
     local_fac(:,:) = local_fac(:,:)/(dffts%nr1*dffts%nr2*dffts%nr3)
  ELSE
     local_fac(:,:) = local_fac(:,:)*fpi*dr/omega
  ENDIF
  ipr(:,:) = ipr(:,:)/(dffts%nr1*dffts%nr2*dffts%nr3)/omega
  !
  ! sum up results
  !
  CALL mp_sum(local_fac,intra_bgrp_comm)
  CALL mp_sum(ipr,intra_bgrp_comm)
  CALL mp_sum(local_fac,inter_image_comm)
  CALL mp_sum(ipr,inter_image_comm)
  !
  ! root writes JSON file
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     !
     ! output localization factor and IPR
     !
     DO iks = 1,k_grid%nps
        !
        WRITE(label_k,'(I5.5)') iks
        !
        CALL json%add('output.L.K'//label_k//'.local_factor',local_fac(:,iks))
        CALL json%add('output.L.K'//label_k//'.ipr',ipr(:,iks))
        !
     ENDDO
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
     !
  ENDIF
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  IF(ALLOCATED(filter)) DEALLOCATE(filter)
  IF(ALLOCATED(filter_loc)) DEALLOCATE(filter_loc)
  IF(ALLOCATED(local_fac)) DEALLOCATE(local_fac)
  IF(ALLOCATED(ipr)) DEALLOCATE(ipr)
  IF(ALLOCATED(auxc)) DEALLOCATE(auxc)
  IF(ALLOCATED(auxg)) DEALLOCATE(auxg)
  IF(ALLOCATED(spav)) DEALLOCATE(spav)
  !
END SUBROUTINE
