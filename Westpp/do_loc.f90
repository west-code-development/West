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
  USE pwcom,                 ONLY : igk_k,npw,npwx,current_k,ngk,nspin
  USE noncollin_module,      ONLY : npol
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
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : aband
  USE class_idistribute,     ONLY : idistribute
  USE control_flags,         ONLY : gamma_only
  USE types_bz_grid,         ONLY : k_grid
  USE cell_base,             ONLY : alat,at,omega
  USE json_module,           ONLY : json_file,json_core,json_value
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
  LOGICAL :: l_box
  INTEGER :: ig, ir, i, iks, ib, ib_g, ib2_g, jb, ir1, ir2, ir3, n_points, nbnd_, iunit
  INTEGER :: dffts_nnr
  REAL(DP), ALLOCATABLE :: local_fac(:,:), ipr(:,:)
  REAL(DP), ALLOCATABLE :: filter(:), filter_loc(:)
  REAL(DP), ALLOCATABLE :: spav(:)
  REAL(DP), ALLOCATABLE :: ovlp_ab(:,:)
  COMPLEX(DP), ALLOCATABLE :: auxc(:), auxg(:)
  !$acc declare device_resident(auxc)
  COMPLEX(DP), ALLOCATABLE :: evc_tmp(:,:,:)
  REAL(DP) :: rho, r_vec(3)
  REAL(DP) :: r, dr
  REAL(DP) :: reduce, reduce2
  CHARACTER(LEN=5) :: label_k
  CHARACTER(LEN=9) :: label_b
  TYPE(bar_type) :: barra
  TYPE(json_file) :: json
  TYPE(json_core) :: jcor
  TYPE(json_value), POINTER :: jval
  !
  nbnd_ = westpp_range(2) - westpp_range(1) + 1
  aband = idistribute()
  CALL aband%init(nbnd_,'i','westpp_range',.TRUE.)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  dffts_nnr = dffts%nnr
  !
  l_box = .TRUE.
  DO i = 1, 7
     IF(westpp_format(i:i) == 's' .OR. westpp_format(i:i) == 'S') l_box = .FALSE.
  ENDDO
  !
  IF(l_box) THEN
     ALLOCATE(filter(dffts%nr1x*dffts%nr2x*dffts%nr3x))
     ALLOCATE(filter_loc(dffts%nnr))
     !$acc enter data create(filter_loc)
  ELSE
     ALLOCATE(auxc(dffts%nnr))
     ALLOCATE(auxg(ngm))
     ALLOCATE(spav(westpp_nr+1))
     !$acc enter data create(auxg,spav)
  ENDIF
  ALLOCATE(local_fac(nbnd_,k_grid%nps))
  ALLOCATE(ipr(nbnd_,k_grid%nps))
  IF(gamma_only .AND. nspin == 2) THEN
     ALLOCATE(evc_tmp(npwx,nbnd_,nspin))
     ALLOCATE(ovlp_ab(nbnd_,nbnd_))
     !$acc enter data create(evc_tmp,ovlp_ab)
  ENDIF
  !
  CALL io_push_title('(L)ocalization Factor')
  !
  CALL start_bar_type(barra, 'westpp', k_grid%nps * MAX(aband%nloc,1))
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
                    r_vec(i) = alat * (REAL(ir1-1,KIND=DP)/REAL(dffts%nr1,KIND=DP)*at(i,1) &
                                   & + REAL(ir2-1,KIND=DP)/REAL(dffts%nr2,KIND=DP)*at(i,2) &
                                   & + REAL(ir3-1,KIND=DP)/REAL(dffts%nr3,KIND=DP)*at(i,3))
                 ENDDO
                 !
                 ! check point is in box
                 !
                 IF((r_vec(1) > westpp_box(1)) .AND. (r_vec(1) < westpp_box(2)) .AND. &
                  & (r_vec(2) > westpp_box(3)) .AND. (r_vec(2) < westpp_box(4)) .AND. &
                  & (r_vec(3) > westpp_box(5)) .AND. (r_vec(3) < westpp_box(6))) THEN
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
     !$acc update device(filter_loc)
     !
     ! broadcast the number of points in box to all FFT processes
     !
     CALL mp_bcast(n_points, root_bgrp, intra_bgrp_comm)
     !
     IF(n_points == 0) CALL errore('do_loc','no point found in integration volume',1)
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
     IF(gamma_only .AND. nspin == 2) THEN
        !
        !$acc parallel loop collapse(2) present(evc_tmp)
        DO ib_g = westpp_range(1), westpp_range(2)
           DO ig = 1, npwx
              evc_tmp(ig,ib_g,iks) = evc_work(ig,ib_g)
           ENDDO
        ENDDO
        !$acc end parallel
        !
     ENDIF
     !
     DO ib = 1,aband%nloc
        !
        ! local -> global
        !
        ib_g = aband%l2g(ib)+westpp_range(1)-1
        ib2_g = aband%l2g(ib)
        !
        reduce = 0._DP
        reduce2 = 0._DP
        !
        IF(gamma_only) THEN
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ib_g),psic,'Wave')
           !
           IF(l_box) THEN
              !
              !$acc parallel loop reduction(+:reduce,reduce2) present(filter_loc) copy(reduce,reduce2)
              DO ir = 1, dffts_nnr
                 rho = REAL(psic(ir),KIND=DP)**2
                 reduce = reduce+filter_loc(ir)*rho
                 reduce2 = reduce2+rho**2
              ENDDO
              !$acc end parallel
              !
              local_fac(ib2_g,iks) = reduce
              ipr(ib2_g,iks) = reduce2
              !
           ELSE
              !
              !$acc parallel loop reduction(+:reduce2) present(auxc) copy(reduce2)
              DO ir = 1, dffts_nnr
                 rho = REAL(psic(ir),KIND=DP)**2
                 auxc(ir) = CMPLX(rho,KIND=DP)
                 reduce2 = reduce2+rho**2
              ENDDO
              !$acc end parallel
              !
              ipr(ib2_g,iks) = reduce2
              !
              !$acc host_data use_device(auxc,auxg)
              CALL single_fwfft_gamma(dffts,ngm,ngm,auxc,auxg,'Rho')
              !$acc end host_data
              !
              !$acc update host(auxg)
              !
              CALL wfc_spav(auxg,westpp_r0,westpp_nr,westpp_rmax,spav)
              !
              !$acc update device(spav)
              !
              !$acc parallel loop reduction(+:reduce) present(spav) copy(reduce)
              DO ir = 1, westpp_nr+1
                 r = REAL(ir-1,KIND=DP) * dr
                 reduce = reduce+(r**2)*spav(ir)
              ENDDO
              !$acc end parallel
              !
              local_fac(ib2_g,iks) = reduce
              !
           ENDIF
           !
        ELSE
           !
           CALL single_invfft_k(dffts,npw,npwx,evc_work(:,ib_g),psic,'Wave',igk_k(:,current_k))
           !
           !$acc parallel loop reduction(+:reduce,reduce2) present(filter_loc) copy(reduce,reduce2)
           DO ir = 1, dffts_nnr
              rho = REAL(CONJG(psic(ir))*psic(ir),KIND=DP)
              reduce = reduce+filter_loc(ir)*rho
              reduce2 = reduce2+rho**2
           ENDDO
           !$acc end parallel
           !
           local_fac(ib2_g,iks) = reduce
           ipr(ib2_g,iks) = reduce2
           !
        ENDIF
        !
        CALL update_bar_type(barra,'westpp', 1)
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
  IF(gamma_only .AND. nspin == 2) THEN
     !
     !$acc host_data use_device(evc_tmp,ovlp_ab)
     CALL glbrak_gamma(evc_tmp(:,:,2),evc_tmp(:,:,1),ovlp_ab,npw,npwx,nbnd_,nbnd_,nbnd_,npol)
     !$acc end host_data
     !
     !$acc update host(ovlp_ab)
     !
     CALL mp_sum(ovlp_ab,intra_bgrp_comm)
     !
  ENDIF
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
     DO iks = 1, k_grid%nps
        WRITE(label_k,'(I5.5)') iks
        CALL json%add('output.L.K'//label_k//'.local_factor',local_fac(:,iks))
        CALL json%add('output.L.K'//label_k//'.ipr',ipr(:,iks))
     ENDDO
     !
     IF(gamma_only .AND. nspin == 2) THEN
        !
        CALL jcor%create_array(jval,'overlap_ab')
        CALL json%add('output.L.overlap_ab',jval)
        !
        DO ib = 1, nbnd_
           !
           jb = MAXLOC(ABS(ovlp_ab(:,ib)),DIM=1)
           !
           WRITE(label_b,'(I9)') ib
           !
           CALL json%add('output.L.overlap_ab('//label_b//').ib',ib+westpp_range(1)-1)
           CALL json%add('output.L.overlap_ab('//label_b//').jb',jb+westpp_range(1)-1)
           CALL json%add('output.L.overlap_ab('//label_b//').val',ABS(ovlp_ab(jb,ib)))
           !
        ENDDO
        !
     ENDIF
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
     !
  ENDIF
  !
  CALL stop_bar_type(barra, 'westpp')
  !
  IF(l_box) THEN
     DEALLOCATE(filter)
     !$acc exit data delete(filter_loc)
     DEALLOCATE(filter_loc)
  ELSE
     DEALLOCATE(auxc)
     !$acc exit data delete(auxg,spav)
     DEALLOCATE(auxg)
     DEALLOCATE(spav)
  ENDIF
  DEALLOCATE(local_fac)
  DEALLOCATE(ipr)
  IF(gamma_only .AND. nspin == 2) THEN
     !$acc exit data delete(evc_tmp,ovlp_ab)
     DEALLOCATE(evc_tmp)
     DEALLOCATE(ovlp_ab)
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
