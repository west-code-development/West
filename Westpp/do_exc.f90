!
! Copyright (C) 2015-2024 M. Govoni
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
SUBROUTINE do_exc()
  !
  USE kinds,                 ONLY : DP
  USE io_push,               ONLY : io_push_title
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE cell_base,             ONLY : omega,at,alat
  USE fft_base,              ONLY : dffts
  USE pwcom,                 ONLY : npw,npwx,igk_k,current_k,ngk,wg
  USE control_flags,         ONLY : gamma_only
  USE mp,                    ONLY : mp_sum,mp_bcast,mp_min
  USE mp_global,             ONLY : me_bgrp,intra_bgrp_comm,my_image_id,inter_image_comm
  USE buffers,               ONLY : get_buffer
  USE westcom,               ONLY : iuwfc,lrwfc,nbndval0x,nbnd_occ,dvg_exc,westpp_range,&
                                  & westpp_n_liouville_to_use,westpp_r0,westpp_save_dir
  USE fft_at_gamma,          ONLY : double_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE plep_db,               ONLY : plep_db_read
  USE distribution_center,   ONLY : pert,kpt_pool,band_group
  USE class_idistribute,     ONLY : idistribute
  USE types_bz_grid,         ONLY : k_grid
  USE wavefunctions,         ONLY : evc,psic
#if defined(__CUDA)
  USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ibnd,nbndval,i,j,k,ir,ip,iks,index0,ir_end,iexc,lexc,drmin_id,dffts_nnr
  INTEGER :: barra_load
  REAL(DP) :: drmin_g,tmp,rcoeff,w1,summ0,inv_nr1,inv_nr2,inv_nr3
  REAL(DP) :: r0(3)
  COMPLEX(DP) :: zcoeff
  REAL(DP), ALLOCATABLE :: r(:,:),dr(:),rho(:)
  COMPLEX(DP), ALLOCATABLE :: psic_aux(:),rho_aux(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: n_ipol = 3
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  dffts_nnr = dffts%nnr
  !
  r0(1) = westpp_r0(1)/alat
  r0(2) = westpp_r0(2)/alat
  r0(3) = westpp_r0(3)/alat
  !
  DO ip = 1, n_ipol
     !
     DO WHILE(r0(ip) < 0)
        r0(ip) = r0(ip) + at(ip,ip)
     ENDDO
     !
     DO WHILE(r0(ip) >= at(ip,ip))
        r0(ip) = r0(ip) - at(ip,ip)
     ENDDO
     !
  ENDDO
  !
  ! ... DISTRIBUTE
  !
  pert = idistribute()
  CALL pert%init(westpp_n_liouville_to_use,'i','nvec',.TRUE.)
  kpt_pool = idistribute()
  CALL kpt_pool%init(k_grid%nps,'p','kpt',.FALSE.)
  band_group = idistribute()
  CALL band_group%init(nbndval0x,'b','nbndval',.FALSE.)
  !
  ! READ EIGENVALUES AND VECTORS FROM OUTPUT
  !
  CALL plep_db_read(westpp_n_liouville_to_use)
  !
  ALLOCATE(r(dffts%nnr,n_ipol))
  ALLOCATE(dr(dffts%nnr))
  ALLOCATE(rho(dffts%nnr))
  ALLOCATE(rho_aux(dffts%nnr))
  !$acc enter data create(rho,rho_aux) copyin(dvg_exc)
  IF(.NOT. gamma_only) THEN
     ALLOCATE(psic_aux(dffts%nnr))
     !$acc enter data create(psic_aux)
  ENDIF
  !
  r(:,:) = 0._DP
  !
  ! Calculate r
  !
  inv_nr1 = 1._DP / REAL(dffts%nr1,KIND=DP)
  inv_nr2 = 1._DP / REAL(dffts%nr2,KIND=DP)
  inv_nr3 = 1._DP / REAL(dffts%nr3,KIND=DP)
  !
  index0 = dffts%nr1x*dffts%nr2x*SUM(dffts%nr3p(1:me_bgrp))
  ir_end = MIN(dffts%nnr,dffts%nr1x*dffts%nr2x*dffts%nr3p(me_bgrp+1))
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dffts%nr1x*dffts%nr2x)
     i = i - (dffts%nr1x*dffts%nr2x)*k
     j = i / dffts%nr1x
     i = i - dffts%nr1x*j
     !
     DO ip = 1, n_ipol
        r(ir,ip) = REAL(i,KIND=DP)*inv_nr1*at(ip,1) + &
                   REAL(j,KIND=DP)*inv_nr2*at(ip,2) + &
                   REAL(k,KIND=DP)*inv_nr3*at(ip,3)
     ENDDO
     !
  ENDDO
  !
  DO ir = 1, dffts%nnr
     dr(ir) = SQRT((r(ir,1) - r0(1))**2 + &
                   (r(ir,2) - r0(2))**2 + &
                   (r(ir,3) - r0(3))**2)
  ENDDO
  !
  drmin_g = MINVAL(dr)
  !
  CALL mp_min(drmin_g,intra_bgrp_comm)
  !
  drmin_id = -1
  DO ir = 1, dffts%nnr
     IF(dr(ir) == drmin_g) THEN
        drmin_id = ir
     ENDIF
  ENDDO
  !
  CALL io_push_title('E(X)citon State')
  !
  barra_load = 0
  DO lexc = 1,pert%nloc
     iexc = pert%l2g(lexc)
     IF(iexc < westpp_range(1) .OR. iexc > westpp_range(2)) CYCLE
     barra_load = barra_load+1
  ENDDO
  !
  CALL start_bar_type(barra,'westpp',barra_load*k_grid%nps*SUM(nbnd_occ))
  !
  DO lexc = 1,pert%nlocx
     !
     ! local -> global
     !
     iexc = pert%l2g(lexc)
     !
     DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
        !
        ! ... Set k-point, spin, kinetic energy, needed by Hpsi
        !
        current_k = iks
        npw = ngk(iks)
        nbndval = nbnd_occ(iks)
        !
        ! ... read in wavefunctions from the previous iteration
        !
        IF(k_grid%nps > 1) THEN
           IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc,0,inter_image_comm)
           !$acc update device(evc)
        ENDIF
        !
        ! CYCLE here because of mp_bcast above
        !
        IF(iexc < westpp_range(1) .OR. iexc > westpp_range(2)) CYCLE
        !
        !$acc kernels present(rho_aux)
        rho_aux(:) = 0._DP
        !$acc end kernels
        !
        DO ibnd = 1, nbndval
           !
           w1 = wg(ibnd,iks)
           !
           IF(gamma_only) THEN
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),dvg_exc(:,ibnd,iks,lexc),psic,'Wave')
              !
              rcoeff = 0._DP
              IF(drmin_id > 0) THEN
                 tmp = psic(drmin_id)
                 rcoeff = REAL(tmp,KIND=DP)
              ENDIF
              CALL mp_sum(rcoeff,intra_bgrp_comm)
              !
              !$acc parallel loop present(rho_aux)
              DO ir = 1, dffts_nnr
                 rho_aux(ir) = rho_aux(ir) + w1 * CMPLX(rcoeff*AIMAG(psic(ir)),KIND=DP)
              ENDDO
              !$acc end parallel
              !
           ELSE
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,current_k))
              CALL single_invfft_k(dffts,npw,npwx,dvg_exc(:,ibnd,iks,lexc),psic_aux,'Wave',igk_k(:,current_k))
              !
              zcoeff = 0._DP
              IF(drmin_id > 0) THEN
                 zcoeff = psic(drmin_id)
              ENDIF
              CALL mp_sum(zcoeff,intra_bgrp_comm)
              !
              !$acc parallel loop present(rho_aux,psic_aux)
              DO ir = 1, dffts_nnr
                 rho_aux(ir) = rho_aux(ir) + w1 * zcoeff * psic_aux(ir)
              ENDDO
              !$acc end parallel
              !
           ENDIF
           !
           CALL update_bar_type(barra,'westpp',1)
           !
        ENDDO
        !
        !$acc parallel loop present(rho,rho_aux)
        DO ir = 1, dffts_nnr
           rho(ir) = (REAL(rho_aux(ir),KIND=DP)**2 + AIMAG(rho_aux(ir))**2) / omega
        ENDDO
        !$acc end parallel
        !
        summ0 = 0._DP
        !$acc parallel loop reduction(+:summ0) present(rho) copy(summ0)
        DO ir = 1, dffts_nnr
           summ0 = summ0 + rho(ir)
        ENDDO
        !$acc end parallel
        !
        summ0 = summ0*omega/(dffts%nr1*dffts%nr2*dffts%nr3)
        !
        CALL mp_sum(summ0,intra_bgrp_comm)
        !
        !$acc parallel loop present(rho)
        DO ir = 1, dffts_nnr
           rho(ir) = rho(ir) / summ0
        ENDDO
        !$acc end parallel
        !
        WRITE(fname,'(a,i6.6,a,i6.6)') TRIM(westpp_save_dir)//'/excK',iks,'E',iexc
        !$acc update host(rho)
        CALL dump_r(rho,TRIM(fname))
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'westpp')
  !
  DEALLOCATE(r)
  DEALLOCATE(dr)
  !$acc exit data delete(rho,rho_aux,dvg_exc)
  DEALLOCATE(rho)
  DEALLOCATE(rho_aux)
  IF(.NOT. gamma_only) THEN
     !$acc exit data delete(psic_aux)
     DEALLOCATE(psic_aux)
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
