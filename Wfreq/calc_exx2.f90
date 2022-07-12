
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
!-----------------------------------------------------------------------
SUBROUTINE calc_exx2(sigma_exx,nb1,nb2)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_exx(n,iks) = < n,iks | V_exx | n,iks >     n = nb1, nb2
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout
  USE gvect,                ONLY : gstart,ngm
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : npw,npwx,nbnd,igk_k,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE class_idistribute,    ONLY : idistribute
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE types_coulomb,        ONLY : pot3D
#if defined(__CUDA)
  USE wavefunctions,        ONLY : evc_host=>evc
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d,psic_nc=>psic_nc_d
  USE west_gpu,             ONLY : allocate_exx_gpu,deallocate_exx_gpu,sqvc_d
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic,psic_nc
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nb1,nb2
  REAL(DP), INTENT(OUT) :: sigma_exx(nb1:nb2,k_grid%nps)
  !
  ! Workspace
  !
  LOGICAL :: l_gammaq
  INTEGER :: barra_load
  INTEGER :: ib,iv,ir,iks,ik,is,ig,ivloc,ibloc,iq,ikqs,ikq,iks_g
  INTEGER :: nbndval
  INTEGER :: npwkq
  INTEGER :: dffts_nnr
  REAL(DP) :: g0(3),peso
  REAL(DP) :: dot_tmp
  COMPLEX(DP), ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  !$acc declare device_resident(pertg,pertr,pertr_nc)
  COMPLEX(DP), ALLOCATABLE :: evckmq(:,:),phase(:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: evckmq,phase
#endif
  TYPE(idistribute) :: vband
  TYPE(bar_type) :: barra
  !
  WRITE(stdout,'(5x,a)') ''
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') '(X)-Sigma'
  CALL io_push_bar()
  !
  IF(gamma_only) THEN
     peso = 2._DP
  ELSE
     peso = 1._DP
     ALLOCATE(phase(dffts%nnr))
     ALLOCATE(evckmq(npwx*npol,nbnd))
     !$acc enter data create(phase,evckmq)
  ENDIF
  ALLOCATE(pertg(ngm))
  IF(noncolin) THEN
     ALLOCATE(pertr_nc(dffts%nnr,npol))
  ELSE
     ALLOCATE(pertr(dffts%nnr))
  ENDIF
  !
  ! Set to zero
  !
  sigma_exx = 0._DP
  !
  CALL band_group%init(nb2-nb1+1,'b','band_group',.FALSE.)
  !
  barra_load = kpt_pool%nloc*band_group%nloc
  CALL start_bar_type(barra,'sigmax',barra_load)
  !
#if defined(__CUDA)
  CALL allocate_exx_gpu()
#endif
  !
  dffts_nnr = dffts%nnr
  !
  ! LOOP
  !
  DO iks = 1,kpt_pool%nloc ! KPOINT-SPIN
     !
     iks_g = kpt_pool%l2g(iks)
     ik = k_grid%ip(iks_g)
     is = k_grid%is(iks_g)
     current_k = iks
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_host,0,inter_image_comm)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
#if defined(__CUDA)
     CALL using_evc(2)
     CALL using_evc_d(0)
#endif
     !
     DO ibloc = 1,band_group%nloc
        !
        ib = band_group%l2g(ibloc)+nb1-1
        !
        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ib),psic,'Wave')
        ELSEIF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc_work(1:npwx,ib),psic_nc(:,1),'Wave',igk_k(:,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc_work(npwx+1:npwx*2,ib),psic_nc(:,2),'Wave',igk_k(:,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc_work(:,ib),psic,'Wave',igk_k(:,current_k))
        ENDIF
        !
        DO iq = 1,q_grid%np
           !
           IF(gamma_only) THEN
              l_gammaq = .TRUE.
              CALL pot3D%init('Rho',.FALSE.,'gb')
              nbndval = nbnd_occ(iks)
           ELSE
              l_gammaq = q_grid%l_pIsGamma(iq)
              CALL pot3D%init('Rho',.FALSE.,'gb',iq)
              !
              CALL k_grid%find(k_grid%p_cart(:,ik)-q_grid%p_cart(:,iq),'cart',ikq,g0)
              ikqs = k_grid%ipis2ips(ikq,is)
              CALL compute_phase(g0,'cart',phase)
              !
              nbndval = nbnd_occ(ikqs)
              npwkq = ngk(ikqs)
              IF(my_image_id == 0) CALL get_buffer(evckmq,lrwfc,iuwfc,ikqs)
              CALL mp_bcast(evckmq,0,inter_image_comm)
              !
              !$acc update device(evckmq,phase)
           ENDIF
           !
#if defined(__CUDA)
           sqvc_d = pot3D%sqvc
#endif
           !
           vband = idistribute()
           CALL vband%init(nbndval,'i','nbndval',.FALSE.)
           !
           DO ivloc = 1,vband%nloc
              iv = vband%l2g(ivloc)
              !
              IF(gamma_only) THEN
                 !$acc host_data use_device(pertr)
                 CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,iv),pertr,'Wave')
                 !$acc end host_data
                 !$acc parallel loop
                 DO ir = 1,dffts_nnr
                    pertr(ir) = psic(ir)*pertr(ir)
                 ENDDO
                 !$acc end parallel
                 !$acc host_data use_device(pertr,pertg)
                 CALL single_fwfft_gamma(dffts,ngm,ngm,pertr,pertg,'Rho')
                 !$acc end host_data
              ELSEIF(noncolin) THEN
                 !$acc host_data use_device(evckmq,pertr_nc)
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1:npwx,iv),pertr_nc(:,1),'Wave',igk_k(:,ikqs))
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1+npwx:npwx*2,iv),pertr_nc(:,2),'Wave',igk_k(:,ikqs))
                 !$acc end host_data
                 !$acc parallel loop present(phase)
                 DO ir = 1,dffts_nnr
                    pertr_nc(ir,1) = CONJG(pertr_nc(ir,1)*phase(ir))*psic_nc(ir,1) &
                    & +CONJG(pertr_nc(ir,2)*phase(ir))*psic_nc(ir,2)
                 ENDDO
                 !$acc end parallel
                 !$acc host_data use_device(pertr_nc,pertg)
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr_nc(:,1),pertg,'Rho') ! no igk
                 !$acc end host_data
              ELSE
                 !$acc host_data use_device(evckmq,pertr)
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(:,iv),pertr,'Wave',igk_k(:,ikqs))
                 !$acc end host_data
                 !$acc parallel loop present(phase)
                 DO ir = 1,dffts_nnr
                    pertr(ir) = CONJG(pertr(ir)*phase(ir))*psic(ir)
                 ENDDO
                 !$acc end parallel
                 !$acc host_data use_device(pertr,pertg)
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr,pertg,'Rho') ! no igk
                 !$acc end host_data
              ENDIF
              !
              !$acc parallel loop
              DO ig = 1,ngm
#if defined(__CUDA)
                 pertg(ig) = pertg(ig)*sqvc_d(ig)
#else
                 pertg(ig) = pertg(ig)*pot3D%sqvc(ig)
#endif
              ENDDO
              !$acc end parallel
              !
              dot_tmp = 0._DP
              !$acc parallel loop reduction(+:dot_tmp) copy(dot_tmp)
              DO ig = 1,ngm
                 dot_tmp = dot_tmp+REAL(pertg(ig),KIND=DP)**2+AIMAG(pertg(ig))**2
              ENDDO
              !$acc end parallel
              !
              sigma_exx(ib,iks_g) = sigma_exx(ib,iks_g)-peso*dot_tmp/omega*q_grid%weight(iq)
              IF(ib == iv .AND. gstart == 2 .AND. l_gammaq) sigma_exx(ib,iks_g) = sigma_exx(ib,iks_g)-pot3D%div
           ENDDO ! ivloc
           !
        ENDDO ! iq
        !
        CALL update_bar_type(barra,'sigmax',1)
        !
     ENDDO ! ibloc
     !
  ENDDO ! iks
  !
#if defined(__CUDA)
  CALL deallocate_exx_gpu()
#endif
  !
  CALL stop_bar_type(barra,'sigmax')
  !
  CALL mp_sum(sigma_exx,intra_bgrp_comm)
  CALL mp_sum(sigma_exx,inter_bgrp_comm)
  CALL mp_sum(sigma_exx,inter_pool_comm)
  CALL mp_sum(sigma_exx,inter_image_comm)
  !
  DEALLOCATE(pertg)
  IF(noncolin) THEN
    DEALLOCATE(pertr_nc)
  ELSE
    DEALLOCATE(pertr)
  ENDIF
  IF(.NOT. gamma_only) THEN
     !$acc exit data delete(phase,evckmq)
     DEALLOCATE(phase)
     DEALLOCATE(evckmq)
  ENDIF
  !
END SUBROUTINE
