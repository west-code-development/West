
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
!-----------------------------------------------------------------------
SUBROUTINE calc_exx2(sigma_exx, l_QDET)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_exx(n,iks) = < qp_bands(n),iks | V_exx | qp_bands(n),iks >     n = 1,n_bands
  !
  ! IF (l_enable_off_diagonal .AND. l_full) store in
  ! sigma_exx_full(ijpmap(m,n),iks) = < qp_bands(m),iks | V_exx | qp_bands(n),iks >     n,m = 1,n_bands & m <= n
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
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,occupation,qp_bands,n_bands,l_enable_off_diagonal,&
                                 & sigma_exx_full,ijpmap
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE class_idistribute,    ONLY : idistribute
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE types_coulomb,        ONLY : pot3D
  USE wavefunctions,        ONLY : evc,psic,psic_nc
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(OUT) :: sigma_exx(n_bands,k_grid%nps)
  LOGICAL, INTENT(IN) :: l_QDET ! True if QDET double-counting term is calculated
  !
  ! Workspace
  !
  LOGICAL :: l_gammaq
  INTEGER :: barra_load
  INTEGER :: ib,iv,ir,iks,ik,is,ig,ivloc,ibloc,iq,ikqs,ikq,iks_g,ib_index,jb,jb_index,ipair
  INTEGER :: nbndval
  INTEGER :: npwkq
  INTEGER :: dffts_nnr
  REAL(DP) :: g0(3),peso
  REAL(DP) :: dot_tmp
  COMPLEX(DP) :: braket
  COMPLEX(DP), ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic1(:),pertg1(:),pertr1(:)
  COMPLEX(DP), ALLOCATABLE :: evckmq(:,:),phase(:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: evckmq
#endif
  TYPE(idistribute) :: vband
  TYPE(bar_type) :: barra
  !
  WRITE(stdout,*)
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') '(X)-Sigma'
  CALL io_push_bar()
  !
  IF(gamma_only) THEN
     peso = 2._DP
     IF(l_enable_off_diagonal) THEN
        ALLOCATE(psic1(dffts%nnr))
        ALLOCATE(pertr1(dffts%nnr))
        ALLOCATE(pertg1(ngm))
        !$acc enter data create(psic1,pertr1,pertg1)
     ENDIF
  ELSE
     peso = 1._DP
     ALLOCATE(phase(dffts%nnr))
     ALLOCATE(evckmq(npwx*npol,nbnd))
     !$acc enter data create(phase,evckmq)
  ENDIF
  ALLOCATE(pertg(ngm))
  !$acc enter data create(pertg)
  IF(noncolin) THEN
     ALLOCATE(pertr_nc(dffts%nnr,npol))
     !$acc enter data create(pertr_nc)
  ELSE
     ALLOCATE(pertr(dffts%nnr))
     !$acc enter data create(pertr)
  ENDIF
  !
  ! Set to zero
  !
  sigma_exx = 0._DP
  IF (l_enable_off_diagonal) sigma_exx_full = 0._DP
  !
  CALL band_group%init(n_bands,'b','band_group',.FALSE.)
  !
  barra_load = 0
  DO ibloc = 1,band_group%nloc
     ib_index = band_group%l2g(ibloc)
     !
     IF(l_enable_off_diagonal) THEN
        barra_load = barra_load+ib_index
     ELSE
        barra_load = barra_load+1
     ENDIF
  ENDDO
  barra_load = barra_load*kpt_pool%nloc
  !
  CALL start_bar_type(barra,'sigmax',barra_load)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
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
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     DO ibloc = 1,band_group%nloc
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index,is)
        !
        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ib),psic,'Wave')
        ELSEIF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1:npwx,ib),psic_nc(:,1),'Wave',igk_k(:,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1:npwx*2,ib),psic_nc(:,2),'Wave',igk_k(:,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(:,ib),psic,'Wave',igk_k(:,current_k))
        ENDIF
        !
        DO jb_index = 1,n_bands
           !
           jb = qp_bands(jb_index,is)
           !
           IF(l_enable_off_diagonal) THEN
              IF(jb > ib) CYCLE
              ipair = ijpmap(jb_index,ib_index)
              !
              IF(jb < ib) THEN
                  IF(gamma_only) THEN
                     CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jb),psic1,'Wave')
                  ENDIF
               ENDIF
           ELSE
              IF(jb /= ib) CYCLE
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
              !$acc enter data copyin(pot3D)
              !$acc enter data copyin(pot3D%sqvc)
              !
              vband = idistribute()
              CALL vband%init(nbndval,'i','nbndval',.FALSE.)
              !
              DO ivloc = 1,vband%nloc
                 !
                 iv = vband%l2g(ivloc)
                 !
                 ! for QDET double counting term, all states need to be within qp_bands
                 !
                 IF(l_QDET) THEN
                    IF(ALL(qp_bands(:,is) /= iv)) CYCLE
                 ENDIF
                 !
                 IF(gamma_only) THEN
                    CALL single_invfft_gamma(dffts,npw,npwx,evc(:,iv),pertr,'Wave')
                    !$acc parallel loop present(pertr1,psic1,pertr)
                    DO ir = 1,dffts_nnr
                       IF(l_enable_off_diagonal .AND. jb < ib) THEN
                          pertr1(ir) = psic1(ir)*pertr(ir)
                       ENDIF
                       pertr(ir) = psic(ir)*pertr(ir)
                    ENDDO
                    !$acc end parallel
                    CALL single_fwfft_gamma(dffts,ngm,ngm,pertr,pertg,'Rho')
                    IF(l_enable_off_diagonal .AND. jb < ib) THEN
                       CALL single_fwfft_gamma(dffts,ngm,ngm,pertr1,pertg1,'Rho')
                    ENDIF
                 ELSEIF(noncolin) THEN
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1:npwx,iv),pertr_nc(:,1),'Wave',igk_k(:,ikqs))
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1+npwx:npwx*2,iv),pertr_nc(:,2),'Wave',igk_k(:,ikqs))
                    !$acc parallel loop present(pertr_nc,phase)
                    DO ir = 1,dffts_nnr
                       pertr_nc(ir,1) = CONJG(pertr_nc(ir,1)*phase(ir))*psic_nc(ir,1) &
                       & +CONJG(pertr_nc(ir,2)*phase(ir))*psic_nc(ir,2)
                    ENDDO
                    !$acc end parallel
                    CALL single_fwfft_k(dffts,ngm,ngm,pertr_nc(:,1),pertg,'Rho') ! no igk
                 ELSE
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(:,iv),pertr,'Wave',igk_k(:,ikqs))
                    !$acc parallel loop present(pertr,phase)
                    DO ir = 1,dffts_nnr
                       pertr(ir) = CONJG(pertr(ir)*phase(ir))*psic(ir)
                    ENDDO
                    !$acc end parallel
                    CALL single_fwfft_k(dffts,ngm,ngm,pertr,pertg,'Rho') ! no igk
                 ENDIF
                 !
                 !$acc parallel loop present(pertg1,pertg,pot3D,pot3D%sqvc)
                 DO ig = 1,ngm
                    IF(l_enable_off_diagonal .AND. jb < ib) THEN
                       pertg1(ig) = pertg1(ig)*pot3D%sqvc(ig)
                    ENDIF
                    pertg(ig) = pertg(ig)*pot3D%sqvc(ig)
                 ENDDO
                 !$acc end parallel
                 !
                 IF(l_enable_off_diagonal .AND. jb < ib) THEN
                    dot_tmp = 0._DP
                    !$acc parallel loop reduction(+:dot_tmp) present(pertg,pertg1) copy(dot_tmp)
                    DO ig = 1,ngm
                       dot_tmp = dot_tmp+REAL(pertg(ig)*CONJG(pertg1(ig)),KIND=DP)
                    ENDDO
                    !$acc end parallel
                    !
                    braket = peso*dot_tmp/omega*q_grid%weight(iq)
                    sigma_exx_full(ipair,iks_g) = sigma_exx_full(ipair,iks_g)-occupation(iv,iks)*braket
                 ELSEIF(jb == ib) THEN
                    dot_tmp = 0._DP
                    !$acc parallel loop reduction(+:dot_tmp) present(pertg) copy(dot_tmp)
                    DO ig = 1,ngm
                       dot_tmp = dot_tmp+REAL(pertg(ig),KIND=DP)**2+AIMAG(pertg(ig))**2
                    ENDDO
                    !$acc end parallel
                    !
                    braket = peso*dot_tmp/omega*q_grid%weight(iq)
                    sigma_exx(ib_index,iks_g) = sigma_exx(ib_index,iks_g)-occupation(iv,iks)*braket
                    !
                    IF(l_enable_off_diagonal) &
                    & sigma_exx_full(ipair,iks_g) = sigma_exx_full(ipair,iks_g)-occupation(iv,iks)*braket
                    !
                    IF(ib == iv .AND. gstart == 2 .AND. l_gammaq) THEN
                       sigma_exx(ib_index,iks_g) = sigma_exx(ib_index,iks_g)-occupation(iv,iks)*pot3D%div
                       IF(l_enable_off_diagonal) &
                       & sigma_exx_full(ipair,iks_g) = sigma_exx_full(ipair,iks_g)-occupation(iv,iks)*pot3D%div
                    ENDIF
                 ENDIF
                 !
              ENDDO ! ivloc
              !
              !$acc exit data delete(pot3D%sqvc)
              !$acc exit data delete(pot3D)
              !
           ENDDO ! iq
           !
           CALL update_bar_type(barra,'sigmax',1)
           !
        ENDDO ! jb
        !
     ENDDO ! ibloc
     !
  ENDDO ! iks
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
  CALL stop_bar_type(barra,'sigmax')
  !
  CALL mp_sum(sigma_exx,intra_bgrp_comm)
  CALL mp_sum(sigma_exx,inter_bgrp_comm)
  CALL mp_sum(sigma_exx,inter_pool_comm)
  CALL mp_sum(sigma_exx,inter_image_comm)
  IF(l_enable_off_diagonal) THEN
     CALL mp_sum(sigma_exx_full,intra_bgrp_comm)
     CALL mp_sum(sigma_exx_full,inter_bgrp_comm)
     CALL mp_sum(sigma_exx_full,inter_pool_comm)
     CALL mp_sum(sigma_exx_full,inter_image_comm)
  ENDIF
  !
  IF(gamma_only) THEN
     IF(l_enable_off_diagonal) THEN
        !$acc exit data delete(psic1,pertr1,pertg1)
        DEALLOCATE(psic1)
        DEALLOCATE(pertr1)
        DEALLOCATE(pertg1)
     ENDIF
  ELSE
     !$acc exit data delete(phase,evckmq)
     DEALLOCATE(phase)
     DEALLOCATE(evckmq)
  ENDIF
  !$acc exit data delete(pertg)
  DEALLOCATE(pertg)
  IF(noncolin) THEN
     !$acc exit data delete(pertr_nc)
     DEALLOCATE(pertr_nc)
  ELSE
     !$acc exit data delete(pertr)
     DEALLOCATE(pertr)
  ENDIF
  !
END SUBROUTINE
