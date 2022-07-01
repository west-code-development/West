
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
SUBROUTINE calc_exx2(sigma_exx)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_exx(n,iks) = < n,iks | V_exx | n,iks >     n = qp_bands(1):qp_bands(SIZE(qp_bands))
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
  USE wavefunctions,        ONLY : evc,psic,psic_nc
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,occupation,qp_bands
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE class_idistribute,    ONLY : idistribute
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(OUT) :: sigma_exx(SIZE(qp_bands),k_grid%nps)
  !
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  COMPLEX(DP), ALLOCATABLE :: evckmq(:,:),phase(:)
  REAL(DP), EXTERNAL :: DDOT
  INTEGER :: ib,iv,ir,iks,ik,is,ig,ivloc,ibloc,iq,ikqs,ikq,iks_g,ib_index
  INTEGER :: nbndval
  INTEGER :: npwkq
  TYPE(idistribute) :: vband
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  LOGICAL :: l_gammaq
  REAL(DP) :: g0(3),peso
  !
  WRITE(stdout,'(5x,a)') ''
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') '(X)-Sigma'
  CALL io_push_bar()
  !
  ALLOCATE(pertg(ngm))
  IF (gamma_only) THEN
     peso = 2._DP
  ELSE
     peso = 1._DP
     ALLOCATE(phase(dffts%nnr))
     ALLOCATE(evckmq(npwx*npol,nbnd))
  ENDIF
  IF(noncolin) THEN
     ALLOCATE(pertr_nc(dffts%nnr,npol))
  ELSE
     ALLOCATE(pertr(dffts%nnr))
  ENDIF
  !
  !
  ! Set to zero
  !
  sigma_exx = 0._DP
  !
  CALL band_group%init(SIZE(qp_bands),'b','band_group',.FALSE.)
  !
  barra_load = kpt_pool%nloc*band_group%nloc
  CALL start_bar_type(barra,'sigmax',barra_load)
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
     ENDIF
     !
     DO ibloc = 1,band_group%nloc
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index)
        !
        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ib),psic,'Wave')
        ELSEIF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic_nc(1,1),'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ib),psic_nc(1,2),'Wave',igk_k(1,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic,'Wave',igk_k(1,current_k))
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
           ENDIF
           !
           vband = idistribute()
           CALL vband%init(nbndval,'i','nbndval',.FALSE.)
           !
           DO ivloc = 1,vband%nloc
              !
              iv = vband%l2g(ivloc)
              !
              ! Bring it to R-space
              IF(gamma_only) THEN
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(1,iv),pertr,'Wave')
                 DO ir = 1,dffts%nnr
                    pertr(ir) = psic(ir)*pertr(ir)
                 ENDDO
                 CALL single_fwfft_gamma(dffts,ngm,ngm,pertr,pertg,'Rho')
              ELSEIF(noncolin) THEN
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1     ,iv),pertr_nc(1,1),'Wave',igk_k(1,ikqs))
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1+npwx,iv),pertr_nc(1,2),'Wave',igk_k(1,ikqs))
                 DO ir = 1,dffts%nnr
                    pertr_nc(ir,1) = CONJG(pertr_nc(ir,1)*phase(ir))*psic_nc(ir,1)+CONJG(pertr_nc(ir,2)*phase(ir))*psic_nc(ir,2)
                 ENDDO
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr_nc(1,1),pertg,'Rho') ! no igk
              ELSE
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1,iv),pertr,'Wave',igk_k(1,ikqs))
                 DO ir = 1,dffts%nnr
                    pertr(ir) = CONJG(pertr(ir)*phase(ir))*psic(ir)
                 ENDDO
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr,pertg,'Rho') ! no igk
              ENDIF
              !
              DO ig = 1,ngm
                 pertg(ig) = pertg(ig)*pot3D%sqvc(ig)
              ENDDO
              sigma_exx(ib_index,iks_g) = sigma_exx(ib_index,iks_g) - &
              & occupation(iv,iks)*peso*DDOT(2*ngm,pertg(1),1,pertg(1),1)/omega*q_grid%weight(iq)
              IF(ib == iv .AND. gstart == 2 .AND. l_gammaq) &
              sigma_exx(ib_index,iks_g) = sigma_exx(ib_index,iks_g) - occupation(iv,iks)*pot3D%div
              !
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
     DEALLOCATE(phase)
     DEALLOCATE(evckmq)
  ENDIF
  !
END SUBROUTINE
