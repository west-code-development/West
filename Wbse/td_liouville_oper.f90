!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE west_apply_liouvillian(evc1,evc1_new)
  !
  ! Applies the linear response operator to response wavefunctions
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE gvect,                ONLY : gstart
  USE uspp,                 ONLY : vkb,nkb
  USE lsda_mod,             ONLY : nspin
  USE pwcom,                ONLY : npw,npwx,et,current_k,current_spin,isk,lsda,nks,xk,ngk,igk_k,nbnd
  USE control_flags,        ONLY : gamma_only
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE mp_global,            ONLY : nimage,my_image_id,inter_image_comm,nbgrp,my_bgrp_id,&
                                 & inter_bgrp_comm
  USE noncollin_module,     ONLY : npol
  USE buffers,              ONLY : get_buffer
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE westcom,              ONLY : lrwfc,iuwfc,nbnd_occ,scissor_ope,nbndval0x,l_bse_calculation,&
                                 & l_qp_correction,l_bse_triplet,l_lanczos,sigma_c_head,&
                                 & sigma_x_head,et_qp
  USE distribution_center,  ONLY : aband
  USE uspp_init,            ONLY : init_us_2
  USE wbse_dv,              ONLY : wbse_dv_of_drho
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,             ONLY : dvrs,hevc1,evc1_loc,reallocate_ps_gpu
  USE cublas
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: evc1(npwx*npol,nbndval0x,nks)
  COMPLEX(DP), INTENT(OUT) :: evc1_new(npwx*npol,nbndval0x,nks)
  !
  ! Local variables
  !
  INTEGER :: ibnd,jbnd,iks,ir,ig,nbndval,nbvalloc,il1
  INTEGER :: dffts_nnr
  COMPLEX(DP) :: factor
#if !defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: hevc1(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc1_loc(:,:)
  !$acc declare device_resident(hevc1,evc1_loc)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('apply_lv')
#else
  CALL start_clock('apply_lv')
#endif
  !
  dffts_nnr = dffts%nnr
  !
  !$acc kernels present(evc1_new)
  evc1_new(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
#if !defined(__CUDA)
  ALLOCATE(hevc1(npwx*npol,aband%nloc))
  ALLOCATE(evc1_loc(npwx*npol,aband%nloc))
  ALLOCATE(dvrs(dffts%nnr,nspin))
#endif
  !
  ! Calculation of the charge density response
  !
  CALL wbse_calc_dens(evc1,dvrs)
  !
  IF(l_bse_calculation) THEN
     CALL wbse_dv_of_drho(dvrs,.TRUE.,.FALSE.)
  ELSE
     CALL wbse_dv_of_drho(dvrs,.FALSE.,.FALSE.)
  ENDIF
  !
  !$acc update device(dvrs)
  !
  DO iks = 1,nks
     !
     nbndval = nbnd_occ(iks)
     nbvalloc = 0
     DO il1 = 1,aband%nloc
        ibnd = aband%l2g(il1)
        IF(ibnd < 1 .OR. ibnd > nbndval) CYCLE
        nbvalloc = nbvalloc+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
#if defined(__CUDA)
     CALL g2_kin_gpu(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read in GS wavefunctions iks
     !
     IF(nks > 1) THEN
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
#if defined(__CUDA)
     !
     ! ... Sync GPU
     !
     CALL using_becp_auto(2)
     CALL using_becp_d_auto(0)
#endif
     !
     IF(.NOT. l_bse_triplet) THEN
        !
        IF(gamma_only) THEN
           !
           ! double bands @ gamma
           !
           DO il1 = 1,nbvalloc-MOD(nbvalloc,2),2
              ibnd = aband%l2g(il1)
              jbnd = aband%l2g(il1+1)
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),evc_work(:,jbnd),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(evc1_new)
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,evc1_new(:,ibnd,iks),evc1_new(:,jbnd,iks),'Wave')
              !$acc end host_data
           ENDDO
           !
           ! single band @ gamma
           !
           IF(MOD(nbvalloc,2) == 1) THEN
              ibnd = aband%l2g(nbvalloc)
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(evc1_new)
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,evc1_new(:,ibnd,iks),'Wave')
              !$acc end host_data
           ENDIF
           !
        ELSE
           !
           ! only single bands
           !
           DO il1 = 1,nbvalloc
              ibnd = aband%l2g(il1)
              !
              CALL single_invfft_k(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave',igk_k(:,current_k))
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*dvrs(ir,current_spin)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(evc1_new)
              CALL single_fwfft_k(dffts,npw,npwx,psic,evc1_new(:,ibnd,iks),'Wave',igk_k(:,current_k))
              !$acc end host_data
           ENDDO
           !
           IF(npol == 2) THEN
              DO il1 = 1,nbvalloc
                 ibnd = aband%l2g(il1)
                 !
                 CALL single_invfft_k(dffts,npw,npwx,evc_work(npwx+1:npwx*2,ibnd),psic,'Wave',igk_k(:,current_k))
                 !
                 !$acc parallel loop present(dvrs)
                 DO ir = 1,dffts_nnr
                    psic(ir) = psic(ir)*dvrs(ir,current_spin)
                 ENDDO
                 !$acc end parallel
                 !
                 !$acc host_data use_device(evc1_new)
                 CALL single_fwfft_k(dffts,npw,npwx,psic,evc1_new(npwx+1:npwx*2,ibnd,iks),'Wave',igk_k(:,current_k))
                 !$acc end host_data
              ENDDO
           ENDIF
           !
        ENDIF
        !
     ENDIF
     !
     !$acc parallel loop collapse(2) present(evc1_loc,evc1)
     DO il1 = 1,nbvalloc
        !
        ! ibnd = aband%l2g(il1)
        !
        DO ig = 1,npw
           IF(l_lanczos) THEN
              ibnd = nimage*(il1-1)+my_image_id+1
           ELSE
              ibnd = nbgrp*(il1-1)+my_bgrp_id+1
           ENDIF
           evc1_loc(ig,il1) = evc1(ig,ibnd,iks)
        ENDDO
        !
     ENDDO
     !$acc end parallel
     !
     ! use h_psi_, i.e. h_psi without band parallelization, as west
     ! handles band parallelization separately
     !
#if defined(__CUDA)
     !$acc host_data use_device(evc1_loc,hevc1)
     CALL h_psi__gpu(npwx,npw,nbvalloc,evc1_loc,hevc1)
     !$acc end host_data
#else
     CALL h_psi_(npwx,npw,nbvalloc,evc1_loc,hevc1)
#endif
     !
     IF(l_qp_correction) THEN
#if defined(__CUDA)
        CALL reallocate_ps_gpu(nbnd,nbvalloc)
#endif
        CALL apply_hqp_to_m_wfcs(iks,nbvalloc,evc1_loc,hevc1)
     ENDIF
     !
     ! Subtract the eigenvalues
     !
     DO il1 = 1,nbvalloc
        !
        ibnd = aband%l2g(il1)
        !
        IF(l_qp_correction) THEN
           IF(l_bse_calculation) THEN
              factor = CMPLX(-(et_qp(ibnd,iks)-scissor_ope+sigma_x_head+sigma_c_head),KIND=DP)
           ELSE
              factor = CMPLX(-(et_qp(ibnd,iks)-scissor_ope),KIND=DP)
           ENDIF
        ELSE
           IF(l_bse_calculation) THEN
              factor = CMPLX(-(et(ibnd,iks)-scissor_ope+sigma_x_head+sigma_c_head),KIND=DP)
           ELSE
              factor = CMPLX(-(et(ibnd,iks)-scissor_ope),KIND=DP)
           ENDIF
        ENDIF
        !
        !$acc host_data use_device(evc1,hevc1)
        CALL ZAXPY(npw,factor,evc1(:,ibnd,iks),1,hevc1(:,il1),1)
        !$acc end host_data
        !
     ENDDO
     !
     !$acc parallel loop collapse(2) present(evc1_new,hevc1)
     DO il1 = 1,nbvalloc
        !
        ! ibnd = aband%l2g(il1)
        !
        DO ig = 1,npw
           IF(l_lanczos) THEN
              ibnd = nimage*(il1-1)+my_image_id+1
           ELSE
              ibnd = nbgrp*(il1-1)+my_bgrp_id+1
           ENDIF
           evc1_new(ig,ibnd,iks) = evc1_new(ig,ibnd,iks)+hevc1(ig,il1)
        ENDDO
        !
     ENDDO
     !$acc end parallel
     !
     !$acc update host(evc1_new)
     IF(l_lanczos) THEN
        CALL mp_sum(evc1_new(:,:,iks),inter_image_comm)
     ELSE
        CALL mp_sum(evc1_new(:,:,iks),inter_bgrp_comm)
     ENDIF
     !$acc update device(evc1_new)
     !
     IF(l_bse_calculation) THEN
        CALL wbse_bse_kernel(current_spin,nbndval,evc1,evc1_new(:,:,iks))
     ENDIF
     !
     IF(gamma_only) THEN
        IF(gstart == 2) THEN
           !$acc parallel loop present(evc1_new)
           DO ibnd = 1,nbndval0x
              evc1_new(1,ibnd,iks) = CMPLX(REAL(evc1_new(1,ibnd,iks),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
     ENDIF
     !
     ! Pc[k]*evc1_new(k)
     !
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,evc1_new(:,:,iks),(1._DP,0._DP))
     !
  ENDDO
  !
#if !defined(__CUDA)
  DEALLOCATE(dvrs)
  DEALLOCATE(hevc1)
  DEALLOCATE(evc1_loc)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('apply_lv')
#else
  CALL stop_clock('apply_lv')
#endif
  !
END SUBROUTINE
