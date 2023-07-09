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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE west_apply_liouvillian(evc1,evc1_new,sf)
  !-----------------------------------------------------------------------
  !
  ! Applies the linear response operator to response wavefunctions
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE gvect,                ONLY : gstart
  USE uspp,                 ONLY : vkb,nkb
  USE lsda_mod,             ONLY : nspin
  USE pwcom,                ONLY : npw,npwx,current_k,current_spin,isk,lsda,xk,ngk,igk_k,nbnd
  USE control_flags,        ONLY : gamma_only
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : inter_image_comm,my_image_id,nbgrp,my_bgrp_id
  USE noncollin_module,     ONLY : npol
  USE buffers,              ONLY : get_buffer
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE westcom,              ONLY : l_bse,l_qp_correction,l_bse_triplet,sigma_c_head,sigma_x_head,&
                                 & nbnd_occ,scissor_ope,n_trunc_bands,et_qp,lrwfc,iuwfc,&
                                 & l_hybrid_tddft,l_spin_flip_kernel
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE uspp_init,            ONLY : init_us_2
  USE exx,                  ONLY : exxalfa
  USE wbse_dv,              ONLY : wbse_dv_of_drho,wbse_dv_of_drho_sf
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE wvfct_gpum,           ONLY : et=>et_d
  USE west_gpu,             ONLY : factors,dvrs,hevc1,reallocate_ps_gpu
  USE cublas
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
  USE wvfct,                ONLY : et
#endif
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: evc1(npwx*npol,band_group%nlocx,kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: evc1_new(npwx*npol,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: sf
  !
  ! Workspace
  !
  LOGICAL :: lrpa,do_k1e
  INTEGER :: ibnd,jbnd,iks,iks_do,ir,ig,nbndval,flnbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  REAL(DP) :: factor
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: factors(:)
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: hevc1(:,:)
#endif
  INTEGER, PARAMETER :: flks(2) = [2,1]
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
  ALLOCATE(factors(band_group%nlocx))
  ALLOCATE(hevc1(npwx*npol,band_group%nlocx))
  ALLOCATE(dvrs(dffts%nnr,nspin))
#endif
  !
  ! Calculation of the charge density response
  !
  CALL wbse_calc_dens(evc1,dvrs,sf)
  !
#if !defined(__NCCL)
  !$acc update device(dvrs)
#endif
  !
  lrpa = l_bse
  !
  If(sf .AND. l_spin_flip_kernel) THEN
     CALL wbse_dv_of_drho_sf(dvrs)
  ELSE
     CALL wbse_dv_of_drho(dvrs,lrpa,.FALSE.)
  ENDIF
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(sf) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     flnbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
#if defined(__CUDA)
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     IF(l_bse_triplet) THEN
        do_k1e = .FALSE.
     ELSEIF(sf .AND. (.NOT. l_spin_flip_kernel)) THEN
        do_k1e = .FALSE.
     ELSEIF(sf .AND. l_spin_flip_kernel) THEN
        do_k1e = .TRUE.
     ELSE
        do_k1e = .TRUE.
     ENDIF
     !
     IF(do_k1e) THEN
        !
        IF(gamma_only) THEN
           !
           ! double bands @ gamma
           !
           DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
              !
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
              jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
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
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,evc1_new(:,lbnd,iks),evc1_new(:,lbnd+1,iks),'Wave')
              !$acc end host_data
              !
           ENDDO
           !
           ! single band @ gamma
           !
           IF(MOD(nbnd_do,2) == 1) THEN
              !
              lbnd = nbnd_do
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
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
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,evc1_new(:,lbnd,iks),'Wave')
              !$acc end host_data
              !
           ENDIF
           !
        ELSE
           !
           ! only single bands
           !
           DO lbnd = 1,nbnd_do
              !
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
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
              CALL single_fwfft_k(dffts,npw,npwx,psic,evc1_new(:,lbnd,iks),'Wave',igk_k(:,current_k))
              !$acc end host_data
              !
           ENDDO
           !
           IF(npol == 2) THEN
              DO lbnd = 1,nbnd_do
                 !
                 ibnd = band_group%l2g(lbnd)+n_trunc_bands
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
                 CALL single_fwfft_k(dffts,npw,npwx,psic,evc1_new(npwx+1:npwx*2,lbnd,iks),'Wave',igk_k(:,current_k))
                 !$acc end host_data
                 !
              ENDDO
           ENDIF
           !
        ENDIF
        !
     ENDIF
     !
     ! use h_psi_, i.e. h_psi without band parallelization, as west
     ! handles band parallelization separately
     !
#if defined(__CUDA)
     !$acc host_data use_device(evc1,hevc1)
     CALL h_psi__gpu(npwx,npw,nbnd_do,evc1(:,:,iks),hevc1)
     !$acc end host_data
#else
     CALL h_psi_(npwx,npw,nbnd_do,evc1(:,:,iks),hevc1)
#endif
     !
     IF(l_qp_correction) THEN
#if defined(__CUDA)
        CALL reallocate_ps_gpu(nbnd,nbnd_do)
#endif
        CALL apply_hqp_to_m_wfcs(iks,nbnd_do,evc1(:,:,iks),hevc1)
     ENDIF
     !
     ! Subtract the eigenvalues
     !
     IF(l_bse) THEN
        factor = -scissor_ope+sigma_x_head+sigma_c_head
     ELSEIF(l_hybrid_tddft) THEN
        factor = -scissor_ope+sigma_x_head*exxalfa
     ELSE
        factor = -scissor_ope
     ENDIF
     !
     IF(l_qp_correction) THEN
        !
        !$acc parallel loop present(factors,et_qp)
        DO lbnd = 1,nbnd_do
           !
           ! ibnd = band_group%l2g(lbnd)+n_trunc_bands
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1+n_trunc_bands
           !
           factors(lbnd) = et_qp(ibnd,iks_do)+factor
           !
        ENDDO
        !$acc end parallel
        !
     ELSE
        !
        !$acc parallel loop present(factors)
        DO lbnd = 1,nbnd_do
           !
           ! ibnd = band_group%l2g(lbnd)+n_trunc_bands
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1+n_trunc_bands
           !
           factors(lbnd) = et(ibnd,iks_do)+factor
           !
        ENDDO
        !$acc end parallel
        !
     ENDIF
     !
     !$acc parallel loop collapse(2) present(evc1_new,hevc1,factors,evc1)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           evc1_new(ig,lbnd,iks) = evc1_new(ig,lbnd,iks)+hevc1(ig,lbnd)-factors(lbnd)*evc1(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(l_bse .OR. l_hybrid_tddft) CALL bse_kernel_gamma(current_spin,evc1,evc1_new(:,:,iks),sf)
     !
     IF(gamma_only) THEN
        IF(gstart == 2) THEN
           !$acc parallel loop present(evc1_new)
           DO lbnd = 1,nbnd_do
              evc1_new(1,lbnd,iks) = CMPLX(REAL(evc1_new(1,lbnd,iks),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
     ENDIF
     !
     ! Pc[k]*evc1_new(k)
     !
     ! load evc from iks to apply Pc of the current spin channel
     !
     IF(kpt_pool%nloc > 1) THEN
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
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,evc1_new(:,:,iks),(1._DP,0._DP))
     !
  ENDDO
  !
#if !defined(__CUDA)
  DEALLOCATE(factors)
  DEALLOCATE(dvrs)
  DEALLOCATE(hevc1)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('apply_lv')
#else
  CALL stop_clock('apply_lv')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE west_apply_liouvillian_btda(evc1,evc1_new,sf)
  !-----------------------------------------------------------------------
  !
  ! Applies the linear response operator to response wavefunctions
  ! beyond the Tamm-Dancoff approximation (btda)
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE gvect,                ONLY : gstart
  USE lsda_mod,             ONLY : nspin
  USE pwcom,                ONLY : npw,npwx,current_k,current_spin,isk,lsda,ngk
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE noncollin_module,     ONLY : npol
  USE buffers,              ONLY : get_buffer
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE westcom,              ONLY : l_bse,l_bse_triplet,nbnd_occ,n_trunc_bands,lrwfc,iuwfc,&
                                 & l_hybrid_tddft,l_spin_flip_kernel
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE wbse_dv,              ONLY : wbse_dv_of_drho,wbse_dv_of_drho_sf
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : dvrs,evc2_new=>hevc1,reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: evc1(npwx*npol,band_group%nlocx,kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: evc1_new(npwx*npol,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: sf
  !
  ! Workspace
  !
  LOGICAL :: lrpa,do_k2e
  INTEGER :: ibnd,jbnd,iks,iks_do,ir,ig,nbndval,flnbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  INTEGER, PARAMETER :: flks(2) = [2,1]
#if !defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc2_new(:,:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('apply_lv_btda')
#else
  CALL start_clock('apply_lv_btda')
#endif
  !
  dffts_nnr = dffts%nnr
  !
#if !defined(__CUDA)
  ALLOCATE(evc2_new(npwx*npol,band_group%nlocx))
  ALLOCATE(dvrs(dffts%nnr,nspin))
#endif
  !
  ! Calculation of the charge density response
  !
  CALL wbse_calc_dens(evc1,dvrs,sf)
  !
#if !defined(__NCCL)
  !$acc update device(dvrs)
#endif
  !
  lrpa = l_bse
  !
  If(sf .AND. l_spin_flip_kernel) THEN
     CALL wbse_dv_of_drho_sf(dvrs)
  ELSE
     CALL wbse_dv_of_drho(dvrs,lrpa,.FALSE.)
  ENDIF
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(sf) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     flnbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     IF(l_bse_triplet) THEN
        do_k2e = .FALSE.
     ELSEIF(sf .AND. (.NOT. l_spin_flip_kernel)) THEN
        do_k2e = .FALSE.
     ELSEIF(sf .AND. l_spin_flip_kernel) THEN
        do_k2e = .TRUE.
     ELSE
        do_k2e = .TRUE.
     ENDIF
     !
     IF(do_k2e) THEN
        !
        ! double bands @ gamma
        !
        DO lbnd = 1, nbnd_do-MOD(nbnd_do,2),2
           !
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),evc_work(:,jbnd),psic,'Wave')
           !
           !$acc parallel loop present(dvrs)
           DO ir = 1,dffts_nnr
              psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(evc2_new)
           CALL double_fwfft_gamma(dffts,npw,npwx,psic,evc2_new(:,lbnd),evc2_new(:,lbnd+1),'Wave')
           !$acc end host_data
           !
        ENDDO
        !
        ! single band @ gamma
        !
        IF(MOD(nbnd_do,2) == 1) THEN
           !
           lbnd = nbnd_do
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
           !
           !$acc parallel loop present(dvrs)
           DO ir = 1,dffts_nnr
              psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(evc2_new)
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,evc2_new(:,lbnd),'Wave')
           !$acc end host_data
           !
        ENDIF
        !
     ENDIF
     !
     ! The other part beyond TDA. exx_div treatment is not needed for this part.
     !
     IF(l_bse) CALL errore('west_apply_liouvillian_btda', 'BSE forces not implemented', 1)
     !
     ! K2d part. exx_div treatment is not needed for this part.
     !
     IF(l_hybrid_tddft) CALL hybrid_kernel_term2(current_spin,evc1,evc2_new,sf)
     !
     IF(gstart == 2) THEN
        !$acc parallel loop present(evc2_new)
        DO lbnd = 1,nbnd_do
           evc2_new(1,lbnd) = CMPLX(REAL(evc2_new(1,lbnd),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
     ENDIF
     !
     ! Pc[k]*evc2_new(k)
     !
     ! load evc from iks to apply Pc of the current spin channel
     !
     IF(kpt_pool%nloc > 1) THEN
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
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,evc2_new,(1._DP,0._DP))
     !
     !$acc parallel loop collapse(2) present(evc1_new,evc2_new)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           evc1_new(ig,lbnd,iks) = evc1_new(ig,lbnd,iks)+evc2_new(ig,lbnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
#if !defined(__CUDA)
  DEALLOCATE(dvrs)
  DEALLOCATE(evc2_new)
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('apply_lv_btda')
#else
  CALL stop_clock('apply_lv_btda')
#endif
  !
END SUBROUTINE
