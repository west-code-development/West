!
! Copyright (C) 2015-2022 M. Govoni
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
!---------------------------------------------------------------------
SUBROUTINE wbse_calc_dens(devc, drho)
  !---------------------------------------------------------------------
  !
  ! This subroutine calculates the response charge density
  ! from linear response orbitals and ground state orbitals.
  !
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE fft_base,               ONLY : dffts
  USE lsda_mod,               ONLY : nspin,lsda
  USE noncollin_module,       ONLY : npol
  USE pwcom,                  ONLY : npw,npwx,igk_k,current_k,nks,current_spin,isk,wg,ngk
  USE control_flags,          ONLY : gamma_only
  USE mp,                     ONLY : mp_sum,mp_bcast
  USE mp_global,              ONLY : my_image_id,inter_image_comm,inter_bgrp_comm
  USE buffers,                ONLY : get_buffer
  USE westcom,                ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_lanczos
  USE fft_at_gamma,           ONLY : double_invfft_gamma
  USE fft_at_k,               ONLY : single_invfft_k
  USE distribution_center,    ONLY : aband
#if defined(__CUDA)
  USE wavefunctions_gpum,     ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,          ONLY : evc_host=>evc
  USE west_gpu,               ONLY : tmp_r,tmp_c,psic2
#else
  USE wavefunctions,          ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: devc(npwx*npol,nbndval0x-n_trunc_bands,nks)
  COMPLEX(DP), INTENT(OUT) :: drho(dffts%nnr,nspin)
  !
  ! Workspace
  !
  INTEGER :: ir, ibnd, ibnd_g, iks, nbndval, lbnd, dffts_nnr
  REAL(DP) :: w1
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: tmp_r(:)
  COMPLEX(DP), ALLOCATABLE :: tmp_c(:)
  COMPLEX(DP), ALLOCATABLE :: psic2(:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('calc_dens')
#else
  CALL start_clock('calc_dens')
#endif
  !
  dffts_nnr = dffts%nnr
  !
#if !defined(__CUDA)
  IF(gamma_only) THEN
     ALLOCATE(tmp_r(dffts%nnr))
  ELSE
     ALLOCATE(tmp_c(dffts%nnr))
     ALLOCATE(psic2(dffts%nnr))
  ENDIF
#endif
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     nbndval = nbnd_occ(iks)
     !
     ! ... Set k-point and spin
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read GS wavefunctions
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
     IF(gamma_only) THEN
        !
        !$acc kernels present(tmp_r)
        tmp_r(:) = 0._DP
        !$acc end kernels
        !
        ! double bands @ gamma
        !
        DO lbnd = 1, aband%nloc
           !
           ibnd = aband%l2g(lbnd)
           ibnd_g = ibnd+n_trunc_bands
           IF(ibnd_g < 1 .OR. ibnd_g > nbndval) CYCLE
           !
           w1 = wg(ibnd_g,iks)/omega
           !
           !$acc host_data use_device(devc)
           CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd_g),devc(:,ibnd,iks),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(tmp_r)
           DO ir = 1, dffts_nnr
              tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        !$acc update host(tmp_r)
        !
        drho(:,current_spin) = CMPLX(tmp_r,KIND=DP)
        !
     ELSE
        !
        !$acc kernels present(tmp_c)
        tmp_c(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        ! only single bands
        !
        DO lbnd = 1, aband%nloc
           !
           ibnd = aband%l2g(lbnd)
           ibnd_g = ibnd+n_trunc_bands
           IF(ibnd_g < 1 .OR. ibnd_g > nbndval) CYCLE
           !
           w1 = wg(ibnd_g,iks)/omega
           !
           CALL single_invfft_k(dffts,npw,npwx,evc_work(:,ibnd_g),psic,'Wave',igk_k(:,current_k))
           !$acc host_data use_device(devc,psic2)
           CALL single_invfft_k(dffts,npw,npwx,devc(:,ibnd,iks),psic2,'Wave',igk_k(:,current_k))
           !$acc end host_data
           !
           !$acc parallel loop present(tmp_c,psic2)
           DO ir = 1, dffts_nnr
              tmp_c(ir) = tmp_c(ir) + w1*CONJG(psic(ir))*psic2(ir)
           ENDDO
           !$acc end parallel
           !
           IF(npol == 2) THEN
              !
              CALL single_invfft_k(dffts,npw,npwx,evc_work(npwx+1:npwx*2,ibnd_g),psic,'Wave',igk_k(:,current_k))
              !$acc host_data use_device(devc,psic2)
              CALL single_invfft_k(dffts,npw,npwx,devc(npwx+1:npwx*2,ibnd,iks),psic2,'Wave',igk_k(:,current_k))
              !$acc end host_data
              !
              !$acc parallel loop present(tmp_c,psic2)
              DO ir = 1, dffts_nnr
                 tmp_c(ir) = tmp_c(ir) + w1*CONJG(psic(ir))*psic2(ir)
              ENDDO
              !$acc end parallel
              !
           ENDIF
           !
        ENDDO
        !
        !$acc update host(tmp_c)
        !
        drho(:,current_spin) = tmp_c
        !
     ENDIF
     !
  ENDDO
  !
  IF(l_lanczos) THEN
     CALL mp_sum(drho,inter_image_comm)
  ELSE
     CALL mp_sum(drho,inter_bgrp_comm)
  ENDIF
  !
  !$acc update device(drho)
  !
#if !defined(__CUDA)
  IF(gamma_only) THEN
     DEALLOCATE(tmp_r)
  ELSE
     DEALLOCATE(tmp_c)
     DEALLOCATE(psic2)
  ENDIF
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('calc_dens')
#else
  CALL stop_clock('calc_dens')
#endif
  !
END SUBROUTINE
