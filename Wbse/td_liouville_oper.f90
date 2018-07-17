!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE west_apply_liouvillian(evc1, evc1_new)
  !
  ! Applies the linear response operator to response wavefunctions
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts, dfftp
  USE gvect,                ONLY : gstart
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE lsda_mod,             ONLY : nspin
  USE wavefunctions_module, ONLY : psic, evc
  USE pwcom,                ONLY : npw, current_k, current_spin, isk, lsda
  USE wvfct,                ONLY : nbnd, npwx, g2kin, et
  USE control_flags,        ONLY : gamma_only
  USE mp,                   ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_global,            ONLY : intra_bgrp_comm,my_image_id,inter_image_comm,&
                                   inter_bgrp_comm
  USE noncollin_module,     ONLY : npol
  USE buffers,              ONLY : get_buffer
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,&
                                   double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE wbsecom,              ONLY : l_diag_term_only,scissor_ope,nbndval0x
  USE wbsecom,              ONLY : l_qp_correction
  USE bse_module,           ONLY : bse_calc,et_qp
  USE westcom,              ONLY : isz,lrwfc,iuwfc,nbnd_occ
  USE wbsecom,              ONLY : l_bse_triplet, l_lanzcos, mac_isz
  USE distribution_center,  ONLY : aband
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: evc1(npwx*npol, nbndval0x, nks)
  COMPLEX(DP), INTENT(OUT) :: evc1_new(npwx*npol, nbndval0x, nks)
  !
  ! Local variables
  !
  INTEGER  :: ibnd, ibnd_1, iks, ir, nbndval
  INTEGER  :: nbvalloc, il1
  REAL(DP) :: scissor
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux_r(:)
  COMPLEX(DP), ALLOCATABLE :: hevc1(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc1_aux(:,:)
  !
  scissor = scissor_ope 
  !
  CALL start_clock('west_apply_liouvillian')
  !
  evc1_new(:,:,:)  = (0.0d0,0.0d0)
  !
  IF (.NOT.ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  ALLOCATE(dvrs(dfftp%nnr,nspin))
  !
  ! Calculation of the charge density response 
  !
  dvrs(:,:) = (0.0_DP, 0.0_DP)
  CALL wbse_calc_dens(evc1, dvrs)
  !
  IF (bse_calc) THEN
     !
     CALL west_dv_of_drho(dvrs, .true., .false.)
     !  
  ELSE
     !
     CALL west_dv_of_drho(dvrs, .false., .false.)
     !  
  ENDIF
  !
  DO iks = 1, nks
     !
     nbndval = nbnd_occ(iks)
     ! 
     nbvalloc = 0
     DO il1 = 1, aband%nloc
        ibnd = aband%l2g(il1)
        IF( ibnd < 1 .OR. ibnd > nbndval ) CYCLE
        nbvalloc = nbvalloc + 1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     !
     CALL g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
     !
     ! ... read in GS wavefunctions iks 
     !
     IF (nks>1) THEN
        !
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        CALL mp_bcast(evc,0,inter_image_comm)
        !
     ENDIF
     ! 
     IF (l_diag_term_only) GOTO 112
     !
     IF (l_bse_triplet)    GOTO 112
     !
     IF (gamma_only) THEN
        !
        ! double bands @ gamma
        !
        DO il1=1, nbvalloc - MOD(nbvalloc,2), 2
           !
           ibnd   = aband%l2g(il1) 
           ibnd_1 = aband%l2g(il1+1) 
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,ibnd_1),psic,'Wave')
           !
           DO ir=1,dffts%nnr
              ! 
              psic(ir) = psic(ir) * CMPLX(DBLE(dvrs(ir,current_spin)), 0.0_DP)
              !
           ENDDO
           !
           CALL double_fwfft_gamma(dffts,npw,npwx,psic,evc1_new(1,ibnd,iks),evc1_new(1,ibnd_1,iks),'Wave')
           !
        ENDDO
        ! 
        ! single band @ gamma
        ! 
        IF ( MOD(nbvalloc,2) == 1 ) THEN
           !
           ibnd=aband%l2g(nbvalloc)
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),psic,'Wave')
           !
           DO ir=1,dffts%nnr
              !
              psic(ir) = CMPLX( REAL(psic(ir),KIND=DP) * REAL(dvrs(ir,current_spin),KIND=DP), 0._DP, KIND=DP)
              !
           ENDDO
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,evc1_new(1,ibnd,iks),'Wave') 
           !
        ENDIF
        !
     ELSE
        !
        ! only single bands
        !   
        DO il1=1, nbvalloc 
           !
           ibnd = aband%l2g(il1) 
           ! 
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,current_k))
           !
           DO ir=1, dffts%nnr
              !   
              psic(ir) = psic(ir) * dvrs(ir,current_spin) 
              !   
           ENDDO
           !
           CALL single_fwfft_k(dffts,npw,npwx,psic,evc1_new(1,ibnd,iks),'Wave',igk_k(1,current_k)) 
           !
        ENDDO
        !
        IF (npol==2) THEN
           !
           DO il1=1, nbvalloc 
              !
              ibnd = aband%l2g(il1) 
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ibnd),psic,'Wave',igk_k(1,current_k))
              !
              DO ir=1, dffts%nnr
                 ! 
                 psic(ir) = psic(ir) * dvrs(ir,current_spin)
                 !
              ENDDO
              ! 
              CALL single_fwfft_k(dffts,npw,npwx,psic,evc1_new(npwx+1,ibnd,iks),'Wave',igk_k(1,current_k)) 
              !
           ENDDO
           ! 
        ENDIF
        !
     ENDIF   
     !
  112 CONTINUE 
     !
     ALLOCATE(hevc1(npwx*npol, nbvalloc))
     ALLOCATE(evc1_aux(npwx*npol, nbvalloc))
     !
     hevc1(:,:)    = (0.0_DP, 0.0_DP) 
     evc1_aux(:,:) = (0.0_DP, 0.0_DP) 
     !
     DO il1 = 1, nbvalloc
        !
        ibnd = aband%l2g(il1) 
        !
        evc1_aux(:,il1) = evc1(:,ibnd,iks)
        !
     ENDDO
     !
     CALL h_psi(npwx,npw,nbvalloc,evc1_aux(:,:),hevc1(:,:))
     !
     IF (l_qp_correction) THEN
        !  
        CALL bse_hqp_psi(iks, current_spin, nbvalloc, evc1_aux(:,:), hevc1(:,:))
        !
     ENDIF
     !
     DEALLOCATE(evc1_aux)
     !
     ! Subtract the eigenvalues
     ! 
     IF (l_qp_correction) THEN
        !
        DO il1 = 1, nbvalloc
           !
           ibnd = aband%l2g(il1) 
           ! 
           IF (bse_calc) THEN
              !
              CALL zaxpy(npw,CMPLX(-(et_qp(ibnd,iks) - scissor + isz + mac_isz),0.0d0,DP), &
                         evc1(:,ibnd,iks),1,hevc1(:,il1),1)
              !
           ELSE 
              !
              CALL zaxpy(npw,CMPLX(-(et_qp(ibnd,iks) - scissor),0.0d0,DP), &
                         evc1(:,ibnd,iks),1,hevc1(:,il1),1)
              !
           ENDIF
           !
        ENDDO
        !
     ELSE
        !   
        DO il1 = 1, nbvalloc
           !
           ibnd = aband%l2g(il1) 
           !
           IF (bse_calc) THEN
              ! 
              CALL zaxpy(npw,CMPLX(-(et(ibnd,iks) - scissor + isz + mac_isz),0.0d0,DP), &
                         evc1(:,ibnd,iks),1,hevc1(:,il1),1)
              !
           ELSE
              !
              CALL zaxpy(npw,CMPLX(-(et(ibnd,iks) - scissor),0.0d0,DP), &
                         evc1(:,ibnd,iks),1,hevc1(:,il1),1)
              !
           ENDIF
           !
        ENDDO
        !
     ENDIF
     !
     DO il1 = 1, nbvalloc
        !
        ibnd = aband%l2g(il1) 
        !
        evc1_new(:,ibnd,iks) = hevc1(:,il1) + evc1_new(:,ibnd,iks)
        !
     ENDDO
     !
     DEALLOCATE(hevc1)
     !
     IF (l_lanzcos) THEN
        ! 
        CALL mp_sum (evc1_new(:,:,iks),inter_image_comm)     
        !
     ELSE
        !
        CALL mp_sum (evc1_new(:,:,iks),inter_bgrp_comm)
        !
     ENDIF
     !
     IF (l_diag_term_only) GOTO 113
     !
     IF (bse_calc) THEN
        !
        CALL wbse_bse_kernel(iks, current_spin, nbndval, evc1(:,:,:), evc1_new(:,:,iks))
        ! 
     ENDIF
     !
  113 CONTINUE 
     !
     IF (gamma_only) THEN
        !
        IF (gstart==2) THEN
           !
           evc1_new(1,:,iks) = CMPLX(DBLE(evc1_new(1,:,iks)),0.0_DP)
           !
        ENDIF
        !
     ENDIF
     !
     ! Pc[k]*evc1_new(k)
     !
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,evc1_new(:,:,iks),(1.0_DP,0.0_DP))
     !
  ENDDO
  !  
  DEALLOCATE(dvrs)
  IF (ALLOCATED(psic)) DEALLOCATE(psic)
  !
  CALL stop_clock('west_apply_liouvillian')
  !
  RETURN
  !
ENDSUBROUTINE west_apply_liouvillian
