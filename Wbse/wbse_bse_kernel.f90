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
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_bse_kernel(current_spin, nbndval_k, evc1, bse_kd1)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE pwcom,                ONLY : npwx,nks
  USE westcom,              ONLY : nbndval0x
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin, nbndval_k
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,nbndval0x,nks)
  COMPLEX(DP), INTENT(INOUT) :: bse_kd1(npwx,nbndval0x)
  !
  IF(gamma_only) THEN
     CALL bse_kernel_finite_field_gamma(current_spin, nbndval_k, evc1, bse_kd1)
  ELSE
     CALL errore('wbse_bse_kernel','Only Gamma is supported',1)
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE bse_kernel_finite_field_gamma(current_spin, nbndval_k, evc1, bse_kd1)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE fft_base,              ONLY : dffts
  USE wavefunctions,         ONLY : psic
  USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,double_invfft_gamma,&
                                  & double_fwfft_gamma
  USE mp_global,             ONLY : inter_image_comm,inter_bgrp_comm
  USE pwcom,                 ONLY : npw,npwx,nks,isk,ngk
  USE westcom,               ONLY : l_lanczos,nbnd_occ,nbndval0x,l_local_repr,u_matrix,idx_matrix,&
                                  & n_bse_idx
  USE distribution_center,   ONLY : aband,bandpair
  USE wbse_io,               ONLY : read_bse_pots_g2r,read_bse_pots_g2g
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin, nbndval_k
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,nbndval0x,nks)
  COMPLEX(DP), INTENT(INOUT) :: bse_kd1(npwx,nbndval0x)
  !
  ! Workspace
  !
  INTEGER :: ibnd, jbnd, ibnd_1, do_idx, summ_index
  INTEGER :: ibnd_index, jbnd_index, il1, ig1
  INTEGER :: nbndval_q, current_spin_ikq, ikq
  INTEGER :: nbvalloc
  !
  REAL(DP), ALLOCATABLE :: raux1(:), raux2(:)
  COMPLEX(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:), kd1_ij(:,:), gaux(:)
  !
  IF(l_local_repr) THEN
     ALLOCATE(tmp1(npwx,nbndval0x))
  ENDIF
  ALLOCATE(tmp2(npwx,nbndval0x))
  ALLOCATE(raux1(dffts%nnr))
  IF(.NOT. l_lanczos) THEN
     ALLOCATE(raux2(dffts%nnr))
  ENDIF
  ALLOCATE(gaux(npwx))
  !
  DO ikq = 1, nks
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     nbndval_q = nbnd_occ(ikq)
     npw = ngk(ikq)
     !
     IF(l_local_repr) THEN
        tmp1(:,:) = (0._DP,0._DP)
        DO ibnd = 1, nbndval0x
           DO jbnd = 1, nbndval0x
              tmp1(:,jbnd) = tmp1(:,jbnd) + u_matrix(ibnd,jbnd,current_spin) * evc1(:,ibnd,ikq)
           ENDDO
        ENDDO
     ENDIF
     !
     do_idx = n_bse_idx(current_spin)
     !
     IF(.NOT. l_lanczos) THEN
        ALLOCATE(kd1_ij(npwx,do_idx))
        kd1_ij(:,:) = (0._DP,0._DP)
     ENDIF
     !
     tmp2 = (0._DP,0._DP)
     !
     DO il1 = 1, bandpair%nlocx
        !
        ig1 = bandpair%l2g(il1) ! global index of n_total
        !
        IF(ig1 < 1 .OR. ig1 > do_idx) CYCLE
        !
        ibnd_index = idx_matrix(ig1,1,current_spin)
        jbnd_index = idx_matrix(ig1,2,current_spin)
        !
        IF(l_lanczos) THEN
           !
           ! READ response at iq,ik,ispin
           !
           CALL read_bse_pots_g2r(raux1,ibnd_index,jbnd_index,current_spin)
           !
           IF(l_local_repr) THEN
              CALL single_invfft_gamma(dffts,npw,npwx,tmp1(:,jbnd_index),psic,'Wave')
           ELSE
              CALL single_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),psic,'Wave')
           ENDIF
           !
           psic(:) = CMPLX(REAL(psic,KIND=DP)*raux1,KIND=DP)
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,gaux,'Wave')
           !
           tmp2(:,ibnd_index) = tmp2(:,ibnd_index) + gaux(:)
           !
        ELSE
           !
           CALL read_bse_pots_g2g(gaux,ibnd_index,jbnd_index,current_spin)
           !
           kd1_ij(:,ig1) = gaux(:)
           !
        ENDIF
        !
     ENDDO
     !
     IF(.NOT. l_lanczos) THEN
        CALL mp_sum(kd1_ij,inter_image_comm)
     ELSE
        CALL mp_sum(tmp2,inter_image_comm)
     ENDIF
     !
     ! LOOP OVER BANDS AT KPOINT
     !
     IF(.NOT. l_lanczos) THEN
        !
        ! Davidson method
        !
        nbvalloc = aband%nloc
        ibnd = 0
        ibnd_1 = 0
        !
        DO il1 = 1, nbvalloc, 2
           !
           ibnd = aband%l2g(il1)
           ibnd_1 = aband%l2g(il1+1)
           !
           raux1(:) = 0._DP
           raux2(:) = 0._DP
           summ_index = 0
           !
           ! LOOP OVER BANDS AT QPOINT
           !
           DO ig1 = 1, do_idx
              !
              IF(summ_index > MIN(do_idx,(nbndval_q+nbndval_k))) CYCLE
              !
              ibnd_index = idx_matrix(ig1,1,current_spin)
              jbnd_index = idx_matrix(ig1,2,current_spin)
              !
              IF(ibnd_index == ibnd) THEN
                 summ_index = summ_index+1
                 psic(:) = (0._DP,0._DP)
                 !
                 IF(l_local_repr) THEN
                    CALL double_invfft_gamma(dffts,npw,npwx,tmp1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
                 ELSE
                    CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),kd1_ij(:,ig1),psic,'Wave')
                 ENDIF
                 !
                 raux1(:) = raux1 + REAL(psic,KIND=DP) * AIMAG(psic)
              ENDIF
              !
              IF(ibnd_index == ibnd_1) THEN
                 summ_index = summ_index+1
                 psic(:) = (0._DP,0._DP)
                 !
                 IF(l_local_repr) THEN
                    CALL double_invfft_gamma(dffts,npw,npwx,tmp1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
                 ELSE
                    CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),kd1_ij(:,ig1),psic,'Wave')
                 ENDIF
                 !
                 raux2(:) = raux2 + REAL(psic,KIND=DP) * AIMAG(psic)
              ENDIF
              !
           ENDDO
           !
           ! Back to reciprocal space
           !
           IF(il1 < nbvalloc) THEN
              psic(:) = CMPLX(raux1,raux2,KIND=DP)
           ELSE
              psic(:) = CMPLX(raux1,KIND=DP)
           ENDIF
           !
           IF(il1 < nbvalloc) THEN
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,tmp2(:,ibnd),tmp2(:,ibnd_1),'Wave')
           ELSE
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,tmp2(:,ibnd),'Wave')
           ENDIF
           !
        ENDDO
        !
        DEALLOCATE(kd1_ij)
        !
        CALL mp_sum(tmp2,inter_bgrp_comm)
        !
     ENDIF
     !
     IF(l_local_repr) THEN
        !
        ! U^{+}(\xi)
        !
        tmp1 = (0._DP,0._DP)
        !
        DO ibnd = 1, nbndval0x
           DO jbnd = 1, nbndval0x
              tmp1(:,jbnd) = tmp1(:,jbnd) + CONJG(u_matrix(jbnd,ibnd,current_spin)) * tmp2(:,ibnd)
           ENDDO
        ENDDO
        !
        bse_kd1(:,:) =  bse_kd1 - tmp1
        !
     ELSE
        !
        bse_kd1(:,:) =  bse_kd1 - tmp2
        !
     ENDIF
     !
  ENDDO
  !
  IF(l_local_repr) THEN
     DEALLOCATE(tmp1)
  ENDIF
  DEALLOCATE(tmp2)
  DEALLOCATE(raux1)
  IF(.NOT. l_lanczos) THEN
     DEALLOCATE(raux2)
  ENDIF
  DEALLOCATE(gaux)
  !
END SUBROUTINE
