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
SUBROUTINE wbse_bse_kernel(current_spin, nbndval_k, evc1, bse_kd1)
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
SUBROUTINE bse_kernel_finite_field_gamma(current_spin, nbndval_k, evc1, bse_kd1)
  !
  USE kinds,                 ONLY : DP
  USE fft_base,              ONLY : dffts
  USE wavefunctions,         ONLY : psic
  USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,double_invfft_gamma,&
                                  & double_fwfft_gamma
  USE mp_global,             ONLY : inter_image_comm
  USE pwcom,                 ONLY : npw,npwx,nks,isk,ngk
  USE westcom,               ONLY : l_lanczos,nbnd_occ,nbndval0x,l_use_localise_repr,u_matrix,&
                                  & index_matrix_lz,size_index_matrix_lz
  USE distribution_center,   ONLY : aband,bseparal
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
  ! local vars
  !
  INTEGER :: ibnd, jbnd, ibnd_1, size_index_matrix, summ_index
  INTEGER :: ibnd_index, jbnd_index, il1, ig1
  INTEGER :: nbndval_q, current_spin_ikq, ikq
  INTEGER :: nbvalloc
  !
  REAL(DP), ALLOCATABLE :: raux1(:), raux2(:)
  COMPLEX(DP), ALLOCATABLE :: aux_bse1(:,:), aux_bse2(:,:), kd1_ij(:,:), gaux(:)
  !
  IF(l_use_localise_repr) THEN
     ALLOCATE(aux_bse1(npwx,nbndval0x))
  ENDIF
  !
  DO ikq = 1, nks
     !
     current_spin_ikq = isk(ikq)
     !
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     nbndval_q = nbnd_occ(ikq)
     !
     npw = ngk(ikq)
     !
     IF(l_use_localise_repr) THEN
        aux_bse1(:,:) = (0._DP,0._DP)
        !
        DO ibnd = 1, nbndval0x
           DO jbnd = 1, nbndval0x
              aux_bse1(:,jbnd) = aux_bse1(:,jbnd) + u_matrix(ibnd,jbnd,current_spin) * evc1(:,ibnd,ikq)
           ENDDO
        ENDDO
     ENDIF
     !
     size_index_matrix = size_index_matrix_lz(current_spin)
     !
     ALLOCATE(aux_bse2(npwx,nbndval0x))
     aux_bse2 = (0._DP,0._DP)
     !
     IF(.NOT. l_lanczos) THEN
        ALLOCATE(kd1_ij(npwx,size_index_matrix))
        kd1_ij(:,:) = (0._DP,0._DP)
     ENDIF
     !
     ALLOCATE(raux1(dffts%nnr))
     !
     DO il1 = 1, bseparal%nlocx
        !
        ig1 = bseparal%l2g(il1) ! global index of n_total
        !
        IF(ig1 < 1 .OR. ig1 > size_index_matrix) CYCLE
        !
        ibnd_index = index_matrix_lz(ig1,1,current_spin)
        jbnd_index = index_matrix_lz(ig1,2,current_spin)
        !
        IF(l_lanczos) THEN
           !
           ! READ response at iq,ik,ispin
           !
           CALL read_bse_pots_g2r(raux1,ibnd_index,jbnd_index,current_spin)
           !
           IF(l_use_localise_repr) THEN
              CALL single_invfft_gamma(dffts,npw,npwx,aux_bse1(:,jbnd_index),psic,'Wave')
           ELSE
              CALL single_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),psic,'Wave')
           ENDIF
           !
           psic(:) = CMPLX(REAL(psic,KIND=DP)*raux1,KIND=DP)
           !
           ALLOCATE(gaux(npwx))
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,gaux,'Wave')
           !
           aux_bse2(:,ibnd_index) = aux_bse2(:,ibnd_index) + gaux(:)
           !
           DEALLOCATE(gaux)
           !
        ELSE
           !
           ALLOCATE(gaux(npwx))
           !
           CALL read_bse_pots_g2g(gaux, ibnd_index, jbnd_index, current_spin)
           !
           kd1_ij(:,ig1) = gaux(:)
           !
           DEALLOCATE(gaux)
           !
        ENDIF
        !
     ENDDO
     !
     IF(.NOT. l_lanczos) THEN
        CALL mp_sum(kd1_ij, inter_image_comm)
     ELSE
        CALL mp_sum(aux_bse2, inter_image_comm)
     ENDIF
     !
     DEALLOCATE(raux1)
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
           ALLOCATE(raux1(dffts%nnr), raux2(dffts%nnr))
           !
           raux1(:) = 0._DP
           raux2(:) = 0._DP
           summ_index = 0
           !
           ! LOOP OVER BANDS AT QPOINT
           !
           DO ig1 = 1, size_index_matrix
              !
              IF(summ_index > MIN(size_index_matrix, (nbndval_q+nbndval_k))) CYCLE
              !
              ibnd_index = index_matrix_lz(ig1,1,current_spin)
              jbnd_index = index_matrix_lz(ig1,2,current_spin)
              !
              IF(ibnd_index == ibnd) THEN
                 summ_index = summ_index+1
                 psic(:) = (0._DP,0._DP)
                 !
                 IF(l_use_localise_repr) THEN
                    CALL double_invfft_gamma(dffts,npw,npwx,aux_bse1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
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
                 IF(l_use_localise_repr) THEN
                    CALL double_invfft_gamma(dffts,npw,npwx,aux_bse1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
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
           DEALLOCATE(raux1, raux2)
           !
           IF(il1 < nbvalloc) THEN
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_bse2(:,ibnd),aux_bse2(:,ibnd_1),'Wave')
           ELSE
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_bse2(:,ibnd),'Wave')
           ENDIF
           !
        ENDDO
        !
        DEALLOCATE(kd1_ij)
        !
        CALL mp_sum(aux_bse2, inter_image_comm)
        !
     ENDIF
     !
     IF(l_use_localise_repr) THEN
        !
        ! U^{+}(\xi)
        !
        aux_bse1 = (0._DP,0._DP)
        !
        DO ibnd = 1, nbndval0x
           DO jbnd = 1, nbndval0x
              aux_bse1(:,jbnd) = aux_bse1(:,jbnd) + CONJG(u_matrix(jbnd,ibnd,current_spin)) * aux_bse2(:,ibnd)
           ENDDO
        ENDDO
        !
        bse_kd1(:,:) =  bse_kd1 - aux_bse1
        !
     ELSE
        !
        bse_kd1(:,:) =  bse_kd1 - aux_bse2
        !
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE(aux_bse2)
  !
  IF(l_use_localise_repr) THEN
     DEALLOCATE(aux_bse1)
  ENDIF
  !
END SUBROUTINE
