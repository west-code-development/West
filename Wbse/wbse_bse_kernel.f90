!
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
SUBROUTINE wbse_bse_kernel(iks, current_spin, nbndval_k, evc1, bse_kd1)
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE pwcom,                ONLY : npwx,nks
  USE westcom,              ONLY : nbndval0x, l_davidson,l_lanzcos
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN)       :: iks, current_spin, nbndval_k
  COMPLEX(DP),INTENT(IN)    :: evc1(npwx,nbndval0x,nks)
  COMPLEX(DP),INTENT(INOUT) :: bse_kd1(npwx,nbndval0x)
  !
  IF (gamma_only) THEN
     !
     IF (l_lanzcos) THEN
        !
        CALL bse_kernel_finite_field_gamma (iks, current_spin, nbndval_k, evc1, bse_kd1, .true.)
        !
     ELSE
        !
        CALL bse_kernel_finite_field_gamma (iks, current_spin, nbndval_k, evc1, bse_kd1, .false.)
        !
     ENDIF
     !
  ELSE
     !
     STOP
     !
  ENDIF
  !
  RETURN
  !
ENDSUBROUTINE wbse_bse_kernel
!
!
!
SUBROUTINE bse_kernel_finite_field_gamma (iks, current_spin, nbndval_k, evc1, bse_kd1, lz_method)
  !
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dfftp,dffts
  USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,&
                                    double_invfft_gamma,double_fwfft_gamma
  USE mp_global,             ONLY : ibnd_start,ibnd_end,intra_bgrp_comm,&
                                    world_comm,inter_image_comm,inter_bgrp_comm
  USE pwcom,                 ONLY : npw,npwx,igk_k,nks,isk
  USE westcom,               ONLY : nbnd_occ, nbndval0x
  USE bse_module,            ONLY : ovl_matrix,u_matrix,ovl_thr,&
                                    l_wannier_repr,index_matrix_lz,&
                                    size_index_matrix_lz
  USE distribution_center,   ONLY : aband,bseparal
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN)       :: iks, current_spin, nbndval_k
  LOGICAL, INTENT(IN)       :: lz_method
  COMPLEX(DP),INTENT(IN)    :: evc1(npwx,nbndval0x,nks)
  COMPLEX(DP),INTENT(INOUT) :: bse_kd1(npwx,nbndval0x)
  !
  ! local vars
  !
  INTEGER  :: ibnd, jbnd, ibnd_1, ir, size_index_matrix, summ_index
  INTEGER  :: ibnd_index, jbnd_index, il1, ig1
  INTEGER  :: nbndval_q, current_spin_ikq, ikq
  INTEGER  :: nbvalloc
  !
  REAL(DP),ALLOCATABLE    :: raux1(:), raux2(:)
  COMPLEX(DP),ALLOCATABLE :: aux_bse1(:,:), aux_bse2(:,:), kd1_ij(:,:)
  COMPLEX(DP),ALLOCATABLE :: caux(:), gaux(:)
  !
  IF (l_wannier_repr) THEN
     !
     ALLOCATE(aux_bse1(npwx,nbndval0x))
     !
  ENDIF
  !
  DO ikq = 1, nks
     !
     current_spin_ikq = isk(ikq)
     !
     IF (current_spin_ikq .NE. current_spin) CYCLE
     !
     nbndval_q = nbnd_occ(ikq)
     !
     IF (l_wannier_repr) THEN
        !
        aux_bse1(:,:) = (0.0_DP, 0.0_DP)
        !
        DO ibnd = 1, nbndval0x
           !
           DO jbnd = 1, nbndval0x
              !
              aux_bse1(:,jbnd) = aux_bse1(:,jbnd) + u_matrix(ibnd,jbnd,current_spin) * evc1(:,ibnd,ikq)
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
     ! parallel loop
     !
     size_index_matrix = size_index_matrix_lz(current_spin)
     !
     ALLOCATE(aux_bse2(npwx,nbndval0x))
     aux_bse2 = (0.0_DP, 0.0_DP)
     !
     IF (.not.lz_method) THEN
        !
        ALLOCATE(kd1_ij(npwx,size_index_matrix))
        !
        kd1_ij(:,:) = (0.0_DP, 0.0_DP)
        !
     ENDIF
     !
     ALLOCATE(caux(dffts%nnr))
     !
     DO il1 = 1, bseparal%nlocx
        !
        ig1 = bseparal%l2g(il1) ! global index of n_total
        !
        IF ((ig1 < 1) .or. (ig1 > size_index_matrix)) CYCLE
        !
        ibnd_index = INT(index_matrix_lz(ig1,1,current_spin))
        jbnd_index = INT(index_matrix_lz(ig1,2,current_spin))
        !
        IF (lz_method) THEN
           !
           ! READ response at iq,ik,ispin
           !
           caux(:) = (0.0_DP,0.0_DP)
           CALL read_bse_pots_g2r (caux(:), ibnd_index, jbnd_index, current_spin, .true.)
           !
           psic(:) = (0.0_DP,0.0_DP)
           !
           IF (l_wannier_repr) THEN
              !
              CALL single_invfft_gamma(dffts,npw,npwx,aux_bse1(:,jbnd_index),psic(:),'Wave')
              !
           ELSE
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),psic(:),'Wave')
              !
           ENDIF
           !
           psic(:) = CMPLX(REAL(psic(:),KIND=DP)*REAL(caux(:),KIND=DP),KIND=DP)
           !
           ALLOCATE (gaux(npwx))
           !
           gaux(:) = (0.0_DP, 0.0_DP)
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,psic(:),gaux(:),'Wave')
           !
           aux_bse2(:,ibnd_index) = aux_bse2(:,ibnd_index) + gaux(:)
           !
           DEALLOCATE (gaux)
           !
        ELSE
           !
           ALLOCATE (gaux(npwx))
           !
           CALL read_bse_pots_g2g(gaux, ibnd_index, jbnd_index, current_spin, .true.)
           !
           kd1_ij(:,ig1) = gaux(:)
           !
           DEALLOCATE (gaux)
           !
        ENDIF
        !
     ENDDO
     !
     IF (.not.lz_method) THEN
        !
        CALL mp_sum(kd1_ij(:,:), inter_image_comm)
        !
     ELSE
        !
        CALL mp_sum(aux_bse2(:,:), inter_image_comm)
        !
     ENDIF
     !
     DEALLOCATE(caux)
     !
     ! LOOP OVER BANDS AT KPOINT
     !
     IF (.NOT. lz_method) THEN
        !
        !   Davidson method
        !
        nbvalloc = aband%nloc
        ibnd     = 0
        ibnd_1   = 0
        !
        DO il1 = 1, nbvalloc, 2
           !
           ibnd   = aband%l2g(il1)
           ibnd_1 = aband%l2g(il1+1)
           !
           ALLOCATE(raux1(dffts%nnr), raux2(dffts%nnr))
           !
           raux1(:) = 0.0_DP
           raux2(:) = 0.0_DP
           !
           summ_index = 0
           !
           ! LOOP OVER BANDS AT QPOINT
           !
           DO ig1 = 1, size_index_matrix
              !
              IF (summ_index > MIN(size_index_matrix, (nbndval_q+nbndval_k))) CYCLE
              !
              ibnd_index = INT(index_matrix_lz(ig1,1,current_spin))
              jbnd_index = INT(index_matrix_lz(ig1,2,current_spin))
              !
              IF (ibnd_index == ibnd) THEN
                 !
                 summ_index = summ_index+1
                 !
                 psic(:)   = (0.0_DP,0.0_DP)
                 !
                 IF (l_wannier_repr) THEN
                    !
                    CALL double_invfft_gamma(dffts,npw,npwx,aux_bse1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
                    !
                 ELSE
                    !
                    CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),kd1_ij(:,ig1),psic,'Wave')
                    !
                 ENDIF
                 !
                 raux1(:) = raux1(:) + REAL(psic(:),KIND=DP) * AIMAG(psic(:))
                 !
              ENDIF
              !
              IF (ibnd_index == ibnd_1) THEN
                 !
                 summ_index = summ_index+1
                 !
                 psic(:)   = (0.0_DP,0.0_DP)
                 !
                 IF (l_wannier_repr) THEN
                    !
                    CALL double_invfft_gamma(dffts,npw,npwx,aux_bse1(:,jbnd_index),kd1_ij(:,ig1),psic,'Wave')
                    !
                 ELSE
                    !
                    CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,jbnd_index,ikq),kd1_ij(:,ig1),psic,'Wave')
                    !
                 ENDIF
                 !
                 raux2(:) = raux2(:) + REAL(psic(:),KIND=DP) * AIMAG(psic(:))
                 !
              ENDIF
              !
           ENDDO
           !
           !   Back to reciprocal space
           !
           psic(:) = (0.0_DP, 0.0_DP)
           IF (il1 < nbvalloc) THEN
              !
              psic(:) = CMPLX(raux1(:),raux2(:),KIND=DP)
              !
           ELSE
              !
              psic(:) = CMPLX(raux1(:),KIND=DP)
              !
           ENDIF
           !
           DEALLOCATE(raux1, raux2)
           !
           IF (il1 < nbvalloc) THEN
              !
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_bse2(:,ibnd),aux_bse2(:,ibnd_1),'Wave')
              !
           ELSE
              !
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_bse2(:,ibnd),'Wave')
              !
           ENDIF
           !
        ENDDO
        !
        DEALLOCATE(kd1_ij)
        !
        CALL mp_sum (aux_bse2(:,:), inter_bgrp_comm)
        !
     ENDIF
     !
     IF (l_wannier_repr) THEN
        !
        ! U^{+}(\xi)
        !
        aux_bse1 = (0.0_DP, 0.0_DP)
        !
        DO ibnd = 1, nbndval0x
           !
           DO jbnd = 1, nbndval0x
              !
              aux_bse1(:,jbnd) = aux_bse1(:,jbnd) + CONJG(u_matrix(jbnd,ibnd,current_spin)) * aux_bse2(:,ibnd)
              !
           ENDDO
           !
        ENDDO
        !
        bse_kd1(:,:) =  bse_kd1(:,:) - aux_bse1(:,:)
        !
     ELSE
        !
        bse_kd1(:,:) =  bse_kd1(:,:) - aux_bse2(:,:)
        !
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE(aux_bse2)
  !
  IF (l_wannier_repr) THEN
     !
     DEALLOCATE (aux_bse1)
     !
  ENDIF
  !
  RETURN
  !
ENDSUBROUTINE
