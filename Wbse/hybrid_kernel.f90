!
! Copyright (C) 2015-2016 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file: 
! Yu Jin
!
!===================================================================================

SUBROUTINE hybrid_kernel_term1_slow (current_spin, evc1, hybrid_kd1, sf)
  !
  ! \sum_{v'} (\int v_c \phi_{v'} \phi_{v}) a_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,& 
                                    double_invfft_gamma,double_fwfft_gamma
  USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
  USE mp_global,             ONLY : my_image_id,inter_bgrp_comm,inter_image_comm
  USE pwcom,                 ONLY : npw,npwx,igk_k,nks,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE control_flags,         ONLY : gamma_only
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : band_group,kpt_pool
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
  USE becmod_subs_gpum,      ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,              ONLY : dvrs,hevc1,reallocate_ps_gpu
  USE cublas
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  ! 
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN)       :: current_spin
  COMPLEX(DP),INTENT(IN)    :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN)       :: sf
  COMPLEX(DP),INTENT(INOUT) :: hybrid_kd1(npwx,band_group%nlocx)
  !
  ! local vars
  !
  INTEGER  :: lbnd, ibnd, ibndp, jbnd, jbndp, ibnd1, ir, ig, iks_do
  INTEGER  :: nbndval, current_spin_ikq, ikq, flnbndval
  INTEGER  :: nbnd_do
  !
  COMPLEX(DP),ALLOCATABLE :: aux_hybrid1(:,:)
  COMPLEX(DP),ALLOCATABLE :: caux(:), gaux(:), raux1(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  ALLOCATE(aux_hybrid1(npwx, nbndval0x-n_trunc_bands))
  aux_hybrid1 = (0.0_DP, 0.0_DP)
  !
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux1(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     !
     IF (current_spin_ikq .NE. current_spin) CYCLE  
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq
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
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        jbndp = jbnd + n_trunc_bands
        !
        raux1(:) = (0._DP, 0._DP)
        !
        DO lbnd = 1, nbnd_do ! index to be summed 
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd + n_trunc_bands
           !
           caux(:) = (0.0_DP, 0.0_DP)
           gaux(:) = (0.0_DP, 0.0_DP)
           !
           IF (gamma_only) THEN ! Only write the gamma_only case
              !
              ! product of evc and evc
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibndp),evc_work(:,jbndp),psic,'Wave')
              !
              DO ir = 1, dffts%nnr
                 !
                 caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir)),0._DP,KIND=DP)
                 !
              ENDDO
              !
           ENDIF
           !
           caux(:) = caux(:) / omega
           !
           ! Apply the bare Coulomb potential
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,caux(:),gaux(:),'Wave')
           !
           DO ig = 1, npw
              !
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
              !
           ENDDO
           !
           caux(:) = (0.0_DP, 0.0_DP)
           !
           CALL double_invfft_gamma(dffts,npw,npwx,gaux(:),evc1(:,lbnd,ikq),caux(:),'Wave')
           !
           psic(:) = (0.0_DP, 0.0_DP)
           !
           DO ir = 1, dffts%nnr
              !
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),0._DP,KIND=DP)
              !
           ENDDO
           !
           raux1(:) = raux1(:) + psic(:)
           !
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,raux1(:),aux_hybrid1(:,jbnd),'Wave')
        !
     ENDDO
     !
     CALL mp_sum (aux_hybrid1(:,:), inter_bgrp_comm)
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd)
        ibndp = ibnd + n_trunc_bands
        !
        hybrid_kd1(:,lbnd) = hybrid_kd1(:,lbnd) - aux_hybrid1(:,ibnd) * exxalfa
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux1)
  DEALLOCATE(aux_hybrid1)
  !
ENDSUBROUTINE


SUBROUTINE hybrid_kernel_term2 (current_spin, evc1, hybrid_kd2, sf)
  !
  ! \sum_{v'} (\int v_c a_{v'} \phi_{v}) \phi_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,&
                                    double_invfft_gamma,double_fwfft_gamma
  USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
  USE mp_global,             ONLY : my_image_id,inter_bgrp_comm,inter_image_comm
  USE pwcom,                 ONLY : npw,npwx,igk_k,nks,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE control_flags,         ONLY : gamma_only
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : band_group,kpt_pool
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
  USE becmod_subs_gpum,      ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,              ONLY : dvrs,hevc1,reallocate_ps_gpu
  USE cublas
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  ! 
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN)       :: current_spin
  COMPLEX(DP),INTENT(IN)    :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN)       :: sf
  COMPLEX(DP),INTENT(INOUT) :: hybrid_kd2(npwx,band_group%nlocx)
  !
  ! local vars
  !
  INTEGER  :: lbnd, ibnd, jbnd, ibndp, jbndp, ibnd1, ir, ig
  INTEGER  :: nbndval_q, current_spin_ikq, ikq, nbndval
  INTEGER  :: nbnd_do
  !
  COMPLEX(DP),ALLOCATABLE :: aux_hybrid2(:,:)
  COMPLEX(DP),ALLOCATABLE :: caux(:), gaux(:), raux1(:)
  !
  IF(sf) CALL errore('hybrid_kernel_gamma_term2', 'spin-flip is not supported', 1)
  !
  ALLOCATE(aux_hybrid2(npwx, nbndval0x-n_trunc_bands))
  aux_hybrid2 = (0._DP, 0._DP)
  !
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux1(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     !
     IF (current_spin_ikq .NE. current_spin) CYCLE  
     !
     nbndval = nbnd_occ(ikq)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        jbndp = jbnd + n_trunc_bands
        !
        raux1(:) = (0._DP, 0._DP)
        !
        DO lbnd = 1, nbnd_do ! index to be summed 
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd + n_trunc_bands
           !
           IF (ibndp > nbndval) CYCLE
           !
           caux(:) = (0.0_DP, 0.0_DP)
           gaux(:) = (0.0_DP, 0.0_DP)
           !
           IF (gamma_only) THEN
              !
              ! product of evc1 and evc
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,lbnd,ikq),evc_work(:,jbndp),psic,'Wave')
              !
              DO ir = 1, dffts%nnr
                 !
                 caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir)),0._DP,KIND=DP)
                 !
              ENDDO
              !
           ENDIF
           !
           caux(:) = caux(:)/omega
           !
           ! Apply the bare Coulomb potential
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,caux(:),gaux(:),'Wave')
           !
           DO ig = 1, npw
              !
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
              !
           ENDDO
           !
           caux(:) = (0._DP, 0._DP)
           !
           CALL double_invfft_gamma(dffts,npw,npwx,gaux(:),evc_work(:,ibndp),caux(:),'Wave')
           !
           psic(:) = (0._DP, 0._DP)
           !
           DO ir = 1, dffts%nnr
              !
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),0._DP,KIND=DP)
              !
           ENDDO
           !
           raux1(:) = raux1(:) + psic(:)
           !
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,raux1(:),aux_hybrid2(:,jbnd),'Wave')
        !
     ENDDO
     !
     CALL mp_sum (aux_hybrid2(:,:), inter_bgrp_comm)
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd)
        ibndp = ibnd + n_trunc_bands
        !
        hybrid_kd2(:,lbnd) = hybrid_kd2(:,lbnd) - aux_hybrid2(:,ibnd) * exxalfa
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux1)
  DEALLOCATE(aux_hybrid2)
  !
ENDSUBROUTINE


SUBROUTINE hybrid_kernel_term3 (current_spin, evc1, hybrid_kd3, sf)
  !
  ! \sum_{v'} (\int v_c a_{v'} \phi_{v}) a_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,&
                                    double_invfft_gamma,double_fwfft_gamma
  USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
  USE mp_global,             ONLY : my_image_id,inter_bgrp_comm,inter_image_comm
  USE pwcom,                 ONLY : npw,npwx,igk_k,nks,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE control_flags,         ONLY : gamma_only
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : band_group,kpt_pool
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
  USE becmod_subs_gpum,      ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,              ONLY : dvrs,hevc1,reallocate_ps_gpu
  USE cublas
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  ! 
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN)       :: current_spin
  COMPLEX(DP),INTENT(IN)    :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN)       :: sf
  COMPLEX(DP),INTENT(INOUT) :: hybrid_kd3(npwx,band_group%nlocx)
  !
  ! local vars
  !
  INTEGER  :: lbnd, ibnd, ibndp, jbnd, jbndp, ibnd1, ir, ig, iks_do, nbnd_do
  INTEGER  :: nbndval_q, current_spin_ikq, ikq, nbndval
  INTEGER  :: flnbndval
  !
  COMPLEX(DP),ALLOCATABLE :: aux_hybrid3(:,:)
  COMPLEX(DP),ALLOCATABLE :: caux(:), gaux(:), raux1(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  ALLOCATE(aux_hybrid3(npwx, nbndval0x-n_trunc_bands))
  aux_hybrid3 = (0._DP, 0._DP)
  !
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux1(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     !
     IF (current_spin_ikq .NE. current_spin) CYCLE  
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(ikq)
     flnbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq 
     !
     ! Note: evc in the current spin channel should be loaded
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        jbndp = jbnd + n_trunc_bands
        !
        raux1(:) = (0._DP, 0._DP)
        !
        DO lbnd = 1, nbnd_do
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd + n_trunc_bands
           !
           caux(:) = (0.0_DP, 0.0_DP)
           gaux(:) = (0.0_DP, 0.0_DP)
           !
           IF (gamma_only) THEN ! Only write the gamma_only case
              !
              ! product of evc1 and evc
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,lbnd,ikq),evc_work(:,jbndp),psic,'Wave')
              !
              DO ir = 1, dffts%nnr
                 !
                 caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir)),0._DP,KIND=DP)
                 !
              ENDDO
              !
           ENDIF
           !
           caux(:) = caux(:)/omega
           !
           ! Apply the bare Coulomb potential
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,caux(:),gaux(:),'Wave')
           !
           DO ig = 1, npw
              !
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
              !
           ENDDO
           !
           caux(:) = (0.0_DP, 0.0_DP)
           !
           CALL double_invfft_gamma(dffts,npw,npwx,gaux(:),evc1(:,lbnd,ikq),caux(:),'Wave')
           !
           psic(:) = (0.0_DP, 0.0_DP)
           !
           DO ir = 1, dffts%nnr
              !
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),0._DP,KIND=DP)
              !
           ENDDO
           !
           raux1(:) = raux1(:) + psic(:)
           !
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,raux1(:),aux_hybrid3(:,jbnd),'Wave')
        !
     ENDDO
     !
     CALL mp_sum (aux_hybrid3(:,:), inter_bgrp_comm)
     !
     ! Note: nbnd_do needs to be recomputed for the corrent spin channel
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd)
        ibndp = ibnd + n_trunc_bands
        !
        hybrid_kd3(:,lbnd) = hybrid_kd3(:,lbnd) - aux_hybrid3(:,ibnd) * exxalfa
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux1)
  DEALLOCATE(aux_hybrid3)
  !
ENDSUBROUTINE



SUBROUTINE hybrid_kernel_term4 (current_spin, evc1, hybrid_kd4, sf)
  !
  ! \sum_{v'} (\int v_c a_{v'} a_{v}) \phi_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma,&
                                    double_invfft_gamma,double_fwfft_gamma
  USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
  USE mp_global,             ONLY : my_image_id,inter_bgrp_comm,inter_image_comm
  USE pwcom,                 ONLY : npw,npwx,igk_k,nks,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE control_flags,         ONLY : gamma_only
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : band_group,kpt_pool
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
  USE becmod_subs_gpum,      ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,              ONLY : dvrs,hevc1,reallocate_ps_gpu
  USE cublas
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  ! 
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN)       :: current_spin
  COMPLEX(DP),INTENT(IN)    :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN)       :: sf
  COMPLEX(DP),INTENT(INOUT) :: hybrid_kd4(npwx,band_group%nlocx)
  !
  ! local vars
  !
  INTEGER  :: lbnd, ibnd, ibndp, jbnd, jbndp, ibnd1, ir, ig, iks_do
  INTEGER  :: nbndval_q, current_spin_ikq, ikq, nbndval, flnbndval
  INTEGER  :: nbnd_do
  !
  COMPLEX(DP),ALLOCATABLE :: aux_hybrid4(:,:), aux_evc1(:,:)
  COMPLEX(DP),ALLOCATABLE :: caux(:), gaux(:), raux1(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  ALLOCATE(aux_evc1(npwx, nbndval0x-n_trunc_bands))
  !
  ALLOCATE(aux_hybrid4(npwx, nbndval0x-n_trunc_bands))
  aux_hybrid4 = (0._DP, 0._DP)
  !
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux1(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     !
     IF (current_spin_ikq .NE. current_spin) CYCLE  
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(ikq)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     ! Collect fragments of evc1 from all bgrps and broadcast to all bgrps
     !
     aux_evc1(:,:) = (0._DP, 0._DP)
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd)
        ibndp = ibnd + n_trunc_bands
        !
        aux_evc1(:,ibnd) = evc1(:,lbnd,iks_do)
        !
     ENDDO
     !
     CALL mp_sum(aux_evc1(:,:), inter_bgrp_comm)
     !
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        jbndp = jbnd + n_trunc_bands
        !
        raux1(:) = (0._DP, 0._DP)
        !
        DO lbnd = 1, nbnd_do ! index to be summed 
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd + n_trunc_bands
           !
           caux(:) = (0._DP, 0._DP)
           gaux(:) = (0._DP, 0._DP)
           !
           IF (gamma_only) THEN ! Only write the gamma_only case
              !
              ! product of evc1 and evc1
              !
              CALL double_invfft_gamma(dffts,npw,npwx,aux_evc1(:,jbnd),evc1(:,lbnd,iks_do),psic,'Wave')
              !
              DO ir = 1, dffts%nnr
                 !
                 caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir)),0._DP,KIND=DP)
                 !
              ENDDO
              !
           ENDIF
           !
           caux(:) = caux(:)/omega
           !
           ! Apply the bare Coulomb potential
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,caux(:),gaux(:),'Wave')
           !
           DO ig = 1, npw
              !
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
              !
           ENDDO
           !
           caux(:) = (0._DP, 0._DP)
           !
           CALL double_invfft_gamma(dffts,npw,npwx,gaux(:),evc_work(:,ibndp),caux(:),'Wave')
           !
           psic(:) = (0._DP, 0._DP)
           !
           DO ir = 1, dffts%nnr
              !
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),0._DP,KIND=DP)
              !
           ENDDO
           !
           raux1(:) = raux1(:) + psic(:)
           !
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,raux1(:),aux_hybrid4(:,jbnd),'Wave')
        !
     ENDDO
     !
     CALL mp_sum (aux_hybrid4(:,:), inter_bgrp_comm)
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd)
        ibndp = ibnd + n_trunc_bands
        !
        hybrid_kd4(:,lbnd) = hybrid_kd4(:,lbnd) - aux_hybrid4(:,ibnd) * exxalfa
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(aux_hybrid4)
  DEALLOCATE(aux_evc1)
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux1)
  !
ENDSUBROUTINE
