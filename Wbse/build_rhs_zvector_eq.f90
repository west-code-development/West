!
! Copyright (C) 2015-2024 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Yu Jin, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE build_rhs_zvector_eq(dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_push,              ONLY : io_push_title
  USE pwcom,                ONLY : ngk,npw,npwx,nspin
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                 & l_hybrid_tddft,l_bse_triplet
  USE mp,                   ONLY : mp_bcast
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE wavefunctions,        ONLY : evc
  USE distribution_center,  ONLY : kpt_pool,band_group
#if defined(__CUDA)
  USE west_gpu,             ONLY : reallocate_ps_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(IN) :: drhox1(dffts%nnr, nspin), drhox2(dffts%nnr, nspin)
  COMPLEX(DP), INTENT(OUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: iks,lbnd,ibnd,nbndval,nbnd_do
  !
  CALL start_clock('build_zvec')
  !
  CALL io_push_title('Build the RHS of the Z vector equation')
  !
  !$acc kernels present(z_rhs_vec)
  z_rhs_vec = (0._DP,0._DP)
  !$acc end kernels
  !
  ! part1: moved to the end
  !
  ! part2: d < a | K1e | a > / d | v >
  !
  IF(.NOT. l_bse_triplet) CALL rhs_zvector_part2( dvg_exc_tmp, z_rhs_vec )
  !
  ! part3: d^2 vxc / d rho^2 contribution to d < a | K1e | a > / d | v >
  !
  IF((.NOT. l_bse) .AND. (.NOT. l_bse_triplet)) CALL rhs_zvector_part3( dvg_exc_tmp, z_rhs_vec )
  !
  ! part4: d < a | K1d | a > / d | v >
  !
  IF(l_hybrid_tddft) THEN
     CALL rhs_zvector_part4( dvg_exc_tmp, z_rhs_vec )
  ELSEIF(l_bse) THEN
     CALL errore('build_rhs_zvector_eq', 'BSE forces not implemented', 1)
  ENDIF
  !
  ! part1: d < a | D | a > / d | v >
  !
  CALL rhs_zvector_part1( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
  !
  DO iks = 1,kpt_pool%nloc
     !
     nbndval = nbnd_occ(iks)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     npw = ngk(iks)
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     ! Pc[k]*z_rhs_vec
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     !
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec(:,:,iks),(1._DP,0._DP))
     !
  ENDDO
  !
  CALL stop_clock('build_zvec')
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part1( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE io_push,              ONLY : io_push_title
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                 & l_hybrid_tddft,l_spin_flip,evc1_all
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
  USE mp,                   ONLY : mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE wbse_dv,              ONLY : wbse_dv_of_drho
  USE wbse_bgrp,            ONLY : gather_bands
  USE west_mp,              ONLY : west_mp_wait
  USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(IN) :: drhox1(dffts%nnr, nspin), drhox2(dffts%nnr, nspin)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  LOGICAL :: lrpa
  INTEGER :: ibnd,jbnd,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  INTEGER :: req
  COMPLEX(DP), ALLOCATABLE :: dotp(:)
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part1(:,:,:),tmp_vec(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: drhox(:,:)
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute d <a|D|a> / d |v>')
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part1(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  !$acc enter data create(z_rhs_vec_part1)
  !
  !$acc kernels present(z_rhs_vec_part1)
  z_rhs_vec_part1(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  IF(l_bse .OR. l_hybrid_tddft) THEN
     !
     ALLOCATE(tmp_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc))
     !$acc enter data create(tmp_vec)
     !
     !$acc kernels present(tmp_vec)
     tmp_vec(:,:,:) = (0._DP,0._DP)
     !$acc end kernels
     !
  ENDIF
  !
  ! Compute drhox
  !
  ALLOCATE(drhox(dffts%nnr, nspin))
  !
  DO iks = 1,nspin
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     drhox(:,iks) = drhox1(:,iks) + drhox2(:,iks_do)
     !
  ENDDO
  !
  !$acc enter data copyin(drhox)
  !
  lrpa = l_bse
  !
  CALL wbse_dv_of_drho(drhox,lrpa,.FALSE.)
  !
  CALL start_bar_type(barra,'zvec1',kpt_pool%nloc)
  !
  DO iks = 1,kpt_pool%nloc
     !
     nbndval = nbnd_occ(iks)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
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
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     ! ... Apply \Delta V_HXC on the z-vector
     !
     ! double bands @ gamma
     !
     DO lbnd = 1, nbnd_do-MOD(nbnd_do,2), 2
        !
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
        !
        CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
        !
        !$acc parallel loop present(drhox)
        DO ir = 1,dffts_nnr
           psic(ir) = psic(ir)*CMPLX(REAL(drhox(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part1(:,lbnd,iks),z_rhs_vec_part1(:,lbnd+1,iks),'Wave')
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
        CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
        !
        !$acc parallel loop present(drhox)
        DO ir = 1,dffts_nnr
           psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(drhox(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part1(:,lbnd,iks),'Wave')
        !
     ENDIF
     !
     IF(l_bse) CALL errore('build_rhs_zvector_eq', 'BSE forces not implemented', 1)
     !
     IF(l_hybrid_tddft) THEN
        !
        CALL hybrid_kernel_term3(current_spin,dvg_exc_tmp,z_rhs_vec_part1(:,:,iks),l_spin_flip)
        !
        IF(l_spin_flip) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        !
        !$acc host_data use_device(evc,dvgdvg_mat,tmp_vec)
        CALL DGEMM('N','N',2*npwx*npol,nbnd_do,nbndval-n_trunc_bands,-1._DP,evc(1,1+n_trunc_bands),&
        & 2*npwx*npol,dvgdvg_mat(1,1,iks_do),nbndval0x-n_trunc_bands,0._DP,tmp_vec(1,1,iks),2*npwx*npol)
        !$acc end host_data
        !
        CALL gather_bands(tmp_vec(:,:,iks),evc1_all(:,:,iks),req)
        CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
        !$acc update device(evc1_all(:,:,iks))
#endif
        CALL bse_kernel_gamma(current_spin,evc1_all(:,:,iks),z_rhs_vec_part1(:,:,iks),.FALSE.)
        !
     ENDIF
     !
     IF(gstart == 2) THEN
        !$acc parallel loop present(z_rhs_vec_part1)
        DO lbnd = 1,nbnd_do
           z_rhs_vec_part1(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part1(1,lbnd,iks),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part1)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part1(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     CALL update_bar_type(barra,'zvec1',1)
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'zvec1')
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part1,z_rhs_vec_part1,band_group%nlocx,dotp)
  !
  WRITE(stdout,*)
  WRITE(stdout,"(5x,'Norm of z_rhs_vec p1 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
  !
  DEALLOCATE(dotp)
  !$acc exit data delete(z_rhs_vec_part1)
  DEALLOCATE(z_rhs_vec_part1)
  IF(l_bse .OR. l_hybrid_tddft) THEN
     !$acc exit data delete(tmp_vec)
     DEALLOCATE(tmp_vec)
  ENDIF
  !$acc exit data delete(drhox)
  DEALLOCATE(drhox)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part2( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE io_push,              ONLY : io_push_title
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,l_spin_flip,&
                                 & l_spin_flip_kernel
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,intra_bgrp_comm
  USE wbse_dv,              ONLY : wbse_dv_of_drho,wbse_dv_of_drho_sf
  USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  LOGICAL :: lrpa
  INTEGER :: ibnd,ibndp,jbnd,jbndp,kbnd,kbndp,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd,flnbndval
  INTEGER :: dffts_nnr,band_group_myoffset
  COMPLEX(DP), ALLOCATABLE :: dotp(:)
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part2(:,:,:),aux_g(:,:),evc_copy(:,:)
  REAL(DP) :: reduce,reduce2
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: dpcpart(:,:)
  REAL(DP), ALLOCATABLE :: dv_vv_mat(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dpcpart,dv_vv_mat
#endif
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute d <a|K1e|a> / d |v>')
  !
  dffts_nnr = dffts%nnr
  band_group_myoffset = band_group%myoffset
  !
  ALLOCATE(z_rhs_vec_part2(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(aux_g(npwx*npol,2))
  ALLOCATE(dv_vv_mat(nbndval0x-n_trunc_bands, band_group%nlocx))
  ALLOCATE(dpcpart(npwx*npol, nbndval0x-n_trunc_bands))
  ALLOCATE(dvrs(dffts%nnr,nspin))
  !$acc enter data create(z_rhs_vec_part2,aux_g,dv_vv_mat,dpcpart,dvrs)
  !
  !$acc kernels present(z_rhs_vec_part2)
  z_rhs_vec_part2(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  !$acc kernels present(dpcpart)
  dpcpart(:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  dv_vv_mat(:,:) = (0._DP,0._DP)
  !
  IF(.NOT. l_spin_flip) THEN
     !
     ! Calculation of the charge density response
     !
     CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.FALSE.)
     !
     !$acc update device(dvrs)
     !
     lrpa = l_bse
     !
     CALL wbse_dv_of_drho(dvrs,lrpa,.FALSE.)
     !
     CALL start_bar_type(barra,'zvec2',kpt_pool%nloc)
     !
     DO iks = 1,kpt_pool%nloc
        !
        nbndval = nbnd_occ(iks)
        !
        nbnd_do = 0
        DO lbnd = 1,band_group%nloc
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
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
           IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc,0,inter_image_comm)
           !$acc update device(evc)
        ENDIF
        !
        ! ... Apply \Delta V_HXC on vector
        !
        ! double bands @ gamma
        !
        DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
           !
           CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),dvg_exc_tmp(:,lbnd+1,iks),psic,'Wave')
           !
           !$acc parallel loop present(dvrs)
           DO ir = 1,dffts_nnr
              psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),z_rhs_vec_part2(:,lbnd+1,iks),'Wave')
           !
        ENDDO
        !
        ! single band @ gamma
        !
        IF(MOD(nbnd_do,2) == 1) THEN
           !
           lbnd = nbnd_do
           !
           CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave')
           !
           !$acc parallel loop present(dvrs)
           DO ir = 1,dffts_nnr
              psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),'Wave')
           !
        ENDIF
        !
        ! Compute the second part: dv_vv_mat
        !
        DO lbnd = 1,nbnd_do
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd+n_trunc_bands
           !
           ! double band @ gamma
           !
           DO jbnd = 1,(nbndval-n_trunc_bands)-MOD((nbndval-n_trunc_bands),2),2
              !
              kbnd = jbnd+1
              jbndp = jbnd+n_trunc_bands
              kbndp = kbnd+n_trunc_bands
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),evc(:,kbndp),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),aux_g(:,2),'Wave')
              !
              reduce = 0._DP
              reduce2 = 0._DP
              !$acc parallel loop reduction(+:reduce,reduce2) present(evc,aux_g) copy(reduce,reduce2)
              DO ig = 1, npw
                 reduce = reduce + REAL(evc(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                 &               + AIMAG(evc(ig,ibndp)) * AIMAG(aux_g(ig,1))
                 reduce2 = reduce2 + REAL(evc(ig,ibndp),KIND=DP) * REAL(aux_g(ig,2),KIND=DP) &
                 &                 + AIMAG(evc(ig,ibndp)) * AIMAG(aux_g(ig,2))
              ENDDO
              !$acc end parallel
              !
              IF(gstart == 2) THEN
                 !$acc update host(aux_g(1,1:2))
                 reduce = reduce - 0.5_DP * REAL(evc(1,ibndp),KIND=DP) * REAL(aux_g(1,1),KIND=DP)
                 reduce2 = reduce2 - 0.5_DP * REAL(evc(1,ibndp),KIND=DP) * REAL(aux_g(1,2),KIND=DP)
              ENDIF
              !
              dv_vv_mat(jbnd,lbnd) = 2._DP * reduce
              dv_vv_mat(kbnd,lbnd) = 2._DP * reduce2
              !
           ENDDO
           !
           ! single band @ gamma
           !
           IF(MOD((nbndval-n_trunc_bands),2) == 1) THEN
              !
              jbnd = nbndval - n_trunc_bands
              jbndp = jbnd + n_trunc_bands
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),'Wave')
              !
              reduce = 0._DP
              !$acc parallel loop reduction(+:reduce) present(evc,aux_g) copy(reduce)
              DO ig = 1, npw
                 reduce = reduce + REAL(evc(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                 &               + AIMAG(evc(ig,ibndp)) * AIMAG(aux_g(ig,1))
              ENDDO
              !$acc end parallel
              !
              IF(gstart == 2) THEN
                 !$acc update host(aux_g(1,1))
                 reduce = reduce - 0.5_DP * REAL(evc(1,ibndp),KIND=DP) * REAL(aux_g(1,1),KIND=DP)
              ENDIF
              !
              dv_vv_mat(jbnd,lbnd) = 2._DP * reduce
              !
           ENDIF
           !
        ENDDO
        !
        CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
        !
        !$acc update device(dv_vv_mat)
        !
        !$acc host_data use_device(dvg_exc_tmp,dv_vv_mat,dpcpart)
        CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,dvg_exc_tmp(1,1,iks),&
        & 2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,dpcpart,2*npwx*npol)
        !$acc end host_data
        !
        !$acc update host(dpcpart)
        CALL mp_sum(dpcpart,inter_bgrp_comm)
        !$acc update device(dpcpart)
        !
        !$acc parallel loop collapse(2) present(z_rhs_vec_part2,dpcpart)
        DO lbnd = 1,nbnd_do
           DO ig = 1,npw
              !
              ! ibnd = band_group%l2g(lbnd)
              !
              ibnd = band_group_myoffset+lbnd
              !
              z_rhs_vec_part2(ig,lbnd,iks) = z_rhs_vec_part2(ig,lbnd,iks)+dpcpart(ig,ibnd)
              !
           ENDDO
        ENDDO
        !$acc end parallel
        !
        IF(gstart == 2) THEN
           !$acc parallel loop present(z_rhs_vec_part2)
           DO lbnd = 1,nbnd_do
              z_rhs_vec_part2(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part2(1,lbnd,iks),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part2)
        DO lbnd = 1,nbnd_do
           DO ig = 1,npw
              z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part2(ig,lbnd,iks)
           ENDDO
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'zvec2',1)
        !
     ENDDO
     !
  ELSE
     !
     IF(l_spin_flip_kernel) THEN
        !
        ALLOCATE(evc_copy(npwx, nbndval0x-n_trunc_bands))
        !$acc enter data create(evc_copy)
        !
        ! Calculation of the charge density response
        !
        CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.TRUE.)
        !
        !$acc update device(dvrs)
        !
        CALL wbse_dv_of_drho_sf(dvrs)
        !
        CALL start_bar_type(barra,'zvec2',kpt_pool%nloc)
        !
        DO iks = 1,kpt_pool%nloc
           !
           iks_do = flks(iks)
           !
           nbndval = nbnd_occ(iks)
           flnbndval = nbnd_occ(iks_do)
           !
           nbnd_do = 0
           DO lbnd = 1,band_group%nloc
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
              IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
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
              IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
              CALL mp_bcast(evc,0,inter_image_comm)
              !$acc update device(evc)
           ENDIF
           !
           ! ... Apply \Delta V_HXC on vector
           !
           ! double bands @ gamma
           !
           DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
              !
              CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks_do),dvg_exc_tmp(:,lbnd+1,iks_do),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,iks_do),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),z_rhs_vec_part2(:,lbnd+1,iks),'Wave')
              !
           ENDDO
           !
           ! single band @ gamma
           !
           IF(MOD(nbnd_do,2) == 1) THEN
              !
              lbnd = nbnd_do
              !
              CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks_do),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,iks_do),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),'Wave')
              !
           ENDIF
           !
           !$acc parallel loop collapse(2) present(evc_copy,evc)
           DO ibnd = 1,nbndval-n_trunc_bands
              DO ig = 1,npwx
                 evc_copy(ig,ibnd) = evc(ig,ibnd+n_trunc_bands)
              ENDDO
           ENDDO
           !$acc end parallel
           !
           IF(kpt_pool%nloc > 1) THEN
              IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
              CALL mp_bcast(evc,0,inter_image_comm)
              !$acc update device(evc)
           ENDIF
           !
           ! evc_copy -> current spin channel
           ! evc -> opposite spin channel
           !
           ! recompute nbnd_do for the opposite spin channel
           !
           nbnd_do = 0
           DO lbnd = 1,band_group%nloc
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
              IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
           ENDDO
           !
           ! Compute the second part: dv_vv_mat
           !
           DO lbnd = 1,nbnd_do
              !
              ibnd = band_group%l2g(lbnd)
              ibndp = ibnd+n_trunc_bands
              !
              ! double band @ gamma
              !
              DO jbnd = 1,(nbndval-n_trunc_bands)-MOD((nbndval-n_trunc_bands),2),2
                 !
                 kbnd = jbnd+1
                 jbndp = jbnd+n_trunc_bands
                 kbndp = kbnd+n_trunc_bands
                 !
                 CALL double_invfft_gamma(dffts,npw,npwx,evc_copy(:,jbnd),evc_copy(:,kbnd),psic,'Wave')
                 !
                 !$acc parallel loop present(dvrs)
                 DO ir = 1,dffts_nnr
                    psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
                 ENDDO
                 !$acc end parallel
                 !
                 CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),aux_g(:,2),'Wave')
                 !
                 reduce = 0._DP
                 reduce2 = 0._DP
                 !$acc parallel loop reduction(+:reduce,reduce2) present(evc,aux_g) copy(reduce,reduce2)
                 DO ig = 1, npw
                    reduce = reduce + REAL(evc(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                    &               + AIMAG(evc(ig,ibndp)) * AIMAG(aux_g(ig,1))
                    reduce2 = reduce2 + REAL(evc(ig,ibndp),KIND=DP) * REAL(aux_g(ig,2),KIND=DP) &
                    &                 + AIMAG(evc(ig,ibndp)) * AIMAG(aux_g(ig,2))
                 ENDDO
                 !$acc end parallel
                 !
                 IF(gstart == 2) THEN
                    !$acc update host(aux_g(1,1:2))
                    reduce = reduce - 0.5_DP * REAL(evc(1,ibndp),KIND=DP) * REAL(aux_g(1,1),KIND=DP)
                    reduce2 = reduce2 - 0.5_DP * REAL(evc(1,ibndp),KIND=DP) * REAL(aux_g(1,2),KIND=DP)
                 ENDIF
                 !
                 dv_vv_mat(jbnd,lbnd) = 2._DP * reduce
                 dv_vv_mat(kbnd,lbnd) = 2._DP * reduce2
                 !
              ENDDO
              !
              ! single band @ gamma
              !
              IF(MOD((nbndval-n_trunc_bands),2) == 1) THEN
                 !
                 jbnd = nbndval - n_trunc_bands
                 jbndp = jbnd + n_trunc_bands
                 !
                 CALL single_invfft_gamma(dffts,npw,npwx,evc_copy(:,jbnd),psic,'Wave')
                 !
                 !$acc parallel loop present(dvrs)
                 DO ir = 1,dffts_nnr
                    psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
                 ENDDO
                 !$acc end parallel
                 !
                 CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),'Wave')
                 !
                 reduce = 0._DP
                 !$acc parallel loop reduction(+:reduce) present(evc,aux_g) copy(reduce)
                 DO ig = 1, npw
                    reduce = reduce + REAL(evc(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                    &               + AIMAG(evc(ig,ibndp)) * AIMAG(aux_g(ig,1))
                 ENDDO
                 !$acc end parallel
                 !
                 IF(gstart == 2) THEN
                    !$acc update host(aux_g(1,1))
                    reduce = reduce - 0.5_DP * REAL(evc(1,ibndp),KIND=DP) * REAL(aux_g(1,1),KIND=DP)
                 ENDIF
                 !
                 dv_vv_mat(jbnd,lbnd) = 2._DP * reduce
                 !
              ENDIF
              !
           ENDDO
           !
           CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
           !
           !$acc update device(dv_vv_mat)
           !
           !$acc host_data use_device(dvg_exc_tmp,dv_vv_mat,dpcpart)
           CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,dvg_exc_tmp(1,1,iks),&
           & 2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,dpcpart,2*npwx*npol)
           !$acc end host_data
           !
           !$acc update host(dpcpart)
           CALL mp_sum(dpcpart,inter_bgrp_comm)
           !$acc update device(dpcpart)
           !
           ! recompute nbnd_do for the current spin channel
           !
           nbnd_do = 0
           DO lbnd = 1,band_group%nloc
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
              IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
           ENDDO
           !
           !$acc parallel loop collapse(2) present(z_rhs_vec_part2,dpcpart)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 !
                 ! ibnd = band_group%l2g(lbnd)
                 !
                 ibnd = band_group_myoffset+lbnd
                 !
                 z_rhs_vec_part2(ig,lbnd,iks) = z_rhs_vec_part2(ig,lbnd,iks)+dpcpart(ig,ibnd)
                 !
              ENDDO
           ENDDO
           !$acc end parallel
           !
           IF(gstart == 2) THEN
              !$acc parallel loop present(z_rhs_vec_part2)
              DO lbnd = 1,nbnd_do
                 z_rhs_vec_part2(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part2(1,lbnd,iks),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
           !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part2)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part2(ig,lbnd,iks)
              ENDDO
           ENDDO
           !$acc end parallel
           !
           CALL update_bar_type(barra,'zvec2',1)
           !
        ENDDO
        !
        !$acc exit data delete(evc_copy)
        DEALLOCATE(evc_copy)
        !
     ENDIF
     !
  ENDIF
  !
  CALL stop_bar_type(barra,'zvec2')
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part2,z_rhs_vec_part2,band_group%nlocx,dotp)
  !
  WRITE(stdout,*)
  WRITE(stdout,"(5x,'Norm of z_rhs_vec p2 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
  !
  DEALLOCATE(dotp)
  !$acc exit data delete(z_rhs_vec_part2,aux_g,dv_vv_mat,dpcpart,dvrs)
  DEALLOCATE(z_rhs_vec_part2)
  DEALLOCATE(aux_g)
  DEALLOCATE(dv_vv_mat)
  DEALLOCATE(dpcpart)
  DEALLOCATE(dvrs)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part3( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE io_push,              ONLY : io_push_title
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,n_trunc_bands,l_spin_flip
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
  USE mp,                   ONLY : mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,iks,ir,ig,nbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: dotp(:)
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part3(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: ddvxc(:,:)
  TYPE(bar_type) :: barra
  !
  CALL io_push_title('Compute d^2 vxc / d rho^2')
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part3(npwx*npol,band_group%nlocx,kpt_pool%nloc))
  !$acc enter data create(z_rhs_vec_part3)
  !
  !$acc kernels present(z_rhs_vec_part3)
  z_rhs_vec_part3(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  ALLOCATE(ddvxc(dffts%nnr,nspin))
  !
  IF(.NOT. l_spin_flip) THEN
     CALL compute_ddvxc_5p(dvg_exc_tmp, ddvxc)
  ELSE
     CALL compute_ddvxc_sf(dvg_exc_tmp, ddvxc)
  ENDIF
  !
  !$acc enter data copyin(ddvxc)
  !
  CALL start_bar_type(barra,'zvec3',kpt_pool%nloc)
  !
  DO iks = 1,kpt_pool%nloc
     !
     nbndval = nbnd_occ(iks)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
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
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     ! ... Apply \Delta V_HXC
     !
     ! double bands @ gamma
     !
     DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
        !
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
        !
        CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
        !
        !$acc parallel loop present(ddvxc)
        DO ir = 1,dffts_nnr
           psic(ir) = psic(ir)*CMPLX(REAL(ddvxc(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part3(:,lbnd,iks),z_rhs_vec_part3(:,lbnd+1,iks),'Wave')
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
        CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
        !
        !$acc parallel loop present(ddvxc)
        DO ir = 1,dffts_nnr
           psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(ddvxc(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part3(:,lbnd,iks),'Wave')
        !
     ENDIF
     !
     IF(gstart == 2) THEN
        !$acc parallel loop present(z_rhs_vec_part3)
        DO lbnd = 1,nbnd_do
           z_rhs_vec_part3(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part3(1,lbnd,iks),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part3)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part3(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     CALL update_bar_type(barra,'zvec3',1)
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'zvec3')
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part3,z_rhs_vec_part3,band_group%nlocx,dotp)
  !
  WRITE(stdout,*)
  WRITE(stdout,"(5x,'Norm of z_rhs_vec p3 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
  !
  DEALLOCATE(dotp)
  !$acc exit data delete(z_rhs_vec_part3,ddvxc)
  DEALLOCATE(z_rhs_vec_part3)
  DEALLOCATE(ddvxc)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_ddvxc_5p( dvg_exc_tmp, ddvxc )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE scf,                  ONLY : rho,rho_core,rhog_core,scf_type,create_scf_type,destroy_scf_type
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE fft_interfaces,       ONLY : fwfft
  USE westcom,              ONLY : ddvxc_fd_coeff
  USE wavefunctions,        ONLY : psic
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: ddvxc(dffts%nnr, nspin)
  !
  ! Workspace
  !
  INTEGER :: iks,ir,indk
  REAL(DP), ALLOCATABLE :: aux_vxc(:,:,:), vxc(:,:), rdvrs(:,:)
  REAL(DP) :: etxc,vtxc
  TYPE (scf_type) :: a_rho
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('ddvxc_5p')
#else
  CALL start_clock('ddvxc_5p')
#endif
  !
  CALL create_scf_type(a_rho)
  !
  ALLOCATE(aux_vxc(dffts%nnr, nspin, 5))
  ALLOCATE(vxc(dffts%nnr, nspin))
  ALLOCATE(dvrs(dffts%nnr,nspin))
  !$acc enter data create(dvrs)
  ALLOCATE(rdvrs(dffts%nnr, nspin))
  !
  ! Calculation of the charge density response
  !
  CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.FALSE.)
  !
  IF(nspin == 1) THEN
     !
     rdvrs(:,1) = REAL(dvrs(:,1),KIND=DP)
     !
  ELSEIF(nspin == 2) THEN
     !
     rdvrs(:,1) = REAL(dvrs(:,1),KIND=DP) + REAL(dvrs(:,2),KIND=DP)
     rdvrs(:,2) = REAL(dvrs(:,1),KIND=DP) - REAL(dvrs(:,2),KIND=DP)
     !
  ELSEIF(nspin == 4) THEN
     !
     CALL errore('compute_ddvxc_5p', 'nspin == 4 not supported', 1)
     !
  ENDIF
  !
  DO indk = 1, 5
     !
     vxc(:,:) = 0._DP
     !
     a_rho%of_r(:,:) = rho%of_r + REAL((indk-3),KIND=DP) * ddvxc_fd_coeff * rdvrs
     !
     DO iks = 1, nspin
        !
        psic(:) = a_rho%of_r(:,iks)
        CALL fwfft ('Rho', psic, dffts)
        a_rho%of_g(:,iks) = psic(dffts%nl)
        !
     ENDDO
     !
     CALL v_xc(a_rho, rho_core, rhog_core, etxc, vtxc, vxc)
     !
     aux_vxc(:,:,indk) = vxc
     !
  ENDDO
  !
  ! compute ddvxc
  !
  DO iks = 1,nspin
     DO ir = 1,dffts%nnr
        ddvxc(ir,iks) = CMPLX( (-1._DP*aux_vxc(ir,iks,1)+16._DP*aux_vxc(ir,iks,2) &
        &                       -30._DP*aux_vxc(ir,iks,3)+16._DP*aux_vxc(ir,iks,4) &
        &                       -1._DP*aux_vxc(ir,iks,5)), KIND=DP ) / (12._DP*ddvxc_fd_coeff**2)
     ENDDO
  ENDDO
  !
  DEALLOCATE(aux_vxc)
  DEALLOCATE(vxc)
  !$acc exit data delete(dvrs)
  DEALLOCATE(dvrs)
  DEALLOCATE(rdvrs)
  !
  CALL destroy_scf_type(a_rho)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('ddvxc_5p')
#else
  CALL stop_clock('ddvxc_5p')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_ddvxc_sf( dvg_exc_tmp, ddvxc )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx
  USE noncollin_module,     ONLY : npol,nspin_gga
  USE xc_lib,               ONLY : xclib_dft_is
  USE fft_base,             ONLY : dffts
  USE gvect,                ONLY : g
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : nlcc_any
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE westcom,              ONLY : wbse_save_dir,sf_kernel,l_spin_flip_kernel,l_spin_flip_alda0,&
                                 & spin_flip_cut,l_print_spin_flip_kernel
  USE qpoint,               ONLY : xq
  USE gc_lr,                ONLY : grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s
  USE eqv,                  ONLY : dmuxc
  USE cubefile,             ONLY : write_wfc_cube_r
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: ddvxc(dffts%nnr, nspin)
  !
  ! Workspace
  !
  INTEGER :: ir,is,is1
  REAL(DP) :: tmp1,tmp2
  COMPLEX(DP), ALLOCATABLE :: drho_sf(:,:),drho_sf_copy(:,:)
  CHARACTER(LEN=:), ALLOCATABLE :: fname
  !
#if defined(__CUDA)
  CALL start_clock_gpu('ddvxc_sf')
#else
  CALL start_clock('ddvxc_sf')
#endif
  !
  IF(nlcc_any) CALL errore('compute_ddvxc_sf', 'nlcc_any not supported', 1)
  !
  ALLOCATE(drho_sf(dffts%nnr,2))
  !$acc enter data create(drho_sf)
  ALLOCATE(drho_sf_copy(dffts%nnr,2))
  !
  CALL wbse_calc_dens(dvg_exc_tmp,drho_sf,.TRUE.)
  !
  DO ir = 1,dffts%nnr
     tmp1 = REAL(drho_sf(ir,1),KIND=DP)**2
     tmp2 = REAL(drho_sf(ir,2),KIND=DP)**2
     drho_sf_copy(ir,1) = CMPLX(tmp1+tmp2,KIND=DP)
     drho_sf_copy(ir,2) = -CMPLX(tmp1+tmp2,KIND=DP)
  ENDDO
  !
  ! divide rho_diff
  !
  DO ir = 1,dffts%nnr
     IF(ABS(rho%of_r(ir,2)) < spin_flip_cut) THEN
        drho_sf_copy(ir,1) = (0._DP,0._DP)
        drho_sf_copy(ir,2) = (0._DP,0._DP)
     ELSE
        drho_sf_copy(ir,1) = drho_sf_copy(ir,1) / rho%of_r(ir,2)
        drho_sf_copy(ir,2) = drho_sf_copy(ir,2) / rho%of_r(ir,2)
     ENDIF
  ENDDO
  !
  ddvxc(:,:) = (0._DP,0._DP)
  !
  IF(l_spin_flip_kernel) THEN
     !
     ! part 2
     !
     DO is = 1,nspin
        DO is1 = 1,nspin
           ddvxc(:,is) = ddvxc(:,is) + dmuxc(:,is,is1) * drho_sf_copy(:,is1)
        ENDDO
     ENDDO
     !
     IF(.NOT. l_spin_flip_alda0) THEN
        IF(xclib_dft_is('gradient')) THEN
           CALL dgradcorr(dffts, rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
                & drho_sf_copy, nspin, nspin_gga, g, ddvxc)
        ENDIF
     ENDIF
     !
     ! part 1
     !
     !$acc update host(sf_kernel)
     !
     DO ir = 1,dffts%nnr
        ddvxc(ir,1) = ddvxc(ir,1) - sf_kernel(ir) * drho_sf_copy(ir,1)
        ddvxc(ir,2) = ddvxc(ir,2) - sf_kernel(ir) * drho_sf_copy(ir,2)
     ENDDO
     !
  ENDIF
  !
  ! print spin flip kernel (keep it for now)
  !
  IF(l_print_spin_flip_kernel) THEN
     !
     fname=TRIM(wbse_save_dir)//'/rho_diff.cube'
     CALL write_wfc_cube_r(dffts, fname, rho%of_r(ir,2))
     !
     fname=TRIM(wbse_save_dir)//'/drho_sf_up.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf(:,1),KIND=DP))
     !
     fname=TRIM(wbse_save_dir)//'/drho_sf_down.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf(:,2),KIND=DP))
     !
     fname=TRIM(wbse_save_dir)//'/drho_sq_rho_diff_up.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf_copy(:,1),KIND=DP))
     !
     fname=TRIM(wbse_save_dir)//'/drho_sq_rho_diff_down.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf_copy(:,2),KIND=DP))
     !
     fname=TRIM(wbse_save_dir)//'/ddvxc_up.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(ddvxc(:,1),KIND=DP))
     !
     fname=TRIM(wbse_save_dir)//'/ddvxc_down.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(ddvxc(:,2),KIND=DP))
     !
  ENDIF
  !
  !$acc exit data delete(drho_sf)
  DEALLOCATE(drho_sf)
  DEALLOCATE(drho_sf_copy)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('ddvxc_sf')
#else
  CALL stop_clock('ddvxc_sf')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part4( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE io_push,              ONLY : io_push_title
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip,evc1_all
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,intra_bgrp_comm
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: ig,lbnd,ibnd,jbnd,jbndp,iks,iks_do,nbnd_do,nbndval,flnbndval
  INTEGER :: band_group_myoffset
  REAL(DP) :: reduce
  COMPLEX(DP), ALLOCATABLE :: dotp(:)
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part4(:,:,:),tmp_vec(:,:)
  REAL(DP), ALLOCATABLE :: dv_vv_mat(:,:)
  COMPLEX(DP), ALLOCATABLE :: dpcpart(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dpcpart
#endif
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute d <a|K1d|a> / d |v>')
  !
  band_group_myoffset = band_group%myoffset
  !
  ALLOCATE(z_rhs_vec_part4(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(dv_vv_mat(nbndval0x-n_trunc_bands, band_group%nlocx))
  ALLOCATE(tmp_vec(npwx*npol, band_group%nlocx))
  ALLOCATE(dpcpart(npwx*npol, nbndval0x-n_trunc_bands))
  !$acc enter data create(z_rhs_vec_part4,dv_vv_mat,tmp_vec,dpcpart)
  !
  !$acc kernels present(z_rhs_vec_part4)
  z_rhs_vec_part4(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  !$acc kernels present(dpcpart)
  dpcpart(:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  !$acc kernels present(dv_vv_mat)
  dv_vv_mat(:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  CALL start_bar_type(barra,'zvec4',kpt_pool%nloc)
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
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
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     ! Compute the first part
     !
     CALL hybrid_kernel_term4(current_spin,dvg_exc_tmp,z_rhs_vec_part4(:,:,iks),l_spin_flip)
     !
     ! Compute the second part: dv_vv_mat
     !
     !$acc kernels present(tmp_vec)
     tmp_vec(:,:) = (0._DP,0._DP)
     !$acc end kernels
     !
     CALL bse_kernel_gamma(current_spin,evc1_all(:,:,iks),tmp_vec,l_spin_flip)
     !
     !$acc parallel vector_length(1024) present(evc,tmp_vec,dv_vv_mat)
     !$acc loop collapse(2)
     DO jbnd = 1, nbndval - n_trunc_bands
        DO lbnd = 1, nbnd_do
           !
           jbndp = jbnd + n_trunc_bands
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1, npw
              reduce = reduce + REAL(evc(ig,jbndp),KIND=DP) * REAL(tmp_vec(ig,lbnd),KIND=DP) &
              &               + AIMAG(evc(ig,jbndp)) * AIMAG(tmp_vec(ig,lbnd))
           ENDDO
           !
           IF(gstart == 2) THEN
              reduce = reduce - 0.5_DP * REAL(evc(1,jbndp),KIND=DP) * REAL(tmp_vec(1,lbnd),KIND=DP)
           ENDIF
           !
           dv_vv_mat(jbnd,lbnd) = 2._DP * reduce
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(dv_vv_mat)
     CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
     !$acc end host_data
     !
     !$acc host_data use_device(dvg_exc_tmp,dv_vv_mat,dpcpart)
     CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,dvg_exc_tmp(1,1,iks),&
     & 2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,dpcpart,2*npwx*npol)
     !$acc end host_data
     !
     !$acc update host(dpcpart)
     CALL mp_sum(dpcpart,inter_bgrp_comm)
     !$acc update device(dpcpart)
     !
     ! compute nbnd_do for the current spin channel
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     !$acc parallel loop collapse(2) present(z_rhs_vec_part4,dpcpart)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = band_group_myoffset+lbnd
           !
           z_rhs_vec_part4(ig,lbnd,iks) = z_rhs_vec_part4(ig,lbnd,iks)+dpcpart(ig,ibnd)
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gstart == 2) THEN
        !$acc parallel loop present(z_rhs_vec_part4)
        DO lbnd = 1,nbnd_do
           z_rhs_vec_part4(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part4(1,lbnd,iks),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part4)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part4(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     CALL update_bar_type(barra,'zvec4',1)
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'zvec4')
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part4,z_rhs_vec_part4,band_group%nlocx,dotp)
  !
  WRITE(stdout,*)
  WRITE(stdout,"(5x,'Norm of z_rhs_vec p4 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
  !
  DEALLOCATE(dotp)
  !$acc exit data delete(z_rhs_vec_part4,dv_vv_mat,tmp_vec,dpcpart)
  DEALLOCATE(z_rhs_vec_part4)
  DEALLOCATE(dv_vv_mat)
  DEALLOCATE(tmp_vec)
  DEALLOCATE(dpcpart)
  !
END SUBROUTINE
