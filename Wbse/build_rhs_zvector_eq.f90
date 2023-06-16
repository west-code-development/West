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
! Yu Jin, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE build_rhs_zvector_eq(dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_push,              ONLY : io_push_title
  USE pwcom,                ONLY : npwx,nspin
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : l_bse,l_hybrid_tddft,l_bse_triplet,nbndval0x,n_trunc_bands
  USE distribution_center,  ONLY : kpt_pool,band_group
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
  REAL(DP), EXTERNAL :: get_clock
  REAL(DP) :: time_spent(2)
  !
  CALL start_clock( 'rhs_z' )
  time_spent(1) = get_clock( 'rhs_z' )
  !
  CALL io_push_title('Forces : Build the RHS of the Z vector equation')
  !
  z_rhs_vec = (0._DP, 0._DP)
  !
  ! part1: d < a | D | a > / d | v >
  !
  CALL rhs_zvector_part1( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
  !
  time_spent(2) = get_clock( 'rhs_z' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'rhs_z' )
  !
  ! part2: d < a | K1e | a > / d | v >
  !
  IF ( .NOT. l_bse_triplet ) THEN
     !
     CALL rhs_zvector_part2( dvg_exc_tmp, z_rhs_vec )
     !
     time_spent(2) = get_clock( 'rhs_z' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'rhs_z' )
     !
  ENDIF
  !
  ! part3: d^2 vxc / d rho^2 contribution to d < a | K1e | a > / d | v >
  !
  IF(.NOT. l_bse) THEN
     !
     IF(.NOT. l_bse_triplet) THEN
        !
        CALL rhs_zvector_part3( dvg_exc_tmp, z_rhs_vec )
        !
        time_spent(2) = get_clock( 'rhs_z' )
        CALL wbse_forces_time(time_spent)
        time_spent(1) = get_clock( 'rhs_z' )
        !
     ENDIF
     !
  ENDIF
  !
  ! part4: d < a | k1d | a > / d | v >
  !
  IF( l_hybrid_tddft ) THEN
     !
     CALL rhs_zvector_part4( dvg_exc_tmp, z_rhs_vec )
     !
     time_spent(2) = get_clock( 'rhs_z' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'rhs_z' )
     !
  ELSEIF( l_bse ) THEN
     !
     CALL errore('build_rhs_zvector_eq', 'BSE forces have not been implemented', 1)
     !
  ENDIF
  !
  CALL stop_clock( 'rhs_z' )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part1( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                 & l_hybrid_tddft,l_spin_flip,l_slow_tddft_k1d
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE wbse_dv,              ONLY : wbse_dv_of_drho
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
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
  ! WORKSPACE
  !
  LOGICAL :: lrpa
  INTEGER :: ibnd,jbnd,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part1(:,:,:),dotp(:),tmp_vec(:,:,:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  COMPLEX(DP), ALLOCATABLE :: drhox(:,:)
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part1(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  !$acc kernels present(z_rhs_vec_part1)
  z_rhs_vec_part1(:,:,:) = (0._DP, 0._DP)
  !$acc end kernels
  !
  IF(l_bse .OR. l_hybrid_tddft) THEN
     ALLOCATE(tmp_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc))
     tmp_vec(:,:,:) = (0._DP, 0._DP)
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
  lrpa = l_bse
  !
  CALL wbse_dv_of_drho(drhox,lrpa,.FALSE.)
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
     ! ... Apply \Delta V_HXC on the z-vector
     !
     ! double bands @ gamma
     !
     DO lbnd = 1, nbnd_do-MOD(nbnd_do,2), 2
        !
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
        !
        CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),evc_work(:,jbnd),psic,'Wave')
        !
        !$acc parallel loop present(drhox)
        DO ir = 1,dffts_nnr
           psic(ir) = psic(ir)*CMPLX(REAL(drhox(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        !$acc host_data use_device(z_rhs_vec_part1)
        CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part1(:,lbnd,iks),z_rhs_vec_part1(:,lbnd+1,iks),'Wave')
        !$acc end host_data
        !
     ENDDO
     !
     ! single band @ gamma
     !
     IF(MOD(nbnd_do,2) == 1 ) THEN
        !
        lbnd = nbnd_do
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
        !
        !$acc parallel loop present(drhox)
        DO ir = 1,dffts_nnr
           psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(drhox(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        !$acc host_data use_device(z_rhs_vec_part1)
        CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part1(:,lbnd,iks),'Wave')
        !$acc end host_data
        !
     ENDIF
     !
     IF (l_bse) THEN
        !
        CALL errore('build_rhs_zvector_eq', 'BSE forces have not been implemented', 1)
        !
     ENDIF
     !
     IF (l_hybrid_tddft) THEN
        !
        CALL hybrid_kernel_term3(current_spin,dvg_exc_tmp,z_rhs_vec_part1(:,:,iks),l_spin_flip)
        !
        DO lbnd = 1, nbnd_do
           !
           DO jbnd = 1, nbndval - n_trunc_bands
              !
              IF(l_spin_flip) THEN
                 iks_do = flks(iks)
              ELSE
                 iks_do = iks
              ENDIF
              !
              DO ig = 1,npw
                 tmp_vec(ig,lbnd,iks) = tmp_vec(ig,lbnd,iks) &
                 & - dvgdvg_mat(jbnd,lbnd,iks_do) * evc_work(ig,jbnd+n_trunc_bands)
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        IF(l_slow_tddft_k1d) THEN
           CALL hybrid_kernel_term1_slow(current_spin,tmp_vec,z_rhs_vec_part1(:,:,iks),.FALSE.)
        ELSE
           CALL bse_kernel_gamma(current_spin,tmp_vec,z_rhs_vec_part1(:,:,iks),.FALSE.)
        ENDIF
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
     ! Pc[k]*z_rhs_vec_part1
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     !
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec_part1(:,:,iks),(1._DP,0._DP))
     !
     !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part1)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           z_rhs_vec(ig,lbnd,iks) = -z_rhs_vec_part1(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part1,z_rhs_vec_part1,band_group%nlocx,dotp)
  !
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  WRITE(stdout, "( 5x,'# Norm of z_rhs_vec p1  = | ', ES15.8, ' |' ) " ) SUM(REAL(dotp,KIND=DP))
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  !
  DEALLOCATE(dotp)
  DEALLOCATE(drhox)
  DEALLOCATE(z_rhs_vec_part1)
  !
  IF (l_bse .OR. l_hybrid_tddft) DEALLOCATE(tmp_vec)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part2( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
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
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,intra_bgrp_comm
  USE wbse_dv,              ONLY : wbse_dv_of_drho,wbse_dv_of_drho_sf
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! WORKSPACE
  !
  LOGICAL :: lrpa
  INTEGER :: ibnd,ibndp,jbnd,jbndp,kbnd,kbndp,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd,flnbndval
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part2(:,:,:),aux_g(:,:),dotp(:),evc_copy(:,:)
  REAL(DP), ALLOCATABLE :: dv_vv_mat(:,:,:)
  REAL(DP) :: reduce
  COMPLEX(DP) :: factor
  INTEGER, PARAMETER :: flks(2) = [2,1]
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: dpcpart(:,:)
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part2(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(dv_vv_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(aux_g(npwx*npol,2))
  ALLOCATE(dpcpart(npwx*npol, nbndval0x-n_trunc_bands))
  ALLOCATE(dvrs(dffts%nnr,nspin))
  !
  !$acc kernels present(z_rhs_vec_part2)
  z_rhs_vec_part2(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  dv_vv_mat(:,:,:) = (0._DP, 0._DP)
  !
  IF( .NOT. l_spin_flip ) THEN
     !
     ! Calculation of the charge density response
     !
     CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.FALSE.)
     !
     lrpa = l_bse
     !
     CALL wbse_dv_of_drho(dvrs,lrpa,.FALSE.)
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
        ! ... Apply \Delta V_HXC on vector
        !
        ! double bands @ gamma
        !
        DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
           !
           !$acc host_data use_device(dvg_exc_tmp)
           CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),dvg_exc_tmp(:,lbnd+1,iks),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(dvrs)
           DO ir = 1,dffts_nnr
              psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(z_rhs_vec_part2)
           CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),z_rhs_vec_part2(:,lbnd+1,iks),'Wave')
           !$acc end host_data
           !
        ENDDO
        !
        ! single band @ gamma
        !
        IF(MOD(nbnd_do,2) == 1) THEN
           !
           lbnd = nbnd_do
           !
           !$acc host_data use_device(dvg_exc_tmp)
           CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(dvrs)
           DO ir = 1,dffts_nnr
              psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(z_rhs_vec_part2)
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),'Wave')
           !$acc end host_data
           !
        ENDIF
        !
        ! Compute the second part: dv_vv_mat
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
              CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,jbndp),evc_work(:,kbndp),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(aux_g)
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),aux_g(:,2),'Wave')
              !$acc end host_data
              !
              reduce = 0._DP
              !$acc loop reduction(+:reduce)
              DO ig = 1, npw
                 reduce = reduce + REAL(evc_work(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                 &               + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,1))
              ENDDO
              !
              dv_vv_mat(jbnd,lbnd,iks) = 2._DP * reduce
              !
              IF (gstart == 2) THEN
                 dv_vv_mat(jbnd,lbnd,iks) = dv_vv_mat(jbnd,lbnd,iks) &
                 & - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,1),KIND=DP)
              ENDIF
              !
              reduce = 0._DP
              !$acc loop reduction(+:reduce)
              DO ig = 1, npw
                 reduce = reduce + REAL(evc_work(ig,ibndp),KIND=DP) * REAL(aux_g(ig,2),KIND=DP) &
                 &               + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,2))
              ENDDO
              !
              dv_vv_mat(kbnd,lbnd,iks) = 2._DP * reduce
              !
              IF (gstart == 2) THEN
                 dv_vv_mat(kbnd,lbnd,iks) = dv_vv_mat(kbnd,lbnd,iks) &
                 & - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,2),KIND=DP)
              ENDIF
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
              CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,jbndp),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(aux_g)
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),'Wave')
              !$acc end host_data
              !
              reduce = 0._DP
              !$acc loop reduction(+:reduce)
              DO ig = 1, npw
                 reduce = reduce + REAL(evc_work(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                 &               + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,1))
              ENDDO
              !
              dv_vv_mat(jbnd,lbnd,iks) = 2._DP * reduce
              !
              IF (gstart == 2) THEN
                 dv_vv_mat(jbnd,lbnd,iks) = dv_vv_mat(jbnd,lbnd,iks) &
                 & - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,1),KIND=DP)
              ENDIF
              !
           ENDIF
           !
        ENDDO
        !
        CALL mp_sum(dv_vv_mat(:,:,iks),intra_bgrp_comm)
        !
        dpcpart(:,:) = (0._DP, 0._DP)
        !
        DO lbnd = 1,nbnd_do
           !
           DO jbnd = 1,nbndval-n_trunc_bands
              !
              factor = CMPLX(-dv_vv_mat(jbnd,lbnd,iks),KIND=DP)
              !
              !$acc host_data use_device(dvg_exc_tmp,dpcpart)
              CALL ZAXPY(npw,factor,dvg_exc_tmp(:,lbnd,iks),1,dpcpart(:,jbnd),1)
              !$acc end host_data
              !
           ENDDO
           !
        ENDDO
        !
        CALL mp_sum(dpcpart,inter_bgrp_comm)
        !
        !$acc parallel loop collapse(2) present(z_rhs_vec_part2,dpcpart)
        DO lbnd = 1,nbnd_do
           DO ig = 1,npw
              ibnd = band_group%l2g(lbnd)
              z_rhs_vec_part2(ig,lbnd,iks) = z_rhs_vec_part2(ig,lbnd,iks)+dpcpart(ig,ibnd)
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
        ! Pc[k]*z_rhs_vec_part2
        !
#if defined(__CUDA)
        CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec_part2(:,:,iks),(1._DP,0._DP))
        !
        !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part2)
        DO lbnd = 1,nbnd_do
           DO ig = 1,npw
              z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part2(ig,lbnd,iks)
           ENDDO
        ENDDO
        !$acc end parallel
        !
     ENDDO
     !
  ELSE
     !
     IF(l_spin_flip_kernel) THEN
        !
        ALLOCATE(evc_copy(npwx, nbndval0x))
        !
        ! Calculation of the charge density response
        CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.TRUE.)
        !
        CALL wbse_dv_of_drho_sf(dvrs)
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
           ! ... Apply \Delta V_HXC on vector
           !
           ! double bands @ gamma
           !
           DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
              !
              !$acc host_data use_device(dvg_exc_tmp)
              CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks_do),dvg_exc_tmp(:,lbnd+1,iks_do),psic,'Wave')
              !$acc end host_data
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,iks_do),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(z_rhs_vec_part2)
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),z_rhs_vec_part2(:,lbnd+1,iks),'Wave')
              !$acc end host_data
              !
           ENDDO
           !
           ! single band @ gamma
           !
           IF(MOD(nbnd_do,2) == 1) THEN
              !
              lbnd = nbnd_do
              !
              !$acc host_data use_device(dvg_exc_tmp)
              CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks_do),psic,'Wave')
              !$acc end host_data
              !
              !$acc parallel loop present(dvrs)
              DO ir = 1,dffts_nnr
                 psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,iks_do),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(z_rhs_vec_part2)
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),'Wave')
              !$acc end host_data
              !
           ENDIF
           !
           evc_copy(:,:) = evc_work
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
           ! now evc_copy contains the evc of the current spin channel, and
           ! evc_work contains the evc of the opposite spin channel
           !
           ! Compute the second part: dv_vv_mat
           ! recompute nbnd_do for the opposite spin channel
           !
           nbnd_do = 0
           DO lbnd = 1,band_group%nloc
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
              IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
           ENDDO
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
                 !$acc host_data use_device(evc_copy)
                 CALL double_invfft_gamma(dffts,npw,npwx,evc_copy(:,jbndp),evc_copy(:,kbndp),psic,'Wave')
                 !$acc end host_data
                 !
                 !$acc parallel loop present(dvrs)
                 DO ir = 1,dffts_nnr
                    psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
                 ENDDO
                 !$acc end parallel
                 !
                 !$acc host_data use_device(aux_g)
                 CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),aux_g(:,2),'Wave')
                 !$acc end host_data
                 !
                 reduce = 0._DP
                 !$acc loop reduction(+:reduce)
                 DO ig = 1, npw
                    reduce = reduce + REAL(evc_work(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                    &               + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,1))
                 ENDDO
                 !
                 dv_vv_mat(jbnd,lbnd,iks) = 2._DP * reduce
                 !
                 IF (gstart == 2) THEN
                    dv_vv_mat(jbnd,lbnd,iks) = dv_vv_mat(jbnd,lbnd,iks) &
                    & - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,1),KIND=DP)
                 ENDIF
                 !
                 reduce = 0._DP
                 !$acc loop reduction(+:reduce)
                 DO ig = 1, npw
                    reduce = reduce + REAL(evc_work(ig,ibndp),KIND=DP) * REAL(aux_g(ig,2),KIND=DP) &
                    &               + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,2))
                 ENDDO
                 !
                 dv_vv_mat(kbnd,lbnd,iks) = 2._DP * reduce
                 !
                 IF (gstart == 2) THEN
                    dv_vv_mat(kbnd,lbnd,iks) = dv_vv_mat(kbnd,lbnd,iks) &
                    & - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,2),KIND=DP)
                 ENDIF
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
                 !$acc host_data use_device(evc_copy)
                 CALL single_invfft_gamma(dffts,npw,npwx,evc_copy(:,jbndp),psic,'Wave')
                 !$acc end host_data
                 !
                 !$acc parallel loop present(dvrs)
                 DO ir = 1,dffts_nnr
                    psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
                 ENDDO
                 !$acc end parallel
                 !
                 !$acc host_data use_device(aux_g)
                 CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,1),'Wave')
                 !$acc end host_data
                 !
                 reduce = 0._DP
                 !$acc loop reduction(+:reduce)
                 DO ig = 1, npw
                    reduce = reduce + REAL(evc_work(ig,ibndp),KIND=DP) * REAL(aux_g(ig,1),KIND=DP) &
                    &               + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,1))
                 ENDDO
                 !
                 dv_vv_mat(jbnd,lbnd,iks) = 2._DP * reduce
                 !
                 IF (gstart == 2) THEN
                    dv_vv_mat(jbnd,lbnd,iks) = dv_vv_mat(jbnd,lbnd,iks) &
                    & - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,1),KIND=DP)
                 ENDIF
                 !
              ENDIF
              !
           ENDDO
           !
           CALL mp_sum(dv_vv_mat(:,:,iks),intra_bgrp_comm)
           !
           dpcpart(:,:) = (0._DP, 0._DP)
           !
           DO lbnd = 1,nbnd_do
              !
              DO jbnd = 1,nbndval-n_trunc_bands
                 !
                 factor = CMPLX(-dv_vv_mat(jbnd,lbnd,iks),KIND=DP)
                 !
                 !$acc host_data use_device(dvg_exc_tmp,dpcpart)
                 CALL ZAXPY(npw,factor,dvg_exc_tmp(:,lbnd,iks),1,dpcpart(:,jbnd),1)
                 !$acc end host_data
                 !
              ENDDO
              !
           ENDDO
           !
           CALL mp_sum(dpcpart,inter_bgrp_comm)
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
                 ibnd = band_group%l2g(lbnd)
                 z_rhs_vec_part2(ig,lbnd,iks) = z_rhs_vec_part2(ig,lbnd,iks)+dpcpart(ig,ibnd)
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
           ! Pc[k]*z_rhs_vec_part2
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
           CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec_part2(:,:,iks),(1._DP,0._DP))
           !
           !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part2)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part2(ig,lbnd,iks)
              ENDDO
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        DEALLOCATE(evc_copy)
        !
     ENDIF
     !
  ENDIF
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part2,z_rhs_vec_part2,band_group%nlocx,dotp)
  !
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  WRITE(stdout, "( 5x,'# Norm of z_rhs_vec p2  = | ', ES15.8, ' |' ) " ) SUM(REAL(dotp,KIND=DP))
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  !
  DEALLOCATE(dotp)
  DEALLOCATE(aux_g)
  DEALLOCATE(dv_vv_mat)
  DEALLOCATE(z_rhs_vec_part2)
  DEALLOCATE(dvrs)
  DEALLOCATE(dpcpart)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part3( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,n_trunc_bands,l_spin_flip
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! WORKSPACE
  !
  INTEGER :: ibnd,jbnd,iks,ir,ig,nbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part3(:,:,:),dotp(:)
  COMPLEX(DP), ALLOCATABLE :: ddvxc(:,:)
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part3(npwx*npol,band_group%nlocx,kpt_pool%nloc))
  ALLOCATE(ddvxc(dffts%nnr,nspin))
  !
  !$acc kernels present(z_rhs_vec_part3)
  z_rhs_vec_part3(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  IF(.NOT. l_spin_flip) THEN
     CALL compute_ddvxc_5p(dvg_exc_tmp, ddvxc)
  ELSE
     CALL compute_ddvxc_sf(dvg_exc_tmp, ddvxc)
  ENDIF
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
     ! ... Apply \Delta V_HXC
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
        !$acc parallel loop present(ddvxc)
        DO ir = 1,dffts_nnr
           psic(ir) = psic(ir)*CMPLX(REAL(ddvxc(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        !$acc host_data use_device(z_rhs_vec_part3)
        CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part3(:,lbnd,iks),z_rhs_vec_part3(:,lbnd+1,iks),'Wave')
        !$acc end host_data
        !
     ENDDO
     !
     ! single band @ gamma
     !
     IF(MOD(nbnd_do,2) == 1 ) THEN
        !
        lbnd = nbnd_do
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
        !
        !$acc parallel loop present(ddvxc)
        DO ir = 1,dffts_nnr
           psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(ddvxc(ir,current_spin),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        !$acc host_data use_device(z_rhs_vec_part3)
        CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part3(:,lbnd,iks),'Wave')
        !$acc end host_data
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
     ! Pc[k]*z_rhs_vec_part3
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec_part3(:,:,iks),(1._DP,0._DP))
     !
     !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part3)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part3(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
  ENDDO
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part3,z_rhs_vec_part3,band_group%nlocx,dotp)
  !
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  WRITE(stdout, "( 5x,'# Norm of z_rhs_vec p3  = | ', ES15.8, ' |' ) " ) SUM(REAL(dotp,KIND=DP))
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  !
  DEALLOCATE(dotp)
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
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : psic=>psic_d
#else
  USE wavefunctions,        ONLY : psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: ddvxc(dffts%nnr, nspin)
  !
  ! WORKSPACE
  !
  INTEGER :: iks,ir,indk
  REAL(DP), ALLOCATABLE :: aux_vxc(:,:,:), vxc(:,:), rdvrs(:,:)
  REAL(DP) :: etxc,vtxc,coeff
  TYPE (scf_type) :: a_rho
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  !
  CALL create_scf_type(a_rho)
  !
  ALLOCATE(aux_vxc(dffts%nnr, nspin, 5))
  ALLOCATE(vxc(dffts%nnr, nspin))
  ALLOCATE(dvrs(dffts%nnr,nspin))
  ALLOCATE(rdvrs(dffts%nnr, nspin))
  rdvrs(:,:) = 0._DP
  !
  ! Calculation of the charge density response
  !
  CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.FALSE.)
  !
  coeff = ddvxc_fd_coeff
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
     a_rho%of_r(:,:) = rho%of_r + REAL((indk-3),KIND=DP) * coeff * rdvrs
     !
     DO iks = 1, nspin
        !
        psic(:) = a_rho%of_r(:,iks)
        CALL fwfft ('Rho', psic, dffts)
        a_rho%of_g(:,iks) = psic(dffts%nl)
        !
     ENDDO
     !
     CALL v_xc( a_rho, rho_core, rhog_core, etxc, vtxc, vxc)
     !
     aux_vxc(:,:,indk) = vxc
     !
  ENDDO
  !
  ! compute ddvxc
  !
  DO iks = 1,nspin
     !
     DO ir = 1,dffts%nnr
        !
        ddvxc(ir,iks) = CMPLX( (-1._DP*aux_vxc(ir,iks,1)+16._DP*aux_vxc(ir,iks,2) &
        &                       -30._DP*aux_vxc(ir,iks,3)+16._DP*aux_vxc(ir,iks,4) &
        &                       -1._DP*aux_vxc(ir,iks,5)), 0._DP, KIND=DP ) / (12._DP*coeff*coeff)
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(aux_vxc)
  DEALLOCATE(vxc)
  DEALLOCATE(dvrs)
  DEALLOCATE(rdvrs)
  !
  CALL destroy_scf_type(a_rho)
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
  USE mp,                   ONLY : mp_barrier
  USE noncollin_module,     ONLY : npol,nspin_gga
  USE xc_lib,               ONLY : xclib_dft_is
  USE fft_base,             ONLY : dffts
  USE gvect,                ONLY : g
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : nlcc_any
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_world,             ONLY : world_comm
  USE westcom,              ONLY : wbse_save_dir,sf_kernel,l_spin_flip_alda0,spin_flip_cut2,&
                                 & l_print_spin_flip_kernel
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
  ! WORKSPACE
  !
  INTEGER :: ir,is,is1,dffts_nnr
  REAL(DP) :: tmp1,tmp2
  COMPLEX(DP), ALLOCATABLE :: drho_sf_copy(:,:),dvsf(:,:)
  CHARACTER(LEN=:), ALLOCATABLE :: fname
  COMPLEX(DP), ALLOCATABLE :: drho_sf(:,:)
  !
  IF(nlcc_any) CALL errore('compute_ddvxc_sf', 'nlcc_any is not supported', 1)
  !
  ALLOCATE(drho_sf(dffts%nnr,2))
  ALLOCATE(drho_sf_copy(dffts%nnr,2))
  ALLOCATE(dvsf(dffts%nnr,2))
  !
  dffts_nnr = dffts%nnr
  !
  CALL wbse_calc_dens(dvg_exc_tmp,drho_sf,.TRUE.)
  !
  DO ir = 1,dffts_nnr
     tmp1 = REAL(drho_sf(ir,1),KIND=DP)**2
     tmp2 = REAL(drho_sf(ir,2),KIND=DP)**2
     drho_sf_copy(ir,1) = CMPLX(tmp1 + tmp2,0._DP,KIND=DP)
     drho_sf_copy(ir,2) = - CMPLX(tmp1 + tmp2,0._DP,KIND=DP)
  ENDDO
  !
  ! divide rho_diff
  !
  DO ir = 1,dffts_nnr
     !
     drho_sf_copy(ir,1) = drho_sf_copy(ir,1) / rho%of_r(ir,2)
     drho_sf_copy(ir,2) = drho_sf_copy(ir,2) / rho%of_r(ir,2)
     !
     IF(ABS(rho%of_r(ir,2)) < spin_flip_cut2) THEN
        drho_sf_copy(ir,1) = (0._DP, 0._DP)
        drho_sf_copy(ir,2) = (0._DP, 0._DP)
     ENDIF
     !
  ENDDO
  !
  ! part 1
  !
  DO ir = 1,dffts_nnr
     ddvxc(ir,1) = - sf_kernel(ir) * drho_sf_copy(ir,1)
     ddvxc(ir,2) = - sf_kernel(ir) * drho_sf_copy(ir,2)
  ENDDO
  !
  ! part 2
  !
  dvsf(:,:) = (0._DP, 0._DP)
  DO is = 1,nspin
     DO is1 = 1,nspin
        dvsf(:,is) = dvsf(:,is) + dmuxc(:,is,is1) * drho_sf_copy(:,is1)
     ENDDO
  ENDDO
  !
  IF(.NOT. l_spin_flip_alda0) THEN
     IF(xclib_dft_is('gradient')) THEN
        CALL dgradcorr(dffts, rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
             & drho_sf_copy, nspin, nspin_gga, g, dvsf)
     ENDIF
  ENDIF
  !
  ! part 1 and part 2 together
  !
  DO ir = 1,dffts_nnr
     !
     ddvxc(ir,1) = ddvxc(ir,1) + dvsf(ir,1)
     ddvxc(ir,2) = ddvxc(ir,2) + dvsf(ir,2)
     !
  ENDDO
  !
  ! print spin flip kernel (keep it for now)
  !
  IF(l_print_spin_flip_kernel) THEN
     !
     CALL mp_barrier(world_comm)
     !
     !$acc update host(rho%of_r)
     fname=TRIM(wbse_save_dir)//'/rho_diff.cube'
     CALL write_wfc_cube_r(dffts, fname, rho%of_r(ir,2))
     !
     !$acc update host(drho_sf)
     fname=TRIM(wbse_save_dir)//'/drho_sf_up.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf(:,1),KIND=DP))
     !
     !$acc update host(drho_sf)
     fname=TRIM(wbse_save_dir)//'/drho_sf_down.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf(:,2),KIND=DP))
     !
     !$acc update host(drho_sf_copy)
     fname=TRIM(wbse_save_dir)//'/drho_sq_rho_diff_up.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf_copy(:,1),KIND=DP))
     !
     !$acc update host(drho_sf_copy)
     fname=TRIM(wbse_save_dir)//'/drho_sq_rho_diff_down.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(drho_sf_copy(:,2),KIND=DP))
     !
     !$acc update host(ddvxc)
     fname=TRIM(wbse_save_dir)//'/ddvxc_up.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(ddvxc(:,1),KIND=DP))
     !
     !$acc update host(ddvxc)
     fname=TRIM(wbse_save_dir)//'/ddvxc_down.cube'
     CALL write_wfc_cube_r(dffts, fname, REAL(ddvxc(:,2),KIND=DP))
     !
     CALL mp_barrier(world_comm)
     !
  ENDIF
  !
  DEALLOCATE(drho_sf)
  DEALLOCATE(drho_sf_copy)
  DEALLOCATE(dvsf)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part4( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip,&
                                 & l_slow_tddft_k1d
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,intra_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! WORKSPACE
  !
  INTEGER :: ig,lbnd,jbnd,jbndp,iks_do,nbnd_do
  INTEGER :: ibnd,nbndval,flnbndval
  INTEGER :: iks
  REAL(DP) :: reduce
  COMPLEX(DP) :: factor
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part4(:,:,:),tmp_vec(:,:),dotp(:)
  REAL(DP), ALLOCATABLE :: dv_vv_mat(:,:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  COMPLEX(DP), ALLOCATABLE :: dpcpart(:,:)
  !
  ALLOCATE(z_rhs_vec_part4(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(dv_vv_mat(nbndval0x-n_trunc_bands, band_group%nlocx))
  ALLOCATE(tmp_vec(npwx*npol, band_group%nlocx))
  ALLOCATE(dpcpart(npwx*npol, nbndval0x-n_trunc_bands))
  !
  !$acc kernels present(z_rhs_vec_part4)
  z_rhs_vec_part4(:,:,:) = (0._DP, 0._DP)
  !$acc end kernels
  !
  dv_vv_mat(:,:) = (0._DP, 0._DP)
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
     ! Compute the first part
     !
     CALL hybrid_kernel_term4(current_spin, dvg_exc_tmp, z_rhs_vec_part4(:,:,iks), l_spin_flip)
     !
     ! Compute the second part: dv_vv_mat
     !
     tmp_vec(:,:) = (0._DP, 0._DP)
     !
     IF(l_slow_tddft_k1d) THEN
        CALL hybrid_kernel_term1_slow(current_spin, dvg_exc_tmp, tmp_vec, l_spin_flip)
     ELSE
        CALL bse_kernel_gamma(current_spin, dvg_exc_tmp, tmp_vec, l_spin_flip)
     ENDIF
     !
     DO jbnd = 1, nbndval - n_trunc_bands
        !
        jbndp = jbnd + n_trunc_bands
        !
        DO lbnd = 1, nbnd_do
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1, npw
              reduce = reduce + REAL(evc_work(ig,jbndp),KIND=DP) * REAL(tmp_vec(ig,lbnd),KIND=DP) &
              &               + AIMAG(evc_work(ig,jbndp)) * AIMAG(tmp_vec(ig,lbnd))
           ENDDO
           !
           dv_vv_mat(jbnd,lbnd) = 2._DP * reduce
           !
           IF( gstart == 2 ) THEN
              dv_vv_mat(jbnd,lbnd) = dv_vv_mat(jbnd,lbnd) &
              & - REAL(evc_work(1,jbndp),KIND=DP)*REAL(tmp_vec(1,lbnd),KIND=DP)
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
     CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
     !
     dpcpart(:,:) = (0._DP, 0._DP)
     !
     DO jbnd = 1, nbndval - n_trunc_bands
        DO lbnd = 1, nbnd_do
           !
           factor = CMPLX(-dv_vv_mat(jbnd,lbnd),KIND=DP)
           !
           !$acc host_data use_device(dvg_exc_tmp,dpcpart)
           CALL ZAXPY(npw,factor,dvg_exc_tmp(:,lbnd,iks),1,dpcpart(:,jbnd),1)
           !$acc end host_data
           !
        ENDDO
     ENDDO
     !
     CALL mp_sum(dpcpart, inter_bgrp_comm)
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
           ibnd = band_group%l2g(lbnd)
           z_rhs_vec_part4(ig,lbnd,iks) = z_rhs_vec_part4(ig,lbnd,iks)+dpcpart(ig,ibnd)
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
     ! Pc[k]*z_rhs_vec_part4
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec_part4(:,:,iks),(1._DP,0._DP))
     !
     !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part4)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part4(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  ALLOCATE(dotp(nspin))
  !
  CALL wbse_dot(z_rhs_vec_part4,z_rhs_vec_part4,band_group%nlocx,dotp)
  !
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  WRITE(stdout, "( 5x,'# Norm of z_rhs_vec p4  = | ', ES15.8, ' |' ) " ) SUM(REAL(dotp,KIND=DP))
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  !
  DEALLOCATE(dotp)
  DEALLOCATE(z_rhs_vec_part4)
  DEALLOCATE(dv_vv_mat)
  DEALLOCATE(tmp_vec)
  DEALLOCATE(dpcpart)
  !
END SUBROUTINE
