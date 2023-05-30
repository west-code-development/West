!
! Copyright (C) 2015 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file:
! Yu Jin, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_build_rhs_zvector_eq( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec)
  !-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE io_push,              ONLY : io_push_title
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx
  USE noncollin_module,     ONLY : npol
  USE klist,                ONLY : nks
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : l_bse,l_hybrid_tddft,l_bse_triplet,nbndval0x,&
                                 & n_trunc_bands,l_spin_flip,l_spin_flip_kernel
  USE distribution_center,  ONLY : kpt_pool,band_group
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks)
  COMPLEX(DP), INTENT(IN) :: drhox1(dffts%nnr, nspin), drhox2(dffts%nnr, nspin)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  REAL(DP), EXTERNAL :: get_clock
  REAL(DP) :: time_spent(2)
  !
  CALL start_clock( 'rhs_z' )
  time_spent(1) = get_clock( 'rhs_z' )
  !
  CALL io_push_title("Forces : Build the RHS of the Z vector Equation")
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
        ! print the time used for part3
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
!  IF( l_hybrid_tddft ) THEN
!     !
!     CALL rhs_zvector_part4( dvg_exc_tmp, z_rhs_vec )
!     !
!     ! print the time used for part4
!     time_spent(2) = get_clock( 'rhs_z' )
!     CALL wbse_forces_time(time_spent)
!     time_spent(1) = get_clock( 'rhs_z' )
!     !
!  ELSEIF( l_bse ) THEN
!     !
!     STOP
!     !
!  ENDIF
  !
  CALL stop_clock( 'rhs_z' )
  !
END SUBROUTINE

!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part1( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                 & l_hybrid_tddft,l_spin_flip
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : npwx,npw
  USE klist,                ONLY : nks
  USE pwcom,                ONLY : isk,igk_k,lsda,current_k,ngk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : noncolin,npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,&
                                 & inter_bgrp_comm,intra_bgrp_comm
  USE wbse_dv,              ONLY : wbse_dv_of_drho
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks)
  COMPLEX(DP), INTENT(IN) :: drhox1(dffts%nnr, nspin), drhox2(dffts%nnr, nspin)
  COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! WORKSPACE
  !
  LOGICAL :: lrpa
  INTEGER :: ibnd,jbnd,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part1(:,:,:),dot_out(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
#if !defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE :: drhox(:,:)
#endif
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part1(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  !$acc kernels present(z_rhs_vec_part1)
  z_rhs_vec_part1(:,:,:) = (0._DP, 0._DP)
  !$acc end kernels
  !
  ! Compute drhox
  !
  ALLOCATE(drhox(dffts%nnr, nspin))
  drhox(:,:) = (0._DP, 0._DP)
  !
  DO iks = 1,kpt_pool%nloc
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
  CALL mp_sum(drhox,inter_pool_comm)
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
     IF (gamma_only) THEN
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
           lbnd  = nbnd_do
           ibnd  = band_group%l2g(lbnd)+n_trunc_bands
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
           !
           !$acc parallel loop present(drhox)
           DO ir=1,dffts_nnr
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
           !$acc parallel loop present(drhox)
           DO ir = 1,dffts_nnr
              psic(ir) = psic(ir)*drhox(ir,current_spin)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(z_rhs_vec_part1)
           CALL single_fwfft_k(dffts,npw,npwx,psic,z_rhs_vec_part1(:,lbnd,iks),'Wave',igk_k(:,current_k))
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
              !$acc parallel loop present(drhox)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*drhox(ir,current_spin)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(z_rhs_vec_part1)
              CALL single_fwfft_k(dffts,npw,npwx,psic,z_rhs_vec_part1(npwx+1:npwx*2,lbnd,iks),'Wave',igk_k(:,current_k))
              !$acc end host_data
              !
           ENDDO
        ENDIF
        !
     ENDIF
     !
     IF (l_bse) THEN
        !
        STOP
        !
     ENDIF
     !
     IF (l_hybrid_tddft) THEN
        !
        STOP
        !
     ENDIF
     !
     IF (gamma_only) THEN
        IF (gstart==2) THEN
           !$acc parallel loop present(z_rhs_vec_part1)
           DO lbnd = 1,nbnd_do
              z_rhs_vec_part1(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part1(1,lbnd,iks),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
     ENDIF
     !
     ! Pc[k]*z_rhs_vec_part1
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
  ALLOCATE(dot_out(nks))
  dot_out(:) = (0._DP, 0._DP)
  !
  CALL wbse_dot(z_rhs_vec_part1,z_rhs_vec_part1,band_group%nlocx,dot_out)
  !
  WRITE(stdout, "( /,5x,'                          *-------------*' ) "  )
  WRITE(stdout, "(   5x,'# Norm of z_rhs_vec p1  = | ', ES11.4, ' |' ) " ) SUM(REAL(dot_out))
  WRITE(stdout, "(   5x,'                          *-------------*' ) "  )
  !
  DEALLOCATE(dot_out)
  !
#if !defined(__CUDA)
  DEALLOCATE(drhox)
#endif
  DEALLOCATE(z_rhs_vec_part1)
  !
END SUBROUTINE

!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part2( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                 & l_hybrid_tddft,l_spin_flip
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : npwx,npw
  USE klist,                ONLY : nks,igk_k
  USE pwcom,                ONLY : isk,lsda,nkstot,current_k,ngk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : noncolin,npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,nbgrp,my_bgrp_id,&
                                 & inter_bgrp_comm,intra_bgrp_comm,world_comm
  USE wbse_dv,              ONLY : wbse_dv_of_drho
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
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
  INTEGER :: ibnd,ibndp,jbnd,jbndp,kbnd,kbndp,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part2(:,:,:),drhox(:,:),aux_g(:,:),dot_out(:)
  REAL(DP), ALLOCATABLE :: dv_vv_mat(:,:,:)
  REAL(DP) :: reduce
  COMPLEX(DP) :: factor
  INTEGER, PARAMETER :: flks(2) = [2,1]
#if !defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: dpcpart(:,:)
#endif
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part2(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(dv_vv_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks))
  ALLOCATE(aux_g(npwx*npol,2))
  ALLOCATE(dpcpart(npwx*npol, nbndval0x-n_trunc_bands))
  dv_vv_mat(:,:,:) = (0._DP, 0._DP)
  aux_g(:,:) = (0._DP, 0._DP)
  !
  !$acc kernels present(z_rhs_vec_part2)
  z_rhs_vec_part2(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
#if !defined(__CUDA)
  ALLOCATE(dvrs(dffts%nnr,nspin))
#endif
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
        IF(gamma_only) THEN
           !
           ! double bands @ gamma
           !
           DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
              !
              CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),dvg_exc_tmp(:,lbnd+1,iks),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir=1,dffts_nnr
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
              CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave')
              !
              !$acc parallel loop present(dvrs)
              DO ir=1,dffts_nnr
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
                 DO ir=1,dffts_nnr
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
                                    + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,1))
                 ENDDO
                 !
                 dv_vv_mat(jbnd,lbnd,iks) = 2._DP * reduce
                 !
                 IF (gstart == 2) THEN 
                    dv_vv_mat(jbnd,lbnd,iks) = dv_vv_mat(jbnd,lbnd,iks) - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,1),KIND=DP)
                 ENDIF
                 !
                 reduce = 0._DP
                 !$acc loop reduction(+:reduce)
                 DO ig = 1, npw
                    reduce = reduce + REAL(evc_work(ig,ibndp),KIND=DP) * REAL(aux_g(ig,2),KIND=DP) &
                                    + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,2))
                 ENDDO
                 !
                 dv_vv_mat(kbnd,lbnd,iks) = 2._DP * reduce
                 !
                 IF (gstart == 2) THEN 
                    dv_vv_mat(kbnd,lbnd,iks) = dv_vv_mat(kbnd,lbnd,iks) - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,2),KIND=DP)
                 ENDIF
                 !
              ENDDO
              ! 
              ! single band @ gamma
              !
              IF(MOD((nbndval-n_trunc_bands),2) == 1) THEN
                 !
                 jbnd  = nbndval - n_trunc_bands
                 jbndp   = jbnd + n_trunc_bands
                 !
                 CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,jbndp),psic,'Wave')
                 !
                 !$acc parallel loop present(dvrs)
                 DO ir=1,dffts_nnr
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
                                    + AIMAG(evc_work(ig,ibndp)) * AIMAG(aux_g(ig,1))
                 ENDDO
                 !
                 dv_vv_mat(jbnd,lbnd,iks) = 2._DP * reduce
                 !
                 IF (gstart == 2) THEN
                    dv_vv_mat(jbnd,lbnd,iks) = dv_vv_mat(jbnd,lbnd,iks) - REAL(evc_work(1,ibndp),KIND=DP)*REAL(aux_g(1,1),KIND=DP)
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
           CALL mp_sum(dpcpart(:,:),inter_bgrp_comm)
           !
        ENDIF
        !
        !$acc parallel loop collapse(2) present(z_rhs_vec_part2,dpcpart)
        DO lbnd = 1,nbnd_do
           ibnd = band_group%l2g(lbnd)
           DO ig = 1,npw
              z_rhs_vec_part2(ig,lbnd,iks) = z_rhs_vec_part2(ig,lbnd,iks)+dpcpart(ig,ibnd)
           ENDDO
        ENDDO
        !$acc end parallel
        !
        IF(gamma_only) THEN
           IF(gstart == 2) THEN
              !$acc parallel loop present(z_rhs_vec_part2)
              DO lbnd = 1,nbnd_do
                 z_rhs_vec_part2(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part2(1,lbnd,iks),KIND=DP),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
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
  ELSE
     !
     STOP
     !
  ENDIF  
  !
#if !defined(__CUDA)
  DEALLOCATE(dvrs)
  DEALLOCATE(dpcpart)
#endif
  !
  ALLOCATE(dot_out(nks))
  dot_out(:) = (0._DP, 0._DP)
  !
  CALL wbse_dot(z_rhs_vec_part2,z_rhs_vec_part2,band_group%nlocx,dot_out)
  !
  WRITE(stdout, "( /,5x,'                          *-------------*' ) "  )
  WRITE(stdout, "(   5x,'# Norm of z_rhs_vec p2  = | ', ES11.4, ' |' ) " ) SUM(REAL(dot_out))
  WRITE(stdout, "(   5x,'                          *-------------*' ) "  )
  !
  DEALLOCATE(dot_out)
  DEALLOCATE(aux_g)
  DEALLOCATE(dv_vv_mat)
  DEALLOCATE(z_rhs_vec_part2)
  !
END SUBROUTINE

!-----------------------------------------------------------------------
SUBROUTINE rhs_zvector_part3( dvg_exc_tmp, z_rhs_vec )
  !-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                 & l_hybrid_tddft,l_spin_flip
  USE lsda_mod,             ONLY : current_spin,nspin,lsda
  USE wvfct,                ONLY : npwx,npw
  USE klist,                ONLY : xk,nks,igk_k
  USE pwcom,                ONLY : isk,lsda,nkstot,current_k,ngk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : noncolin,npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,nbgrp,my_bgrp_id,&
                                 & inter_bgrp_comm,intra_bgrp_comm,world_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY :
using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
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
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part3(:,:,:),dot_out(:)
#if !defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE :: ddvxc(:,:)
#endif
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(z_rhs_vec_part3(npwx*npol,band_group%nlocx,kpt_pool%nloc))
  !
  !$acc kernels present(z_rhs_vec_part3)
  z_rhs_vec_part3(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
#if !defined(__CUDA)
  ALLOCATE(ddvxc(dffts%nnr,nspin))
#endif
  !
  IF(.NOT. l_spin_flip) THEN
     !
     CALL compute_ddvxc_5p( dvg_exc_tmp, ddvxc )
     !
  ELSE
     !
     STOP
     !
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
     IF (gamma_only) THEN
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
           DO ir=1,dffts_nnr
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
           !
           ibnd  = band_group%l2g(lbnd)+n_trunc_bands
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
           !
           !$acc parallel loop present(ddvxc)
           DO ir=1,dffts_nnr
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
           !$acc parallel loop present(ddvxc)
           DO ir = 1,dffts_nnr
              psic(ir) = psic(ir)*ddvxc(ir,current_spin)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(z_rhs_vec_part3)
           CALL single_fwfft_k(dffts,npw,npwx,psic,z_rhs_vec_part3(:,lbnd,iks),'Wave',igk_k(:,current_k))
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
              !$acc parallel loop present(ddvxc)
              DO ir = 1,dffts_nnr
                 psic(ir) = psic(ir)*ddvxc(ir,current_spin)
              ENDDO
              !$acc end parallel
              !
              !$acc host_data use_device(z_rhs_vec_part3)
              CALL single_fwfft_k(dffts,npw,npwx,psic,z_rhs_vec_part3(npwx+1:npwx*2,lbnd,iks),'Wave',igk_k(:,current_k))
              !$acc end host_data
              !
           ENDDO
        ENDIF
        !
     ENDIF
     !
     IF(gamma_only) THEN
        IF(gstart == 2) THEN
           !$acc parallel loop present(z_rhs_vec_part3)
           DO lbnd = 1,nbnd_do
              z_rhs_vec_part3(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part3(1,lbnd,iks),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
     ENDIF
     !
     ! Pc[k]*z_rhs_vec_part3
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
  ALLOCATE(dot_out(nks))
  dot_out = (0._DP, 0._DP)
  !
  CALL wbse_dot(z_rhs_vec_part3,z_rhs_vec_part3,band_group%nlocx,dot_out)
  !
  WRITE(stdout, "( /,5x,'                          *-------------*' ) "  )
  WRITE(stdout, "(   5x,'# Norm of z_rhs_vec p3  = | ', ES11.4, ' |' ) " ) SUM(REAL(dot_out))
  WRITE(stdout, "(   5x,'                          *-------------*' ) "  )
  !
  DEALLOCATE(dot_out)
  DEALLOCATE(ddvxc)
  DEALLOCATE(z_rhs_vec_part3)
  !
END SUBROUTINE


!-----------------------------------------------------------------------
SUBROUTINE compute_ddvxc_5p( dvg_exc_tmp, ddvxc )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE lsda_mod,             ONLY : current_spin,nspin,lsda
  USE wvfct,                ONLY : npwx
  USE klist,                ONLY : nks,wk,igk_k
  USE pwcom,                ONLY : isk,lsda,current_k,ngk
  USE mp,                   ONLY : mp_barrier
  USE noncollin_module,     ONLY : noncolin,npol
  USE fft_base,             ONLY : dffts
  USE scf,                  ONLY : rho,rho_core,rhog_core,scf_type,&
                                 & create_scf_type,destroy_scf_type
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE fft_interfaces,       ONLY : fwfft
  USE mp_world,             ONLY : world_comm
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
  COMPLEX(DP), INTENT(INOUT) :: ddvxc(dffts%nnr, nspin)
  !
  ! WORKSPACE
  !
  INTEGER  :: iks,ir,indk
  REAL(DP), ALLOCATABLE :: aux_vxc(:,:,:), vxc(:,:), rdvrs(:,:)
  REAL(DP) :: etxc,vtxc,coeff
  TYPE (scf_type) :: a_rho
#if !defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
#endif
  !
  CALL mp_barrier(world_comm)
  !
  CALL create_scf_type(a_rho) 
  !
  ALLOCATE(aux_vxc(dffts%nnr, nspin, 5))
  aux_vxc(:,:,:) = 0._DP
  ALLOCATE(vxc(dffts%nnr, nspin))
  vxc(:,:) = 0._DP
  !
  ddvxc(:,:) = (0._DP, 0._DP)
  !
#if !defined(__CUDA)
  ALLOCATE(dvrs(dffts%nnr,nspin))
#endif
  !
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
     a_rho%of_r(:,:) = rho%of_r(:,:) + REAL((indk-3),KIND=DP) * coeff * rdvrs(:,:)
     !
     DO iks = 1, nspin
        !
        psic(:) = a_rho%of_r(:,iks)
        CALL fwfft ('Rho', psic, dffts)
        a_rho%of_g(:,iks) = psic(dffts%nl(:))
        !
     ENDDO
     !
     CALL v_xc( a_rho, rho_core, rhog_core, etxc, vtxc, vxc)
     !
     aux_vxc(:,:,indk) = vxc(:,:)
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
                                -30._DP*aux_vxc(ir,iks,3)+16._DP*aux_vxc(ir,iks,4) &
                                -1._DP*aux_vxc(ir,iks,5)), 0._DP, KIND=DP ) / (12._DP*coeff*coeff)
        !
     ENDDO
     !
  ENDDO
  !
#if !defined(__CUDA)
  DEALLOCATE(dvrs)
#endif
  DEALLOCATE(rdvrs)
  DEALLOCATE(aux_vxc)
  DEALLOCATE(vxc)
  CALL destroy_scf_type(a_rho)
  !
  CALL mp_barrier(world_comm)
  !
END SUBROUTINE
