!
! Copyright (C) 2015-2023 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
!-----------------------------------------------------------------------
MODULE west_gpu_data
   !-----------------------------------------------------------------------
   !
   USE kinds,         ONLY : DP,i8b
#if defined(__CUDA)
   USE becmod_gpum,   ONLY : bec_type_d
   USE cublas
   USE cusolverdn
   !
   IMPLICIT NONE
   !
   ! Linsolve
   !
   LOGICAL, ALLOCATABLE :: is_conv(:)
   INTEGER, ALLOCATABLE :: ibnd_todo(:)
   REAL(DP), ALLOCATABLE :: eu(:)
   REAL(DP), ALLOCATABLE :: a(:)
   REAL(DP), ALLOCATABLE :: c(:)
   REAL(DP), ALLOCATABLE :: rho(:)
   REAL(DP), ALLOCATABLE :: rhoold(:)
   COMPLEX(DP), ALLOCATABLE :: g(:,:)
   COMPLEX(DP), ALLOCATABLE :: t(:,:)
   COMPLEX(DP), ALLOCATABLE :: h(:,:)
   !$acc declare device_resident(is_conv,ibnd_todo,eu,a,c,rho,rhoold,g,t,h)
   !
   ! GW
   !
   COMPLEX(DP), PINNED, ALLOCATABLE :: dvpsi_h(:,:)
   !
   ! Macropol
   !
   REAL(DP), ALLOCATABLE :: gk(:,:)
   REAL(DP), ALLOCATABLE :: deff(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
   COMPLEX(DP), ALLOCATABLE :: ps2(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: psc(:,:,:,:)
   COMPLEX(DP), ALLOCATABLE :: dvkb(:,:)
   COMPLEX(DP), ALLOCATABLE :: work(:,:)
   !$acc declare device_resident(gk,deff,deff_nc,ps2,psc,dvkb,work)
   REAL(DP), DEVICE, POINTER :: becp1_d_r_d(:,:)
   REAL(DP), DEVICE, POINTER :: becp2_d_r_d(:,:)
   COMPLEX(DP), DEVICE, POINTER :: becp1_d_k_d(:,:)
   COMPLEX(DP), DEVICE, POINTER :: becp2_d_k_d(:,:)
   COMPLEX(DP), DEVICE, POINTER :: becp1_d_nc_d(:,:,:)
   COMPLEX(DP), DEVICE, POINTER :: becp2_d_nc_d(:,:,:)
   TYPE(bec_type_d), TARGET :: becp1_d ! (nkb,m)
   TYPE(bec_type_d), TARGET :: becp2_d ! (nkb,m)
   !
   ! Chi invert
   !
   INTEGER(i8b) :: l_inv
   INTEGER(i8b) :: l_inv_h
#if CUDA_VERSION > 11010 && (__PGIC__ > 21 || (__PGIC__ == 21 && __PGIC_MINOR > 2))
   INTEGER(i8b), ALLOCATABLE :: piv(:)
#else
   INTEGER, ALLOCATABLE :: piv(:)
#endif
   !$acc declare device_resident(piv)
   REAL(DP), ALLOCATABLE :: work_r_h(:)
   REAL(DP), ALLOCATABLE :: work_r(:)
   REAL(DP), ALLOCATABLE :: x_r(:,:)
   REAL(DP), ALLOCATABLE :: wh_r(:,:)
   REAL(DP), ALLOCATABLE :: wl_r(:,:)
   REAL(DP), ALLOCATABLE :: temph_r(:,:)
   REAL(DP), ALLOCATABLE :: templ_r(:,:)
   REAL(DP), ALLOCATABLE :: tempt_r(:,:)
   !$acc declare device_resident(work_r,x_r,wh_r,wl_r,temph_r,templ_r,tempt_r)
   COMPLEX(DP), ALLOCATABLE :: work_c_h(:)
   COMPLEX(DP), ALLOCATABLE :: work_c(:)
   COMPLEX(DP), ALLOCATABLE :: x_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: wh_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: wl_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: temph_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: templ_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: tempt_c(:,:)
   !$acc declare device_resident(work_c,x_c,wh_c,wl_c,temph_c,templ_c,tempt_c)
   !
   ! Lanczos
   !
   COMPLEX(DP), ALLOCATABLE :: r(:,:)
   !$acc declare device_resident(r)
   !
   ! BSE
   !
   REAL(DP), ALLOCATABLE :: factors(:)
   REAL(DP), ALLOCATABLE :: raux1(:)
   REAL(DP), ALLOCATABLE :: raux2(:)
   COMPLEX(DP), ALLOCATABLE :: caux1(:,:)
   COMPLEX(DP), ALLOCATABLE :: caux2(:,:)
   COMPLEX(DP), ALLOCATABLE :: caux3(:,:)
   COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
   COMPLEX(DP), ALLOCATABLE :: hevc1(:,:)
   COMPLEX(DP), ALLOCATABLE :: psic2(:)
   COMPLEX(DP), PINNED, ALLOCATABLE :: dvaux(:,:)
   COMPLEX(DP), ALLOCATABLE :: dvhart(:)
   !$acc declare device_resident(factors,raux1,raux2,caux3,hevc1,psic2,dvhart)
   COMPLEX(DP), PINNED, ALLOCATABLE :: gaux(:)
   !
   ! Workspace
   !
   REAL(DP), ALLOCATABLE :: tmp_r(:)
   REAL(DP), ALLOCATABLE :: tmp_r3(:,:,:)
   REAL(DP), ALLOCATABLE :: ps_r(:,:)
   COMPLEX(DP), ALLOCATABLE :: tmp_c(:)
   COMPLEX(DP), ALLOCATABLE :: tmp_c3(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
   !$acc declare device_resident(tmp_r3,ps_r,tmp_c3,ps_c)
   INTEGER, DEVICE, POINTER :: dfft_nl_d(:)
   INTEGER, DEVICE, POINTER :: dfft_nlm_d(:)
   TYPE(cusolverDnHandle) :: cusolv_h
   TYPE(cusolverDnParams) :: cusolv_p
   !
   CONTAINS
   !
   !-----------------------------------------------------------------------
   SUBROUTINE west_gpu_start()
   !-----------------------------------------------------------------------
   !
   USE control_flags,         ONLY : use_gpu
   USE io_global,             ONLY : stdout
   !
   IMPLICIT NONE
   !
   ! Workspace
   !
   INTEGER :: istat
   !
   LOGICAL, EXTERNAL :: check_gpu_support
   !
   use_gpu = check_gpu_support
   IF(.NOT. use_gpu) CALL errore('gpu_start','use_gpu .FALSE.',1)
   !
   istat = cusolverDnCreate(cusolv_h)
   IF(istat /= 0) CALL errore('gpu_start','coSOLVER init failed',istat)
   !
   istat = cusolverDnCreateParams(cusolv_p)
   IF(istat /= 0) CALL errore('gpu_start','coSOLVER params init failed',istat)
   !
   WRITE(stdout,'(/,5X,A)') 'GPU acceleration enabled'
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE west_gpu_end()
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   ! Workspace
   !
   INTEGER :: istat
   !
   istat = cusolverDnDestroyParams(cusolv_p)
   istat = cusolverDnDestroy(cusolv_h)
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_gpu()
   !-----------------------------------------------------------------------
   !
   USE fft_base,              ONLY : dffts
   USE pwcom,                 ONLY : wg,ngk,nks
   USE westcom,               ONLY : igq_q
   USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d
   USE wvfct_gpum,            ONLY : using_et,using_et_d
   USE becmod_subs_gpum,      ONLY : using_becp_auto,using_becp_d_auto
   !
   IMPLICIT NONE
   !
   dfft_nl_d => dffts%nl_d
   dfft_nlm_d => dffts%nlm_d
   !
   CALL using_et(2)
   CALL using_et_d(0)
   CALL using_becp_auto(2)
   CALL using_becp_d_auto(0)
   !
   IF(nks == 1) THEN
      CALL using_evc(2)
      CALL using_evc_d(0)
   ENDIF
   !
   !$acc enter data copyin(wg,ngk,igq_q)
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_gpu()
   !-----------------------------------------------------------------------
   !
   USE pwcom,                 ONLY : wg,ngk
   USE westcom,               ONLY : igq_q
   !
   IMPLICIT NONE
   !
   IF(ASSOCIATED(dfft_nl_d)) THEN
      NULLIFY(dfft_nl_d)
   ENDIF
   IF(ASSOCIATED(dfft_nlm_d)) THEN
      NULLIFY(dfft_nlm_d)
   ENDIF
   !
   !$acc exit data delete(wg,ngk,igq_q)
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_linsolve_gpu(nbndloc)
   !-----------------------------------------------------------------------
   !
   USE noncollin_module,      ONLY : npol
   USE pwcom,                 ONLY : npwx
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER, INTENT(IN) :: nbndloc
   !
   ALLOCATE(is_conv(nbndloc))
   ALLOCATE(ibnd_todo(nbndloc))
   ALLOCATE(eu(nbndloc))
   ALLOCATE(a(nbndloc))
   ALLOCATE(c(nbndloc))
   ALLOCATE(rho(nbndloc))
   ALLOCATE(rhoold(nbndloc))
   ALLOCATE(g(npwx*npol,nbndloc))
   ALLOCATE(t(npwx*npol,nbndloc))
   ALLOCATE(h(npwx*npol,nbndloc))
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_linsolve_gpu()
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   IF(ALLOCATED(is_conv)) THEN
      DEALLOCATE(is_conv)
   ENDIF
   IF(ALLOCATED(ibnd_todo)) THEN
      DEALLOCATE(ibnd_todo)
   ENDIF
   IF(ALLOCATED(eu)) THEN
      DEALLOCATE(eu)
   ENDIF
   IF(ALLOCATED(a)) THEN
      DEALLOCATE(a)
   ENDIF
   IF(ALLOCATED(c)) THEN
      DEALLOCATE(c)
   ENDIF
   IF(ALLOCATED(rho)) THEN
      DEALLOCATE(rho)
   ENDIF
   IF(ALLOCATED(rhoold)) THEN
      DEALLOCATE(rhoold)
   ENDIF
   IF(ALLOCATED(g)) THEN
      DEALLOCATE(g)
   ENDIF
   IF(ALLOCATED(t)) THEN
      DEALLOCATE(t)
   ENDIF
   IF(ALLOCATED(h)) THEN
      DEALLOCATE(h)
   ENDIF
   IF(ALLOCATED(ps_r)) THEN
      DEALLOCATE(ps_r)
   ENDIF
   IF(ALLOCATED(ps_c)) THEN
      DEALLOCATE(ps_c)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_gw_gpu(nlocx,nloc)
   !-----------------------------------------------------------------------
   !
   USE control_flags,         ONLY : gamma_only
   USE noncollin_module,      ONLY : npol
   USE pwcom,                 ONLY : npwx
   USE westcom,               ONLY : n_lanczos
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER, INTENT(IN) :: nlocx
   INTEGER, INTENT(IN) :: nloc
   !
   IF(gamma_only) THEN
      ALLOCATE(tmp_r3(nlocx,nloc,n_lanczos))
   ELSE
      ALLOCATE(tmp_c3(nlocx,nloc,n_lanczos))
   ENDIF
   ALLOCATE(dvpsi_h(npwx*npol,nlocx))
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_gw_gpu()
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   IF(ALLOCATED(tmp_r3)) THEN
      DEALLOCATE(tmp_r3)
   ENDIF
   IF(ALLOCATED(tmp_c3)) THEN
      DEALLOCATE(tmp_c3)
   ENDIF
   IF(ALLOCATED(dvpsi_h)) THEN
      DEALLOCATE(dvpsi_h)
   ENDIF
   IF(ALLOCATED(ps_r)) THEN
      DEALLOCATE(ps_r)
   ENDIF
   IF(ALLOCATED(ps_c)) THEN
      DEALLOCATE(ps_c)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_lanczos_gpu(nloc)
   !-----------------------------------------------------------------------
   !
   USE control_flags,         ONLY : gamma_only
   USE noncollin_module,      ONLY : npol
   USE pwcom,                 ONLY : npwx
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER, INTENT(IN) :: nloc
   !
   ALLOCATE(r(npwx*npol,nloc))
   ALLOCATE(tmp_r(nloc))
   !$acc enter data create(tmp_r)
   IF(.NOT. gamma_only) THEN
      ALLOCATE(tmp_c(nloc))
      !$acc enter data create(tmp_c)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_lanczos_gpu()
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   IF(ALLOCATED(r)) THEN
      DEALLOCATE(r)
   ENDIF
   IF(ALLOCATED(tmp_r)) THEN
      !$acc exit data delete(tmp_r)
      DEALLOCATE(tmp_r)
   ENDIF
   IF(ALLOCATED(tmp_c)) THEN
      !$acc exit data delete(tmp_c)
      DEALLOCATE(tmp_c)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE reallocate_ps_gpu(nbndval,m)
   !-----------------------------------------------------------------------
   !
   USE control_flags,         ONLY : gamma_only
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER, INTENT(IN) :: nbndval
   INTEGER, INTENT(IN) :: m
   !
   IF(gamma_only) THEN
      IF(ALLOCATED(ps_r)) THEN
         IF(SIZE(ps_r,DIM=1) /= nbndval .OR. SIZE(ps_r,DIM=2) /= m) THEN
            DEALLOCATE(ps_r)
         ENDIF
      ENDIF
      IF(.NOT. ALLOCATED(ps_r)) THEN
         ALLOCATE(ps_r(nbndval,m))
      ENDIF
   ELSE
      IF(ALLOCATED(ps_c)) THEN
         IF(SIZE(ps_c,DIM=1) /= nbndval .OR. SIZE(ps_c,DIM=2) /= m) THEN
            DEALLOCATE(ps_c)
         ENDIF
      ENDIF
      IF(.NOT. ALLOCATED(ps_c)) THEN
         ALLOCATE(ps_c(nbndval,m))
      ENDIF
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_macropol_gpu(m)
   !-----------------------------------------------------------------------
   !
   USE control_flags,         ONLY : gamma_only
   USE ions_base,             ONLY : nat
   USE lsda_mod,              ONLY : nspin
   USE noncollin_module,      ONLY : noncolin,npol
   USE uspp,                  ONLY : nkb
   USE uspp_param,            ONLY : nhm
   USE wvfct,                 ONLY : npwx
   USE becmod_subs_gpum,      ONLY : allocate_bec_type_gpu
   USE westcom,               ONLY : l_skip_nl_part_of_hcomr
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER, INTENT(IN) :: m
   !
   CALL allocate_linsolve_gpu(3)
   !
   ALLOCATE(gk(3,npwx))
   IF(.NOT. l_skip_nl_part_of_hcomr .AND. nkb > 0) THEN
      ALLOCATE(dvkb(npwx,nkb))
      ALLOCATE(work(npwx,nkb))
      IF(noncolin) THEN
         ALLOCATE(deff_nc(nhm,nhm,nat,nspin))
         ALLOCATE(psc(nkb,npol,m,2))
      ELSE
         ALLOCATE(deff(nhm,nhm,nat))
         ALLOCATE(ps2(nkb,m,2))
      ENDIF
      !
      CALL allocate_bec_type_gpu(nkb,m,becp1_d)
      CALL allocate_bec_type_gpu(nkb,m,becp2_d)
      !
      IF(noncolin) THEN
         becp1_d_nc_d => becp1_d%nc_d
         becp2_d_nc_d => becp2_d%nc_d
      ELSEIF(gamma_only) THEN
         becp1_d_r_d => becp1_d%r_d
         becp2_d_r_d => becp2_d%r_d
      ELSE
         becp1_d_k_d => becp1_d%k_d
         becp2_d_k_d => becp2_d%k_d
      ENDIF
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_macropol_gpu()
   !-----------------------------------------------------------------------
   !
   USE becmod_subs_gpum,      ONLY : deallocate_bec_type_gpu
   USE uspp,                  ONLY : nkb
   USE westcom,               ONLY : l_skip_nl_part_of_hcomr

   IMPLICIT NONE
   !
   CALL deallocate_linsolve_gpu()
   !
   IF(ALLOCATED(gk)) THEN
      DEALLOCATE(gk)
   ENDIF
   IF(ALLOCATED(dvkb)) THEN
      DEALLOCATE(dvkb)
   ENDIF
   IF(ALLOCATED(work)) THEN
      DEALLOCATE(work)
   ENDIF
   IF(ALLOCATED(deff_nc)) THEN
      DEALLOCATE(deff_nc)
   ENDIF
   IF(ALLOCATED(psc)) THEN
      DEALLOCATE(psc)
   ENDIF
   IF(ALLOCATED(deff)) THEN
      DEALLOCATE(deff)
   ENDIF
   IF(ALLOCATED(ps2)) THEN
      DEALLOCATE(ps2)
   ENDIF
   IF(ASSOCIATED(becp1_d_r_d)) THEN
      NULLIFY(becp1_d_r_d)
   ENDIF
   IF(ASSOCIATED(becp2_d_r_d)) THEN
      NULLIFY(becp2_d_r_d)
   ENDIF
   IF(ASSOCIATED(becp1_d_k_d)) THEN
      NULLIFY(becp1_d_k_d)
   ENDIF
   IF(ASSOCIATED(becp2_d_k_d)) THEN
      NULLIFY(becp2_d_k_d)
   ENDIF
   IF(ASSOCIATED(becp1_d_nc_d)) THEN
      NULLIFY(becp1_d_nc_d)
   ENDIF
   IF(ASSOCIATED(becp2_d_nc_d)) THEN
      NULLIFY(becp2_d_nc_d)
   ENDIF
   !
   IF(.NOT. l_skip_nl_part_of_hcomr .AND. nkb > 0) THEN
      CALL deallocate_bec_type_gpu(becp1_d)
      CALL deallocate_bec_type_gpu(becp2_d)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_chi_gpu(l_real)
   !-----------------------------------------------------------------------
   !
   USE westcom,               ONLY : n_pdep_eigen_to_use,l_macropol
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   LOGICAL, INTENT(IN) :: l_real
   !
   ! Workspace
   !
   INTEGER :: istat
   INTEGER :: lwork
   INTEGER(i8b) :: n8
   !
   ALLOCATE(piv(n_pdep_eigen_to_use))
   IF(l_real) THEN
      ALLOCATE(x_r(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      IF(l_macropol) THEN
         ALLOCATE(wh_r(n_pdep_eigen_to_use,3))
         ALLOCATE(wl_r(3,n_pdep_eigen_to_use))
         ALLOCATE(temph_r(n_pdep_eigen_to_use,3))
         ALLOCATE(templ_r(3,n_pdep_eigen_to_use))
         ALLOCATE(tempt_r(3,3))
      ENDIF
      !
#if CUDA_VERSION > 11010 && (__PGIC__ > 21 || (__PGIC__ == 21 && __PGIC_MINOR > 2))
      n8 = INT(n_pdep_eigen_to_use,KIND=i8b)
      !
      !$acc host_data use_device(x_r)
      istat = cusolverDnXgetrf_bufferSize(cusolv_h,cusolv_p,n8,n8,cudaDataType(CUDA_R_64F),x_r,n8,&
            & cudaDataType(CUDA_R_64F),l_inv,l_inv_h)
      !$acc end host_data
      !
      l_inv = MAX(l_inv,(n8**2)*8)
      !
      ALLOCATE(work_r(l_inv/8))
      ALLOCATE(work_r_h(l_inv_h/8))
#else
      !$acc host_data use_device(x_r)
      istat = cusolverDnDgetrf_bufferSize(cusolv_h,n_pdep_eigen_to_use,n_pdep_eigen_to_use,x_r,&
            & n_pdep_eigen_to_use,lwork)
      !$acc end host_data
      !
      lwork = MAX(lwork,n_pdep_eigen_to_use**2)
      !
      ALLOCATE(work_r(lwork))
#endif
   ELSE
      ALLOCATE(x_c(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      IF(l_macropol) THEN
         ALLOCATE(wh_c(n_pdep_eigen_to_use,3))
         ALLOCATE(wl_c(3,n_pdep_eigen_to_use))
         ALLOCATE(temph_c(n_pdep_eigen_to_use,3))
         ALLOCATE(templ_c(3,n_pdep_eigen_to_use))
         ALLOCATE(tempt_c(3,3))
      ENDIF
      !
#if CUDA_VERSION > 11010 && (__PGIC__ > 21 || (__PGIC__ == 21 && __PGIC_MINOR > 2))
      n8 = INT(n_pdep_eigen_to_use,KIND=i8b)
      !
      !$acc host_data use_device(x_c)
      istat = cusolverDnXgetrf_bufferSize(cusolv_h,cusolv_p,n8,n8,cudaDataType(CUDA_C_64F),x_c,n8,&
            & cudaDataType(CUDA_C_64F),l_inv,l_inv_h)
      !$acc end host_data
      !
      l_inv = MAX(l_inv,(n8**2)*16)
      !
      ALLOCATE(work_c(l_inv/16))
      ALLOCATE(work_c_h(l_inv_h/16))
#else
      !$acc host_data use_device(x_c)
      istat = cusolverDnZgetrf_bufferSize(cusolv_h,n_pdep_eigen_to_use,n_pdep_eigen_to_use,x_c,&
            & n_pdep_eigen_to_use,lwork)
      !$acc end host_data
      !
      lwork = MAX(lwork,n_pdep_eigen_to_use**2)
      !
      ALLOCATE(work_c(lwork))
#endif
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_chi_gpu()
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   IF(ALLOCATED(piv)) THEN
      DEALLOCATE(piv)
   ENDIF
   IF(ALLOCATED(x_r)) THEN
      DEALLOCATE(x_r)
   ENDIF
   IF(ALLOCATED(wh_r)) THEN
      DEALLOCATE(wh_r)
   ENDIF
   IF(ALLOCATED(wl_r)) THEN
      DEALLOCATE(wl_r)
   ENDIF
   IF(ALLOCATED(temph_r)) THEN
      DEALLOCATE(temph_r)
   ENDIF
   IF(ALLOCATED(templ_r)) THEN
      DEALLOCATE(templ_r)
   ENDIF
   IF(ALLOCATED(tempt_r)) THEN
      DEALLOCATE(tempt_r)
   ENDIF
   IF(ALLOCATED(work_r_h)) THEN
      DEALLOCATE(work_r_h)
   ENDIF
   IF(ALLOCATED(work_r)) THEN
      DEALLOCATE(work_r)
   ENDIF
   IF(ALLOCATED(x_c)) THEN
      DEALLOCATE(x_c)
   ENDIF
   IF(ALLOCATED(wh_c)) THEN
      DEALLOCATE(wh_c)
   ENDIF
   IF(ALLOCATED(wl_c)) THEN
      DEALLOCATE(wl_c)
   ENDIF
   IF(ALLOCATED(temph_c)) THEN
      DEALLOCATE(temph_c)
   ENDIF
   IF(ALLOCATED(templ_c)) THEN
      DEALLOCATE(templ_c)
   ENDIF
   IF(ALLOCATED(tempt_c)) THEN
      DEALLOCATE(tempt_c)
   ENDIF
   IF(ALLOCATED(work_c_h)) THEN
      DEALLOCATE(work_c_h)
   ENDIF
   IF(ALLOCATED(work_c)) THEN
      DEALLOCATE(work_c)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_bse_gpu(nbndlocx)
   !-----------------------------------------------------------------------
   !
   USE control_flags,         ONLY : gamma_only
   USE lsda_mod,              ONLY : nspin
   USE wvfct,                 ONLY : npwx
   USE noncollin_module,      ONLY : npol
   USE fft_base,              ONLY : dffts
   USE westcom,               ONLY : nbndval0x,n_trunc_bands,l_bse,l_hybrid_tddft,l_local_repr,&
                                   & et_qp,u_matrix
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER, INTENT(IN) :: nbndlocx
   !
   ALLOCATE(factors(nbndlocx))
   IF(l_bse .OR. l_hybrid_tddft) THEN
      ALLOCATE(raux1(dffts%nnr))
      ALLOCATE(raux2(dffts%nnr))
      ALLOCATE(caux1(npwx,nbndval0x-n_trunc_bands))
      ALLOCATE(caux2(npwx,nbndlocx))
      IF(l_local_repr) THEN
         ALLOCATE(caux3(npwx,nbndval0x-n_trunc_bands))
      ENDIF
      ALLOCATE(gaux(npwx))
      !$acc enter data create(caux1,caux2,gaux)
   ENDIF
   ALLOCATE(hevc1(npwx*npol,nbndlocx))
   ALLOCATE(dvrs(dffts%nnr,nspin))
   !$acc enter data create(dvrs)
   IF(gamma_only) THEN
      ALLOCATE(tmp_r(dffts%nnr))
      !$acc enter data create(tmp_r)
   ENDIF
   ALLOCATE(tmp_c(dffts%nnr))
   !$acc enter data create(tmp_c)
   IF(.NOT. gamma_only) THEN
      ALLOCATE(psic2(dffts%nnr))
   ENDIF
   ALLOCATE(dvaux(dffts%nnr,nspin))
   !$acc enter data create(dvaux)
   ALLOCATE(dvhart(dffts%nnr))
   !
   !$acc enter data copyin(et_qp,u_matrix)
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_bse_gpu()
   !-----------------------------------------------------------------------
   !
   USE westcom,               ONLY : et_qp,u_matrix
   !
   IMPLICIT NONE
   !
   IF(ALLOCATED(factors)) THEN
      DEALLOCATE(factors)
   ENDIF
   IF(ALLOCATED(raux1)) THEN
      DEALLOCATE(raux1)
   ENDIF
   IF(ALLOCATED(raux2)) THEN
      DEALLOCATE(raux2)
   ENDIF
   IF(ALLOCATED(caux1)) THEN
      !$acc exit data delete(caux1)
      DEALLOCATE(caux1)
   ENDIF
   IF(ALLOCATED(caux2)) THEN
      !$acc exit data delete(caux2)
      DEALLOCATE(caux2)
   ENDIF
   IF(ALLOCATED(caux3)) THEN
      DEALLOCATE(caux3)
   ENDIF
   IF(ALLOCATED(gaux)) THEN
      !$acc exit data delete(gaux)
      DEALLOCATE(gaux)
   ENDIF
   IF(ALLOCATED(hevc1)) THEN
      DEALLOCATE(hevc1)
   ENDIF
   IF(ALLOCATED(dvrs)) THEN
      !$acc exit data delete(dvrs)
      DEALLOCATE(dvrs)
   ENDIF
   IF(ALLOCATED(tmp_r)) THEN
      !$acc exit data delete(tmp_r)
      DEALLOCATE(tmp_r)
   ENDIF
   IF(ALLOCATED(tmp_c)) THEN
      !$acc exit data delete(tmp_c)
      DEALLOCATE(tmp_c)
   ENDIF
   IF(ALLOCATED(psic2)) THEN
      DEALLOCATE(psic2)
   ENDIF
   IF(ALLOCATED(dvaux)) THEN
      !$acc exit data delete(dvaux)
      DEALLOCATE(dvaux)
   ENDIF
   IF(ALLOCATED(dvhart)) THEN
      DEALLOCATE(dvhart)
   ENDIF
   IF(ALLOCATED(ps_r)) THEN
      DEALLOCATE(ps_r)
   ENDIF
   IF(ALLOCATED(ps_c)) THEN
      DEALLOCATE(ps_c)
   ENDIF
   !
   !$acc exit data delete(et_qp,u_matrix)
   !
   END SUBROUTINE
#endif
END MODULE
