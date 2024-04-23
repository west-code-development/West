!
! Copyright (C) 2015-2024 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
!-----------------------------------------------------------------------
MODULE west_gpu
   !-----------------------------------------------------------------------
   !
   USE kinds,         ONLY : DP,i8b
#if defined(__CUDA)
   USE becmod,        ONLY : bec_type
   USE cublas
   USE cusolverdn
   !
   IMPLICIT NONE
   !
   ! Linsolve
   !
   LOGICAL, ALLOCATABLE :: is_conv(:)
   INTEGER, ALLOCATABLE :: l2i_map(:)
   REAL(DP), ALLOCATABLE :: eu(:)
   REAL(DP), ALLOCATABLE :: a(:)
   REAL(DP), ALLOCATABLE :: c(:)
   REAL(DP), ALLOCATABLE :: rho(:)
   REAL(DP), ALLOCATABLE :: rhoold(:)
   COMPLEX(DP), ALLOCATABLE :: g(:,:)
   COMPLEX(DP), ALLOCATABLE :: t(:,:)
   COMPLEX(DP), ALLOCATABLE :: h(:,:)
   !
   ! GW
   !
   COMPLEX(DP), ALLOCATABLE :: dvpsi2(:,:)
   ATTRIBUTES(PINNED) :: dvpsi2
   !
   ! Macropol
   !
   REAL(DP), ALLOCATABLE :: ep_pol(:)
   REAL(DP), ALLOCATABLE :: e_pol(:)
   COMPLEX(DP), ALLOCATABLE :: phi(:,:)
   REAL(DP), ALLOCATABLE :: gk(:,:)
   COMPLEX(DP), ALLOCATABLE :: ps2(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: psc(:,:,:,:)
   COMPLEX(DP), ALLOCATABLE :: dvkb(:,:)
   COMPLEX(DP), ALLOCATABLE :: work(:,:)
   TYPE(bec_type) :: becp1 ! (nkb,m)
   TYPE(bec_type) :: becp2 ! (nkb,m)
   !
   ! Chi invert
   !
   INTEGER(i8b) :: l_inv
   INTEGER(i8b) :: l_inv_h
   INTEGER(i8b), ALLOCATABLE :: piv(:)
   REAL(DP), ALLOCATABLE :: work_r_h(:)
   REAL(DP), ALLOCATABLE :: work_r(:)
   REAL(DP), ALLOCATABLE :: x_r(:,:)
   REAL(DP), ALLOCATABLE :: wh_r(:,:)
   REAL(DP), ALLOCATABLE :: wl_r(:,:)
   REAL(DP), ALLOCATABLE :: temph_r(:,:)
   REAL(DP), ALLOCATABLE :: templ_r(:,:)
   REAL(DP), ALLOCATABLE :: tempt_r(:,:)
   COMPLEX(DP), ALLOCATABLE :: work_c_h(:)
   COMPLEX(DP), ALLOCATABLE :: work_c(:)
   COMPLEX(DP), ALLOCATABLE :: x_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: wh_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: wl_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: temph_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: templ_c(:,:)
   COMPLEX(DP), ALLOCATABLE :: tempt_c(:,:)
   !
   ! Lanczos
   !
   REAL(DP), ALLOCATABLE :: alpha(:)
   ATTRIBUTES(PINNED) :: alpha
   REAL(DP), ALLOCATABLE :: beta(:)
   ATTRIBUTES(PINNED) :: beta
   COMPLEX(DP), ALLOCATABLE :: r(:,:)
   !
   ! BSE
   !
   REAL(DP), ALLOCATABLE :: factors(:)
   REAL(DP), ALLOCATABLE :: raux1(:)
   REAL(DP), ALLOCATABLE :: raux2(:)
   COMPLEX(DP), ALLOCATABLE :: caux1(:,:)
   COMPLEX(DP), ALLOCATABLE :: caux2(:,:)
   COMPLEX(DP), ALLOCATABLE :: caux3(:,:)
   COMPLEX(DP), ALLOCATABLE :: caux4(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
   COMPLEX(DP), ALLOCATABLE :: hevc1(:,:)
   COMPLEX(DP), ALLOCATABLE :: psic2(:)
   COMPLEX(DP), ALLOCATABLE :: dvaux(:,:)
   COMPLEX(DP), ALLOCATABLE :: gdrho(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: dvhart(:)
   COMPLEX(DP), ALLOCATABLE :: gaux(:)
   ATTRIBUTES(PINNED) :: gaux
   !
   ! Workspace
   !
   REAL(DP), ALLOCATABLE :: tmp_r(:)
   REAL(DP), ALLOCATABLE :: tmp_r3(:,:,:)
   REAL(DP), ALLOCATABLE :: ps_r(:,:)
   COMPLEX(DP), ALLOCATABLE :: tmp_c(:)
   COMPLEX(DP), ALLOCATABLE :: tmp_c3(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
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
   USE fft_base,              ONLY : dffts,dfftp
   USE pwcom,                 ONLY : wg,ngk,nks
   USE cell_base,             ONLY : bg
   USE westcom,               ONLY : igq_q
   USE wavefunctions,         ONLY : evc,psic,psic_nc
   USE wvfct,                 ONLY : et
   !
   IMPLICIT NONE
   !
   !$acc enter data create(evc,psic) copyin(wg,ngk,bg,igq_q,et)
   IF(nks == 1) THEN
      !$acc update device(evc)
   ENDIF
   IF(ALLOCATED(psic_nc)) THEN
      !$acc enter data create(psic_nc)
   ENDIF
   !$acc enter data copyin(dffts)
   IF(ALLOCATED(dffts%nl)) THEN
      !$acc enter data copyin(dffts%nl)
   ENDIF
   IF(ALLOCATED(dffts%nlm)) THEN
      !$acc enter data copyin(dffts%nlm)
   ENDIF
   !$acc enter data copyin(dfftp)
   IF(ALLOCATED(dfftp%nl)) THEN
      !$acc enter data copyin(dfftp%nl)
   ENDIF
   IF(ALLOCATED(dfftp%nlm)) THEN
      !$acc enter data copyin(dfftp%nlm)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_gpu()
   !-----------------------------------------------------------------------
   !
   USE fft_base,              ONLY : dffts,dfftp
   USE pwcom,                 ONLY : wg,ngk
   USE cell_base,             ONLY : bg
   USE westcom,               ONLY : igq_q
   USE wavefunctions,         ONLY : evc,psic,psic_nc
   USE wvfct,                 ONLY : et
   !
   IMPLICIT NONE
   !
   !$acc exit data delete(wg,ngk,bg,igq_q,et,evc,psic)
   IF(ALLOCATED(psic_nc)) THEN
      !$acc exit data delete(psic_nc)
   ENDIF
   IF(ALLOCATED(dffts%nl)) THEN
      !$acc exit data delete(dffts%nl)
   ENDIF
   IF(ALLOCATED(dffts%nlm)) THEN
      !$acc exit data delete(dffts%nlm)
   ENDIF
   !$acc exit data delete(dffts)
   IF(ALLOCATED(dfftp%nl)) THEN
      !$acc exit data delete(dfftp%nl)
   ENDIF
   IF(ALLOCATED(dfftp%nlm)) THEN
      !$acc exit data delete(dfftp%nlm)
   ENDIF
   !$acc exit data delete(dfftp)
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
   !$acc enter data create(is_conv)
   ALLOCATE(l2i_map(nbndloc))
   !$acc enter data create(l2i_map)
   ALLOCATE(eu(nbndloc))
   !$acc enter data create(eu)
   ALLOCATE(a(nbndloc))
   !$acc enter data create(a)
   ALLOCATE(c(nbndloc))
   !$acc enter data create(c)
   ALLOCATE(rho(nbndloc))
   !$acc enter data create(rho)
   ALLOCATE(rhoold(nbndloc))
   !$acc enter data create(rhoold)
   ALLOCATE(g(npwx*npol,nbndloc))
   !$acc enter data create(g)
   ALLOCATE(t(npwx*npol,nbndloc))
   !$acc enter data create(t)
   ALLOCATE(h(npwx*npol,nbndloc))
   !$acc enter data create(h)
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
      !$acc exit data delete(is_conv)
      DEALLOCATE(is_conv)
   ENDIF
   IF(ALLOCATED(l2i_map)) THEN
      !$acc exit data delete(l2i_map)
      DEALLOCATE(l2i_map)
   ENDIF
   IF(ALLOCATED(eu)) THEN
      !$acc exit data delete(eu)
      DEALLOCATE(eu)
   ENDIF
   IF(ALLOCATED(a)) THEN
      !$acc exit data delete(a)
      DEALLOCATE(a)
   ENDIF
   IF(ALLOCATED(c)) THEN
      !$acc exit data delete(c)
      DEALLOCATE(c)
   ENDIF
   IF(ALLOCATED(rho)) THEN
      !$acc exit data delete(rho)
      DEALLOCATE(rho)
   ENDIF
   IF(ALLOCATED(rhoold)) THEN
      !$acc exit data delete(rhoold)
      DEALLOCATE(rhoold)
   ENDIF
   IF(ALLOCATED(g)) THEN
      !$acc exit data delete(g)
      DEALLOCATE(g)
   ENDIF
   IF(ALLOCATED(t)) THEN
      !$acc exit data delete(t)
      DEALLOCATE(t)
   ENDIF
   IF(ALLOCATED(h)) THEN
      !$acc exit data delete(h)
      DEALLOCATE(h)
   ENDIF
   IF(ALLOCATED(ps_r)) THEN
      !$acc exit data delete(ps_r)
      DEALLOCATE(ps_r)
   ENDIF
   IF(ALLOCATED(ps_c)) THEN
      !$acc exit data delete(ps_c)
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
      !$acc enter data create(tmp_r3)
   ELSE
      ALLOCATE(tmp_c3(nlocx,nloc,n_lanczos))
      !$acc enter data create(tmp_c3)
   ENDIF
   ALLOCATE(dvpsi2(npwx*npol,nlocx))
   !$acc enter data create(dvpsi2)
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
      !$acc exit data delete(tmp_r3)
      DEALLOCATE(tmp_r3)
   ENDIF
   IF(ALLOCATED(tmp_c3)) THEN
      !$acc exit data delete(tmp_c3)
      DEALLOCATE(tmp_c3)
   ENDIF
   IF(ALLOCATED(dvpsi2)) THEN
      !$acc exit data delete(dvpsi2)
      DEALLOCATE(dvpsi2)
   ENDIF
   IF(ALLOCATED(ps_r)) THEN
      !$acc exit data delete(ps_r)
      DEALLOCATE(ps_r)
   ENDIF
   IF(ALLOCATED(ps_c)) THEN
      !$acc exit data delete(ps_c)
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
   ALLOCATE(alpha(nloc))
   !$acc enter data create(alpha)
   ALLOCATE(beta(nloc))
   !$acc enter data create(beta)
   ALLOCATE(r(npwx*npol,nloc))
   !$acc enter data create(r)
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
   IF(ALLOCATED(alpha)) THEN
      !$acc exit data delete(alpha)
      DEALLOCATE(alpha)
   ENDIF
   IF(ALLOCATED(beta)) THEN
      !$acc exit data delete(beta)
      DEALLOCATE(beta)
   ENDIF
   IF(ALLOCATED(r)) THEN
      !$acc exit data delete(r)
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
         IF(SIZE(ps_r,1) /= nbndval .OR. SIZE(ps_r,2) /= m) THEN
            !$acc exit data delete(ps_r)
            DEALLOCATE(ps_r)
         ENDIF
      ENDIF
      IF(.NOT. ALLOCATED(ps_r)) THEN
         ALLOCATE(ps_r(nbndval,m))
         !$acc enter data create(ps_r)
      ENDIF
   ELSE
      IF(ALLOCATED(ps_c)) THEN
         IF(SIZE(ps_c,1) /= nbndval .OR. SIZE(ps_c,2) /= m) THEN
            !$acc exit data delete(ps_c)
            DEALLOCATE(ps_c)
         ENDIF
      ENDIF
      IF(.NOT. ALLOCATED(ps_c)) THEN
         ALLOCATE(ps_c(nbndval,m))
         !$acc enter data create(ps_c)
      ENDIF
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_macropol_gpu(m)
   !-----------------------------------------------------------------------
   !
   USE ions_base,             ONLY : nat
   USE lsda_mod,              ONLY : nspin
   USE noncollin_module,      ONLY : noncolin,npol
   USE uspp,                  ONLY : nkb
   USE uspp_param,            ONLY : nhm
   USE wvfct,                 ONLY : npwx
   USE becmod,                ONLY : allocate_bec_type_acc
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
   !$acc enter data create(gk)
   ALLOCATE(ep_pol(3))
   !$acc enter data create(ep_pol)
   ALLOCATE(e_pol(3))
   !$acc enter data create(e_pol)
   ALLOCATE(phi(npwx*npol,3))
   !$acc enter data create(phi)
   IF(.NOT. l_skip_nl_part_of_hcomr .AND. nkb > 0) THEN
      ALLOCATE(dvkb(npwx,nkb))
      !$acc enter data create(dvkb)
      ALLOCATE(work(npwx,nkb))
      !$acc enter data create(work)
      IF(noncolin) THEN
         ALLOCATE(psc(nkb,npol,m,2))
         !$acc enter data create(psc)
      ELSE
         ALLOCATE(ps2(nkb,m,2))
         !$acc enter data create(ps2)
      ENDIF
      !
      CALL allocate_bec_type_acc(nkb,m,becp1)
      CALL allocate_bec_type_acc(nkb,m,becp2)
   ENDIF
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_macropol_gpu()
   !-----------------------------------------------------------------------
   !
   USE becmod,                ONLY : deallocate_bec_type_acc
   USE uspp,                  ONLY : nkb
   USE westcom,               ONLY : l_skip_nl_part_of_hcomr

   IMPLICIT NONE
   !
   CALL deallocate_linsolve_gpu()
   !
   IF(ALLOCATED(gk)) THEN
      !$acc exit data delete(gk)
      DEALLOCATE(gk)
   ENDIF
   IF(ALLOCATED(ep_pol)) THEN
      !$acc exit data delete(ep_pol)
      DEALLOCATE(ep_pol)
   ENDIF
   IF(ALLOCATED(e_pol)) THEN
      !$acc exit data delete(e_pol)
      DEALLOCATE(e_pol)
   ENDIF
   IF(ALLOCATED(phi)) THEN
      !$acc exit data delete(phi)
      DEALLOCATE(phi)
   ENDIF
   !
   IF(ALLOCATED(dvkb)) THEN
      !$acc exit data delete(dvkb)
      DEALLOCATE(dvkb)
   ENDIF
   IF(ALLOCATED(work)) THEN
      !$acc exit data delete(work)
      DEALLOCATE(work)
   ENDIF
   IF(ALLOCATED(psc)) THEN
      !$acc exit data delete(psc)
      DEALLOCATE(psc)
   ENDIF
   IF(ALLOCATED(ps2)) THEN
      !$acc exit data delete(ps2)
      DEALLOCATE(ps2)
   ENDIF
   !
   IF(.NOT. l_skip_nl_part_of_hcomr .AND. nkb > 0) THEN
      CALL deallocate_bec_type_acc(becp1)
      CALL deallocate_bec_type_acc(becp2)
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
   INTEGER(i8b) :: n8
   !
   ALLOCATE(piv(n_pdep_eigen_to_use))
   !$acc enter data create(piv)
   IF(l_real) THEN
      ALLOCATE(x_r(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      !$acc enter data create(x_r)
      IF(l_macropol) THEN
         ALLOCATE(wh_r(n_pdep_eigen_to_use,3))
         !$acc enter data create(wh_r)
         ALLOCATE(wl_r(3,n_pdep_eigen_to_use))
         !$acc enter data create(wl_r)
         ALLOCATE(temph_r(n_pdep_eigen_to_use,3))
         !$acc enter data create(temph_r)
         ALLOCATE(templ_r(3,n_pdep_eigen_to_use))
         !$acc enter data create(templ_r)
         ALLOCATE(tempt_r(3,3))
         !$acc enter data create(tempt_r)
      ENDIF
      !
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
      !$acc enter data create(work_r)
      ALLOCATE(work_r_h(l_inv_h/8))
   ELSE
      ALLOCATE(x_c(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
      !$acc enter data create(x_c)
      IF(l_macropol) THEN
         ALLOCATE(wh_c(n_pdep_eigen_to_use,3))
         !$acc enter data create(wh_c)
         ALLOCATE(wl_c(3,n_pdep_eigen_to_use))
         !$acc enter data create(wl_c)
         ALLOCATE(temph_c(n_pdep_eigen_to_use,3))
         !$acc enter data create(temph_c)
         ALLOCATE(templ_c(3,n_pdep_eigen_to_use))
         !$acc enter data create(templ_c)
         ALLOCATE(tempt_c(3,3))
         !$acc enter data create(tempt_c)
      ENDIF
      !
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
      !$acc enter data create(work_c)
      ALLOCATE(work_c_h(l_inv_h/16))
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
      !$acc exit data delete(piv)
      DEALLOCATE(piv)
   ENDIF
   IF(ALLOCATED(x_r)) THEN
      !$acc exit data delete(x_r)
      DEALLOCATE(x_r)
   ENDIF
   IF(ALLOCATED(wh_r)) THEN
      !$acc exit data delete(wh_r)
      DEALLOCATE(wh_r)
   ENDIF
   IF(ALLOCATED(wl_r)) THEN
      !$acc exit data delete(wl_r)
      DEALLOCATE(wl_r)
   ENDIF
   IF(ALLOCATED(temph_r)) THEN
      !$acc exit data delete(temph_r)
      DEALLOCATE(temph_r)
   ENDIF
   IF(ALLOCATED(templ_r)) THEN
      !$acc exit data delete(templ_r)
      DEALLOCATE(templ_r)
   ENDIF
   IF(ALLOCATED(tempt_r)) THEN
      !$acc exit data delete(tempt_r)
      DEALLOCATE(tempt_r)
   ENDIF
   IF(ALLOCATED(work_r_h)) THEN
      DEALLOCATE(work_r_h)
   ENDIF
   IF(ALLOCATED(work_r)) THEN
      !$acc exit data delete(work_r)
      DEALLOCATE(work_r)
   ENDIF
   IF(ALLOCATED(x_c)) THEN
      !$acc exit data delete(x_c)
      DEALLOCATE(x_c)
   ENDIF
   IF(ALLOCATED(wh_c)) THEN
      !$acc exit data delete(wh_c)
      DEALLOCATE(wh_c)
   ENDIF
   IF(ALLOCATED(wl_c)) THEN
      !$acc exit data delete(wl_c)
      DEALLOCATE(wl_c)
   ENDIF
   IF(ALLOCATED(temph_c)) THEN
      !$acc exit data delete(temph_c)
      DEALLOCATE(temph_c)
   ENDIF
   IF(ALLOCATED(templ_c)) THEN
      !$acc exit data delete(templ_c)
      DEALLOCATE(templ_c)
   ENDIF
   IF(ALLOCATED(tempt_c)) THEN
      !$acc exit data delete(tempt_c)
      DEALLOCATE(tempt_c)
   ENDIF
   IF(ALLOCATED(work_c_h)) THEN
      DEALLOCATE(work_c_h)
   ENDIF
   IF(ALLOCATED(work_c)) THEN
      !$acc exit data delete(work_c)
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
   USE noncollin_module,      ONLY : npol,nspin_gga
   USE fft_base,              ONLY : dffts,dfftp
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
   !$acc enter data create(factors)
   IF(l_bse .OR. l_hybrid_tddft) THEN
      ALLOCATE(raux1(dffts%nnr))
      !$acc enter data create(raux1)
      ALLOCATE(raux2(dffts%nnr))
      !$acc enter data create(raux2)
      ALLOCATE(caux1(npwx,nbndval0x-n_trunc_bands))
      !$acc enter data create(caux1)
      ALLOCATE(caux2(npwx,nbndlocx))
      !$acc enter data create(caux2)
      IF(l_local_repr) THEN
         ALLOCATE(caux3(npwx,nbndval0x-n_trunc_bands))
         !$acc enter data create(caux3)
      ENDIF
      ALLOCATE(gaux(npwx))
      !$acc enter data create(gaux)
   ENDIF
   IF(.NOT. l_bse) THEN
      ALLOCATE(caux4(3,dfftp%nnr,nspin_gga))
      !$acc enter data create(caux4)
      ALLOCATE(gdrho(3,dfftp%nnr,nspin_gga))
      !$acc enter data create(gdrho)
   ENDIF
   ALLOCATE(hevc1(npwx*npol,nbndlocx))
   !$acc enter data create(hevc1)
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
      !$acc enter data create(psic2)
   ENDIF
   ALLOCATE(dvaux(dfftp%nnr,nspin))
   !$acc enter data create(dvaux)
   ALLOCATE(dvhart(dfftp%nnr))
   !$acc enter data create(dvhart)
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
      !$acc exit data delete(factors)
      DEALLOCATE(factors)
   ENDIF
   IF(ALLOCATED(raux1)) THEN
      !$acc exit data delete(raux1)
      DEALLOCATE(raux1)
   ENDIF
   IF(ALLOCATED(raux2)) THEN
      !$acc exit data delete(raux2)
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
      !$acc exit data delete(caux3)
      DEALLOCATE(caux3)
   ENDIF
   IF(ALLOCATED(caux4)) THEN
      !$acc exit data delete(caux4)
      DEALLOCATE(caux4)
   ENDIF
   IF(ALLOCATED(gaux)) THEN
      !$acc exit data delete(gaux)
      DEALLOCATE(gaux)
   ENDIF
   IF(ALLOCATED(gdrho)) THEN
      !$acc exit data delete(gdrho)
      DEALLOCATE(gdrho)
   ENDIF
   IF(ALLOCATED(hevc1)) THEN
      !$acc exit data delete(hevc1)
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
      !$acc exit data delete(psic2)
      DEALLOCATE(psic2)
   ENDIF
   IF(ALLOCATED(dvaux)) THEN
      !$acc exit data delete(dvaux)
      DEALLOCATE(dvaux)
   ENDIF
   IF(ALLOCATED(dvhart)) THEN
      !$acc exit data delete(dvhart)
      DEALLOCATE(dvhart)
   ENDIF
   IF(ALLOCATED(ps_r)) THEN
      !$acc exit data delete(ps_r)
      DEALLOCATE(ps_r)
   ENDIF
   IF(ALLOCATED(ps_c)) THEN
      !$acc exit data delete(ps_c)
      DEALLOCATE(ps_c)
   ENDIF
   !
   !$acc exit data delete(et_qp,u_matrix)
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE allocate_forces_gpu()
   !-----------------------------------------------------------------------
   !
   USE fft_base,              ONLY : dffts
   !
   IMPLICIT NONE
   !
   ALLOCATE(tmp_r(dffts%nnr))
   !$acc enter data create(tmp_r)
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE deallocate_forces_gpu()
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   IF(ALLOCATED(tmp_r)) THEN
      !$acc exit data delete(tmp_r)
      DEALLOCATE(tmp_r)
   ENDIF
   !
   END SUBROUTINE
#endif
   !
END MODULE
