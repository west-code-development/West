!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE wbse_dv
  !
  IMPLICIT NONE
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE wbse_dv_setup(l_skip)
    !-----------------------------------------------------------------------
    !
    !  This subroutine prepares some variables which are needed for derivatives
    !  1) Set the nonlinear core correction
    !  2) Computes dmuxc (derivative of the XC potential)
    !  3) Set gradient correction (GGA) if needed
    !
    USE kinds,                 ONLY : DP
    USE ions_base,             ONLY : ntyp => nsp
    USE fft_base,              ONLY : dfftp
    USE uspp_param,            ONLY : upf
    USE uspp,                  ONLY : nlcc_any
    USE eqv,                   ONLY : dmuxc
    USE xc_lib,                ONLY : xclib_dft_is
    USE lsda_mod,              ONLY : nspin
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: l_skip
    !
    CALL start_clock('dv_setup')
    !
    ! 0) allocate dmuxc
    !
    ALLOCATE(dmuxc(dfftp%nnr,nspin,nspin))
    !
    IF(l_skip) THEN
       !
       dmuxc = (0._DP,0._DP)
       !
    ELSE
       !
       ! 1) Set the nonlinear core correction
       !
       nlcc_any = ANY(upf(1:ntyp)%nlcc)
       !
       ! 2) Compute the derivative of the XC potential
       !
       CALL setup_dmuxc()
       !
       ! 3) Setup gradient correction
       !
       IF(xclib_dft_is('gradient')) THEN
          !
          CALL setup_dgc()
          !
       ENDIF
       !
    ENDIF
    !
    CALL stop_clock('dv_setup')
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE wbse_dv_of_drho(dvscf, lrpa, add_nlcc, drhoc)
    !-----------------------------------------------------------------------
    !
    !  This routine computes the change of the self consistent potential
    !  (Hartree and XC) due to the perturbation.
    !  Note: gamma_only is disregarded for PHonon calculations,
    !  TDDFPT purposes only.
    !
    USE kinds,             ONLY : DP
    USE constants,         ONLY : e2,fpi
    USE fft_base,          ONLY : dfftp
    USE fft_interfaces,    ONLY : fwfft,invfft
    USE gvect,             ONLY : ngm,g
    USE cell_base,         ONLY : tpiba2,omega
    USE noncollin_module,  ONLY : nspin_gga
    USE lsda_mod,          ONLY : nspin
    USE xc_lib,            ONLY : xclib_dft_is
    USE funct,             ONLY : dft_is_nonlocc
    USE scf,               ONLY : rho,rho_core
    USE uspp,              ONLY : nlcc_any
    USE control_flags,     ONLY : gamma_only
    USE martyna_tuckerman, ONLY : wg_corr_h,do_comp_mt
    USE qpoint,            ONLY : xq
    USE gc_lr,             ONLY : grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s
    USE eqv,               ONLY : dmuxc
#if defined(__CUDA)
    USE west_gpu,          ONLY : dvaux,dvhart,dfft_nl_d,dfft_nlm_d
#endif
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin)
    ! input:  response charge density
    ! output: response Hartree-and-XC potential
    !
    LOGICAL, INTENT(IN) :: lrpa
    ! input: if true, only hatree term is computed
    !
    LOGICAL, INTENT(IN) :: add_nlcc
    ! input: if true, add core charge density
    !
    COMPLEX(DP), INTENT(IN), OPTIONAL :: drhoc(dfftp%nnr)
    ! input: response core charge density
    ! (needed only for PHonon when add_nlcc=.true.)
    !
    INTEGER :: is, is1, ig, ir, dfftp_nnr
    ! counter on r vectors
    ! counter on spin polarizations
    ! counter on g vectors
    !
    REAL(DP) :: xq1, xq2, xq3
    REAL(DP) :: qg2, fac, eh_corr
    ! qg2: the modulus of (q+G)^2
    ! fac: the structure factor
    ! eh_corr: the correction to response Hartree energy due
    ! to Martyna-Tuckerman correction (calculated, but not used).
    !
#if !defined(__CUDA)
    COMPLEX(DP), ALLOCATABLE :: dvaux(:,:), dvhart(:), dvaux_mt(:)
#endif
    ! dvaux: response XC potential
    ! dvhart: response Hartree potential
    ! dvaux_mt: auxiliary array for Martyna-Tuckerman correction
    !
#if defined(__CUDA)
    CALL start_clock_gpu('dv_drho')
#else
    CALL start_clock('dv_drho')
#endif
    !
    dfftp_nnr = dfftp%nnr
    xq1 = xq(1)
    xq2 = xq(2)
    xq3 = xq(3)
    !
#if defined(__CUDA)
    IF(add_nlcc) CALL errore('wbse_dv_of_drho', 'add_nlcc not supported on GPUs', 1)
    IF(do_comp_mt) CALL errore('wbse_dv_of_drho', 'do_comp_mt not supported on GPUs', 1)
#endif
    !
    IF(add_nlcc .AND. .NOT. PRESENT(drhoc)) &
    & CALL errore('wbse_dv_of_drho', 'drhoc is not present in the input of the routine', 1)
    !
    IF(lrpa) THEN
#if !defined(__CUDA)
       ALLOCATE(dvaux(1,1))
#endif
    ELSE
#if !defined(__CUDA)
       ALLOCATE(dvaux(dfftp%nnr, nspin))
#endif
       !
       ! 1) The exchange-correlation contribution is computed in real space
       !
       fac = 1._DP / REAL(nspin,KIND=DP)
       !
       IF(nlcc_any .AND. add_nlcc) THEN
          DO is = 1, nspin
             rho%of_r(:,is) = rho%of_r(:,is) + fac * rho_core(:)
             dvscf(:,is) = dvscf(:,is) + fac * drhoc(:)
          ENDDO
       ENDIF
       !
       dvaux(:,:) = (0._DP,0._DP)
       DO is = 1, nspin
          DO is1 = 1, nspin
             dvaux(:,is) = dvaux(:,is) + dmuxc(:,is,is1) * dvscf(:,is1)
          ENDDO
       ENDDO
       !
       ! Add gradient correction to the response XC potential.
       ! NB: If nlcc=.true. we need to add here its contribution.
       ! grho contains already the core charge
       !
       IF(xclib_dft_is('gradient')) THEN
          CALL dgradcorr(dfftp, rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
               & dvscf, nspin, nspin_gga, g, dvaux)
       ENDIF
       !
       IF(dft_is_nonlocc()) THEN
          CALL dnonloccorr(rho%of_r, dvscf, xq, dvaux)
       ENDIF
       !
       IF(nlcc_any .AND. add_nlcc) THEN
          DO is = 1, nspin
             rho%of_r(:,is) = rho%of_r(:,is) - fac * rho_core(:)
             dvscf(:,is) = dvscf(:,is) - fac * drhoc(:)
          ENDDO
       ENDIF
       !
       !$acc update device(dvaux)
       !
    ENDIF
    !
    ! 2) Hartree contribution is computed in reciprocal space
    !
    IF(nspin == 2) THEN
       !$acc parallel loop present(dvscf)
       DO ir = 1, dfftp_nnr
          dvscf(ir,1) = dvscf(ir,1) + dvscf(ir,2)
       ENDDO
       !$acc end parallel
    ENDIF
    !
    !$acc host_data use_device(dvscf)
    CALL fwfft('Rho', dvscf(:,1), dfftp)
    !$acc end host_data
    !
#if defined(__CUDA)
    !$acc kernels present(dvhart)
    dvhart(:) = (0._DP,0._DP)
    !$acc end kernels
    !
    !$acc parallel loop present(g,dvhart,dvscf)
    DO ig = 1, ngm
       !
       qg2 = (g(1,ig)+xq1)**2 + (g(2,ig)+xq2)**2 + (g(3,ig)+xq3)**2
       !
       IF(qg2 > 1.E-8_DP) THEN
          dvhart(dfft_nl_d(ig)) = e2 * fpi * dvscf(dfft_nl_d(ig),1) / (tpiba2 * qg2)
       ENDIF
       !
    ENDDO
    !$acc end parallel
#else
    ALLOCATE(dvhart(dfftp%nnr))
    dvhart(:) = (0._DP,0._DP)
    !
    DO ig = 1, ngm
       !
       qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
       !
       IF(qg2 > 1.E-8_DP) THEN
          dvhart(dfftp%nl(ig)) = e2 * fpi * dvscf(dfftp%nl(ig),1) / (tpiba2 * qg2)
       ENDIF
       !
    ENDDO
    !
    IF(do_comp_mt) THEN
       !
       ALLOCATE(dvaux_mt(ngm))
       !
       CALL wg_corr_h(omega, ngm, dvscf(dfftp%nl(:),1), dvaux_mt, eh_corr)
       !
       DO ig = 1, ngm
          dvhart(dfftp%nl(ig)) = dvhart(dfftp%nl(ig)) + dvaux_mt(ig)
       ENDDO
       !
       DEALLOCATE(dvaux_mt)
       !
    ENDIF
#endif
    !
    IF(gamma_only) THEN
#if defined(__CUDA)
       !$acc parallel loop present(dvhart)
       DO ig = 1, ngm
          dvhart(dfft_nlm_d(ig)) = CONJG(dvhart(dfft_nl_d(ig)))
       ENDDO
       !$acc end parallel
#else
       DO ig = 1, ngm
          dvhart(dfftp%nlm(ig)) = CONJG(dvhart(dfftp%nl(ig)))
       ENDDO
#endif
    ENDIF
    !
    ! Transformed back to real space
    !
    !$acc host_data use_device(dvhart)
    CALL invfft('Rho', dvhart, dfftp)
    !$acc end host_data
    !
    IF(lrpa) THEN
       !$acc parallel loop collapse(2) present(dvscf,dvhart)
       DO is = 1, nspin
          DO ir = 1, dfftp_nnr
             dvscf(ir,is) = dvhart(ir)
          ENDDO
       ENDDO
       !$acc end parallel
    ELSE
       !$acc parallel loop collapse(2) present(dvscf,dvhart,dvaux)
       DO is = 1, nspin
          DO ir = 1, dfftp_nnr
             dvscf(ir,is) = dvhart(ir) + dvaux(ir,is)
          ENDDO
       ENDDO
       !$acc end parallel
    ENDIF
    !
#if !defined(__CUDA)
    DEALLOCATE(dvaux)
    DEALLOCATE(dvhart)
#endif
    !
#if defined(__CUDA)
    CALL stop_clock_gpu('dv_drho')
#else
    CALL stop_clock('dv_drho')
#endif
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE wbse_sf_kernel_setup()
    !-----------------------------------------------------------------------
    !
    !  This subroutine computes the spin-flip kernel for spin-flip TDDFT
    !
    USE kinds,                 ONLY : DP
    USE ions_base,             ONLY : ntyp => nsp
    USE fft_base,              ONLY : dfftp
    USE uspp_param,            ONLY : upf
    USE uspp,                  ONLY : nlcc_any
    USE xc_lib,                ONLY : xclib_dft_is
    USE lsda_mod,              ONLY : nspin
    USE westcom,               ONLY : wbse_save_dir,sf_kernel,l_spin_flip_alda0,spin_flip_cut1,&
                                    & l_print_spin_flip_kernel
    USE scf,                   ONLY : rho,rho_core,rhog_core
    USE xc_lib,                ONLY : xc
    USE martyna_tuckerman,     ONLY : do_comp_mt
    USE funct,                 ONLY : nlc,dft_is_nonlocc
    USE mp,                    ONLY : mp_barrier
    USE mp_world,              ONLY : world_comm
    USE constants,             ONLY : e2
    USE cubefile,              ONLY : write_wfc_cube_r
    !
    IMPLICIT NONE
    !
    ! local variables
    !
    INTEGER :: ir, is
    INTEGER :: dfftp_nnr
    REAL(DP), ALLOCATABLE :: vxc(:,:), ex(:), ec(:), vx(:,:), vc(:,:)
    REAL(DP) :: etxc, vtxc
    CHARACTER(LEN=:), ALLOCATABLE :: fname
    !
    CALL start_clock('sf_kernel')
    !
    nlcc_any = ANY(upf(1:ntyp)%nlcc)
    !
    IF(nlcc_any) CALL errore('wbse_sf_kernel_setup', 'add_nlcc not supported', 1)
    IF(do_comp_mt) CALL errore('wbse_sf_kernel_setup', 'do_comp_mt not supported', 1)
    IF(dft_is_nonlocc()) CALL errore('wbse_sf_kernel_setup', 'dft_is_nonlocc not supported', 1)
    !
    ALLOCATE(sf_kernel(dfftp%nnr))
    ALLOCATE(vxc(dfftp%nnr,nspin))
    ALLOCATE(ex(dfftp%nnr))
    ALLOCATE(vx(dfftp%nnr,nspin))
    ALLOCATE(ec(dfftp%nnr))
    ALLOCATE(vc(dfftp%nnr,nspin))
    ex(:) = 0._DP
    ec(:) = 0._DP
    vx(:,:) = 0._DP
    vc(:,:) = 0._DP
    !
    !$acc enter data create(sf_kernel,vxc,ex,ec,vx,vc) copyin(rho%of_r,rho%of_g)
    !
    ! ... spin-polarized case
    !
    CALL xc( dfftp%nnr, 2, 2, rho%of_r, ex, ec, vx, vc, gpu_args_=.TRUE. )
    !
    dfftp_nnr = dfftp%nnr
    !
    !$acc parallel loop collapse(2) present(vxc,vx,vc)
    DO is = 1, nspin
       DO ir = 1, dfftp_nnr
          vxc(ir,is) = e2*(vx(ir,is)+vc(ir,is))
       ENDDO
    ENDDO
    !$acc end parallel
    !
    !$acc exit data delete(ex,ec,vx,vc)
    DEALLOCATE(ex)
    DEALLOCATE(vx)
    DEALLOCATE(ec)
    DEALLOCATE(vc)
    !
    ! ... add gradient corrections (if any)
    !
    IF(.NOT. l_spin_flip_alda0 .AND. xclib_dft_is('gradient')) THEN
       etxc = 0._DP
       vtxc = 0._DP
       !
       !$acc enter data copyin(rho_core,rhog_core)
       CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, vxc )
       !$acc exit data delete(rho_core,rhog_core)
    ENDIF
    !
    ! IF nspin == 2, rho%of_r(:,1) and rho%of_r(:,2) stores up+down and up-down respectively
    !
    !$acc parallel loop present(sf_kernel,vxc,rho%of_r)
    DO ir = 1, dfftp_nnr
       sf_kernel(ir) = (vxc(ir,1) - vxc(ir,2)) / rho%of_r(ir,2)
    ENDDO
    !$acc end parallel
    !
    IF(.NOT. l_spin_flip_alda0 .AND. xclib_dft_is('gradient')) THEN
       !$acc parallel loop present(sf_kernel)
       DO ir = 1, dfftp_nnr
          IF(ABS(sf_kernel(ir)) > spin_flip_cut1) sf_kernel(ir) = 0._DP
       ENDDO
       !$acc end parallel
    ENDIF
    !
    !$acc exit data delete(vxc,rho%of_r,rho%of_g)
    DEALLOCATE(vxc)
    !
    ! print spin flip kernel
    !
    IF(l_print_spin_flip_kernel) THEN
       CALL mp_barrier( world_comm )
       !
       !$acc update host(sf_kernel)
       fname = TRIM(wbse_save_dir)//'/sf_kernel.cube'
       CALL write_wfc_cube_r(dfftp, fname, sf_kernel)
       !
       CALL mp_barrier( world_comm )
    ENDIF
    !
    CALL stop_clock('sf_kernel')
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE wbse_dv_of_drho_sf(dvscf)
    !-----------------------------------------------------------------------
    !
    !  This routine computes the change of the self consistent potential
    !  (XC) due to the perturbation in spin-flip calculations
    !  Note: gamma_only is disregarded for PHonon calculations,
    !  TDDFPT purposes only.
    !
    USE kinds,             ONLY : DP
    USE fft_base,          ONLY : dfftp
    USE lsda_mod,          ONLY : nspin
    USE uspp,              ONLY : nlcc_any
    USE martyna_tuckerman, ONLY : do_comp_mt
    USE westcom,           ONLY : sf_kernel
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin)
    ! input:  response charge density
    ! output: response XC potential
    !
    INTEGER :: is, ir, dfftp_nnr
    !
#if defined(__CUDA)
    CALL start_clock_gpu('dv_drho_sf')
#else
    CALL start_clock('dv_drho_sf')
#endif
    !
    IF(nlcc_any) CALL errore('wbse_dv_of_drho_sf', 'nlcc_any not supported', 1)
    IF(do_comp_mt) CALL errore('wbse_dv_of_drho_sf', 'do_comp_mt not supported', 1)
    !
    dfftp_nnr = dfftp%nnr
    !
    !$acc parallel loop collapse(2) present(dvscf,sf_kernel)
    DO is = 1, nspin
       DO ir = 1, dfftp_nnr
          dvscf(ir,is) = sf_kernel(ir) * dvscf(ir,is)
       ENDDO
    ENDDO
    !$acc end parallel
    !
#if defined(__CUDA)
    CALL stop_clock_gpu('dv_drho_sf')
#else
    CALL stop_clock('dv_drho_sf')
#endif
    !
  END SUBROUTINE
  !
END MODULE
