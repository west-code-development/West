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
    USE constants,         ONLY : e2, fpi
    USE fft_base,          ONLY : dfftp
    USE fft_interfaces,    ONLY : fwfft, invfft
    USE gvect,             ONLY : ngm, g
    USE cell_base,         ONLY : tpiba2, omega
    USE noncollin_module,  ONLY : nspin_gga
    USE lsda_mod,          ONLY : nspin
    USE xc_lib,            ONLY : xclib_dft_is
    USE funct,             ONLY : dft_is_nonlocc
    USE scf,               ONLY : rho, rho_core
    USE uspp,              ONLY : nlcc_any
    USE control_flags,     ONLY : gamma_only
    USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
    USE qpoint,            ONLY : xq
    USE gc_lr,             ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
    USE eqv,               ONLY : dmuxc
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
    INTEGER :: is, is1, ig
    ! counter on r vectors
    ! counter on spin polarizations
    ! counter on g vectors
    !
    REAL(DP) :: qg2, fac, eh_corr
    ! qg2: the modulus of (q+G)^2
    ! fac: the structure factor
    ! eh_corr: the correction to response Hartree energy due
    ! to Martyna-Tuckerman correction (calculated, but not used).
    !
    COMPLEX(DP), ALLOCATABLE :: dvaux(:,:), dvhart(:,:), dvaux_mt(:)
    ! dvaux: response XC potential
    ! dvhart: response Hartree potential
    ! dvaux_mt: auxiliary array for Martyna-Tuckerman correction
    !
    CALL start_clock('dv_of_drho')
    !
    IF(add_nlcc .AND. .NOT. PRESENT(drhoc)) &
    & CALL errore('wbse_dv_of_drho', 'drhoc is not present in the input of the routine', 1)
    !
    IF(lrpa) THEN
       ALLOCATE(dvaux(1,1))
    ELSE
       ALLOCATE(dvaux(dfftp%nnr, nspin))
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
    ENDIF
    !
    ! 2) Hartree contribution is computed in reciprocal space
    !
    IF(nspin == 2) THEN
       dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
    ENDIF
    !
    CALL fwfft('Rho', dvscf(:,1), dfftp)
    !
    ALLOCATE(dvhart(dfftp%nnr, nspin))
    !
    dvhart(:,:) = (0._DP,0._DP)
    !
    DO ig = 1, ngm
       !
       qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
       !
       IF(qg2 > 1.E-8_DP) THEN
          dvhart(dfftp%nl(ig),:) = e2 * fpi * dvscf(dfftp%nl(ig),1) / (tpiba2 * qg2)
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
          dvhart(dfftp%nl(ig),:) = dvhart(dfftp%nl(ig),1) + dvaux_mt(ig)
       ENDDO
       !
       DEALLOCATE(dvaux_mt)
       !
    ENDIF
    !
    IF(gamma_only) THEN
       DO ig = 1, ngm
          dvhart(dfftp%nlm(ig),:) = CONJG(dvhart(dfftp%nl(ig),:))
       ENDDO
    ENDIF
    !
    ! Transformed back to real space
    !
    DO is = 1, nspin
       CALL invfft('Rho', dvhart(:,is), dfftp)
    ENDDO
    !
    IF(lrpa) THEN
       dvscf(:,:) = dvhart
    ELSE
       dvscf(:,:) = dvaux + dvhart
    ENDIF
    !
    DEALLOCATE(dvaux)
    DEALLOCATE(dvhart)
    !
    CALL stop_clock('dv_of_drho')
    !
  END SUBROUTINE
  !
#if defined(__CUDA)
  !-----------------------------------------------------------------------
  SUBROUTINE wbse_dv_of_drho_gpu(dvscf, lrpa, add_nlcc)
    !-----------------------------------------------------------------------
    !
    USE kinds,             ONLY : DP
    USE constants,         ONLY : e2, fpi
    USE fft_base,          ONLY : dfftp
    USE fft_interfaces,    ONLY : fwfft, invfft
    USE gvect,             ONLY : ngm, g
    USE cell_base,         ONLY : tpiba2
    USE noncollin_module,  ONLY : nspin_gga
    USE lsda_mod,          ONLY : nspin
    USE xc_lib,            ONLY : xclib_dft_is
    USE scf,               ONLY : rho
    USE control_flags,     ONLY : gamma_only
    USE martyna_tuckerman, ONLY : do_comp_mt
    USE qpoint,            ONLY : xq
    USE gc_lr,             ONLY : grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s
    USE eqv,               ONLY : dmuxc
    USE west_gpu,          ONLY : dvaux,dvhart,dfft_nl_d,dfft_nlm_d
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin)
    LOGICAL, INTENT(IN) :: lrpa
    LOGICAL, INTENT(IN) :: add_nlcc
    !
    ! Workspace
    !
    INTEGER :: is, is1, ig, ir, dfftp_nnr
    REAL(DP) :: xq1, xq2, xq3, qg2
    !
    CALL start_clock_gpu('dv_of_drho')
    !
    IF(add_nlcc) CALL errore('wbse_dv_of_drho_gpu', 'add_nlcc is not supported', 1)
    IF(do_comp_mt) CALL errore('wbse_dv_of_drho_gpu', 'do_comp_mt is not supported', 1)
    !
    IF(.NOT. lrpa) THEN
       !
       ! Compute exchange-correlation contribution in real space
       !
       dvaux(:,:) = (0._DP,0._DP)
       DO is = 1, nspin
          DO is1 = 1, nspin
             dvaux(:,is) = dvaux(:,is) + dmuxc(:,is,is1) * dvscf(:,is1)
          ENDDO
       ENDDO
       !
       ! Add gradient correction to the response XC potential
       !
       IF(xclib_dft_is('gradient')) THEN
          CALL dgradcorr(dfftp, rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
               & dvscf, nspin, nspin_gga, g, dvaux)
       ENDIF
       !
       !$acc update device(dvaux)
       !
    ENDIF
    !
    ! Compute Hartree contribution in reciprocal space
    !
    dfftp_nnr = dfftp%nnr
    xq1 = xq(1)
    xq2 = xq(2)
    xq3 = xq(3)
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
    !
    IF(gamma_only) THEN
       !$acc parallel loop present(dvhart)
       DO ig = 1, ngm
          dvhart(dfft_nlm_d(ig)) = CONJG(dvhart(dfft_nl_d(ig)))
       ENDDO
       !$acc end parallel
    ENDIF
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
    CALL stop_clock_gpu('dv_of_drho')
    !
  END SUBROUTINE
#endif
  !
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
    USE eqv,                   ONLY : dmuxc
    USE xc_lib,                ONLY : xclib_dft_is
    USE lsda_mod,              ONLY : nspin
    USE westcom,               ONLY : sf_kernel, l_sf_alda0, l_sf_cut1, l_print_sf_kernel
    USE scf,                   ONLY : rho, rho_core, rhog_core
    USE xc_lib,                ONLY : xc
    USE martyna_tuckerman,     ONLY : do_comp_mt
    USE cell_base,             ONLY : omega
    USE funct,                 ONLY : nlc, dft_is_nonlocc
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum, mp_barrier
    USE mp_world,              ONLY : world_comm
    USE constants,             ONLY : e2, eps8
    USE cubefile,              ONLY : write_wfc_cube_r
    !
    IMPLICIT NONE
    !
    ! local variables
    !
    INTEGER :: ir
    REAL(DP), ALLOCATABLE :: vxc(:,:), ex(:), ec(:), vx(:,:), vc(:,:)
    REAL(DP) :: etxc, vtxc
    REAL(DP) :: rho_up, rho_down, rho_diff
    REAL(DP), PARAMETER :: small = 1.E-30_DP, medium = 1.E-15_DP, large = 1.E-10_DP
    CHARACTER(LEN=20) :: prefix, ind, filename
    !
    CALL start_clock('sf_kernel_setup')
    !
    nlcc_any = ANY(upf(1:ntyp)%nlcc)
    !
    IF(nlcc_any) CALL errore('wbse_sf_kernel_setup', 'add_nlcc is not supported', 1)
    IF(do_comp_mt) CALL errore('wbse_sf_kernel_setup', 'do_comp_mt is not supported', 1)
    !
    ALLOCATE(sf_kernel(dfftp%nnr))
    sf_kernel(:) = 0._DP
    !
    ALLOCATE(vxc(dfftp%nnr,nspin))
    vxc(:,:) = 0._DP
    !
    etxc = 0.D0
    vtxc = 0.D0
    !
    !$acc data copyin( rho%of_r, rho%of_g )
    !
    ALLOCATE( ex(dfftp%nnr), vx(dfftp%nnr,nspin) )
    ALLOCATE( ec(dfftp%nnr), vc(dfftp%nnr,nspin) )
    ex(:) = 0._DP
    ec(:) = 0._DP
    vx(:,:) = 0._DP
    vc(:,:) = 0._DP
    !
    !$acc data create( ex, ec, vx, vc )
    !
    ! ... spin-polarized case
    !
    CALL xc( dfftp%nnr, 2, 2, rho%of_r, ex, ec, vx, vc, gpu_args_=.TRUE. )
    !
    !$acc parallel loop reduction(+:etxc,vtxc) reduction(-:rhoneg1,rhoneg2) &
    !$acc&              present(rho)
    DO ir = 1, dfftp%nnr
       vxc(ir,1) = e2*( vx(ir,1) + vc(ir,1) )
       vxc(ir,2) = e2*( vx(ir,2) + vc(ir,2) )
       etxc = etxc + e2*( (ex(ir) + ec(ir))*rho%of_r(ir,1) )
       vtxc = vtxc + ( ( vxc(ir,1) + vxc(ir,2) )*rho%of_r(ir,1) + &
                       ( vxc(ir,1) - vxc(ir,2) )*rho%of_r(ir,2) )*0.5d0
    ENDDO
    !
    !$acc end data
    DEALLOCATE( ex, vx )
    DEALLOCATE( ec, vc )
    !
    ! ... energy terms, local-density contribution
    !
    vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
    etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
    !
    ! ... add gradient corrections (if any)
    !
    IF(xclib_dft_is('gradient')) THEN
       !
       IF(.NOT. l_sf_alda0) THEN
          !
          CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, vxc )
          !
       ENDIF
       !
    ENDIF
    !
    !$acc end data
    !$acc end data
    !
    ! ... add non local corrections (if any)
    ! ... should not work in principle
    IF(dft_is_nonlocc()) CALL errore('wbse_sf_kernel_setup', 'dft_is_nonlocc is not supported', 1)
    !IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
    !
    CALL mp_sum(  vtxc , intra_bgrp_comm )
    CALL mp_sum(  etxc , intra_bgrp_comm )
    !
    vtxc = omega * vtxc / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
    etxc = omega * etxc / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
    !
    DO ir = 1, dfftp%nnr
       !
       !rho_up   = rho%of_r (ir, 1)
       !rho_down = rho%of_r (ir, 2)
       rho_up   = 0.5_DP * (rho%of_r(ir,1) + rho%of_r(ir,2))
       rho_down = 0.5_DP * (rho%of_r(ir,1) - rho%of_r(ir,2))
       !
       rho_diff = 0.0_DP
       rho_diff = rho_up - rho_down
       !
       sf_kernel(ir) = (vxc(ir,1) - vxc(ir,2)) / rho_diff
       !
       IF (.NOT. l_sf_alda0) THEN
          !
          IF (xclib_dft_is('gradient') .AND. ABS(sf_kernel(ir)) > l_sf_cut1) THEN
             ! remove the divergent part
             sf_kernel(ir) = 0.0_DP
             !
          ENDIF
          !
       ENDIF
       !
    ENDDO
    !
    DEALLOCATE(vxc)
    !
    ! TEST print sflip_kernel
    !
    CALL mp_barrier( world_comm )
    !
    IF (l_print_sf_kernel) THEN
       !
       prefix='sf_kernel'
       !
       filename=TRIM(prefix)//'.cube'
       !
       CALL write_wfc_cube_r(dfftp, 1002, filename, sf_kernel)
       !
    ENDIF
    !
    CALL mp_barrier( world_comm )
    !
    CALL stop_clock('sf_kernel_setup')
    !
  END SUBROUTINE
  !
!  !-----------------------------------------------------------------------
!  SUBROUTINE wbse_sf_kernel_setup()
!    !-----------------------------------------------------------------------
!    !
!    !  This subroutine computes the spin-flip kernel for spin-flip TDDFT
!    !
!    USE kinds,                 ONLY : DP
!    USE ions_base,             ONLY : ntyp => nsp
!    USE fft_base,              ONLY : dfftp
!    USE uspp_param,            ONLY : upf
!    USE uspp,                  ONLY : nlcc_any
!    USE eqv,                   ONLY : dmuxc
!    USE xc_lib,                ONLY : xclib_dft_is
!    USE lsda_mod,              ONLY : nspin
!    USE westcom,               ONLY : sf_kernel, l_sf_alda0, l_sf_cut1, l_print_sf_kernel
!    USE scf,                   ONLY : rho, rho_core, rhog_core
!    USE xc_lib,                ONLY : xc
!    USE martyna_tuckerman,     ONLY : do_comp_mt
!    USE cell_base,             ONLY : omega
!    USE funct,                 ONLY : nlc, dft_is_nonlocc
!    USE mp_bands,              ONLY : intra_bgrp_comm
!    USE mp,                    ONLY : mp_sum, mp_barrier
!    USE mp_world,              ONLY : world_comm
!    USE constants,             ONLY : e2, eps8
!    USE cubefile,              ONLY : write_wfc_cube_r
!    !
!    IMPLICIT NONE
!    !
!    ! local variables
!    !
!    INTEGER :: ir
!    REAL(DP), ALLOCATABLE :: vxc(:,:), ex(:), ec(:), vx(:,:), vc(:,:)
!    REAL(DP) :: etxc, vtxc
!    REAL(DP) :: rho_up, rho_down, rho_diff
!    REAL(DP), PARAMETER :: small = 1.E-30_DP, medium = 1.E-15_DP, large = 1.E-10_DP
!    CHARACTER(LEN=20) :: prefix, ind, filename
!    !
!    CALL start_clock('sf_kernel_setup')
!    !
!    nlcc_any = ANY(upf(1:ntyp)%nlcc)
!    !
!    IF(nlcc_any) CALL errore('wbse_sf_kernel_setup', 'add_nlcc is not supported', 1)
!    IF(do_comp_mt) CALL errore('wbse_sf_kernel_setup', 'do_comp_mt is not supported', 1)
!    !
!    ALLOCATE(sf_kernel(dfftp%nnr))
!    sf_kernel(:) = 0._DP
!    !
!    ALLOCATE (vxc(dfftp%nnr,nspin))
!    CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
!    !
!    DO ir = 1, dfftp%nnr
!       !
!       rho_up   = rho%of_r (ir, 1)
!       rho_down = rho%of_r (ir, 2)
!       !
!       rho_diff = 0.0_DP
!       rho_diff = rho_up - rho_down
!       !
!       sf_kernel(ir) = (vxc(ir,1) - vxc(ir,2)) / rho_diff
!       !
!       IF (.NOT. l_sf_alda0) THEN
!          !
!          IF (xclib_dft_is('gradient') .AND. ABS(sf_kernel(ir)) > l_sf_cut1) THEN
!             ! remove the divergent part
!             sf_kernel(ir) = 0.0_DP
!             !
!          ENDIF
!          !
!       ENDIF
!       !
!    ENDDO
!    !
!    DEALLOCATE(vxc)
!    !
!    ! TEST print sflip_kernel
!    !
!    CALL mp_barrier( world_comm )
!    !
!    IF (l_print_sf_kernel) THEN
!       !
!       prefix='sf_kernel'
!       !
!       filename=TRIM(prefix)//'.cube'
!       !
!       CALL write_wfc_cube_r(dfftp, 1002, filename, sf_kernel)
!       !
!    ENDIF
!    !
!    CALL mp_barrier( world_comm )
!    !
!    CALL stop_clock('sf_kernel_setup')
!    !
!  END SUBROUTINE

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
    USE constants,         ONLY : e2, fpi
    USE fft_base,          ONLY : dfftp
    USE fft_interfaces,    ONLY : fwfft, invfft
    USE gvect,             ONLY : ngm, g
    USE cell_base,         ONLY : tpiba2, omega
    USE noncollin_module,  ONLY : nspin_gga
    USE lsda_mod,          ONLY : nspin
    USE xc_lib,            ONLY : xclib_dft_is
    USE funct,             ONLY : dft_is_nonlocc
    USE scf,               ONLY : rho, rho_core
    USE uspp,              ONLY : nlcc_any
    USE control_flags,     ONLY : gamma_only
    USE martyna_tuckerman, ONLY : do_comp_mt
    USE qpoint,            ONLY : xq
    USE gc_lr,             ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
    USE westcom,           ONLY : sf_kernel
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin)
    ! input:  response charge density
    ! output: response XC potential
    !
    INTEGER :: is
    ! counter on r vectors
    ! counter on spin polarizations
    ! counter on g vectors
    !
    COMPLEX(DP), ALLOCATABLE :: dvaux(:,:)
    ! dvaux: response XC potential
    !
    CALL start_clock('dv_of_drho_sf')
    !
    IF(nlcc_any) CALL errore('wbse_dv_of_drho_sf', 'nlcc_any is not supported', 1)
    IF(do_comp_mt) CALL errore('wbse_dv_of_drho_sf', 'do_comp_mt is not supported', 1)
    !
    ALLOCATE(dvaux(dfftp%nnr, nspin))
    !
    dvaux(:,:) = (0._DP,0._DP)
    DO is = 1, nspin
       dvaux(:,is) = sf_kernel(:) * dvscf(:,is)
    ENDDO
    !
    dvscf(:,:) = dvaux(:,:)
    !
    DEALLOCATE(dvaux)
    !
    CALL stop_clock('dv_of_drho_sf')
    !
  END SUBROUTINE
  !
#if defined(__CUDA)
  !-----------------------------------------------------------------------
  SUBROUTINE wbse_dv_of_drho_sf_gpu(dvscf)
    !-----------------------------------------------------------------------
    !
    USE kinds,             ONLY : DP
    USE constants,         ONLY : e2, fpi
    USE fft_base,          ONLY : dfftp
    USE fft_interfaces,    ONLY : fwfft, invfft
    USE gvect,             ONLY : ngm, g
    USE cell_base,         ONLY : tpiba2
    USE noncollin_module,  ONLY : nspin_gga
    USE lsda_mod,          ONLY : nspin
    USE xc_lib,            ONLY : xclib_dft_is
    USE scf,               ONLY : rho
    USE control_flags,     ONLY : gamma_only
    USE martyna_tuckerman, ONLY : do_comp_mt
    USE qpoint,            ONLY : xq
    USE gc_lr,             ONLY : grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s
    USE eqv,               ONLY : dmuxc
    USE west_gpu,          ONLY : dvaux,dvhart,dfft_nl_d,dfft_nlm_d
    USE westcom,           ONLY : sf_kernel
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin)
    LOGICAL, INTENT(IN) :: lrpa
    LOGICAL, INTENT(IN) :: add_nlcc
    !
    ! Workspace
    !
    INTEGER :: is, is1, ig, ir, dfftp_nnr
    REAL(DP) :: xq1, xq2, xq3, qg2
    !
    CALL start_clock_gpu('dv_of_drho_sf')
    !
    IF(nlcc_any) CALL errore('wbse_dv_of_drho_sf_gpu', 'nlcc_any is not supported', 1)
    IF(do_comp_mt) CALL errore('wbse_dv_of_drho_sf_gpu', 'do_comp_mt is not supported', 1)
    !
    ! Compute exchange-correlation contribution in real space
    !
    dvaux(:,:) = (0._DP,0._DP)
    DO is = 1, nspin
       dvaux(:,is) = dvaux(:,is) + sf_kernel(:) * dvscf(:,is)
    ENDDO
    !
    !$acc parallel loop collapse(2) present(dvscf,dvhart,dvaux)
    DO is = 1, nspin
       DO ir = 1, dfftp_nnr
          dvscf(ir,is) = dvaux(ir,is)
       ENDDO
    ENDDO
    !$acc end parallel
    ENDIF
    !
    CALL stop_clock_gpu('dv_of_drho_sf')
    !
  END SUBROUTINE
#endif
  !
END MODULE
