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
  SUBROUTINE wbse_dv_setup(l_bse_calc)
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
    USE noncollin_module,      ONLY : noncolin,domag
    USE eqv,                   ONLY : dmuxc
    USE xc_lib,                ONLY : xclib_dft_is
    USE wavefunctions,         ONLY : psic
    USE lsda_mod,              ONLY : nspin
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: l_bse_calc
    !
    CALL start_clock('dv_setup')
    !
    ! 0) allocate dmuxc
    !
    ALLOCATE(dmuxc(dfftp%nnr,nspin,nspin))
    !
    IF(l_bse_calc) THEN
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
          IF(noncolin .AND. domag) THEN
             IF(.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
             psic(:) = (0._DP,0._DP)
          ENDIF
          !
          CALL setup_dgc()
          !
          IF(ALLOCATED(psic)) DEALLOCATE(psic)
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
       IF (qg2 > 1.E-8_DP) THEN
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
       CALL invfft ('Rho', dvhart (:,is), dfftp)
    ENDDO
    !
    IF(lrpa) THEN
       dvscf(:,:) = dvhart
    ELSE
       dvscf(:,:) = dvaux + dvhart
    ENDIF
    !
    DEALLOCATE(dvaux)
    !
    CALL stop_clock('dv_of_drho')
    !
  END SUBROUTINE
  !
END MODULE
