!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine west_dv_of_drho (dvscf, lrpa, add_nlcc, drhoc)
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
  USE gvect,             ONLY : nl, ngm, g,nlm, gstart
  USE cell_base,         ONLY : alat, tpiba2, omega
  USE noncollin_module,  ONLY : nspin_gga
  USE lsda_mod,          ONLY : nspin
  USE funct,             ONLY : dft_is_gradient, dft_is_nonlocc
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
  INTEGER :: ir, is, is1, ig
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
  CALL start_clock ('wbse_dv_of_drho')
  !
  if (add_nlcc .and. .not.present(drhoc)) &
     & CALL errore( 'wbse_dv_of_drho', 'drhoc is not present in the input of the routine', 1 )
  !
  if (lrpa) then
     !
     allocate(dvaux(1,1))
     goto 111
     !
  else
     !
     allocate(dvaux(dfftp%nnr, nspin))
     !
  endif
  !
  ! 1) The exchange-correlation contribution is computed in real space
  !
  fac = 1.d0 / DBLE (nspin)
  !
  if (nlcc_any.and.add_nlcc) then
     !
     do is = 1, nspin
        !
        rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core (:)
        !
        dvscf(:, is) = dvscf(:, is) + fac * drhoc (:)
        !
     enddo
     !
  endif
  !
  dvaux (:,:) = (0.d0, 0.d0)
  do is = 1, nspin
     !
     do is1 = 1, nspin
        !
        dvaux(:,is) = dvaux(:,is) + dmuxc(:,is,is1) * dvscf(:,is1)
        !
     enddo
     !
  enddo
  !
  ! Add gradient correction to the response XC potential.
  ! NB: If nlcc=.true. we need to add here its contribution.
  ! grho contains already the core charge
  !
  if ( dft_is_gradient() ) then
     !
     call dgradcorr (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
            dvscf, dfftp%nnr, nspin, nspin_gga, nl, ngm, g, alat, dvaux )
     !
  endif
  !
  if (dft_is_nonlocc()) then
     !
     call dnonloccorr(rho%of_r, dvscf, xq, dvaux)
     !
  endif
  !
  if (nlcc_any.and.add_nlcc) then
     !
     do is = 1, nspin
        !
        rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) - fac * drhoc (:)
        !
     enddo
     !
  endif
  !
111 continue
  !
  ! 2) Hartree contribution is computed in reciprocal space
  !
  if (nspin == 2) then
     !
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
     !
  endif
  !
  call fwfft ('Dense', dvscf(:,1), dfftp)
  !
  allocate (dvhart(dfftp%nnr, nspin))
  !
  dvhart(:,:) = (0.d0,0.d0)
  !
  do ig = 1, ngm
     !
     qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
     !
     if (qg2 > 1.d-8) then
        !
        dvhart(nl(ig),:) = e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
        !
     endif
     !
  enddo
  !
  if (do_comp_mt) then
     !
     allocate(dvaux_mt(ngm))
     !
     call wg_corr_h (omega, ngm, dvscf(nl(:),1), dvaux_mt, eh_corr)
     !
     do ig = 1, ngm
        !
        dvhart(nl(ig),:) = dvhart(nl(ig),1) + dvaux_mt(ig)
        !
     enddo
     !
     deallocate(dvaux_mt)
     !
  endif
  !
  if (gamma_only) then
     !
     do ig = 1, ngm
        !
        dvhart(nlm(ig),:) = conjg(dvhart(nl(ig),:))
        !
     enddo
     !
  endif
  !
  ! Transformed back to real space
  !
  do is = 1, nspin
     !
     call invfft ('Dense', dvhart (:,is), dfftp)
     !
  enddo
  !
  dvscf(:,:) = (0.0_DP,0.0_DP)
  !
  if (lrpa) then
     !
     dvscf(:,:) = dvhart(:,:)
     !
  else
     !
     dvscf(:,:) = dvaux(:,:) + dvhart(:,:)
     !
  endif
  !
  deallocate (dvaux)
  !
  CALL stop_clock ('wbse_dv_of_drho')
  !
  RETURN
  !
end subroutine west_dv_of_drho
