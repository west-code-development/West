!! Copyright (C) 2015-2016 M. Govoni
!! This file is distributed under the terms of the
!! GNU General Public License. See the file `License'
!! in the root directory of the present distribution,
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! This file is part of WEST.
!!
!! Contributors to this file:
!! Marco Govoni
!!
!#define ZERO ( 0.D0, 0.D0 )
!#define ONE  ( 1.D0, 0.D0 )
!!
!SUBROUTINE wbse_init_fock_energy()
!  !
!  USE kinds,                ONLY : DP
!  USE io_global,            ONLY : stdout
!  USE constants,            ONLY : e2, fpi
!  USE cell_base,            ONLY : alat, tpiba2, omega
!  USE io_push,              ONLY : io_push_title
!  USE types_coulomb,         ONLY : pot3D
!  USE westcom,              ONLY :  exx_etot, nbndval0x
!  !USE westcom,              ONLY : sqvc,isz, exx_etot
!  USE control_flags,        ONLY : gamma_only
!  USE wavefunctions_module, ONLY : evc,psic
!  USE fft_base,             ONLY : dfftp,dffts
!  USE fft_at_gamma,         ONLY : double_invfft_gamma, single_fwfft_gamma
!  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
!  USE pwcom,                ONLY : wk,nks,nelup,neldw,isk,g,igk_k,tpiba2,xk,npw,npwx,lsda,nkstot,&
!                                 & current_k,ngk,nbnd,wg,gstart,ngm,ngms
!  USE mp_bands,             ONLY : intra_bgrp_comm
!  USE mp_global,            ONLY : inter_image_comm
!  USE mp,                   ONLY : mp_sum
!  USE bse_module,           ONLY : l_wannier_repr
!  !wbsecom combined into westcom
!  !USE wbsecom,              ONLY : nbndval0x
!  USE class_idistribute,    ONLY : idistribute
!  USE distribution_center,  ONLY : aband
!  !
!  USE wvfct,                ONLY : g2kin
!  !
!  IMPLICIT NONE
!  !
!  INTEGER :: ibnd, jbnd, ig, itr, num_intergration, nbndval
!  REAL(DP):: ovl_value, etmp0, etmp1, fock_energy, ovl_thr, ekin
!  !REAL(DP), PARAMETER :: list_ovl(8) = (/0.0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.2, 0.99/)
!  !REAL(DP), PARAMETER :: list_ovl(1) = (/0.001/)
!  !REAL(DP), PARAMETER :: list_ovl(1) = (/0.1/)
!  REAL(DP), PARAMETER :: list_ovl(1) = (/0.001/)
!  !REAL(DP), PARAMETER :: list_ovl(1) = (/0.00015/)
!  REAL(DP),    ALLOCATABLE  :: rho_aux(:)
!  REAL(DP),    ALLOCATABLE  :: ovl_matrix(:,:)
!  COMPLEX(DP), ALLOCATABLE  :: evc_loc(:,:), aux1_g(:)
!  COMPLEX(DP), ALLOCATABLE  :: dvg(:), psic_aux(:)
!  !
!  REAL(kind=dp), EXTERNAL   :: DDOT
!  !
!  CALL wbse_init_memory_report()
!  !
!  nbndval = nbndval0x
!  !
!  aband = idistribute()
!  CALL aband%init(nbndval,'i','bse_nbndval',.TRUE.)
!  !
!  IF (.NOT.ALLOCATED(psic)) ALLOCATE (psic(dffts%nnr))
!  IF (.NOT.gamma_only) ALLOCATE (psic_aux(dffts%nnr))
!  !
!  ALLOCATE (ovl_matrix(nbndval,nbndval))
!  !
!  ovl_matrix(:,:)   = 0.0_DP
!  !
!  ! print test kinetic only test for gamma_only
!  !
!  call g2_kin( 1 )
!  !
!  ekin = 0.0_DP
!  DO ibnd = 1, nbndval
!     DO ig = 1, npw
!        ekin = ekin + wg(ibnd,1)*CONJG(evc(ig,ibnd)) * g2kin (ig) * evc(ig,ibnd)
!     ENDDO
!     !
!     ! Note in gamma-only case, so |k+G|=0 for G=0
!     !ekin = ekin - DBLE(evc(1,ibnd)) * g2kin (1) * DBLE(evc(1,ibnd))
!  ENDDO
!  !
!  IF (gamma_only) THEN
!     !
!     ekin = 2.0_DP*ekin
!     !
!  ENDIF
!  !
!  CALL mp_sum(ekin, intra_bgrp_comm)
!  !
!  WRITE(stdout,*) 'Kinetic energy:', ekin
!  !
!  !
!  IF (l_wannier_repr) THEN
!     !
!     ALLOCATE (evc_loc(npwx,nbndval))
!     !
!     CALL bse_do_localization (1, nbndval, evc_loc, ovl_matrix, .false.)
!     !
!  ENDIF
!  !
!  ! Here is the main part of this code
!  !
!  ALLOCATE (rho_aux(dffts%nnr))
!  ALLOCATE (dvg(ngm))
!  ALLOCATE (aux1_g(ngms))
!  !
!  DO itr = 1, size(list_ovl)
!     !
!     ovl_thr = list_ovl(itr)
!     !
!     num_intergration = 0
!     fock_energy = 0.0_DP
!     !
!     DO ibnd = 1, nbndval
!        !
!        DO jbnd = 1, nbndval
!           !
!           ovl_value = ovl_matrix(ibnd,jbnd)
!           !
!           IF (ovl_value .GE. ovl_thr) THEN
!              !
!              IF (gamma_only) THEN
!                 !
!                 IF (l_wannier_repr) THEN
!                    !
!                    CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(1,ibnd),evc_loc(1,jbnd), psic,'Wave')
!                    !
!                 ELSE
!                    !
!                    CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,jbnd), psic,'Wave')
!                    !
!                 ENDIF
!                 !
!                 rho_aux(:) = DBLE(psic(:)) * AIMAG(psic(:))
!                 !
!              ELSE
!                 !
!                 IF (l_wannier_repr) THEN
!                    !
!                    CALL single_invfft_k(dffts,npw,npwx,evc_loc(1,ibnd),psic,'Wave',igk_k(1,1)) !only 1 kpoint
!                    CALL single_invfft_k(dffts,npw,npwx,evc_loc(1,jbnd),psic_aux,'Wave',igk_k(1,1)) !only 1 kpoint
!                    !
!                 ELSE
!                    !
!                    CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,1)) !only 1 kpoint
!                    CALL single_invfft_k(dffts,npw,npwx,evc(1,jbnd),psic_aux,'Wave',igk_k(1,1)) !only 1 kpoint
!                    !
!                 ENDIF
!                 !
!                 rho_aux(:) = DBLE(CONJG(psic(:)) * psic_aux(:))/omega
!                 !
!              ENDIF
!              !
!              psic(:) = CMPLX(rho_aux(:), 0.0_DP)
!              !
!              ! aux_r -> aux1_g
!              !
!              IF (gamma_only) THEN
!                 !
!                 CALL single_fwfft_gamma(dffts,ngm,ngms,psic,aux1_g,'Dense' )
!                 !
!              ELSE
!                 !
!                 CALL single_fwfft_k(dffts,ngm,ngms,psic,aux1_g,'Dense')
!                 !
!              ENDIF
!              !
!              ! vc in fock like term
!              !
!              dvg(:) = (0.0_DP, 0.0_DP)
!              DO ig = 1, ngm
!                 !
!                 dvg(ig) = aux1_g(ig) * pot3D%sqvc(ig)
!                 !
!              ENDDO
!              !
!              IF (gamma_only) THEN
!                 !
!                 etmp0 = - 2.0 * DDOT( 2*ngm, dvg(1), 1, dvg(1), 1) / omega
!                 !
!              ELSE
!                 !
!                 etmp0 = - 1.0 * DDOT( 2*ngm, dvg(1), 1, dvg(1), 1) / omega
!                 !
!              ENDIF
!              !
!              IF( ibnd == jbnd .AND. gstart == 2 ) etmp0 = etmp0 - pot3D%div
!              !IF( ibnd == jbnd .AND. gstart == 2 ) etmp0 = etmp0 - isz
!              !
!              CALL mp_sum(etmp0, intra_bgrp_comm)
!              !
!              WRITE(stdout, *) " <n_ij|n_ij> :", ibnd, jbnd, etmp0
!              !
!              num_intergration = num_intergration + 1
!              !
!              fock_energy = fock_energy + etmp0
!              !
!           ENDIF
!           !
!        ENDDO
!        !
!     ENDDO
!     !
!     WRITE(stdout, *) " Overlap threshold : ", ovl_thr
!     WRITE(stdout, *) " The number of integrals : ", num_intergration
!     WRITE(stdout, *) " Fock energy =", fock_energy, " Ry"
!     !
!  ENDDO
!  !
!  DEALLOCATE(rho_aux, dvg, aux1_g)
!  !
!  DEALLOCATE (ovl_matrix)
!  !
!  IF (ALLOCATED(psic))    DEALLOCATE(psic)
!  IF (ALLOCATED(psic_aux))DEALLOCATE(psic_aux)
!  IF (ALLOCATED(evc_loc)) DEALLOCATE(evc_loc)
!  !
!  RETURN
!  !
!ENDSUBROUTINE
