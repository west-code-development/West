!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE solve_e_psi()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE mp,                   ONLY : mp_barrier
  USE uspp,                 ONLY : okvan
  USE gvect,                ONLY : gstart
  USE control_flags,        ONLY : gamma_only
  USE westcom,              ONLY : d0psi,l_macropol
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: n_ipol = 3
  !
  IF(okvan) CALL errore('solve_e_psi','Real space dipole + USPP not supported',1)
  !
  CALL start_clock('solve_e_psi')
  !
  ! Compute dipole in the R space. This option can be used
  ! only for finite systems (e.g. molecules).
  !
  IF(l_macropol) THEN
     CALL compute_d0psi_dfpt()
  ELSE
     CALL compute_d0psi_rs()
  ENDIF
  !
  IF(gstart == 2 .AND. gamma_only) &
  & d0psi(1,:,:,:) = CMPLX(REAL(d0psi(1,:,:,:),KIND=DP),KIND=DP)
  !
  CALL stop_clock('solve_e_psi')
  !
END SUBROUTINE
!
SUBROUTINE compute_d0psi_rs()
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at,alat
  USE fft_base,             ONLY : dfftp,dffts
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : my_image_id,inter_image_comm,me_bgrp
  USE mp,                   ONLY : mp_sum,mp_barrier,mp_bcast
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE buffers,              ONLY : get_buffer
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc,psic
  USE noncollin_module,     ONLY : npol
  USE pwcom,                ONLY : isk,igk_k,ngk,lsda,npw,npwx,nks
  USE westcom,              ONLY : nbnd_occ,iuwfc,lrwfc,d0psi
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER :: i, j, k, ip, ir, ir_end, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  INTEGER :: ibnd, nbndval
  INTEGER :: iks, current_k, current_spin
  INTEGER, PARAMETER :: n_ipol = 3
  REAL(DP), ALLOCATABLE :: r(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux_r(:)
  !
  WRITE(stdout,'(/,5X,"Calculation of the dipole in real space")')
  !
  ALLOCATE(aux_r(dfftp%nnr))
  ALLOCATE(r(dfftp%nnr,n_ipol))
  !
  r(:,:) = 0._DP
  !
  ! Calculate r
  !
  inv_nr1 = 1._DP / REAL(dfftp%nr1,KIND=DP)
  inv_nr2 = 1._DP / REAL(dfftp%nr2,KIND=DP)
  inv_nr3 = 1._DP / REAL(dfftp%nr3,KIND=DP)
  !
  index0 = dfftp%nr1x*dfftp%nr2x*SUM(dfftp%nr3p(1:me_bgrp))
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%nr3p(me_bgrp+1))
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     !
     DO ip = 1, n_ipol
        r(ir,ip) = REAL(i,KIND=DP)*inv_nr1*at(ip,1) &
               & + REAL(j,KIND=DP)*inv_nr2*at(ip,2) &
               & + REAL(k,KIND=DP)*inv_nr3*at(ip,3)
     ENDDO
     !
  ENDDO
  !
  CALL shift_d0psi(r,n_ipol)
  !
  ! Calculate the product r * psi(r)
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, number pw
     !
     current_k = iks
     !
     IF(lsda) current_spin = isk(iks)
     !
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     !
     IF(gamma_only) THEN
        !
        DO ibnd = 1,nbndval,2
           !
           IF(ibnd < nbndval) THEN
              !
              ! double bands @ gamma
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,ibnd+1),psic,'Wave')
              !
              aux_r(:) = psic
              !
              DO ip = 1, n_ipol
                 !
                 psic(:) = aux_r(:) * r(:,ip) * alat
                 !
                 CALL double_fwfft_gamma(dffts,npw,npwx,psic,d0psi(:,ibnd,iks,ip),d0psi(:,ibnd+1,iks,ip),'Wave')
                 !
              ENDDO
              !
           ELSE
              !
              ! single band @ gamma
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
              !
              aux_r(:) = psic
              !
              DO ip = 1, n_ipol
                 !
                 psic(:) = CMPLX(REAL(aux_r(:),KIND=DP) * r(:,ip) * alat, KIND=DP)
                 !
                 CALL single_fwfft_gamma(dffts,npw,npwx,psic,d0psi(:,ibnd,iks,ip),'Wave')
                 !
              ENDDO
              !
           ENDIF
           !
        ENDDO
        !
     ELSE
        !
        ! only single bands
        !
        DO ibnd = 1,nbndval
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,current_k))
           !
           aux_r(:) = psic
           !
           DO ip = 1, n_ipol
              !
              psic(:) = aux_r(:) * r(:,ip) * alat
              !
              CALL single_fwfft_k(dffts,npw,npwx,psic,d0psi(:,ibnd,iks,ip),'Wave',igk_k(:,current_k))
              !
           ENDDO
           !
        ENDDO
        !
        IF(npol == 2) THEN
           !
           DO ibnd = 1,nbndval
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1:npwx*2,ibnd),psic,'Wave',igk_k(:,current_k))
              !
              aux_r(:) = psic
              !
              DO ip = 1, n_ipol
                 !
                 psic(:) = aux_r(:) * r(:,ip) * alat
                 !
                 CALL single_fwfft_k(dffts,npw,npwx,psic,d0psi(npwx+1:npwx*2,ibnd,iks,ip),'Wave',igk_k(:,current_k))
                 !
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
     ENDIF
     !
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,d0psi(:,:,iks,ip),(1._DP,0._DP))
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(aux_r)
  DEALLOCATE(r)
  !
END SUBROUTINE
!
SUBROUTINE shift_d0psi(r, n_ipol)
  !
  ! Shift a position operator r to the center of the molecule
  ! for a proper calculation of d0psi in R-space.
  !
  USE fft_base,         ONLY : dfftp
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat,tau
  USE io_global,        ONLY : stdout
  USE cell_base,        ONLY : at
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n_ipol
  REAL(DP), INTENT(INOUT) :: r(dfftp%nnr,n_ipol)
  !
  ! local vars
  !
  REAL(DP) :: mmin(3), mmax(3), center(3), origin(3), check_cell
  INTEGER :: ip, iatm, ir, ip1, ip2
  !
  WRITE(stdout,'(/,5X,"Dipole is shifted to the center of cell for the calculation of d0psi")')
  !
  check_cell = 0._DP
  !
  DO ip1 = 1,3
     DO ip2 = 1,3
        IF(ip1 /= ip2) check_cell = check_cell + at(ip1,ip2)**2
     ENDDO
  ENDDO
  !
  IF(check_cell > 1.E-5_DP) CALL errore('shift_d0psi','This type of the supercell is not supported',1)
  !
  mmin(:) = 2000._DP
  mmax(:) = -2000._DP
  !
  DO ip = 1, n_ipol
     DO iatm = 1, nat
        mmin(ip) = MIN(mmin(ip), tau(ip,iatm))
        mmax(ip) = MAX(mmax(ip), tau(ip,iatm))
     ENDDO
  ENDDO
  !
  center(:) = 0.5_DP*(mmin(:)+mmax(:))
  !
  DO ip = 1, n_ipol
     origin(ip)= center(ip)-0.5_DP*at(ip,ip)
  ENDDO
  !
  DO ir = 1, dfftp%nnr
     !
     DO ip = 1, n_ipol
        r(ir,ip)= r(ir,ip) - origin(ip)
     ENDDO
     !
     DO ip = 1, n_ipol
        IF(r(ir,ip) < 0) r(ir,ip) = r(ir,ip)+at(ip,ip)
        IF(r(ir,ip) > at(ip,ip)) r(ir,ip) = r(ir,ip)-at(ip,ip)
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE
!
SUBROUTINE compute_d0psi_dfpt()
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_sum,mp_barrier,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : npol
  USE uspp,                 ONLY : vkb,nkb
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,lsda,igk_k,current_k,ngk
  USE westcom,              ONLY : nbnd_occ,iuwfc,lrwfc,tr2_dfpt,n_dfpt_maxiter,l_kinetic_only,&
                                 & d0psi,l_lanczos
  USE distribution_center,  ONLY : aband
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER :: il1, iv, ip
  INTEGER :: nbndval
  INTEGER :: iks, ierr
  INTEGER, PARAMETER :: n_ipol = 3
  LOGICAL, PARAMETER :: l_skip_nl_part_of_hcomr = .FALSE.
  REAL(DP), ALLOCATABLE :: eprec(:), e(:)
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), phi_tmp(:,:)
  !
  WRITE(stdout,'(/,5X,"Calculation of the dipole using DFPT method")')
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     !
     IF(lsda) current_spin = isk(iks)
     !
     npw = ngk(iks)
     !
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF(nkb > 0) CALL init_us_2(ngk(iks), igk_k(1,iks), xk(1,iks), vkb)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     !
     ! LOOP over band states
     !
     d0psi(:,:,iks,:) = (0._DP,0._DP)
     !
     DO il1 = 1, aband%nloc
        !
        iv = aband%l2g(il1)
        !
        ALLOCATE(phi(npwx*npol,3))
        ALLOCATE(phi_tmp(npwx*npol,3))
        !
        CALL commut_Hx_psi(iks, 1, 1, evc(1,iv), phi_tmp(1,1), l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi(iks, 1, 2, evc(1,iv), phi_tmp(1,2), l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi(iks, 1, 3, evc(1,iv), phi_tmp(1,3), l_skip_nl_part_of_hcomr)
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,3,phi_tmp,(1._DP,0._DP))
        !
        ALLOCATE(eprec(3))
        ALLOCATE(e(3))
        !
        CALL set_eprec(1, evc(1,iv), eprec(1))
        eprec(2) = eprec(1)
        eprec(3) = eprec(1)
        !
        e(1) = et(iv,iks)
        e(2) = et(iv,iks)
        e(3) = et(iv,iks)
        !
        CALL precondition_m_wfcts(3, phi_tmp, phi, eprec)
        !
        tr2_dfpt = 1.E-12_DP
        n_dfpt_maxiter = 250
        l_kinetic_only = .FALSE.
        !
        CALL linsolve_sternheimer_m_wfcts(nbndval, 3, phi_tmp, phi, e, eprec, tr2_dfpt, ierr)
        !
        IF(ierr /= 0) WRITE(stdout,'(7X,"** WARNING : MACROPOL not converged, ierr = ",i8)') ierr
        !
        d0psi(:,iv,iks,:) = phi
        !
        DEALLOCATE(eprec)
        DEALLOCATE(e)
        DEALLOCATE(phi_tmp)
        DEALLOCATE(phi)
        !
     ENDDO
     !
     IF(l_lanczos) THEN
        CALL mp_sum(d0psi(:,:,iks,:),inter_image_comm)
     ENDIF
     !
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,d0psi(1,1,iks,ip),(1._DP,0._DP))
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE
