!
! Copyright (C) 2015-2024 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Ngoc Linh Nguyen, Victor Yu
!
!-------------------------------------------------------------------------
SUBROUTINE solve_e_psi()
  !-----------------------------------------------------------------------
  !
  USE uspp,                 ONLY : okvan
  USE westcom,              ONLY : l_dipole_realspace
  !
  IMPLICIT NONE
  !
  IF(okvan) CALL errore('solve_e_psi','Real space dipole + USPP not supported',1)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('solve_e_psi')
#else
  CALL start_clock('solve_e_psi')
#endif
  !
  ! Compute dipole in the R space. This option can be used
  ! only for finite systems (e.g. molecules).
  !
  IF(l_dipole_realspace) THEN
     CALL compute_d0psi_rs()
  ELSE
     CALL compute_d0psi_dfpt()
  ENDIF
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('solve_e_psi')
#else
  CALL stop_clock('solve_e_psi')
#endif
  !
END SUBROUTINE
!
SUBROUTINE compute_d0psi_rs()
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at,alat
  USE fft_base,             ONLY : dffts
  USE mp_global,            ONLY : my_image_id,inter_image_comm,me_bgrp
  USE mp,                   ONLY : mp_bcast
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE buffers,              ONLY : get_buffer
  USE gvect,                ONLY : gstart
  USE io_push,              ONLY : io_push_title
  USE pwcom,                ONLY : isk,ngk,lsda,npw,npwx
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands,iuwfc,lrwfc,d0psi
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
  USE west_gpu,             ONLY : reallocate_ps_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER :: i, j, k, ip, ir, ir_end, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  INTEGER :: ibnd, jbnd, lbnd, nbndval, nbnd_do
  INTEGER :: iks, current_k, current_spin
  INTEGER :: dffts_nnr, kpt_pool_nloc, band_group_nloc
  INTEGER, PARAMETER :: n_ipol = 3
  REAL(DP), ALLOCATABLE :: r(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux_r(:)
  !$acc declare device_resident(aux_r)
  !
  CALL io_push_title('Calculation of the dipole in real space')
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(aux_r(dffts%nnr))
  ALLOCATE(r(dffts%nnr,n_ipol))
  !$acc enter data create(r)
  !
  r(:,:) = 0._DP
  !
  ! Calculate r
  !
  inv_nr1 = 1._DP / REAL(dffts%nr1,KIND=DP)
  inv_nr2 = 1._DP / REAL(dffts%nr2,KIND=DP)
  inv_nr3 = 1._DP / REAL(dffts%nr3,KIND=DP)
  !
  index0 = dffts%nr1x*dffts%nr2x*SUM(dffts%nr3p(1:me_bgrp))
  ir_end = MIN(dffts%nnr,dffts%nr1x*dffts%nr2x*dffts%nr3p(me_bgrp+1))
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dffts%nr1x*dffts%nr2x)
     i = i - (dffts%nr1x*dffts%nr2x)*k
     j = i / dffts%nr1x
     i = i - dffts%nr1x*j
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
  !$acc update device(r)
  !
  ! Calculate the product r * psi(r)
  !
  DO iks = 1, kpt_pool%nloc  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, number pw
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     !
     nbnd_do = 0
     DO lbnd = 1, band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     DO lbnd = 1, nbnd_do, 2
        !
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
        !
        IF(lbnd < nbnd_do) THEN
           !
           ! double bands @ gamma
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
           !
           !$acc kernels present(aux_r)
           aux_r(:) = psic
           !$acc end kernels
           !
           DO ip = 1, n_ipol
              !
              !$acc parallel loop present(aux_r,r)
              DO ir = 1, dffts_nnr
                 psic(ir) = aux_r(ir)*r(ir,ip)*alat
              ENDDO
              !$acc end parallel
              !
              CALL double_fwfft_gamma(dffts,npw,npwx,psic,d0psi(:,lbnd,iks,ip),d0psi(:,lbnd+1,iks,ip),'Wave')
              !
           ENDDO
           !
        ELSE
           !
           ! single band @ gamma
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
           !
           !$acc kernels present(aux_r)
           aux_r(:) = psic
           !$acc end kernels
           !
           DO ip = 1, n_ipol
              !
              !$acc parallel loop present(aux_r,r)
              DO ir = 1, dffts_nnr
                 psic(ir) = CMPLX(REAL(aux_r(ir),KIND=DP)*r(ir,ip)*alat,KIND=DP)
              ENDDO
              !$acc end parallel
              !
              CALL single_fwfft_gamma(dffts,npw,npwx,psic,d0psi(:,lbnd,iks,ip),'Wave')
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     !
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,d0psi(:,:,iks,ip),(1._DP,0._DP))
     ENDDO
     !
  ENDDO
  !
  IF(gstart == 2) THEN
     kpt_pool_nloc = kpt_pool%nloc
     band_group_nloc = band_group%nloc
     !
     !$acc parallel loop collapse(3) present(d0psi)
     DO ip = 1, n_ipol
        DO iks = 1, kpt_pool_nloc
           DO lbnd = 1, band_group_nloc
              d0psi(1,lbnd,iks,ip) = CMPLX(REAL(d0psi(1,lbnd,iks,ip),KIND=DP),KIND=DP)
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
  ENDIF
  !
  DEALLOCATE(aux_r)
  !$acc exit data delete(r)
  DEALLOCATE(r)
  !
END SUBROUTINE
!
SUBROUTINE shift_d0psi(r,n_ipol)
  !
  ! Shift a position operator r to the center of the molecule
  ! for a proper calculation of d0psi in R-space.
  !
  USE fft_base,             ONLY : dffts
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps6
  USE ions_base,            ONLY : nat,tau
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : at
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n_ipol
  REAL(DP), INTENT(INOUT) :: r(dffts%nnr,n_ipol)
  !
  ! local vars
  !
  REAL(DP) :: mmin(3), mmax(3), center(3), origin(3), check_cell
  INTEGER :: ip, iatm, ir, ip1, ip2
  !
  WRITE(stdout,'(5X,"Dipole is shifted to the center of cell for the calculation of d0psi")')
  !
  check_cell = 0._DP
  !
  DO ip1 = 1,3
     DO ip2 = 1,3
        IF(ip1 /= ip2) check_cell = check_cell + at(ip1,ip2)**2
     ENDDO
  ENDDO
  !
  IF(check_cell > eps6) CALL errore('shift_d0psi','This type of the supercell is not supported',1)
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
  DO ir = 1, dffts%nnr
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
  USE mp,                   ONLY : mp_bcast
  USE buffers,              ONLY : get_buffer
  USE uspp_init,            ONLY : init_us_2
  USE gvect,                ONLY : gstart
  USE noncollin_module,     ONLY : npol
  USE io_push,              ONLY : io_push_title
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,lsda,igk_k,current_k,ngk
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands,iuwfc,lrwfc,tr2_dfpt,n_dfpt_maxiter,&
                                 & l_kinetic_only,d0psi,l_skip_nl_part_of_hcomr
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE cell_base,            ONLY : bg
  USE uspp,                 ONLY : vkb,nkb
  USE wavefunctions,        ONLY : evc
  USE wvfct,                ONLY : et
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_macropol_gpu,deallocate_macropol_gpu,reallocate_ps_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER :: iv, ip, ibnd, lbnd, iks, ig, ierr
  INTEGER :: nbndval, nbnd_do
  INTEGER :: kpt_pool_nloc, band_group_nloc
  INTEGER, PARAMETER :: n_ipol = 3
  REAL(DP), ALLOCATABLE :: eprec(:), e(:)
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), phi_tmp(:,:)
  !$acc declare device_resident(eprec,e,phi)
  !
  CALL io_push_title('Calculation of the dipole using DFPT method')
  !
  DO iks = 1, kpt_pool%nloc  ! KPOINT-SPIN
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     !
     nbnd_do = 0
     DO lbnd = 1, band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
#if defined(__CUDA)
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
#if defined(__CUDA)
     CALL allocate_macropol_gpu(1)
     CALL reallocate_ps_gpu(nbndval,n_ipol)
#endif
     !
     ALLOCATE(eprec(n_ipol))
     ALLOCATE(e(n_ipol))
     ALLOCATE(phi_tmp(npwx*npol,n_ipol))
     ALLOCATE(phi(npwx*npol,n_ipol))
     !$acc enter data create(phi_tmp)
     !
     ! LOOP over band states
     !
     d0psi(:,:,iks,:) = (0._DP,0._DP)
     !
     DO lbnd = 1, band_group%nloc
        !
        iv = band_group%l2g(lbnd)+n_trunc_bands
        !
        CALL commut_Hx_psi(iks,1,1,evc(1,iv),phi_tmp(1,1),l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi(iks,1,2,evc(1,iv),phi_tmp(1,2),l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi(iks,1,3,evc(1,iv),phi_tmp(1,3),l_skip_nl_part_of_hcomr)
        !
        !$acc parallel loop collapse(2) present(phi,phi_tmp,bg)
        DO ip = 1, n_ipol
           DO ig = 1, npwx*npol
              phi(ig,ip) = phi_tmp(ig,1)*bg(ip,1)+phi_tmp(ig,2)*bg(ip,2)+phi_tmp(ig,3)*bg(ip,3)
           ENDDO
        ENDDO
        !$acc end parallel
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,n_ipol,phi,(1._DP,0._DP))
        !
        CALL set_eprec(1,evc(:,iv),eprec(1))
        !
        !$acc parallel loop present(eprec,e,et)
        DO ip = 1, n_ipol
           eprec(ip) = eprec(1)
           e(ip) = et(iv,iks)
        ENDDO
        !$acc end parallel
        !
        CALL precondition_m_wfcts(n_ipol,phi,phi_tmp,eprec)
        !
        tr2_dfpt = 1.E-12_DP
        n_dfpt_maxiter = 250
        l_kinetic_only = .FALSE.
        !
#if defined(__CUDA)
        CALL linsolve_sternheimer_m_wfcts_gpu(nbndval,n_ipol,phi,phi_tmp,e,eprec,tr2_dfpt,ierr)
#else
        CALL linsolve_sternheimer_m_wfcts(nbndval,n_ipol,phi,phi_tmp,e,eprec,tr2_dfpt,ierr)
#endif
        !
        IF(ierr /= 0) WRITE(stdout,'(7X,"** WARNING : MACROPOL not converged, ierr = ",i8)') ierr
        !
        !$acc update host(phi_tmp)
        !
        d0psi(:,lbnd,iks,:) = phi_tmp
        !
     ENDDO
     !
     !$acc update device(d0psi)
     !
     DEALLOCATE(eprec)
     DEALLOCATE(e)
     !$acc exit data delete(phi_tmp)
     DEALLOCATE(phi_tmp)
     DEALLOCATE(phi)
     !
#if defined(__CUDA)
     CALL deallocate_macropol_gpu()
#endif
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     !
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,d0psi(:,:,iks,ip),(1._DP,0._DP))
     ENDDO
     !
  ENDDO
  !
  IF(gstart == 2) THEN
     kpt_pool_nloc = kpt_pool%nloc
     band_group_nloc = band_group%nloc
     !
     !$acc parallel loop collapse(3) present(d0psi)
     DO ip = 1, n_ipol
        DO iks = 1, kpt_pool_nloc
           DO lbnd = 1, band_group_nloc
              d0psi(1,lbnd,iks,ip) = CMPLX(REAL(d0psi(1,lbnd,iks,ip),KIND=DP),KIND=DP)
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
  ENDIF
  !
END SUBROUTINE
