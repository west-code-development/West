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
  USE uspp,                 ONLY : okvan
  USE westcom,              ONLY : l_macropol
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: n_ipol = 3
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
  IF(l_macropol) THEN
     CALL compute_d0psi_dfpt()
  ELSE
     CALL compute_d0psi_rs()
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
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE buffers,              ONLY : get_buffer
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart
  USE noncollin_module,     ONLY : npol
  USE io_push,              ONLY : io_push_title
  USE pwcom,                ONLY : isk,igk_k,ngk,lsda,npw,npwx,nks
  USE westcom,              ONLY : nbndval0x,nbnd_occ,iuwfc,lrwfc,d0psi
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER :: i, j, k, ip, ir, ir_end, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  INTEGER :: ibnd, nbndval
  INTEGER :: iks, current_k, current_spin
  INTEGER :: dffts_nnr
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
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, number pw
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks > 1) THEN
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
     IF(gamma_only) THEN
        !
        DO ibnd = 1, nbndval, 2
           !
           IF(ibnd < nbndval) THEN
              !
              ! double bands @ gamma
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),evc_work(:,ibnd+1),psic,'Wave')
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
                 !$acc host_data use_device(d0psi)
                 CALL double_fwfft_gamma(dffts,npw,npwx,psic,d0psi(:,ibnd,iks,ip),d0psi(:,ibnd+1,iks,ip),'Wave')
                 !$acc end host_data
                 !
              ENDDO
              !
           ELSE
              !
              ! single band @ gamma
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
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
                 !$acc host_data use_device(d0psi)
                 CALL single_fwfft_gamma(dffts,npw,npwx,psic,d0psi(:,ibnd,iks,ip),'Wave')
                 !$acc end host_data
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
        DO ibnd = 1, nbndval
           !
           CALL single_invfft_k(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave',igk_k(:,current_k))
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
              !$acc host_data use_device(d0psi)
              CALL single_fwfft_k(dffts,npw,npwx,psic,d0psi(:,ibnd,iks,ip),'Wave',igk_k(:,current_k))
              !$acc end host_data
              !
           ENDDO
           !
        ENDDO
        !
        IF(npol == 2) THEN
           !
           DO ibnd = 1, nbndval
              !
              CALL single_invfft_k(dffts,npw,npwx,evc_work(npwx+1:npwx*2,ibnd),psic,'Wave',igk_k(:,current_k))
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
                 !$acc host_data use_device(d0psi)
                 CALL single_fwfft_k(dffts,npw,npwx,psic,d0psi(npwx+1:npwx*2,ibnd,iks,ip),'Wave',igk_k(:,current_k))
                 !$acc end host_data
                 !
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
     ENDIF
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbndval)
#endif
     !
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,d0psi(:,:,iks,ip),(1._DP,0._DP))
     ENDDO
     !
  ENDDO
  !
  IF(gstart == 2 .AND. gamma_only) THEN
     !$acc parallel loop collapse(3) present(d0psi)
     DO ip = 1, n_ipol
        DO iks = 1, nks
           DO ibnd = 1, nbndval0x
              d0psi(1,ibnd,iks,ip) = CMPLX(REAL(d0psi(1,ibnd,iks,ip),KIND=DP),KIND=DP)
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
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE uspp_init,            ONLY : init_us_2
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart
  USE noncollin_module,     ONLY : npol
  USE io_push,              ONLY : io_push_title
  USE pwcom,                ONLY : npw,npwx,nks,current_spin,isk,xk,lsda,igk_k,current_k,ngk
  USE westcom,              ONLY : nbndval0x,nbnd_occ,iuwfc,lrwfc,tr2_dfpt,n_dfpt_maxiter,&
                                 & l_kinetic_only,d0psi,l_lanczos,l_skip_nl_part_of_hcomr
  USE distribution_center,  ONLY : aband
#if defined(__CUDA)
  USE uspp,                 ONLY : vkb,nkb,deeq,deeq_d,qq_at,qq_at_d
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE wvfct_gpum,           ONLY : using_et,using_et_d,et=>et_d
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,             ONLY : allocate_macropol_gpu,deallocate_macropol_gpu,reallocate_ps_gpu
#else
  USE uspp,                 ONLY : vkb,nkb
  USE wavefunctions,        ONLY : evc_work=>evc
  USE wvfct,                ONLY : et
#endif
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER :: il1, iv, ip, ibnd, iks, ie, ierr
  INTEGER :: nbndval
  INTEGER, PARAMETER :: n_ipol = 3
  REAL(DP), ALLOCATABLE :: eprec(:), e(:)
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), phi_tmp(:,:)
  !$acc declare device_resident(eprec,e,phi,phi_tmp)
  !
  CALL io_push_title('Calculation of the dipole using DFPT method')
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     !
#if defined(__CUDA)
     CALL g2_kin_gpu(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks > 1) THEN
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
     !
     ! ... Sync GPU
     !
     CALL using_becp_auto(2)
     CALL using_becp_d_auto(0)
     CALL using_et(2)
     CALL using_et_d(0)
     !
     deeq_d(:,:,:,:) = deeq
     qq_at_d(:,:,:) = qq_at
     !
     CALL allocate_macropol_gpu()
     CALL reallocate_ps_gpu(nbndval,3)
#endif
     !
     ALLOCATE(eprec(3))
     ALLOCATE(e(3))
     ALLOCATE(phi_tmp(npwx*npol,3))
     ALLOCATE(phi(npwx*npol,3))
     !
     ! LOOP over band states
     !
     !$acc kernels present(d0psi)
     d0psi(:,:,iks,:) = (0._DP,0._DP)
     !$acc end kernels
     !
     DO il1 = 1, aband%nloc
        !
        iv = aband%l2g(il1)
        !
#if defined(__CUDA)
        !$acc host_data use_device(phi_tmp)
        CALL commut_Hx_psi_gpu(iks,1,1,evc_work(1,iv),phi_tmp(1,1),l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi_gpu(iks,1,2,evc_work(1,iv),phi_tmp(1,2),l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi_gpu(iks,1,3,evc_work(1,iv),phi_tmp(1,3),l_skip_nl_part_of_hcomr)
        !$acc end host_data
#else
        CALL commut_Hx_psi(iks,1,1,evc_work(1,iv),phi_tmp(1,1),l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi(iks,1,2,evc_work(1,iv),phi_tmp(1,2),l_skip_nl_part_of_hcomr)
        CALL commut_Hx_psi(iks,1,3,evc_work(1,iv),phi_tmp(1,3),l_skip_nl_part_of_hcomr)
#endif
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,3,phi_tmp,(1._DP,0._DP))
        !
        CALL set_eprec(1,evc_work(:,iv),eprec(1))
        !
        !$acc parallel loop present(eprec,e)
        DO ie = 1,3
           eprec(ie) = eprec(1)
           e(ie) = et(iv,iks)
        ENDDO
        !$acc end parallel
        !
        CALL precondition_m_wfcts(3,phi_tmp,phi,eprec)
        !
        tr2_dfpt = 1.E-12_DP
        n_dfpt_maxiter = 250
        l_kinetic_only = .FALSE.
        !
#if defined(__CUDA)
        CALL linsolve_sternheimer_m_wfcts_gpu(nbndval,3,phi_tmp,phi,e,eprec,tr2_dfpt,ierr)
#else
        CALL linsolve_sternheimer_m_wfcts(nbndval,3,phi_tmp,phi,e,eprec,tr2_dfpt,ierr)
#endif
        !
        IF(ierr /= 0) WRITE(stdout,'(7X,"** WARNING : MACROPOL not converged, ierr = ",i8)') ierr
        !
        d0psi(:,iv,iks,:) = phi
        !
     ENDDO
     !
     DEALLOCATE(eprec)
     DEALLOCATE(e)
     DEALLOCATE(phi_tmp)
     DEALLOCATE(phi)
     !
#if defined(__CUDA)
     CALL deallocate_macropol_gpu()
#endif
     !
     IF(l_lanczos) THEN
        !$acc update host(d0psi)
        CALL mp_sum(d0psi(:,:,iks,:),inter_image_comm)
        !$acc update device(d0psi)
     ENDIF
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbndval)
#endif
     !
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,d0psi(:,:,iks,ip),(1._DP,0._DP))
     ENDDO
     !
  ENDDO
  !
  IF(gstart == 2 .AND. gamma_only) THEN
     !$acc parallel loop collapse(3) present(d0psi)
     DO ip = 1, n_ipol
        DO iks = 1, nks
           DO ibnd = 1, nbndval0x
              d0psi(1,ibnd,iks,ip) = CMPLX(REAL(d0psi(1,ibnd,iks,ip),KIND=DP),KIND=DP)
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
  ENDIF
  !
END SUBROUTINE
