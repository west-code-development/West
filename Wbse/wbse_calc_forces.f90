!
! Copyright (C) 2015 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_forces( dvg_exc_tmp )
!-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx
  USE klist,                ONLY : nks
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : nbndval0x,n_trunc_bands
  USE distribution_center,  ONLY : band_group,kpt_pool
  !
  IMPLICIT NONE
  !
  ! !/O
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp( npwx*npol, band_group%nlocx, kpt_pool%nloc )
  !
  ! Workspace
  !
  INTEGER :: n, na, ipol
  REAL(DP),ALLOCATABLE :: forces(:), dvgdvg_mat(:,:,:)
  REAL(DP) :: sumforces, time_spent(2)
  REAL(DP), EXTERNAL :: get_clock
  COMPLEX(DP),ALLOCATABLE :: z_rhs_vec(:,:,:), zvector(:,:,:), drhox1(:,:), drhox2(:,:)
  LOGICAL :: do_zvector, poor_of_ram_drhox2
  !
  CALL start_clock( 'l_forces' )
  time_spent(1) = get_clock( 'l_forces' )
  !
  CALL io_push_bar()
  CALL io_push_title("Calculating forces")
  CALL io_push_bar()
  !
  n = 3 * nat
  !
  ALLOCATE ( forces(n) )
  forces = 0._DP
  ALLOCATE( dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks) )
  dvgdvg_mat = 0._DP
  !
  ALLOCATE( drhox1(dffts%nnr, nspin) )
  drhox1(:,:) = (0._DP, 0._DP)
  !
  !
  !
  CALL io_push_title("Compute drhox1")
  !
  CALL wbse_calc_drhox1( dvg_exc_tmp, drhox1 )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title("Compute forces of drhox1")
  !
  CALL wbse_forces_drhox1( n, dvg_exc_tmp, drhox1, forces )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title("Compute < dvg | dvg >")
  !
  CALL wbse_calc_dvgdvg_mat( dvg_exc_tmp, dvgdvg_mat )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title("Compute drhox2")
  !
  ALLOCATE( drhox2(dffts%nnr, nspin) )
  drhox2(:,:) = (0._DP, 0._DP)
  !
  poor_of_ram_drhox2 = .FALSE.
  !
  IF(poor_of_ram_drhox2) THEN
     !
     CALL wbse_calc_drhox2_slow( dvgdvg_mat, drhox2 )
     !
  ELSE
     !
     CALL wbse_calc_drhox2( dvgdvg_mat, drhox2 )
     !
  ENDIF
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title("Compute forces of drhox2")
  !
  CALL wbse_forces_drhox2( n, dvg_exc_tmp, dvgdvg_mat, drhox2, forces )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title("Build and solve the Z-vector equations")
  !
  ! Zvector part
  do_zvector = .FALSE.
  !
  IF( do_zvector ) THEN
     !
     ALLOCATE( z_rhs_vec( npwx, band_group%nlocx, kpt_pool%nloc ) )
     z_rhs_vec = ( 0.0_DP, 0.0_DP )
     !
!     CALL wbse_build_rhs_zvector_eq( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
     !
     time_spent(2) = get_clock( 'l_forces' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'l_forces' )
     !
     !
     !
     ALLOCATE( zvector( npwx, band_group%nlocx, kpt_pool%nloc ) )
     zvector = ( 0.0_DP, 0.0_DP )
     !
!     CALL wbse_solve_zvector_eq_cg( z_rhs_vec, zvector )
     !
     time_spent(2) = get_clock( 'l_forces' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'l_forces' )
     !
     !
     !
     CALL io_push_title("Compute forces of Z vector")
     !
     CALL wbse_forces_drhoz( n, zvector, forces )
     !
     time_spent(2) = get_clock( 'l_forces' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'l_forces' )
     !
     DEALLOCATE( z_rhs_vec )
     DEALLOCATE( zvector )
     !
  ENDIF
  !
  CALL io_push_title("Forces Total")
  !
  DO na = 1, nat
     !
     ! Note: forces are the negative of the gradients
     WRITE( stdout, 9035) na, ityp(na), ( - forces( 3 * na - 3 + ipol ), ipol = 1, 3 )
     !
  ENDDO
  !
  ! enforce total forces to be 0 in each direction
  !
  DO ipol = 1, 3
     !
     sumforces = 0._DP
     !
     DO na = 1, nat
        !
        sumforces = sumforces + forces( 3 * na - 3 + ipol )
        !
     ENDDO
     !
     DO na = 1, nat
        !
        forces( 3 * na - 3 + ipol ) = forces( 3 * na - 3 + ipol ) - sumforces / DBLE(nat)
        !
     ENDDO
     !
  ENDDO
  !
  CALL io_push_title("Forces Corrected")
  !
  DO na = 1, nat
     !
     ! Note: forces are the negative of the gradients
     WRITE( stdout, 9035) na, ityp(na), ( - forces( 3 * na - 3 + ipol ), ipol = 1, 3 )
     !
  ENDDO
  !
  WRITE(stdout,'(7X,a)') ' '
  !
  DEALLOCATE( forces )
  DEALLOCATE( dvgdvg_mat )
  DEALLOCATE( drhox1 )
  DEALLOCATE( drhox2 )
  !
  CALL stop_clock( 'l_forces' )
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE



!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_time(time)
!-----------------------------------------------------------------------
   !
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   !
   IMPLICIT NONE
   !
   REAL(DP),INTENT(IN) :: time(2)
   !
   INTEGER :: i,j
   CHARACTER(20),EXTERNAL :: human_readable_time
   !
   WRITE(stdout, "(5x,'Tot. elapsed time ',a,',  time spent in last step ',a) ") &
   TRIM(human_readable_time(time(2))), TRIM(human_readable_time(time(2)-time(1)))
   !
END SUBROUTINE



!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox1( dvg_exc_tmp, drhox1 )
!-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : npwx,npw
  USE klist,                ONLY : nks
  USE pwcom,                ONLY : isk,igk_k,lsda,current_k,wg,ngk
  USE mp,                   ONLY : mp_sum,mp_barrier,mp_bcast
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_invfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands,nbndval0x,l_spin_flip
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : my_image_id,inter_pool_comm,inter_bgrp_comm
#if defined(__CUDA)
  USE west_gpu,             ONLY : tmp_r,tmp_c
#else
  USE wavefunctions,        ONLY : psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: drhox1(dffts%nnr, nspin)
  !
  ! Workspace
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ir, lbnd, ibnd, jbnd, dffts_nnr
  REAL(DP) :: prod, w1, w2
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: tmp_r(:)
  COMPLEX(DP), ALLOCATABLE :: tmp_c(:)
#endif
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
#if defined(__CUDA)
  CALL start_clock_gpu('calc_drhox1')
#else
  CALL start_clock('calc_drhox1')
#endif
  !
  dffts_nnr = dffts%nnr
  drhox1(:,:) = (0._DP,0._DP)
  !
#if !defined(__CUDA)
  IF(gamma_only) THEN
     ALLOCATE(tmp_r(dffts%nnr))
  ELSE
     ALLOCATE(tmp_c(dffts%nnr))
  ENDIF
#endif
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point and spin
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     IF (gamma_only) THEN
        !
        !$acc kernels present(tmp_r)
        tmp_r(:) = 0._DP
        !$acc end kernels
        !
        ! double bands @ gamma
        !
        DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
           !
           ibnd = band_group%l2g(lbnd) + n_trunc_bands
           jbnd = band_group%l2g(lbnd+1) + n_trunc_bands
           !
           w1 = wg(ibnd,iks_do)/omega
           w2 = wg(jbnd,iks_do)/omega
           !
           psic(:) = (0._DP, 0._DP)
           !$acc host_data use_device(dvg_exc_tmp)
           CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),dvg_exc_tmp(:,lbnd+1,iks),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(tmp_r)
           DO ir = 1, dffts_nnr
              !
              tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)*REAL(psic(ir),KIND=DP)
              tmp_r(ir) = tmp_r(ir) + w2*AIMAG(psic(ir))*AIMAG(psic(ir))
              !
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        ! single band @ gamma
        !
        IF(MOD(nbnd_do,2) == 1) THEN
           !
           lbnd = nbnd_do
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           !
           w1 = wg(ibnd,iks_do)/omega
           !
           psic(:) = (0._DP, 0._DP)
           !$acc host_data use_device(dvg_exc_tmp)
           CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(tmp_r)
           DO ir=1, dffts_nnr
              !
              tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)*REAL(psic(ir),KIND=DP)
              !
           ENDDO
           !$acc end parallel
           !
        ENDIF
        !
        !$acc update host(tmp_r)
        !
        drhox1(:,current_spin) = CMPLX(tmp_r,KIND=DP)
        !
     ELSE
        !
        !$acc kernels present(tmp_c)
        tmp_c(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        ! only single bands
        !
        DO lbnd = 1, nbnd_do
           !
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           !
           w1 = wg(ibnd,iks_do)/omega
           !
           !$acc host_data use_device(dvg_exc_tmp)
           CALL single_invfft_k(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave',igk_k(:,current_k))
           !$acc end host_data
           !
           !$acc parallel loop present(tmp_c)
           DO ir = 1, dffts_nnr
              tmp_c(ir) = tmp_c(ir) + w1*CONJG(psic(ir))*psic(ir)
           ENDDO
           !$acc end parallel
           IF(npol == 2) THEN
              !
              !$acc host_data use_device(dvg_exc_tmp)
              CALL single_invfft_k(dffts,npw,npwx,dvg_exc_tmp(npwx+1:npwx*2,lbnd,iks),psic,'Wave',igk_k(:,current_k))
              !$acc end host_data
              !
              !$acc parallel loop present(tmp_c)
              DO ir = 1, dffts_nnr
                 tmp_c(ir) = tmp_c(ir) + w1*CONJG(psic(ir))*psic(ir)
              ENDDO
              !$acc end parallel
              !
           ENDIF
           !
        ENDDO
        !
        !$acc update host(tmp_c)
        !
        drhox1(:,current_spin) = tmp_c
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_sum(drhox1,inter_bgrp_comm)
  CALL mp_sum(drhox1,inter_pool_comm)
  !
  !$acc update device(drhox1)
  !
#if !defined(__CUDA)
  IF(gamma_only) THEN
     DEALLOCATE(tmp_r)
  ELSE
     DEALLOCATE(tmp_c)
  ENDIF
#endif
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('calc_drhox1')
#else
  CALL stop_clock('calc_drhox1')
#endif
  !
END SUBROUTINE


!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhox1( n, dvg_exc_tmp, drhox1, forces )
!-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : npwx,npw,et
  USE klist,                ONLY : xk,nks,wk
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_k,ngk
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : nbnd_occ,nbndval0x,n_trunc_bands, l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,             ONLY : dvrs,hevc1,reallocate_ps_gpu
  USE cublas
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc), drhox1(dffts%nnr, nspin)
  REAL(DP), INTENT(INOUT) :: forces(n)
  !
  ! WORKSPACE
  !
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, ir, lbnd, ibnd
  REAL(DP),ALLOCATABLE :: forces_aux(:), forces_drhox1(:), forcelc(:,:), rdrhox1(:,:)
  REAL(DP), external :: DDOT
  REAL(DP) :: prod, w1
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title("Forces : drhox1 part")
  !
  ALLOCATE( forces_drhox1(n) )
  ALLOCATE( forces_aux(n) )
  ALLOCATE( forcelc(3, nat) )
  ALLOCATE( dvpsi( npwx, n) )
  ALLOCATE( rdrhox1(dffts%nnr, nspin) )
  !
  forces_aux = 0._DP
  forces_drhox1 = 0._DP
  rdrhox1(:,:) = 0._DP
  !
  ! nonlocal part
  !
  DO iks = 1, kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
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
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
#if defined(__CUDA)
     !
     ! ... Sync GPU
     !
     CALL using_becp_auto(2)
     CALL using_becp_d_auto(0)
#endif
     !
     DO lbnd = 1, nbnd_do
        !
        ! 1) | dvpsi_i >
        !
        dvpsi = (0._DP, 0._DP)
        !
        CALL wbse_get_dvpsi_per_state_gamma_nonlocal(n, dvg_exc_tmp(:,lbnd,iks), dvpsi)
        !
        ! 2) forces_aux_i = < dvg | dvpsi_i > 
        !
        forces_aux = 0._DP
        !
        DO ia = 1, n
           !
           forces_aux(ia) = 2._DP * DDOT(2*npw, dvg_exc_tmp(:,lbnd,iks), 1, dvpsi(:,ia), 1)
           !
           IF (gstart == 2) forces_aux(ia) = forces_aux(ia) - DBLE(dvg_exc_tmp(1,lbnd,iks))*DBLE(dvpsi(1,ia))
           !
        ENDDO
        !
        IF (nspin == 2) THEN
           !
           forces_drhox1(:) = forces_drhox1(:) + 1.0_DP * wk(iks) * forces_aux(:)
           !
        ELSEIF ( nspin == 1 ) THEN
           !
           forces_drhox1(:) = forces_drhox1(:) + 0.5_DP * wk(iks) * forces_aux(:)
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(forces_drhox1, intra_bgrp_comm)
  CALL mp_sum(forces_drhox1, inter_bgrp_comm)
  CALL mp_sum(forces_drhox1, inter_pool_comm)
  !
  ! local part
  !
  rdrhox1(:,:) = DBLE(drhox1(:,:))
  !
  IF (nspin == 2) THEN
     !
     rdrhox1(:,1) = rdrhox1(:,1) + rdrhox1(:,2)
     !
  ENDIF
  !
  CALL force_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, g, rdrhox1(:,1), dffts%nl, gstart, gamma_only, vloc, forcelc)
  !
  IF (nspin == 2) THEN
     !
     forcelc(:,:) = -1.0_DP * forcelc(:,:)
     !
  ELSEIF (nspin == 1) THEN
     !
     forcelc(:,:) = -0.5_DP * forcelc(:,:)
     !
  ENDIF
  !
  DO na = 1, nat
     !
     DO ipol = 1, 3
        !
        forces_drhox1( 3 * na - 3 + ipol ) = forces_drhox1(3 * na - 3 + ipol)&
                                           + forcelc(ipol, na)
        !
     ENDDO
     !
  ENDDO
  !
  forces(:) = forces(:) + forces_drhox1(:)
  !
  DO na = 1, nat
     !
     ! Note: forces are the negative of gradients
     WRITE( stdout, 9035) na, ityp(na), (- forces_drhox1(3 * na - 3 + ipol), ipol = 1, 3)
     !
  ENDDO
  !
  WRITE(stdout,'(7X,a)') ' '
  !
  DEALLOCATE(dvpsi)
  DEALLOCATE(forces_aux)
  DEALLOCATE(forces_drhox1)
  DEALLOCATE(forcelc)
  DEALLOCATE(rdrhox1)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE



!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_dvgdvg_mat( dvg_exc_tmp, dvgdvg_mat )
  !-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npwx,npw
  USE gvect,                ONLY : gstart
  USE pwcom,                ONLY : ngk
  USE klist,                ONLY : nks
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(INOUT) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks)
  !
  ! WORKSPACE
  !
  INTEGER :: iks, iks_do, iv, iv2, nbndval, nbnd_do, lbnd, ibnd, ibndp, lbnd2, jbnd, jbndp
  COMPLEX(DP), ALLOCATABLE :: dvg_exc_tmp_copy(:,:)
  REAL(DP), external :: DDOT
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  ALLOCATE( dvg_exc_tmp_copy(npwx*npol, nbndval0x-n_trunc_bands) )
  !
  DO iks = 1, kpt_pool%nloc
     !
     npw = ngk(iks)
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     dvg_exc_tmp_copy(:,:) = (0._DP, 0._DP)
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd)
        !
        dvg_exc_tmp_copy(:,ibnd) = dvg_exc_tmp(:,lbnd,iks)
        !
     ENDDO
     !
     CALL mp_sum(dvg_exc_tmp_copy(:,:), inter_bgrp_comm)
     !
     DO ibnd = 1, nbndval - n_trunc_bands
        !
        DO lbnd = 1, nbnd_do
           !
           dvgdvg_mat(ibnd, lbnd, iks) = 2._DP * DDOT( 2*npw, dvg_exc_tmp_copy(:,ibnd), 1, dvg_exc_tmp(:,lbnd,iks), 1 )
           !
           IF ( gstart == 2 ) dvgdvg_mat(ibnd, lbnd, iks) = dvgdvg_mat(ibnd, lbnd, iks) &
                              - DBLE(dvg_exc_tmp_copy(1,ibnd)) * DBLE(dvg_exc_tmp(1,lbnd,iks))
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum( dvgdvg_mat, intra_bgrp_comm )
  CALL mp_sum( dvgdvg_mat, inter_pool_comm )
  !
  DEALLOCATE( dvg_exc_tmp_copy )
  !
END SUBROUTINE



!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox2( dvgdvg_mat, drhox2 )
  !-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : g2kin,npwx, npw
  USE klist,                ONLY : xk,nks,wk
  USE pwcom,                ONLY : isk,lsda,current_k,wg,ngk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_invfft_gamma
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,&
                                 & inter_bgrp_comm,intra_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : tmp_r,tmp_c,psic2
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks)
  COMPLEX(DP), INTENT(INOUT) :: drhox2(dffts%nnr, nspin)
  !
  ! WORKSPACE
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, ir,&
             &lbnd, ibnd, ibndp, jbnd, jbndp
  REAL(DP), external :: DDOT
  REAL(DP) :: prod, w1
  COMPLEX(DP), ALLOCATABLE :: aux_all_r(:,:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  ALLOCATE( aux_all_r(dffts%nnr, nbndval0x-n_trunc_bands) )
  !
  aux_all_r(:,:) = (0._DP, 0._DP)
  !
  drhox2(:,:) = (0._DP, 0._DP)
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     npw = ngk(iks)
     !
     ! ... read GS wavefunctions
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     ! INVFFT all evc into aux_all_r
     !
     DO ibnd = 1, nbndval - n_trunc_bands
        !
        ibndp = ibnd + n_trunc_bands
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibndp),psic,'Wave')
        !
        aux_all_r(:,ibnd) = psic(:)
        !
     ENDDO
     !
     IF (gamma_only) THEN
        !
        DO lbnd = 1, nbnd_do
           !
           ibnd = band_group%l2g(lbnd)
           !
           w1 = wg(ibnd + n_trunc_bands, iks_do)/omega
           !
           DO jbnd = 1, nbndval - n_trunc_bands
              !
              DO ir=1, dffts%nnr
                 !
                 prod = REAL(aux_all_r(ir,ibnd),KIND=DP) * REAL(aux_all_r(ir,jbnd),KIND=DP) * dvgdvg_mat(jbnd, lbnd, iks)
                 !
                 drhox2(ir,iks) = drhox2(ir,iks) - w1 * CMPLX(prod, 0._DP, KIND=DP)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_sum(drhox2(:,:), inter_bgrp_comm)
  CALL mp_sum(drhox2(:,:), inter_pool_comm)
  !
  DEALLOCATE( aux_all_r )
  !
END SUBROUTINE



!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox2_slow( dvgdvg_mat, drhox2 )
  !-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : g2kin,npwx,npw
  USE klist,                ONLY : xk,nks,wk
  USE pwcom,                ONLY : isk,lsda,current_k,wg,ngk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_invfft_gamma
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,&
                                 & inter_bgrp_comm,intra_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : tmp_r,tmp_c,psic2
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif

  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks)
  COMPLEX(DP), INTENT(INOUT) :: drhox2(dffts%nnr, nspin)
  !
  ! WORKSPACE
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, ir,&
             &lbnd, ibnd, ibndp, jbnd, jbndp
  REAL(DP), external :: DDOT
  REAL(DP) :: prod, w1
  COMPLEX(DP), ALLOCATABLE :: aux_r(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  ALLOCATE( aux_r(dffts%nnr) )
  !
  aux_r(:) = (0._DP, 0._DP)
  !
  drhox2(:,:) = (0._DP, 0._DP)
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read GS wavefunctions
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     IF (gamma_only) THEN
        !
        DO lbnd = 1, nbnd_do
           !
           ibnd = band_group%l2g(lbnd) + n_trunc_bands
           !
           w1 = wg(ibnd, iks_do)/omega
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
           !
           DO jbnd = 1, nbndval - n_trunc_bands
              !
              jbndp = jbnd + n_trunc_bands ! index for evc
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,jbndp),aux_r,'Wave')
              !
              DO ir=1, dffts%nnr
                 !
                 prod = REAL(psic(ir),KIND=DP) * REAL(aux_r(ir),KIND=DP) * dvgdvg_mat(jbnd, lbnd, iks)
                 !
                 drhox2(ir,iks) = drhox2(ir,iks) - w1 * CMPLX(prod, 0._DP, KIND=DP)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_sum(drhox2(:,:), inter_bgrp_comm)
  CALL mp_sum(drhox2(:,:), inter_pool_comm)
  !
  DEALLOCATE( aux_r )
  !
END SUBROUTINE




!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhox2( n, dvg_exc_tmp, dvgdvg_mat, drhox2, forces )
!-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : g2kin,npwx,npw
  USE klist,                ONLY : xk,nks,wk
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_k,wg,ngk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,&
                                 & inter_bgrp_comm,intra_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,             ONLY : tmp_r,tmp_c,psic2
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif

  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, nks)
  COMPLEX(DP), INTENT(IN) :: drhox2(dffts%nnr, nspin)
  REAL(DP), INTENT(INOUT) :: forces(n)
  !
  ! WORKSPACE
  !
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:), aux_g(:)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, na,&
             &lbnd, ibnd, ibndp, jbnd, jbndp
  REAL(DP),ALLOCATABLE :: forces_aux(:), forces_drhox2(:), rdrhox2(:,:), forcelc(:,:)
  REAL(DP), external :: DDOT
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title("Forces : drhox2 part")
  !
  ALLOCATE( forces_aux(n) )
  ALLOCATE( forces_drhox2(n) )
  ALLOCATE( forcelc(3, nat) )
  ALLOCATE( dvpsi( npwx, n) )
  ALLOCATE( rdrhox2( dffts%nnr, nspin) )
  ALLOCATE( aux_g(npwx*npol) )
  !
  aux_g(:) = (0._DP, 0._DP)
  forces_aux = 0._DP
  forces_drhox2 = 0._DP
  forcelc(:,:) = 0._DP
  rdrhox2(:,:) = 0._DP
  !
  ! nonlocal part
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
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
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks_do)
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
#endif
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        !
        ! 1) | dvpsi_i >
        !
        dvpsi = (0._DP, 0._DP)
        !
        CALL wbse_get_dvpsi_per_state_gamma_nonlocal(n, evc_work(:,ibnd), dvpsi)
        !
        ! 2) force_aux = < evc_iv2 | dvpsi_ia_iv > 
        !
        forces_aux = 0._DP
        !
        aux_g(:) = (0._DP, 0._DP)
        !
        DO jbnd = 1, nbndval - n_trunc_bands
           !
           jbndp = jbnd + n_trunc_bands
           !
           CALL ZAXPY(npw, CMPLX(dvgdvg_mat(jbnd, lbnd, iks), 0.0d0, DP), evc_work(:,jbndp), 1, aux_g(:), 1)
           !
        ENDDO
        !
        DO ia = 1, n
           !
           forces_aux(ia) = 2._DP * DDOT( 2*npw, aux_g(:), 1, dvpsi(:,ia), 1 )
           !
           IF ( gstart == 2 ) forces_aux(ia) = forces_aux(ia) - DBLE(aux_g(1))*DBLE(dvpsi(1,ia))
           !
        ENDDO
        !
        IF ( nspin == 2) THEN
           !
           forces_drhox2(:) = forces_drhox2(:) - 1.0_DP * wk(iks) * forces_aux(:)
           !
        ELSEIF ( nspin == 1 ) THEN
           !
           forces_drhox2(:) = forces_drhox2(:) - 0.5_DP * wk(iks) * forces_aux(:)
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum( forces_drhox2, intra_bgrp_comm )
  CALL mp_sum( forces_drhox2, inter_bgrp_comm )
  CALL mp_sum( forces_drhox2, inter_pool_comm )
  !
  ! local part
  !
  rdrhox2(:,:) = DBLE(drhox2(:,:))
  !
  IF (nspin == 2) THEN
     !
     rdrhox2(:,1) = rdrhox2(:,1) + rdrhox2(:,2)
     !
  ENDIF
  !
  CALL force_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, g, rdrhox2(:,1), dffts%nl, gstart, gamma_only, vloc, forcelc)
  !
  IF ( nspin == 2 ) THEN
     !
     forcelc(:,:) = -1.0_DP * forcelc(:,:)
     !
  ELSEIF ( nspin == 1 ) THEN
     !
     forcelc(:,:) = -0.5_DP * forcelc(:,:)
     !
  ENDIF
  !
  DO na = 1, nat
     !
     DO ipol = 1, 3
        !
        forces_drhox2( 3 * na - 3 + ipol ) = forces_drhox2( 3 * na - 3 + ipol )&
                                           + forcelc(ipol, na)
        !
     ENDDO
     !
  ENDDO
  !
  forces(:) = forces(:) + forces_drhox2(:)
  !
  DO na = 1, nat
     !
     ! Note: forces are the negative gradients
     WRITE( stdout, 9035) na, ityp(na), ( - forces_drhox2( 3 * na - 3 + ipol ), ipol = 1, 3 )
     !
  ENDDO
  !
  WRITE(stdout,'(7X,a)') ' '
  !
  DEALLOCATE( dvpsi )
  DEALLOCATE( rdrhox2 )
  DEALLOCATE( aux_g )
  !
  DEALLOCATE( forces_aux )
  DEALLOCATE( forces_drhox2 )
  DEALLOCATE( forcelc )
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE




!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhoz( n, zvector, forces )
  !-----------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE lsda_mod,             ONLY : current_spin,nspin
  USE wvfct,                ONLY : g2kin,npwx,npw
  USE klist,                ONLY : xk,nks,wk
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_k,wg,ngk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,&
                                 & inter_bgrp_comm,intra_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE west_gpu,             ONLY : tmp_r,tmp_c,psic2
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: zvector(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(INOUT) :: forces(n)
  !
  ! WORKSPACE
  !
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:), drhoz(:,:)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, lbnd, ibnd, ibndp
  REAL(DP),ALLOCATABLE :: forces_aux(:), forces_drhoz(:), forcelc(:,:), rdrhoz(:,:)
  REAL(DP), external :: DDOT
  !
  CALL io_push_title("Forces : Z-Vector Part")
  !
  ALLOCATE( forces_drhoz(n) )
  ALLOCATE( forces_aux(n) )
  ALLOCATE( forcelc(3,nat) )
  ALLOCATE( dvpsi( npwx, n) )
  ALLOCATE( drhoz(dffts%nnr, nspin) )
  ALLOCATE( rdrhoz(dffts%nnr, nspin) )
  !
  forces_aux = 0._DP
  forces_drhoz = 0._DP
  drhoz(:,:) = (0._DP, 0._DP)
  rdrhoz(:,:) = 0._DP
  !
  ! nonlocal
  !
  DO iks = 1, kpt_pool%nloc
     !
     ! Note: the Z vector part is always spin-conserving
     IF(l_spin_flip) THEN
        iks_do = iks
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
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
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks_do)
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
#endif
     !
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        !
        ! 1) | dvpsi_i >
        !
        dvpsi = (0._DP, 0._DP)
        !
        CALL wbse_get_dvpsi_per_state_gamma_nonlocal( n, evc_work(:,ibnd), dvpsi )
        !
        ! 2) force_aux_i = < z_vector | dvpsi_i > 
        !
        forces_aux = 0._DP
        !
        DO ia = 1, n
           !
           forces_aux(ia) = 2._DP * DDOT( 2*npw, zvector(:,lbnd,iks), 1, dvpsi(:,ia), 1 )
           !
           IF ( gstart == 2 ) forces_aux(ia) = forces_aux(ia) - DBLE(zvector(1,lbnd,iks))*DBLE(dvpsi(1,ia))
           !
        ENDDO
        !
        IF ( nspin == 2 ) THEN
           !
           forces_drhoz(:) = forces_drhoz(:) + 2.0_DP * wk(iks) * forces_aux(:) ! complex conjugate
           !
        ELSEIF ( nspin == 1 ) THEN
           !
           forces_drhoz(:) = forces_drhoz(:) + 1.0_DP * wk(iks) * forces_aux(:)
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum( forces_drhoz, intra_bgrp_comm )
  CALL mp_sum( forces_drhoz, inter_bgrp_comm )
  CALL mp_sum( forces_drhoz, inter_pool_comm )
  !
  ! local part
  !
  CALL wbse_calc_dens( zvector, drhoz, .FALSE. )
  !
  drhoz(:,:) = 2._DP * drhoz(:,:)
  !
  rdrhoz(:,:) = DBLE(drhoz(:,:))
  !
  IF (nspin == 2) THEN
     !
     rdrhoz(:,1) = rdrhoz(:,1) + rdrhoz(:,2)
     !
  ENDIF
  !
  CALL force_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, g, rdrhoz(:,1), dffts%nl, gstart, gamma_only, vloc, forcelc)
  !
  IF ( nspin == 2 ) THEN
     !
     forcelc(:,:) = -1.0_DP * forcelc(:,:)
     !
  ELSEIF ( nspin == 1 ) THEN
     !
     forcelc(:,:) = -0.5_DP * forcelc(:,:)
     !
  ENDIF
  !
  DO na = 1, nat
     !
     DO ipol = 1, 3
        !
        forces_drhoz( 3 * na - 3 + ipol ) = forces_drhoz( 3 * na - 3 + ipol )&
                                          + forcelc(ipol, na)
        !
     ENDDO
     !
  ENDDO
  !
  forces(:) = forces(:) + forces_drhoz(:)
  !
  DO na = 1, nat
     !
     ! Note: forces are the negative of the gradients
     WRITE( stdout, 9035) na, ityp(na), ( - forces_drhoz( 3 * na - 3 + ipol ), ipol = 1, 3 )
     !
  ENDDO
  !
  WRITE(stdout,'(7X,a)') ' '
  !
  DEALLOCATE( dvpsi )
  DEALLOCATE( drhoz )
  DEALLOCATE( rdrhoz )
  DEALLOCATE( forcelc )
  !
  DEALLOCATE( forces_aux )
  DEALLOCATE( forces_drhoz )
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE



!-----------------------------------------------------------------------
SUBROUTINE wbse_get_dvpsi_per_state_gamma_nonlocal (n, dvg_tmp, dvpsi)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector u.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  ! It implements Eq. B29 of PRB 64, 235118 (2001). The contribution
  ! of the local pseudopotential is calculated here, that of the nonlocal
  ! pseudopotential in dvqpsi_us_only.
  !
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat,ityp,ntyp => nsp
  USE cell_base,             ONLY : tpiba
  USE fft_base,              ONLY : dffts
  USE fft_interfaces,        ONLY : fwfft,invfft
  USE gvect,                 ONLY : g,gstart
  USE noncollin_module,      ONLY : npol
  use uspp_param,            ONLY : upf,nh,nhm
  USE wvfct,                 ONLY : npw,npwx
  USE uspp,                  ONLY : dvan,nkb,vkb
  USE vlocal,                ONLY : vloc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE klist,                 ONLY : xk,nks,wk
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: dvg_tmp(npwx*npol)
  COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx,n)
  !
  ! Workspace
  !
  INTEGER :: ia, mu, ikk, ig, nt, ir, is, ip, ih,ik,jkb,ic,it,glob_mu,new_ia
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! counter on G vectors
  ! the type of atom
  ! counter on bands
  ! counter on real mesh
  REAL(DP) :: gu
  real(DP), POINTER :: bec1(:,:), bec2(:,:)
  COMPLEX(DP) :: fact, gtau,exc
  COMPLEX(DP),POINTER :: work(:,:)
  !TYPE(bar_type) :: barra
  !INTEGER :: barra_load
  !LOGICAL, PARAMETER :: show_barra =.FALSE.
  !
  jkb=0
  fact = tpiba * cmplx(0._DP,-1._DP,kind=DP) 
  DO nt = 1,ntyp
     ALLOCATE (work( npwx, nh(nt)))
     ALLOCATE (bec1( nh(nt),1))
     ALLOCATE (bec2( nh(nt),1))
     DO ia = 1, nat
        IF (ityp(ia) == nt) THEN
           IF( nh(nt) > 0) THEN
              DO ic = 1,3
                 mu =3*(ia-1) +ic
                 !  first term: sum_l sum_G' [ i V_l(G) V^*_l(G') (G'*u) psi(G')
                 !  second term: sum_l sum_G' [-i (G*u) V_l(G) V^*_l(G') psi(G')
                 !
                 DO ih = 1,nh(nt)
                    DO ig = 1,npw
                       work(ig,ih) = vkb(ig,jkb+ih) * g(ic,ig)*fact
                    ENDDO
                 ENDDO
                 !
                 !CALL calbec ( npw, work, evc(1), bec1 , 1)
                 !Calbec can't be called on an arbitrary band of evc so this is the equivalent
                 CALL DGEMV( 'C', 2*npw, nh(nt), 2.0_DP, work, 2*npwx, dvg_tmp, 1, 0.0_DP, &
                        bec1, 1 )
                 IF ( gstart == 2 ) bec1(:,1) = bec1(:,1) - work(1,:)*dvg_tmp(1)
                 CALL mp_sum( bec1( :, 1 ), intra_bgrp_comm )
                 !
                 !CALL calbec ( npw, vkb(:,jkb+1:jkb+nh(nt)), evc(1), bec2, 1)
                 !
                 CALL DGEMV( 'C', 2*npw, nh(nt), 2.0_DP, vkb(:,jkb+1:jkb+nh(nt)), 2*npwx, dvg_tmp, 1, 0.0_DP, &
                        bec2, 1 )
                 IF ( gstart == 2 ) bec2(:,1) = bec2(:,1) - vkb(1,jkb+1:jkb+nh(nt))*dvg_tmp(1)
                 CALL mp_sum( bec2( :, 1 ), intra_bgrp_comm )
                 !
                 DO ih = 1,nh(nt)
                    bec1(ih,1) = dvan(ih,ih,nt) * bec1(ih,1)
                    bec2(ih,1) = dvan(ih,ih,nt) * bec2(ih,1)
                 ENDDO
                 !
                 !
                 CALL DGEMM ('N', 'N', 2*npw, 1, nh(nt), 1._DP, vkb(1,jkb+1), &
                      2*npwx, bec1, max(nh(nt),1), 1._DP, dvpsi(1,mu), 2*npwx)
                 CALL DGEMM ('N', 'N', 2*npw, 1, nh(nt), 1._DP, work, &
                      2*npwx, bec2, max(nh(nt),1), 1._DP, dvpsi(1,mu), 2*npwx)
                 !
              ENDDO
              jkb = jkb + nh(nt)
           ENDIF
           !IF (show_barra) CALL update_bar_type( barra, 'dvpnl', 1 )
        ENDIF
     ENDDO
     DEALLOCATE(work)
     DEALLOCATE(bec2)
     DEALLOCATE(bec1)
  ENDDO
  !IF (jkb/=nkb) CALL errore('dvpsi_kb','unexpected error',1)
  !
  !IF (show_barra) CALL stop_bar_type( barra, 'dvpnl' )
  !
END SUBROUTINE
