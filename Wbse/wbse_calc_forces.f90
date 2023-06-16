!
! Copyright (C) 2015-2023 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Yu Jin, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_forces( dvg_exc_tmp )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp
  USE io_push,              ONLY : io_push_title
  USE pwcom,                ONLY : nspin,npwx
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : logfile,nbndval0x,n_trunc_bands
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  !
  IMPLICIT NONE
  !
  ! !/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp( npwx*npol, band_group%nlocx, kpt_pool%nloc )
  !
  ! Workspace
  !
  INTEGER :: n, na, ipol
  REAL(DP), ALLOCATABLE :: forces(:), dvgdvg_mat(:,:,:)
  !$acc declare device_resident(dvgdvg_mat)
  REAL(DP) :: sumforces, time_spent(2)
  REAL(DP), EXTERNAL :: get_clock
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec(:,:,:), zvector(:,:,:), drhox1(:,:), drhox2(:,:)
  LOGICAL :: do_zvector, poor_of_ram_drhox2
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  CALL start_clock( 'l_forces' )
  time_spent(1) = get_clock( 'l_forces' )
  !
  CALL io_push_title('Calculating forces')
  !
  n = 3 * nat
  !
  ALLOCATE( forces(n) )
  forces(:) = 0._DP
  ALLOCATE( dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc) )
  ALLOCATE( drhox1(dffts%nnr, nspin) )
  !
  !
  !
  CALL io_push_title('Compute drhox1')
  !
  CALL wbse_calc_drhox1( dvg_exc_tmp, drhox1 )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title('Compute forces of drhox1')
  !
  CALL wbse_forces_drhox1( n, dvg_exc_tmp, drhox1, forces )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title('Compute < dvg | dvg >')
  !
  CALL wbse_calc_dvgdvg_mat( dvg_exc_tmp, dvgdvg_mat )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title('Compute drhox2')
  !
  ALLOCATE( drhox2(dffts%nnr, nspin) )
  !
  poor_of_ram_drhox2 = .TRUE.
  !
  IF(poor_of_ram_drhox2) THEN
     CALL wbse_calc_drhox2_slow( dvgdvg_mat, drhox2 )
  ELSE
     CALL wbse_calc_drhox2( dvgdvg_mat, drhox2 )
  ENDIF
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title('Compute forces of drhox2')
  !
  CALL wbse_forces_drhox2( n, dvgdvg_mat, drhox2, forces )
  !
  time_spent(2) = get_clock( 'l_forces' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'l_forces' )
  !
  !
  !
  CALL io_push_title('Build and solve the Z-vector equations')
  !
  ! Zvector part
  !
  do_zvector = .TRUE.
  !
  IF( do_zvector ) THEN
     !
     ALLOCATE( z_rhs_vec( npwx, band_group%nlocx, kpt_pool%nloc ) )
     ALLOCATE( zvector( npwx, band_group%nlocx, kpt_pool%nloc ) )
     !
     CALL build_rhs_zvector_eq( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
     !
     time_spent(2) = get_clock( 'l_forces' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'l_forces' )
     !
     !
     !
     CALL solve_zvector_eq_cg( z_rhs_vec, zvector )
     !
     time_spent(2) = get_clock( 'l_forces' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'l_forces' )
     !
     !
     !
     CALL io_push_title('Compute forces of Z vector')
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
  CALL io_push_title('Forces Total')
  !
  DO na = 1, nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) na, ityp(na), ( - forces( 3 * na - 3 + ipol ), ipol = 1, 3 )
     !
  ENDDO
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('exec.forces.forces_total', - forces(1:n))
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
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
        forces( 3 * na - 3 + ipol ) = forces( 3 * na - 3 + ipol ) - sumforces / REAL(nat,KIND=DP)
        !
     ENDDO
     !
  ENDDO
  !
  CALL io_push_title('Forces Corrected')
  !
  DO na = 1,nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) na, ityp(na), (-forces(3*na-3+ipol), ipol=1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('exec.forces.forces_corrected', -forces(1:n))
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
  DEALLOCATE(forces)
  DEALLOCATE(dvgdvg_mat)
  DEALLOCATE(drhox1)
  DEALLOCATE(drhox2)
  !
  CALL stop_clock( 'l_forces' )
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_time(time)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: time(2)
  !
  CHARACTER(20), EXTERNAL :: human_readable_time
  !
  WRITE(stdout, "(5x,'Tot. elapsed time ',a,',  time spent in last step ',a) ") &
  TRIM(human_readable_time(time(2))), TRIM(human_readable_time(time(2)-time(1)))
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox1( dvg_exc_tmp, drhox1 )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,wg,ngk,npwx,npw
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_invfft_gamma,double_invfft_gamma
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : psic=>psic_d
#else
  USE wavefunctions,        ONLY : psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: drhox1(dffts%nnr, nspin)
  !
  ! Workspace
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ir, lbnd, ibnd, jbnd, dffts_nnr
  REAL(DP) :: w1, w2
  REAL(DP), ALLOCATABLE :: tmp_r(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  dffts_nnr = dffts%nnr
  drhox1(:,:) = (0._DP,0._DP)
  !
  ALLOCATE(tmp_r(dffts%nnr))
  !$acc enter data create(tmp_r)
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
        !$acc host_data use_device(dvg_exc_tmp)
        CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),dvg_exc_tmp(:,lbnd+1,iks),psic,'Wave')
        !$acc end host_data
        !
        !$acc parallel loop present(tmp_r)
        DO ir = 1, dffts_nnr
           tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)**2 + w2*AIMAG(psic(ir))**2
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
        !$acc host_data use_device(dvg_exc_tmp)
        CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave')
        !$acc end host_data
        !
        !$acc parallel loop present(tmp_r)
        DO ir = 1, dffts_nnr
           tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)**2
        ENDDO
        !$acc end parallel
        !
     ENDIF
     !
     !$acc update host(tmp_r)
     !
     drhox1(:,current_spin) = CMPLX(tmp_r,KIND=DP)
     !
  ENDDO
  !
  CALL mp_sum(drhox1,inter_bgrp_comm)
  CALL mp_sum(drhox1,inter_pool_comm)
  !
  !$acc exit data delete(tmp_r)
  DEALLOCATE(tmp_r)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhox1( n, dvg_exc_tmp, drhox1, forces )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ntyp=>nsp,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE io_push,              ONLY : io_push_title
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,nspin,current_spin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : logfile,nbnd_occ,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
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
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:)
  !$acc declare device_resident(dvpsi)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, lbnd, ibnd, ig
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_aux(:), forces_drhox1(:), forcelc(:,:), rdrhox1(:,:)
  !$acc declare device_resident(forces_aux)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  CALL io_push_title('Forces : drhox1 part')
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
  ALLOCATE( forces_aux(n) )
  ALLOCATE( forces_drhox1(n) )
  !$acc enter data create(forces_drhox1)
  ALLOCATE( forcelc(3, nat) )
  ALLOCATE( dvpsi(npwx, n) )
  ALLOCATE( rdrhox1(dffts%nnr, nspin) )
  !
  !$acc kernels present(forces_drhox1)
  forces_drhox1(:) = 0._DP
  !$acc end kernels
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
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     this_wk = wk(iks)
     !
     DO lbnd = 1, nbnd_do
        !
        ! 1) | dvpsi_i >
        !
        !$acc kernels present(dvpsi)
        dvpsi(:,:) = (0._DP,0._DP)
        !$acc end kernels
        !
        !$acc host_data use_device(dvg_exc_tmp,dvpsi)
        CALL wbse_get_dvpsi_per_state_gamma_nonlocal(n, dvg_exc_tmp(:,lbnd,iks), dvpsi)
        !$acc end host_data
        !
        ! 2) forces_aux_i = < dvg | dvpsi_i >
        !
        !$acc parallel vector_length(1024) present(dvg_exc_tmp,dvpsi,forces_aux)
        !$acc loop
        DO ia = 1, n
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1, npw
              reduce = reduce + REAL(dvg_exc_tmp(ig,lbnd,iks),KIND=DP) * REAL(dvpsi(ig,ia),KIND=DP) &
              &               + AIMAG(dvg_exc_tmp(ig,lbnd,iks)) * AIMAG(dvpsi(ig,ia))
           ENDDO
           !
           forces_aux(ia) = 2._DP * reduce
           !
        ENDDO
        !$acc end parallel
        !
        IF (gstart == 2) THEN
           !$acc parallel loop present(forces_aux,dvg_exc_tmp,dvpsi)
           DO ia = 1, n
              forces_aux(ia) = forces_aux(ia) - REAL(dvg_exc_tmp(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,ia),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc parallel loop present(forces_drhox1,forces_aux)
        DO ia = 1, n
           forces_drhox1(ia) = forces_drhox1(ia) + factor * this_wk * forces_aux(ia)
        ENDDO
        !$acc end parallel
        !
     ENDDO
     !
  ENDDO
  !
  !$acc update host(forces_drhox1)
  !
  CALL mp_sum(forces_drhox1,intra_bgrp_comm)
  CALL mp_sum(forces_drhox1,inter_bgrp_comm)
  CALL mp_sum(forces_drhox1,inter_pool_comm)
  !
  ! local part
  !
  rdrhox1(:,:) = REAL(drhox1,KIND=DP)
  !
  IF(nspin == 2) THEN
     rdrhox1(:,1) = rdrhox1(:,1) + rdrhox1(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhox1(:,1), gstart, &
  & gamma_only, vloc, forcelc)
  !
  forcelc(:,:) = - factor * forcelc
  !
  DO na = 1, nat
     DO ipol = 1, 3
        forces_drhox1(3 * na - 3 + ipol) = forces_drhox1(3 * na - 3 + ipol) + forcelc(ipol,na)
     ENDDO
  ENDDO
  !
  forces(:) = forces + forces_drhox1
  !
  DO na = 1, nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) na, ityp(na), ( - forces_drhox1(3 * na - 3 + ipol), ipol = 1, 3 )
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('exec.forces.forces_drhox1', - forces_drhox1(1:n))
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
  DEALLOCATE(forces_aux)
  !$acc exit data delete(forces_drhox1)
  DEALLOCATE(forces_drhox1)
  DEALLOCATE(forcelc)
  DEALLOCATE(dvpsi)
  DEALLOCATE(rdrhox1)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_dvgdvg_mat( dvg_exc_tmp, dvgdvg_mat )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE pwcom,                ONLY : ngk,npwx,npw
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE westcom,              ONLY : nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_bgrp_comm,nbgrp,my_bgrp_id,intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(OUT) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  !
  ! WORKSPACE
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, lbnd, ibnd, ig
  REAL(DP) :: reduce
  COMPLEX(DP), ALLOCATABLE :: dvg_exc_tmp_copy(:,:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  !$acc kernels present(dvgdvg_mat)
  dvgdvg_mat(:,:,:) = 0._DP
  !$acc end kernels
  !
  ALLOCATE( dvg_exc_tmp_copy(npwx*npol, nbndval0x-n_trunc_bands) )
  !$acc enter data create(dvg_exc_tmp_copy)
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
     !$acc kernels present(dvg_exc_tmp_copy)
     dvg_exc_tmp_copy(:,:) = (0._DP,0._DP)
     !$acc end kernels
     !
     !$acc parallel loop collapse(2) present(dvg_exc_tmp_copy,dvg_exc_tmp)
     DO lbnd = 1, nbnd_do
        DO ig = 1, npwx
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
           !
           dvg_exc_tmp_copy(ig,ibnd) = dvg_exc_tmp(ig,lbnd,iks)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc update host(dvg_exc_tmp_copy)
     CALL mp_sum(dvg_exc_tmp_copy,inter_bgrp_comm)
     !$acc update device(dvg_exc_tmp_copy)
     !
     !$acc parallel present(dvgdvg_mat,dvg_exc_tmp_copy,dvg_exc_tmp)
     !$acc loop collapse(2)
     DO lbnd = 1, nbnd_do
        DO ibnd = 1, nbndval - n_trunc_bands
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1, npw
              reduce = reduce + REAL(dvg_exc_tmp_copy(ig,ibnd),KIND=DP)*REAL(dvg_exc_tmp(ig,lbnd,iks),KIND=DP) &
              &               + AIMAG(dvg_exc_tmp_copy(ig,ibnd))*AIMAG(dvg_exc_tmp(ig,lbnd,iks))
           ENDDO
           !
           dvgdvg_mat(ibnd,lbnd,iks) = 2._DP * reduce
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gstart == 2) THEN
        !$acc parallel present(dvgdvg_mat,dvg_exc_tmp_copy,dvg_exc_tmp)
        !$acc loop collapse(2)
        DO lbnd = 1, nbnd_do
           DO ibnd = 1, nbndval - n_trunc_bands
              dvgdvg_mat(ibnd,lbnd,iks) = dvgdvg_mat(ibnd,lbnd,iks) &
              & - REAL(dvg_exc_tmp_copy(1,ibnd),KIND=DP) * REAL(dvg_exc_tmp(1,lbnd,iks),KIND=DP)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
  ENDDO
  !
  !$acc host_data use_device(dvgdvg_mat)
  CALL mp_sum(dvgdvg_mat,intra_bgrp_comm)
  !$acc end host_data
  !
  !$acc exit data delete(dvg_exc_tmp_copy)
  DEALLOCATE( dvg_exc_tmp_copy )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox2( dvgdvg_mat, drhox2 )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE pwcom,                ONLY : wg,ngk,nspin,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_invfft_gamma
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: drhox2(dffts%nnr, nspin)
  !
  ! WORKSPACE
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ir, lbnd, ibnd, jbnd
  REAL(DP) :: prod, w1
  COMPLEX(DP), ALLOCATABLE :: aux_all_r(:,:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  ALLOCATE( aux_all_r(dffts%nnr, nbndval0x-n_trunc_bands) )
  !
  aux_all_r(:,:) = (0._DP,0._DP)
  !
  drhox2(:,:) = (0._DP,0._DP)
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
     DO ibnd = n_trunc_bands+1, nbndval
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
        !
        aux_all_r(:,ibnd) = psic(:)
        !
     ENDDO
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
  ENDDO
  !
  CALL mp_sum(drhox2,inter_bgrp_comm)
  CALL mp_sum(drhox2,inter_pool_comm)
  !
  DEALLOCATE( aux_all_r )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox2_slow( dvgdvg_mat, drhox2 )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE pwcom,                ONLY : isk,lsda,wg,ngk,current_spin,nspin,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma,single_invfft_gamma
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,        ONLY : evc_host=>evc
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: drhox2(dffts%nnr, nspin)
  !
  ! WORKSPACE
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ir, lbnd, ibnd, jbnd, jbndp, dffts_nnr
  REAL(DP) :: prod, w1
  REAL(DP), ALLOCATABLE :: aux_r(:)
  !$acc declare device_resident(aux_r)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( aux_r(dffts%nnr) )
  !
  !$acc enter data create(drhox2)
  !
  !$acc kernels present(drhox2)
  drhox2(:,:) = (0._DP,0._DP)
  !$acc end kernels
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
     IF(lsda) current_spin = isk(iks)
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
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        !
        w1 = wg(ibnd,iks_do)/omega
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
        !
        !$acc parallel loop present(aux_r)
        DO ir = 1, dffts_nnr
           aux_r(ir) = REAL(psic(ir),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        DO jbnd = 1, nbndval - n_trunc_bands, 2
           !
           jbndp = jbnd + n_trunc_bands ! index for evc
           !
           IF(jbnd < nbndval - n_trunc_bands) THEN
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,jbndp),evc_work(:,jbndp+1),psic,'Wave')
              !
              !$acc parallel loop present(aux_r,dvgdvg_mat,drhox2)
              DO ir = 1, dffts_nnr
                 prod = aux_r(ir) * ( REAL(psic(ir),KIND=DP) * dvgdvg_mat(jbnd,lbnd,iks) &
                 &                  + AIMAG(psic(ir)) * dvgdvg_mat(jbnd+1,lbnd,iks) )
                 drhox2(ir,current_spin) = drhox2(ir,current_spin) - w1 * CMPLX(prod,KIND=DP)
              ENDDO
              !$acc end parallel
              !
           ELSE
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,jbndp),psic,'Wave')
              !
              !$acc parallel loop present(aux_r,dvgdvg_mat,drhox2)
              DO ir = 1, dffts_nnr
                 prod = aux_r(ir) * REAL(psic(ir),KIND=DP) * dvgdvg_mat(jbnd,lbnd,iks)
                 drhox2(ir,current_spin) = drhox2(ir,current_spin) - w1 * CMPLX(prod,KIND=DP)
              ENDDO
              !$acc end parallel
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  !$acc exit data copyout(drhox2)
  !
  CALL mp_sum(drhox2,inter_bgrp_comm)
  CALL mp_sum(drhox2,inter_pool_comm)
  !
  DEALLOCATE( aux_r )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhox2( n, dvgdvg_mat, drhox2, forces )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ntyp=>nsp,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE io_push,              ONLY : io_push_title
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_spin,nspin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : iuwfc,lrwfc,logfile,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,&
                                 & intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE cublas
#else
  USE wavefunctions,        ONLY : evc_work=>evc
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(IN) :: drhox2(dffts%nnr, nspin)
  REAL(DP), INTENT(INOUT) :: forces(n)
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  ! WORKSPACE
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:), aux_g(:), aux_z(:)
  !$acc declare device_resident(dvpsi,aux_g,aux_z)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, na, lbnd, ibnd, jbnd, ig
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_aux(:), forces_drhox2(:), rdrhox2(:,:), forcelc(:,:)
  !$acc declare device_resident(forces_aux)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Forces : drhox2 part')
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
  ALLOCATE( forces_aux(n) )
  ALLOCATE( forces_drhox2(n) )
  !$acc enter data create(forces_drhox2)
  ALLOCATE( forcelc(3, nat) )
  ALLOCATE( dvpsi(npwx, n) )
  ALLOCATE( rdrhox2(dffts%nnr, nspin) )
  ALLOCATE( aux_g(npwx*npol) )
  ALLOCATE( aux_z(nbndval0x-n_trunc_bands) )
  !
  !$acc kernels present(forces_drhox2)
  forces_drhox2(:) = 0._DP
  !$acc end kernels
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
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     this_wk = wk(iks)
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
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        !
        ! 1) | dvpsi_i >
        !
        !$acc kernels present(dvpsi)
        dvpsi(:,:) = (0._DP,0._DP)
        !$acc end kernels
        !
        !$acc host_data use_device(dvpsi)
        CALL wbse_get_dvpsi_per_state_gamma_nonlocal(n, evc_work(:,ibnd), dvpsi)
        !$acc end host_data
        !
        ! 2) force_aux = < evc_iv2 | dvpsi_ia_iv >
        !
        !$acc parallel loop present(aux_z,dvgdvg_mat)
        DO jbnd = 1, nbndval - n_trunc_bands
           aux_z(jbnd) = CMPLX(dvgdvg_mat(jbnd,lbnd,iks),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        !$acc host_data use_device(dvgdvg_mat,aux_g)
        CALL ZGEMM('N', 'N', npw, 1, nbndval-n_trunc_bands, (1._DP,0._DP), evc_work(1,n_trunc_bands+1), &
        & npwx, aux_z, nbndval-n_trunc_bands, (0._DP,0._DP), aux_g, npwx)
        !$acc end host_data
        !
        !$acc parallel vector_length(1024) present(aux_g,dvpsi,forces_aux)
        !$acc loop
        DO ia = 1, n
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1, npw
              reduce = reduce + REAL(aux_g(ig),KIND=DP) * REAL(dvpsi(ig,ia),KIND=DP) &
              &               + AIMAG(aux_g(ig)) * AIMAG(dvpsi(ig,ia))
           ENDDO
           !
           forces_aux(ia) = 2._DP * reduce
           !
        ENDDO
        !$acc end parallel
        !
        IF (gstart == 2) THEN
           !$acc parallel loop present(forces_aux,aux_g,dvpsi)
           DO ia = 1, n
              forces_aux(ia) = forces_aux(ia) - REAL(aux_g(1),KIND=DP)*REAL(dvpsi(1,ia),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc parallel loop present(forces_drhox2,forces_aux)
        DO ia = 1, n
           forces_drhox2(ia) = forces_drhox2(ia) - factor * this_wk * forces_aux(ia)
        ENDDO
        !$acc end parallel
        !
     ENDDO
     !
  ENDDO
  !
  !$acc update host(forces_drhox2)
  !
  CALL mp_sum(forces_drhox2,intra_bgrp_comm)
  CALL mp_sum(forces_drhox2,inter_bgrp_comm)
  CALL mp_sum(forces_drhox2,inter_pool_comm)
  !
  ! local part
  !
  rdrhox2(:,:) = REAL(drhox2,KIND=DP)
  !
  IF(nspin == 2) THEN
     rdrhox2(:,1) = rdrhox2(:,1) + rdrhox2(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhox2(:,1), gstart, &
  & gamma_only, vloc, forcelc)
  !
  forcelc(:,:) = - factor * forcelc
  !
  DO na = 1, nat
     DO ipol = 1, 3
        forces_drhox2(3 * na - 3 + ipol) = forces_drhox2(3 * na - 3 + ipol) + forcelc(ipol,na)
     ENDDO
  ENDDO
  !
  forces(:) = forces + forces_drhox2
  !
  DO na = 1, nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) na, ityp(na), ( - forces_drhox2( 3 * na - 3 + ipol ), ipol = 1, 3 )
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('exec.forces.forces_drhox2', - forces_drhox2(1:n))
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
  DEALLOCATE( forces_aux )
  !$acc exit data delete(forces_drhox2)
  DEALLOCATE( forces_drhox2 )
  DEALLOCATE( forcelc )
  DEALLOCATE( dvpsi )
  DEALLOCATE( rdrhox2 )
  DEALLOCATE( aux_g )
  DEALLOCATE( aux_z )
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhoz( n, zvector, forces )
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ntyp=>nsp,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE io_push,              ONLY : io_push_title
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_spin,nspin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : iuwfc,lrwfc,logfile,nbnd_occ,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,&
                                 & intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE wavefunctions,        ONLY : evc_host=>evc
#else
  USE wavefunctions,        ONLY : evc_work=>evc
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
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:), drhoz(:,:)
  !$acc declare device_resident(dvpsi)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, lbnd, ibnd, ig
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_aux(:), forces_drhoz(:), forcelc(:,:), rdrhoz(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  CALL io_push_title('Forces : drhoz part')
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
  ALLOCATE( forces_aux(n) )
  ALLOCATE( forces_drhoz(n) )
  !$acc enter data create(forces_drhoz)
  ALLOCATE( forcelc(3, nat) )
  ALLOCATE( dvpsi(npwx, n) )
  ALLOCATE( rdrhoz(dffts%nnr, nspin) )
  ALLOCATE( drhoz(dffts%nnr, nspin) )
  !
  !$acc kernels present(forces_drhoz)
  forces_drhoz(:) = 0._DP
  !$acc end kernels
  !
  ! nonlocal part
  !
  DO iks = 1, kpt_pool%nloc
     !
     ! Z vector always spin-conserving
     !
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
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     this_wk = wk(iks)
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
     DO lbnd = 1, nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        !
        ! 1) | dvpsi_i >
        !
        !$acc kernels present(dvpsi)
        dvpsi(:,:) = (0._DP,0._DP)
        !$acc end kernels
        !
        !$acc host_data use_device(dvpsi)
        CALL wbse_get_dvpsi_per_state_gamma_nonlocal(n, evc_work(:,ibnd), dvpsi)
        !$acc end host_data
        !
        ! 2) force_aux_i = < z_vector | dvpsi_i >
        !
        !$acc parallel vector_length(1024) present(zvector,dvpsi,forces_aux)
        !$acc loop
        DO ia = 1, n
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1, npw
              reduce = reduce + REAL(zvector(ig,lbnd,iks),KIND=DP) * REAL(dvpsi(ig,ia),KIND=DP) &
              &               + AIMAG(zvector(ig,lbnd,iks)) * AIMAG(dvpsi(ig,ia))
           ENDDO
           !
           forces_aux(ia) = 2._DP * reduce
           !
        ENDDO
        !$acc end parallel
        !
        IF (gstart == 2) THEN
           !$acc parallel loop present(forces_aux,zvector,dvpsi)
           DO ia = 1, n
              forces_aux(ia) = forces_aux(ia) - REAL(zvector(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,ia),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc parallel loop present(forces_drhoz,forces_aux)
        DO ia = 1, n
           forces_drhoz(ia) = forces_drhoz(ia) + 2._DP * factor * this_wk * forces_aux(ia)
        ENDDO
        !$acc end parallel
        !
     ENDDO
     !
  ENDDO
  !
  !$acc update host(forces_drhoz)
  !
  CALL mp_sum(forces_drhoz,intra_bgrp_comm)
  CALL mp_sum(forces_drhoz,inter_bgrp_comm)
  CALL mp_sum(forces_drhoz,inter_pool_comm)
  !
  ! local part
  !
  CALL wbse_calc_dens( zvector, drhoz, .FALSE. )
  !
  drhoz(:,:) = 2._DP * drhoz
  rdrhoz(:,:) = REAL(drhoz,KIND=DP)
  !
  IF(nspin == 2) THEN
     rdrhoz(:,1) = rdrhoz(:,1) + rdrhoz(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhoz(:,1), gstart, &
  & gamma_only, vloc, forcelc)
  !
  forcelc(:,:) = - factor * forcelc
  !
  DO na = 1, nat
     DO ipol = 1, 3
        forces_drhoz(3 * na - 3 + ipol) = forces_drhoz(3 * na - 3 + ipol) + forcelc(ipol,na)
     ENDDO
  ENDDO
  !
  forces(:) = forces + forces_drhoz
  !
  DO na = 1, nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) na, ityp(na), ( - forces_drhoz( 3 * na - 3 + ipol ), ipol = 1, 3 )
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('exec.forces.forces_drhoz', - forces_drhoz(1:n))
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
  DEALLOCATE( forces_aux )
  !$acc exit data delete(forces_drhoz)
  DEALLOCATE( forces_drhoz )
  DEALLOCATE( forcelc )
  DEALLOCATE( dvpsi )
  DEALLOCATE( rdrhoz )
  DEALLOCATE( drhoz )
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_get_dvpsi_per_state_gamma_nonlocal (n, dvg_tmp, dvpsi)
  !-----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector u.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  ! It implements Eq. B29 of PRB 64, 235118 (2001). The contribution
  ! of the local pseudopotential is calculated here, that of the nonlocal
  ! pseudopotential in dvqpsi_us_only.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp,ntyp => nsp
  USE cell_base,            ONLY : tpiba
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE gvect,                ONLY : g,gstart
  USE noncollin_module,     ONLY : npol
  use uspp_param,           ONLY : nh
  USE pwcom,                ONLY : npw,npwx
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
#if defined(__CUDA)
  USE uspp,                 ONLY : dvan_work=>dvan_d,dvan_host=>dvan,vkb
  USE cublas
#else
  USE uspp,                 ONLY : dvan_work=>dvan,vkb
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: dvg_tmp(npwx*npol)
  COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx,n)
#if defined(__CUDA)
  ATTRIBUTES(DEVICE) :: dvg_tmp,dvpsi
#endif
  !
  ! Workspace
  !
  INTEGER :: ia, mu, ig, nt, ih, jkb, ic, nh_nt, nh_nt_max
  REAL(DP), ALLOCATABLE :: bec1(:,:), bec2(:,:)
  !$acc declare device_resident(bec1,bec2)
  COMPLEX(DP) :: fact
  COMPLEX(DP), ALLOCATABLE :: work(:,:)
  !$acc declare device_resident(work)
  !
  nh_nt_max = MAXVAL(nh)
  fact = tpiba * (0._DP,-1._DP)
  jkb = 0
#if defined(__CUDA)
  dvan_work(:,:,:) = dvan_host
#endif
  !
  ALLOCATE(work(npwx,nh_nt_max))
  ALLOCATE(bec1(nh_nt_max,1))
  ALLOCATE(bec2(nh_nt_max,1))
  !
  DO nt = 1,ntyp
     nh_nt = nh(nt)
     DO ia = 1, nat
        IF (ityp(ia) == nt) THEN
           IF( nh_nt > 0) THEN
              DO ic = 1,3
                 mu = 3*(ia-1) +ic
                 !
                 ! first term: sum_l sum_G' [ i V_l(G) V^*_l(G') (G'*u) psi(G')
                 ! second term: sum_l sum_G' [-i (G*u) V_l(G) V^*_l(G') psi(G')
                 !
                 !$acc parallel loop collapse(2) present(work,vkb,g)
                 DO ih = 1,nh_nt
                    DO ig = 1,npw
                       work(ig,ih) = vkb(ig,jkb+ih) * g(ic,ig)*fact
                    ENDDO
                 ENDDO
                 !$acc end parallel
                 !
                 !$acc host_data use_device(work,bec1)
                 CALL cudaDGEMV('C', 2*npw, nh_nt, 2._DP, work, 2*npwx, dvg_tmp, 1, 0._DP, bec1, 1)
                 !$acc end host_data
                 !
                 IF(gstart == 2) THEN
                    !$acc parallel loop present(bec1,work)
                    DO ih = 1,nh_nt
                       bec1(ih,1) = bec1(ih,1) - work(1,ih)*dvg_tmp(1)
                    ENDDO
                    !$acc end parallel
                 ENDIF
                 !
                 !$acc host_data use_device(bec1)
                 CALL mp_sum(bec1,intra_bgrp_comm)
                 !$acc end host_data
                 !
                 !$acc host_data use_device(vkb,bec2)
                 CALL cudaDGEMV('C', 2*npw, nh_nt, 2._DP, vkb(:,jkb+1:jkb+nh_nt), 2*npwx, dvg_tmp, &
                 & 1, 0._DP, bec2, 1)
                 !$acc end host_data
                 !
                 IF(gstart == 2) THEN
                    !$acc parallel loop present(bec2,vkb)
                    DO ih = 1,nh_nt
                       bec2(ih,1) = bec2(ih,1) - vkb(1,jkb+ih)*dvg_tmp(1)
                    ENDDO
                    !$acc end parallel
                 ENDIF
                 !
                 !$acc host_data use_device(bec2)
                 CALL mp_sum(bec2,intra_bgrp_comm)
                 !$acc end host_data
                 !
                 !$acc parallel loop present(bec1,bec2)
                 DO ih = 1,nh_nt
                    bec1(ih,1) = dvan_work(ih,ih,nt) * bec1(ih,1)
                    bec2(ih,1) = dvan_work(ih,ih,nt) * bec2(ih,1)
                 ENDDO
                 !$acc end parallel
                 !
                 !$acc host_data use_device(vkb,bec1,work,bec2)
                 CALL DGEMM('N', 'N', 2*npw, 1, nh_nt, 1._DP, vkb(1,jkb+1), 2*npwx, bec1, MAX(nh_nt,1), 1._DP, dvpsi(1,mu), 2*npwx)
                 CALL DGEMM('N', 'N', 2*npw, 1, nh_nt, 1._DP, work, 2*npwx, bec2, MAX(nh_nt,1), 1._DP, dvpsi(1,mu), 2*npwx)
                 !$acc end host_data
                 !
              ENDDO
              jkb = jkb + nh_nt
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  !
  DEALLOCATE(work)
  DEALLOCATE(bec1)
  DEALLOCATE(bec2)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_zvector_eq_cg(z_rhs, z_out)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE mp,                   ONLY : mp_max
  USE noncollin_module,     ONLY : npol
  USE pwcom,                ONLY : npwx,nspin
  USE westcom,              ONLY : forces_zeq_cg_tr,l_pre_shift
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_world,             ONLY : world_comm
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: z_rhs(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: z_out(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! ... Local variables
  !
  INTEGER :: iter
  INTEGER, PARAMETER :: max_cg_iters=1000
  REAL(DP) :: threshold,residual_sq
  COMPLEX(DP) :: alpha,beta
  COMPLEX(DP), ALLOCATABLE :: residual_old(:),residual_new(:),dotp(:),rz_new(:),rz_old(:)
  COMPLEX(DP), ALLOCATABLE :: r_old(:,:,:),r_new(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: p(:,:,:),Ap(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z(:,:,:)
  LOGICAL :: cg_prec, turn_shift
  !
  REAL(DP) :: time_spent(2)
  REAL(DP), EXTERNAL :: get_clock
  !
  CALL start_clock( 'zvector_cg' )
  time_spent(1) = get_clock( 'zvector_cg' )
  !
  cg_prec = .TRUE.
  turn_shift = l_pre_shift
  !
  ALLOCATE(r_new(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(r_old(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(p    (npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(Ap   (npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(z    (npwx*npol, band_group%nlocx, kpt_pool%nloc))
  !
  ALLOCATE(residual_old(nspin))
  ALLOCATE(residual_new(nspin))
  ALLOCATE(dotp(nspin))
  ALLOCATE(rz_new(nspin))
  ALLOCATE(rz_old(nspin))
  !
  CALL io_push_title('Forces : Solve the Z vector equation using the CG algorithm')
  !
  CALL wbse_dot(z_rhs,z_rhs,band_group%nlocx,dotp)
  !
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  WRITE(stdout, "( 5x,'# Norm of z_rhs_vec     = | ', ES15.8, ' |' ) " ) SUM(REAL(dotp,KIND=DP))
  WRITE(stdout, "( 5x,'                          *-----------------*' ) " )
  !
  threshold = forces_zeq_cg_tr
  !
  ! Initial guess
  z_out(:,:,:) = z_rhs
  !
  ! z-vector equation won't have the spin flip case
  CALL west_apply_liouvillian(z_out, r_new, .FALSE.)
  time_spent(2) = get_clock( 'zvector_cg' )
  CALL wbse_forces_time(time_spent)
  time_spent(1) = get_clock( 'zvector_cg' )
  !
  CALL west_apply_liouvillian_btda(z_out, r_new, .FALSE.)
  time_spent(2) = get_clock( 'zvector_cg' )
  CALL wbse_forces_time(time_spent)
  !
  r_new(:,:,:) = z_rhs - r_new
  !
  CALL wbse_dot(r_new,r_new,band_group%nlocx,residual_new)
  !
  WRITE(stdout,"(5x,'Initial residual = ',E15.8)") REAL(SUM(residual_new),KIND=DP)
  WRITE(stdout,*)
  !
  IF (cg_prec) THEN
     !
     CALL cg_precondition(r_new,z,turn_shift)
     !
     CALL wbse_dot(r_new,z,band_group%nlocx,rz_new)
     !
     p(:,:,:) = z
     !
  ELSE
     !
     p(:,:,:) = r_new
     !
  ENDIF
  !
  r_old(:,:,:) = r_new
  residual_old(:) = residual_new
  !
  IF (cg_prec) THEN
     rz_old(:) = rz_new
  ENDIF
  !
  DO iter = 1, max_cg_iters
     !
     time_spent(1) = get_clock( 'zvector_cg' )
     !
     CALL west_apply_liouvillian(p, Ap, .FALSE.)
     time_spent(2) = get_clock( 'zvector_cg' )
     CALL wbse_forces_time(time_spent)
     time_spent(1) = get_clock( 'zvector_cg' )
     !
     CALL west_apply_liouvillian_btda(p, Ap, .FALSE.)
     time_spent(2) = get_clock( 'zvector_cg' )
     CALL wbse_forces_time(time_spent)
     !
     CALL wbse_dot(p,Ap,band_group%nlocx,dotp)
     !
     IF (cg_prec) THEN
        alpha = SUM(rz_old) / SUM(dotp)
     ELSE
        alpha = SUM(residual_old) / SUM(dotp)
     ENDIF
     !
     z_out(:,:,:) = z_out + alpha * p
     r_new(:,:,:) = r_old - alpha * Ap
     !
     CALL wbse_dot(r_new,r_new,band_group%nlocx,residual_new)
     !
     IF ( MOD(iter,1) == 0 ) THEN
        WRITE(stdout,"(5x,'Residual(',I5.5,') = ',E15.8)") iter, REAL(SUM(residual_new),KIND=DP)
        WRITE(stdout,*)
     ENDIF
     !
     residual_sq = ABS(SUM(residual_new))
     CALL mp_max(residual_sq,world_comm)
     !
     IF (residual_sq < threshold) EXIT
     !
     IF (cg_prec) THEN
        !
        CALL cg_precondition(r_new,z,turn_shift)
        !
        CALL wbse_dot(r_new,z,band_group%nlocx,rz_new)
        !
        beta = SUM(rz_new) / SUM(rz_old)
        !
        p(:,:,:) = z + beta * p
        !
     ELSE
        !
        beta = SUM(residual_new) / SUM(residual_old)
        !
        p(:,:,:) = r_new + beta * p
        !
     ENDIF
     !
     r_old(:,:,:) = r_new
     residual_old(:) = residual_new
     !
     IF (cg_prec) THEN
        rz_old(:) = rz_new
     ENDIF
     !
  ENDDO ! iter
  !
  DEALLOCATE(r_new)
  DEALLOCATE(r_old)
  DEALLOCATE(p)
  DEALLOCATE(Ap)
  DEALLOCATE(z)
  DEALLOCATE(residual_old)
  DEALLOCATE(residual_new)
  DEALLOCATE(dotp)
  DEALLOCATE(rz_old)
  DEALLOCATE(rz_new)
  !
  CALL stop_clock( 'zvector_cg' )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE west_apply_liouvillian_btda(evc1,evc1_new,sf)
  !-----------------------------------------------------------------------
  !
  ! Applies the linear response operator to response wavefunctions, the part
  ! beyond the Tamm-Dancoff approximation (btda)
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE mp,                   ONLY : mp_bcast
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                 & double_invfft_gamma
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,lsda,nspin,current_k,ngk
  USE gvect,                ONLY : gstart
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE westcom,              ONLY : l_bse,l_bse_triplet,nbnd_occ,n_trunc_bands,lrwfc,iuwfc,&
                                 & l_hybrid_tddft,l_spin_flip_kernel
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE wbse_dv,              ONLY : wbse_dv_of_drho,wbse_dv_of_drho_sf
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
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: evc1(npwx*npol,band_group%nlocx,kpt_pool%nloc)
  COMPLEX(DP), INTENT(INOUT) :: evc1_new(npwx*npol,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: sf
  !
  ! ... Local variables
  !
  LOGICAL :: lrpa,do_k2e
  INTEGER :: ibnd,jbnd,iks,iks_do,ir,ig,nbndval,flnbndval,nbnd_do,lbnd
  INTEGER :: dffts_nnr
  INTEGER, PARAMETER :: flks(2) = [2,1]
  COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc2_new(:,:)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('apply_lv_btda')
#else
  CALL start_clock('apply_lv_btda')
#endif
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(dvrs(dffts%nnr,nspin))
  ALLOCATE(evc2_new(npwx*npol,band_group%nlocx))
  !
  ! Calculation of the charge density response
  !
  CALL wbse_calc_dens(evc1,dvrs,sf)
  !
  lrpa = l_bse
  !
  If(sf .AND. l_spin_flip_kernel) THEN
     CALL wbse_dv_of_drho_sf(dvrs)
  ELSE
     CALL wbse_dv_of_drho(dvrs,lrpa,.FALSE.)
  ENDIF
  !
  DO iks = 1,kpt_pool%nloc
     !
     !$acc kernels present(evc2_new)
     evc2_new(:,:) = (0._DP,0._DP)
     !$acc end kernels
     !
     IF(sf) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     flnbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
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
     IF(l_bse_triplet) THEN
        do_k2e = .FALSE.
     ELSEIF(sf .AND. (.NOT. l_spin_flip_kernel)) THEN
        do_k2e = .FALSE.
     ELSEIF(sf .AND. l_spin_flip_kernel) THEN
        do_k2e = .TRUE.
     ELSE
        do_k2e = .TRUE.
     ENDIF
     !
     IF(do_k2e) THEN
        !
        ! double bands @ gamma
        !
        DO lbnd=1, nbnd_do-MOD(nbnd_do,2),2
           !
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),evc_work(:,jbnd),psic,'Wave')
           !
           !$acc parallel loop present(dvrs)
           DO ir=1,dffts_nnr
              psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(evc2_new)
           CALL double_fwfft_gamma(dffts,npw,npwx,psic,evc2_new(:,lbnd),evc2_new(:,lbnd+1),'Wave')
           !$acc end host_data
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
           CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
           !
           !$acc parallel loop present(dvrs)
           DO ir = 1,dffts_nnr
              psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(evc2_new)
           CALL single_fwfft_gamma(dffts,npw,npwx,psic,evc2_new(:,lbnd),'Wave')
           !$acc end host_data
           !
        ENDIF
        !
     ENDIF
     !
     IF (l_bse) THEN
        !
        ! The other part beyond TDA. exx_div treatment is not needed for this part.
        !
        CALL errore('west_apply_liouvillian_btda', 'BSE forces have not been implemented', 1)
        !
     ENDIF
     !
     IF (l_hybrid_tddft) THEN
        !
        ! K2d part. exx_div treatment is not needed for this part.
        !
        CALL hybrid_kernel_term2(current_spin,evc1,evc2_new,sf)
        !
     ENDIF
     !
     IF(gstart == 2) THEN
        !$acc parallel loop present(evc1_new)
        DO lbnd = 1,nbnd_do
           evc2_new(1,lbnd) = CMPLX(REAL(evc2_new(1,lbnd),KIND=DP),KIND=DP)
        ENDDO
        !$acc end parallel
     ENDIF
     !
     ! Pc[k]*evc2_new(k)
     !
     ! load evc from iks to apply Pc of the current spin channel
     !
     IF(kpt_pool%nloc > 1) THEN
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
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,evc2_new,(1._DP,0._DP))
     !
     !$acc parallel loop collapse(2) present(evc1_new,evc2_new)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           evc1_new(ig,lbnd,iks) = evc1_new(ig,lbnd,iks)+evc2_new(ig,lbnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  DEALLOCATE(dvrs)
  DEALLOCATE(evc2_new)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('apply_lv_btda')
#else
  CALL stop_clock('apply_lv_btda')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE cg_precondition(x, px, turn_shift)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_bcast
  USE buffers,              ONLY : get_buffer
  USE pwcom,                ONLY : npwx
  USE westcom,              ONLY : nbnd_occ,lrwfc,iuwfc,n_trunc_bands
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,nbgrp,my_bgrp_id
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE wvfct,                ONLY : g2kin
  USE wvfct_gpum,           ONLY : et=>et_d
  USE west_gpu,             ONLY : reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc
  USE wvfct,                ONLY : g2kin,et
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: x(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: turn_shift
  COMPLEX(DP), INTENT(OUT) :: px(npwx,band_group%nlocx,kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: ig, ibnd, nbndval, nbnd_do, lbnd, iks
  REAL(DP):: tmp,tmp_abs,tmp_sgn
  REAL(DP), ALLOCATABLE :: g2kin_save(:,:)
  REAL(DP), PARAMETER :: minimum = 1._DP
  !
#if defined(__CUDA)
  CALL start_clock_gpu('precd_cg')
#else
  CALL start_clock('precd_cg')
#endif
  !
  px(:,:,:) = (0._DP, 0._DP)
  !
  ALLOCATE(g2kin_save(npwx,kpt_pool%nloc))
  !$acc enter data create(g2kin_save)
  !
  DO iks = 1,kpt_pool%nloc
     !
     CALL g2_kin(iks)
     !
     !$acc kernels present(g2kin_save,g2kin)
     g2kin_save(:,iks) = g2kin
     !$acc end kernels
     !
     nbndval = nbnd_occ(iks)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     !$acc parallel vector_length(1024) present(g2kin_save,x)
     !$acc loop
     DO lbnd = 1,nbnd_do
        DO ig = 1,npwx
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
           !
           IF (turn_shift) THEN
              tmp = g2kin_save(ig,iks)-et(ibnd+n_trunc_bands,iks)
           ELSE
              tmp = g2kin_save(ig,iks)
           ENDIF
           !
           ! Same as the following line but without thread divergence
           ! IF(ABS(tmp) < minimum) tmp = SIGN(minimum,tmp)
           !
           tmp_abs = MAX(ABS(tmp),minimum)
           tmp_sgn = SIGN(1._DP,tmp)
           tmp = tmp_sgn*tmp_abs
           !
           px(ig,lbnd,iks) = x(ig,lbnd,iks)/tmp
           !
        ENDDO
        !
     ENDDO
     !$acc end parallel
     !
     ! Pc[k]*px(k)
     !
     ! load evc from iks to apply Pc of the current spin channel
     !
     IF(kpt_pool%nloc > 1) THEN
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
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     !
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,px(:,:,iks),(1._DP,0._DP))
     !
  ENDDO
  !
  !$acc exit data delete(g2kin_save) copyout(px)
  DEALLOCATE(g2kin_save)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('precd_cg')
#else
  CALL stop_clock('precd_cg')
#endif
  !
END SUBROUTINE
