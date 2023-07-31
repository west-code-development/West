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
  USE pwcom,                ONLY : nspin,npwx
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : logfile,nbndval0x,n_trunc_bands,evc1_all
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
  USE wbse_bgrp,            ONLY : gather_bands
  USE mp,                   ONLY : mp_waitall
  !
  IMPLICIT NONE
  !
  ! !/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp( npwx*npol, band_group%nlocx, kpt_pool%nloc )
  !
  ! Workspace
  !
  INTEGER :: iks, n, na, ipol
  INTEGER, ALLOCATABLE :: reqs(:)
  REAL(DP), ALLOCATABLE :: forces(:), dvgdvg_mat(:,:,:)
  !$acc declare device_resident(dvgdvg_mat)
  REAL(DP) :: sumforces
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec(:,:,:), zvector(:,:,:), drhox1(:,:), drhox2(:,:)
  !$acc declare device_resident(z_rhs_vec,zvector)
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  CALL start_clock('calc_force')
  !
  CALL io_push_title('Compute forces')
  !
  n = 3 * nat
  !
  ALLOCATE( reqs(kpt_pool%nloc) )
  ALLOCATE( forces(n) )
  forces(:) = 0._DP
  ALLOCATE( dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc) )
  ALLOCATE( drhox1(dffts%nnr, nspin) )
  !
  DO iks = 1,kpt_pool%nloc
     CALL gather_bands( dvg_exc_tmp(:,:,iks), evc1_all(:,:,iks), reqs(iks) )
  ENDDO
  !
  ! drhox1
  !
  CALL wbse_calc_drhox1( dvg_exc_tmp, drhox1 )
  !
  CALL wbse_forces_drhox1( n, dvg_exc_tmp, drhox1, forces )
  !
  ! < dvg | dvg >
  !
  CALL mp_waitall(reqs)
#if !defined(__GPU_MPI)
  !$acc update device(evc1_all)
#endif
  !
  CALL wbse_calc_dvgdvg_mat( dvg_exc_tmp, dvgdvg_mat )
  !
  ! drhox2
  !
  ALLOCATE( drhox2(dffts%nnr, nspin) )
  !
  CALL wbse_calc_drhox2( dvgdvg_mat, drhox2 )
  !
  CALL wbse_forces_drhox2( n, dvgdvg_mat, drhox2, forces )
  !
  ! Z vector
  !
  ALLOCATE( z_rhs_vec( npwx, band_group%nlocx, kpt_pool%nloc ) )
  ALLOCATE( zvector( npwx, band_group%nlocx, kpt_pool%nloc ) )
  !
  CALL build_rhs_zvector_eq( dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec )
  !
  CALL solve_zvector_eq_cg( z_rhs_vec, zvector )
  !
  CALL wbse_forces_drhoz( n, zvector, forces )
  !
  DEALLOCATE( z_rhs_vec )
  DEALLOCATE( zvector )
  !
  CALL io_push_title('Forces total')
  !
  DO na = 1,nat
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
  DO ipol = 1,3
     !
     sumforces = 0._DP
     !
     DO na = 1,nat
        sumforces = sumforces + forces( 3 * na - 3 + ipol )
     ENDDO
     !
     DO na = 1,nat
        forces( 3 * na - 3 + ipol ) = forces( 3 * na - 3 + ipol ) - sumforces / REAL(nat,KIND=DP)
     ENDDO
     !
  ENDDO
  !
  CALL io_push_title('Forces corrected')
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
  DEALLOCATE(reqs)
  DEALLOCATE(forces)
  DEALLOCATE(dvgdvg_mat)
  DEALLOCATE(drhox1)
  DEALLOCATE(drhox2)
  !
  CALL stop_clock('calc_force')
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
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
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm
  USE io_push,              ONLY : io_push_title
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
  INTEGER :: barra_load
  REAL(DP) :: w1, w2
  REAL(DP), ALLOCATABLE :: tmp_r(:)
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute drhox1')
  !
  dffts_nnr = dffts%nnr
  drhox1(:,:) = (0._DP,0._DP)
  !
  ALLOCATE(tmp_r(dffts%nnr))
  !$acc enter data create(tmp_r)
  !
  barra_load = 0
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
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) barra_load = barra_load+1
     ENDDO
     !
  ENDDO
  !
  CALL start_bar_type(barra,'drhox1',barra_load)
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
        DO ir = 1,dffts_nnr
           tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)**2 + w2*AIMAG(psic(ir))**2
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'drhox1',2)
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
        DO ir = 1,dffts_nnr
           tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)**2
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'drhox1',1)
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
  CALL stop_bar_type(barra,'drhox1')
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
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,nspin,current_spin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : logfile,nbnd_occ,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc), drhox1(dffts%nnr, nspin)
  REAL(DP), INTENT(INOUT) :: forces(n)
  !
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:)
  !$acc declare device_resident(dvpsi)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, lbnd, ibnd, ig
  INTEGER :: barra_load
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_aux(:), forces_drhox1(:), forcelc(:,:), rdrhox1(:,:)
  !$acc declare device_resident(forces_aux)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute forces of drhox1')
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
  barra_load = 0
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
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) barra_load = barra_load+1
     ENDDO
     !
  ENDDO
  !
  CALL start_bar_type(barra,'f_drhox1',barra_load)
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
     DO lbnd = 1,nbnd_do
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
        DO ia = 1,n
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1,npw
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
           DO ia = 1,n
              forces_aux(ia) = forces_aux(ia) - REAL(dvg_exc_tmp(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,ia),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc parallel loop present(forces_drhox1,forces_aux)
        DO ia = 1,n
           forces_drhox1(ia) = forces_drhox1(ia) + factor * this_wk * forces_aux(ia)
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'f_drhox1',1)
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
  CALL stop_bar_type(barra,'f_drhox1')
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
  DO na = 1,nat
     DO ipol = 1,3
        forces_drhox1(3 * na - 3 + ipol) = forces_drhox1(3 * na - 3 + ipol) + forcelc(ipol,na)
     ENDDO
  ENDDO
  !
  forces(:) = forces + forces_drhox1
  !
  CALL io_push_title('Forces drhox1')
  !
  DO na = 1,nat
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
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip,evc1_all
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : intra_bgrp_comm
  USE io_push,              ONLY : io_push_title
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(OUT) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, lbnd, ibnd, ig
  REAL(DP) :: reduce
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute <dvg|dvg>')
  !
  !$acc kernels present(dvgdvg_mat)
  dvgdvg_mat(:,:,:) = 0._DP
  !$acc end kernels
  !
  CALL start_bar_type(barra,'dvgdvg',kpt_pool%nloc)
  !
  DO iks = 1,kpt_pool%nloc
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
     !$acc parallel present(dvgdvg_mat,evc1_all,dvg_exc_tmp)
     !$acc loop collapse(2)
     DO lbnd = 1,nbnd_do
        DO ibnd = 1,nbndval - n_trunc_bands
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1,npw
              reduce = reduce &
              & + REAL(evc1_all(ig,ibnd,iks),KIND=DP)*REAL(dvg_exc_tmp(ig,lbnd,iks),KIND=DP) &
              & + AIMAG(evc1_all(ig,ibnd,iks))*AIMAG(dvg_exc_tmp(ig,lbnd,iks))
           ENDDO
           !
           dvgdvg_mat(ibnd,lbnd,iks) = 2._DP * reduce
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gstart == 2) THEN
        !$acc parallel present(dvgdvg_mat,evc1_all,dvg_exc_tmp)
        !$acc loop collapse(2)
        DO lbnd = 1,nbnd_do
           DO ibnd = 1,nbndval - n_trunc_bands
              dvgdvg_mat(ibnd,lbnd,iks) = dvgdvg_mat(ibnd,lbnd,iks) &
              & - REAL(evc1_all(1,ibnd,iks),KIND=DP) * REAL(dvg_exc_tmp(1,lbnd,iks),KIND=DP)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     CALL update_bar_type(barra,'dvgdvg',1)
     !
  ENDDO
  !
  !$acc host_data use_device(dvgdvg_mat)
  CALL mp_sum(dvgdvg_mat,intra_bgrp_comm)
  !$acc end host_data
  !
  CALL stop_bar_type(barra,'dvgdvg')
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox2( dvgdvg_mat, drhox2 )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE pwcom,                ONLY : isk,lsda,wg,ngk,current_spin,nspin,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma,single_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm
  USE io_push,              ONLY : io_push_title
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
  ! Workspace
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ir, lbnd, ibnd, jbnd, jbndp, dffts_nnr
  INTEGER :: barra_load
  REAL(DP) :: prod, w1
  REAL(DP), ALLOCATABLE :: aux_r(:)
  !$acc declare device_resident(aux_r)
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute drhox2')
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
  barra_load = 0
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
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) barra_load = barra_load+1
     ENDDO
     !
  ENDDO
  !
  CALL start_bar_type(barra,'drhox2',barra_load)
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
     DO lbnd = 1,nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        !
        w1 = wg(ibnd,iks_do)/omega
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd),psic,'Wave')
        !
        !$acc parallel loop present(aux_r)
        DO ir = 1,dffts_nnr
           aux_r(ir) = REAL(psic(ir),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        DO jbnd = 1,nbndval - n_trunc_bands, 2
           !
           jbndp = jbnd + n_trunc_bands ! index for evc
           !
           IF(jbnd < nbndval - n_trunc_bands) THEN
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,jbndp),evc_work(:,jbndp+1),psic,'Wave')
              !
              !$acc parallel loop present(aux_r,dvgdvg_mat,drhox2)
              DO ir = 1,dffts_nnr
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
              DO ir = 1,dffts_nnr
                 prod = aux_r(ir) * REAL(psic(ir),KIND=DP) * dvgdvg_mat(jbnd,lbnd,iks)
                 drhox2(ir,current_spin) = drhox2(ir,current_spin) - w1 * CMPLX(prod,KIND=DP)
              ENDDO
              !$acc end parallel
              !
           ENDIF
           !
        ENDDO
        !
        CALL update_bar_type(barra,'drhox2',1)
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
  CALL stop_bar_type(barra,'drhox2')
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
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_spin,nspin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : iuwfc,lrwfc,logfile,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,&
                                 & intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
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
  !
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:), aux_g(:)
  !$acc declare device_resident(dvpsi,aux_g)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, na, lbnd, ibnd, ig
  INTEGER :: barra_load
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_aux(:), forces_drhox2(:), rdrhox2(:,:), forcelc(:,:)
  !$acc declare device_resident(forces_aux)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute forces of drhox2')
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
  !
  !$acc kernels present(forces_drhox2)
  forces_drhox2(:) = 0._DP
  !$acc end kernels
  !
  barra_load = 0
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
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) barra_load = barra_load+1
     ENDDO
     !
  ENDDO
  !
  CALL start_bar_type(barra,'f_drhox2',barra_load)
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
     DO lbnd = 1,nbnd_do
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
        !$acc host_data use_device(dvgdvg_mat,aux_g)
        CALL DGEMM('N', 'N', 2*npw, 1, nbndval-n_trunc_bands, 1._DP, evc_work(1,n_trunc_bands+1), &
        & 2*npwx, dvgdvg_mat(1,lbnd,iks), nbndval0x-n_trunc_bands, 0._DP, aux_g, 2*npwx)
        !$acc end host_data
        !
        !$acc parallel vector_length(1024) present(aux_g,dvpsi,forces_aux)
        !$acc loop
        DO ia = 1,n
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1,npw
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
           DO ia = 1,n
              forces_aux(ia) = forces_aux(ia) - REAL(aux_g(1),KIND=DP)*REAL(dvpsi(1,ia),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc parallel loop present(forces_drhox2,forces_aux)
        DO ia = 1,n
           forces_drhox2(ia) = forces_drhox2(ia) - factor * this_wk * forces_aux(ia)
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'f_drhox2',1)
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
  CALL stop_bar_type(barra,'f_drhox2')
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
  DO na = 1,nat
     DO ipol = 1,3
        forces_drhox2(3 * na - 3 + ipol) = forces_drhox2(3 * na - 3 + ipol) + forcelc(ipol,na)
     ENDDO
  ENDDO
  !
  forces(:) = forces + forces_drhox2
  !
  CALL io_push_title('Forces drhox2')
  !
  DO na = 1,nat
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
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_spin,nspin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : iuwfc,lrwfc,logfile,nbnd_occ,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,&
                                 & intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
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
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:), drhoz(:,:)
  !$acc declare device_resident(dvpsi)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, na, ipol, lbnd, ibnd, ig
  INTEGER :: barra_load
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_aux(:), forces_drhoz(:), forcelc(:,:), rdrhoz(:,:)
  !$acc declare device_resident(forces_aux)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  !
  CALL io_push_title('Compute forces of Z vector')
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
  !$acc enter data create(drhoz)
  !
  !$acc kernels present(forces_drhoz)
  forces_drhoz(:) = 0._DP
  !$acc end kernels
  !
  barra_load = 0
  !
  DO iks = 1,kpt_pool%nloc
     !
     nbndval = nbnd_occ(iks)
     !
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) barra_load = barra_load+1
     ENDDO
     !
  ENDDO
  !
  CALL start_bar_type(barra,'f_drhoxz',barra_load)
  !
  ! nonlocal part
  !
  DO iks = 1,kpt_pool%nloc
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
     DO lbnd = 1,nbnd_do
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
        DO ia = 1,n
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1,npw
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
           DO ia = 1,n
              forces_aux(ia) = forces_aux(ia) - REAL(zvector(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,ia),KIND=DP)
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc parallel loop present(forces_drhoz,forces_aux)
        DO ia = 1,n
           forces_drhoz(ia) = forces_drhoz(ia) + 2._DP * factor * this_wk * forces_aux(ia)
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'f_drhoxz',1)
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
  CALL stop_bar_type(barra,'f_drhoxz')
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
  DO na = 1,nat
     DO ipol = 1,3
        forces_drhoz(3 * na - 3 + ipol) = forces_drhoz(3 * na - 3 + ipol) + forcelc(ipol,na)
     ENDDO
  ENDDO
  !
  forces(:) = forces + forces_drhoz
  !
  CALL io_push_title('Forces drhoz')
  !
  DO na = 1,nat
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
  !$acc exit data delete(drhoz)
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
