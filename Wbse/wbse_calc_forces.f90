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
! Yu Jin, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_forces(dvg_exc_tmp)
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
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_bse_gpu,deallocate_bse_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! !/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: iks, n, ia, ipol
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
  ALLOCATE(reqs(kpt_pool%nloc))
  ALLOCATE(forces(n))
  forces(:) = 0._DP
  ALLOCATE(dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(drhox1(dffts%nnr, nspin))
  !
  DO iks = 1,kpt_pool%nloc
     CALL gather_bands(dvg_exc_tmp(:,:,iks), evc1_all(:,:,iks), reqs(iks))
  ENDDO
  !
  ! drhox1
  !
  CALL wbse_calc_drhox1(dvg_exc_tmp, drhox1)
  !
  CALL wbse_forces_drhox1(n, dvg_exc_tmp, drhox1, forces)
  !
  ! < dvg | dvg >
  !
  CALL mp_waitall(reqs)
#if !defined(__GPU_MPI)
  !$acc update device(evc1_all)
#endif
  !
  CALL wbse_calc_dvgdvg_mat(dvg_exc_tmp, dvgdvg_mat)
  !
  ! drhox2
  !
  ALLOCATE(drhox2(dffts%nnr, nspin))
  !
  CALL wbse_calc_drhox2(dvgdvg_mat, drhox2)
  !
  CALL wbse_forces_drhox2(n, dvgdvg_mat, drhox2, forces)
  !
  ! Z vector
  !
  ALLOCATE(z_rhs_vec(npwx, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(zvector(npwx, band_group%nlocx, kpt_pool%nloc))
  !
#if defined(__CUDA)
  CALL allocate_bse_gpu(band_group%nlocx)
#endif
  !
  CALL build_rhs_zvector_eq(dvg_exc_tmp, dvgdvg_mat, drhox1, drhox2, z_rhs_vec)
  !
  CALL solve_zvector_eq_cg(z_rhs_vec, zvector)
  !
#if defined(__CUDA)
  CALL deallocate_bse_gpu()
#endif
  !
  CALL wbse_forces_drhoz(n, zvector, forces)
  !
  DEALLOCATE(z_rhs_vec)
  DEALLOCATE(zvector)
  !
  CALL io_push_title('Forces total')
  !
  DO ia = 1,nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) ia, ityp(ia), (-forces(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.forces.forces_total', -forces(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
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
     DO ia = 1,nat
        sumforces = sumforces + forces(3*ia-3+ipol)
     ENDDO
     !
     DO ia = 1,nat
        forces(3*ia-3+ipol) = forces(3*ia-3+ipol) - sumforces/REAL(nat,KIND=DP)
     ENDDO
     !
  ENDDO
  !
  CALL io_push_title('Forces corrected')
  !
  DO ia = 1,nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) ia, ityp(ia), (-forces(3*ia-3+ipol), ipol=1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.forces.forces_corrected', -forces(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
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
SUBROUTINE wbse_calc_drhox1(dvg_exc_tmp, drhox1)
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
  USE wavefunctions,        ONLY : psic
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
        CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),dvg_exc_tmp(:,lbnd+1,iks),psic,'Wave')
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
        CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave')
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
SUBROUTINE wbse_forces_drhox1(n, dvg_exc_tmp, drhox1, forces)
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
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_forces_gpu,deallocate_forces_gpu
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
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:,:)
  !$acc declare device_resident(dvpsi)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, lbnd, ibnd, ig
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_drhox1(:), forcelc(:,:), rdrhox1(:,:)
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
#if defined(__CUDA)
  CALL allocate_forces_gpu()
#endif
  !
  ALLOCATE(forces_drhox1(n))
  ALLOCATE(forcelc(3, nat))
  ALLOCATE(dvpsi(npwx, band_group%nlocx, 3))
  ALLOCATE(rdrhox1(dffts%nnr, nspin))
  !
  forces_drhox1(:) = 0._DP
  !
  CALL start_bar_type(barra,'f_drhox1',kpt_pool%nloc*nat)
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
     this_wk = wk(iks)*factor
     !
     DO ia = 1,nat
        !
        ! 1) | dvpsi_i >
        !
        CALL wbse_get_dvpsi_gamma_nonlocal(ia, dvg_exc_tmp(:,:,iks), dvpsi)
        !
        ! 2) forces_drhox1_i = < dvg | dvpsi_i >
        !
        DO ipol = 1,3
           !
           reduce = 0._DP
           !
           !$acc parallel loop collapse(2) reduction(+:reduce) present(dvg_exc_tmp,dvpsi) copy(reduce)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 reduce = reduce + REAL(dvg_exc_tmp(ig,lbnd,iks),KIND=DP)*REAL(dvpsi(ig,lbnd,ipol),KIND=DP) &
                 &               + AIMAG(dvg_exc_tmp(ig,lbnd,iks))*AIMAG(dvpsi(ig,lbnd,ipol))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           reduce = 2._DP*reduce
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:reduce) present(dvg_exc_tmp,dvpsi) copy(reduce)
              DO lbnd = 1,nbnd_do
                 reduce = reduce - REAL(dvg_exc_tmp(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,lbnd,ipol),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
           forces_drhox1(3*ia-3+ipol) = forces_drhox1(3*ia-3+ipol) + this_wk*reduce
           !
        ENDDO
        !
        CALL update_bar_type(barra,'f_drhox1',1)
        !
     ENDDO
     !
  ENDDO
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
     rdrhox1(:,1) = rdrhox1(:,1)+rdrhox1(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhox1(:,1), gstart, &
  & gamma_only, vloc, forcelc)
  !
  forcelc(:,:) = -factor*forcelc
  !
  DO ia = 1,nat
     DO ipol = 1,3
        forces_drhox1(3*ia-3+ipol) = forces_drhox1(3*ia-3+ipol) + forcelc(ipol,ia)
     ENDDO
  ENDDO
  !
  forces(:) = forces+forces_drhox1
  !
  CALL io_push_title('Forces drhox1')
  !
  DO ia = 1,nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) ia, ityp(ia), (-forces_drhox1(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.forces.forces_drhox1', -forces_drhox1(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_forces_gpu()
#endif
  !
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
SUBROUTINE wbse_calc_dvgdvg_mat(dvg_exc_tmp, dvgdvg_mat)
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
        DO ibnd = 1,nbndval-n_trunc_bands
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1,npw
              reduce = reduce &
              & + REAL(evc1_all(ig,ibnd,iks),KIND=DP)*REAL(dvg_exc_tmp(ig,lbnd,iks),KIND=DP) &
              & + AIMAG(evc1_all(ig,ibnd,iks))*AIMAG(dvg_exc_tmp(ig,lbnd,iks))
           ENDDO
           !
           dvgdvg_mat(ibnd,lbnd,iks) = 2._DP*reduce
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gstart == 2) THEN
        !$acc parallel loop collapse(2) present(dvgdvg_mat,evc1_all,dvg_exc_tmp)
        DO lbnd = 1,nbnd_do
           DO ibnd = 1,nbndval - n_trunc_bands
              dvgdvg_mat(ibnd,lbnd,iks) = dvgdvg_mat(ibnd,lbnd,iks) &
              & - REAL(evc1_all(1,ibnd,iks),KIND=DP)*REAL(dvg_exc_tmp(1,lbnd,iks),KIND=DP)
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
SUBROUTINE wbse_calc_drhox2(dvgdvg_mat, drhox2)
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
  USE wavefunctions,        ONLY : evc,psic
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
  ALLOCATE(aux_r(dffts%nnr))
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
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     DO lbnd = 1,nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        !
        w1 = wg(ibnd,iks_do)/omega
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
        !
        !$acc parallel loop present(aux_r)
        DO ir = 1,dffts_nnr
           aux_r(ir) = REAL(psic(ir),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        DO jbnd = 1,nbndval-n_trunc_bands,2
           !
           jbndp = jbnd + n_trunc_bands
           !
           IF(jbnd < nbndval-n_trunc_bands) THEN
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),evc(:,jbndp+1),psic,'Wave')
              !
              !$acc parallel loop present(aux_r,dvgdvg_mat,drhox2)
              DO ir = 1,dffts_nnr
                 prod = aux_r(ir) * (REAL(psic(ir),KIND=DP)*dvgdvg_mat(jbnd,lbnd,iks) &
                 &                + AIMAG(psic(ir))*dvgdvg_mat(jbnd+1,lbnd,iks))
                 drhox2(ir,current_spin) = drhox2(ir,current_spin) - w1*CMPLX(prod,KIND=DP)
              ENDDO
              !$acc end parallel
              !
           ELSE
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),psic,'Wave')
              !
              !$acc parallel loop present(aux_r,dvgdvg_mat,drhox2)
              DO ir = 1,dffts_nnr
                 prod = aux_r(ir) * REAL(psic(ir),KIND=DP) * dvgdvg_mat(jbnd,lbnd,iks)
                 drhox2(ir,current_spin) = drhox2(ir,current_spin) - w1*CMPLX(prod,KIND=DP)
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
  DEALLOCATE(aux_r)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhox2(n, dvgdvg_mat, drhox2, forces)
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
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_forces_gpu,deallocate_forces_gpu
  USE cublas
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
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:,:), aux1(:,:), aux2(:,:)
  !$acc declare device_resident(dvpsi,aux1,aux2)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, lbnd, ibnd, ig
  INTEGER :: band_group_myoffset
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_drhox2(:), rdrhox2(:,:), forcelc(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute forces of drhox2')
  !
  band_group_myoffset = band_group%myoffset
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_forces_gpu()
#endif
  !
  ALLOCATE(forces_drhox2(n))
  ALLOCATE(forcelc(3, nat))
  ALLOCATE(dvpsi(npwx, band_group%nlocx, 3))
  ALLOCATE(rdrhox2(dffts%nnr, nspin))
  ALLOCATE(aux1(npwx, band_group%nlocx))
  ALLOCATE(aux2(npwx, band_group%nlocx))
  !
  forces_drhox2(:) = 0._DP
  !
  CALL start_bar_type(barra,'f_drhox2',kpt_pool%nloc*nat)
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
     this_wk = wk(iks)*factor
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     !$acc parallel loop collapse(2) present(aux1,evc)
     DO lbnd = 1,nbnd_do
        !
        ! ibnd = band_group%l2g(lbnd)+n_trunc_bands
        !
        DO ig = 1,npw
           ibnd = band_group_myoffset+lbnd+n_trunc_bands
           aux1(ig,lbnd) = evc(ig,ibnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     DO ia = 1,nat
        !
        ! 1) | dvpsi_i >
        !
        CALL wbse_get_dvpsi_gamma_nonlocal(ia, aux1, dvpsi)
        !
        ! 2) forces_drhox2 = < evc_iv2 | dvpsi_ia_iv >
        !
        !$acc host_data use_device(evc,dvgdvg_mat,aux2)
        CALL DGEMM('N', 'N', 2*npw, band_group%nloc, nbndval-n_trunc_bands, 1._DP, &
        & evc(1,n_trunc_bands+1), 2*npwx, dvgdvg_mat(1,1,iks), nbndval0x-n_trunc_bands, &
        & 0._DP, aux2, 2*npwx)
        !$acc end host_data
        !
        DO ipol = 1,3
           !
           reduce = 0._DP
           !
           !$acc parallel loop collapse(2) reduction(+:reduce) present(aux2,dvpsi) copy(reduce)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 reduce = reduce + REAL(aux2(ig,lbnd),KIND=DP)*REAL(dvpsi(ig,lbnd,ipol),KIND=DP) &
                 &               + AIMAG(aux2(ig,lbnd))*AIMAG(dvpsi(ig,lbnd,ipol))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           reduce = 2._DP*reduce
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:reduce) present(aux2,dvpsi) copy(reduce)
              DO lbnd = 1,nbnd_do
                 reduce = reduce - REAL(aux2(1,lbnd),KIND=DP)*REAL(dvpsi(1,lbnd,ipol),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
           forces_drhox2(3*ia-3+ipol) = forces_drhox2(3*ia-3+ipol) - this_wk*reduce
           !
        ENDDO
        !
        CALL update_bar_type(barra,'f_drhox2',1)
        !
     ENDDO
     !
  ENDDO
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
     rdrhox2(:,1) = rdrhox2(:,1)+rdrhox2(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhox2(:,1), gstart, &
  & gamma_only, vloc, forcelc)
  !
  forcelc(:,:) = -factor*forcelc
  !
  DO ia = 1,nat
     DO ipol = 1,3
        forces_drhox2(3*ia-3+ipol) = forces_drhox2(3*ia-3+ipol)+forcelc(ipol,ia)
     ENDDO
  ENDDO
  !
  forces(:) = forces+forces_drhox2
  !
  CALL io_push_title('Forces drhox2')
  !
  DO ia = 1,nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) ia, ityp(ia), (-forces_drhox2(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.forces.forces_drhox2', -forces_drhox2(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_forces_gpu()
#endif
  !
  DEALLOCATE(forces_drhox2)
  DEALLOCATE(forcelc)
  DEALLOCATE(dvpsi)
  DEALLOCATE(rdrhox2)
  DEALLOCATE(aux1)
  DEALLOCATE(aux2)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_forces_drhoz(n, zvector, forces)
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
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_forces_gpu,deallocate_forces_gpu
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
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:,:), aux1(:,:), drhoz(:,:)
  !$acc declare device_resident(dvpsi,aux1)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, lbnd, ibnd, ig
  INTEGER :: band_group_myoffset
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: forces_drhoz(:), forcelc(:,:), rdrhoz(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  !
  CALL io_push_title('Compute forces of Z vector')
  !
  band_group_myoffset = band_group%myoffset
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_forces_gpu()
#endif
  !
  ALLOCATE(forces_drhoz(n))
  ALLOCATE(forcelc(3, nat))
  ALLOCATE(dvpsi(npwx, band_group%nlocx, 3))
  ALLOCATE(rdrhoz(dffts%nnr, nspin))
  ALLOCATE(drhoz(dffts%nnr, nspin))
  !$acc enter data create(drhoz)
  ALLOCATE(aux1(npwx, band_group%nlocx))
  !
  forces_drhoz(:) = 0._DP
  !
  CALL start_bar_type(barra,'f_drhoxz',kpt_pool%nloc*nat)
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
     this_wk = wk(iks)*factor
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     !$acc parallel loop collapse(2) present(aux1,evc)
     DO lbnd = 1,nbnd_do
        !
        ! ibnd = band_group%l2g(lbnd)+n_trunc_bands
        !
        DO ig = 1,npw
           ibnd = band_group_myoffset+lbnd+n_trunc_bands
           aux1(ig,lbnd) = evc(ig,ibnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     DO ia = 1,nat
        !
        ! 1) | dvpsi_i >
        !
        CALL wbse_get_dvpsi_gamma_nonlocal(ia, aux1, dvpsi)
        !
        ! 2) forces_drhoz_i = < z_vector | dvpsi_i >
        !
        DO ipol = 1,3
           !
           reduce = 0._DP
           !
           !$acc parallel loop collapse(2) reduction(+:reduce) present(zvector,dvpsi) copy(reduce)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 reduce = reduce &
                 & + REAL(zvector(ig,lbnd,iks),KIND=DP)*REAL(dvpsi(ig,lbnd,ipol),KIND=DP) &
                 & + AIMAG(zvector(ig,lbnd,iks))*AIMAG(dvpsi(ig,lbnd,ipol))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           reduce = 2._DP*reduce
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:reduce) present(zvector,dvpsi) copy(reduce)
              DO lbnd = 1,nbnd_do
                 reduce = reduce - REAL(zvector(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,lbnd,ipol),KIND=DP)
              ENDDO
             !$acc end parallel
           ENDIF
           !
           forces_drhoz(3*ia-3+ipol) = forces_drhoz(3*ia-3+ipol) + 2._DP*this_wk*reduce
           !
        ENDDO
        !
        CALL update_bar_type(barra,'f_drhoxz',1)
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(forces_drhoz,intra_bgrp_comm)
  CALL mp_sum(forces_drhoz,inter_bgrp_comm)
  CALL mp_sum(forces_drhoz,inter_pool_comm)
  !
  CALL stop_bar_type(barra,'f_drhoxz')
  !
  ! local part
  !
  CALL wbse_calc_dens(zvector, drhoz, .FALSE.)
  !
  drhoz(:,:) = 2._DP*drhoz
  rdrhoz(:,:) = REAL(drhoz,KIND=DP)
  !
  IF(nspin == 2) THEN
     rdrhoz(:,1) = rdrhoz(:,1)+rdrhoz(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhoz(:,1), gstart, &
  & gamma_only, vloc, forcelc)
  !
  forcelc(:,:) = -factor*forcelc
  !
  DO ia = 1,nat
     DO ipol = 1,3
        forces_drhoz(3*ia-3+ipol) = forces_drhoz(3*ia-3+ipol) + forcelc(ipol,ia)
     ENDDO
  ENDDO
  !
  forces(:) = forces+forces_drhoz
  !
  CALL io_push_title('Forces drhoz')
  !
  DO ia = 1,nat
     !
     ! forces = - gradients
     !
     WRITE(stdout, 9035) ia, ityp(ia), (-forces_drhoz(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.forces.forces_drhoz', -forces_drhoz(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_forces_gpu()
#endif
  !
  DEALLOCATE(forces_drhoz)
  DEALLOCATE(forcelc)
  DEALLOCATE(dvpsi)
  DEALLOCATE(rdrhoz)
  !$acc exit data delete(drhoz)
  DEALLOCATE(drhoz)
  DEALLOCATE(aux1)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_get_dvpsi_gamma_nonlocal(i_at, dvg_tmp, dvpsi)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp,ntyp=>nsp
  USE cell_base,            ONLY : tpiba
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE gvect,                ONLY : g,gstart
  USE noncollin_module,     ONLY : npol
  USE uspp_param,           ONLY : nh
  USE uspp,                 ONLY : dvan,vkb
  USE pwcom,                ONLY : npw,npwx
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE distribution_center,  ONLY : band_group
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: i_at
  COMPLEX(DP), INTENT(IN) :: dvg_tmp(npwx*npol, band_group%nlocx)
  COMPLEX(DP), INTENT(OUT) :: dvpsi(npwx, band_group%nlocx, 3)
  !
  ! Workspace
  !
  INTEGER :: ia, ib, ig, nt, ih, jkb, ic, nh_nt
  INTEGER :: band_group_nloc
  COMPLEX(DP) :: factor
  REAL(DP), ALLOCATABLE :: bec1(:,:), bec2(:,:)
  COMPLEX(DP), ALLOCATABLE :: work(:,:)
  !$acc declare device_resident(bec1,bec2,work)
  !
  !$acc kernels present(dvpsi)
  dvpsi(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  band_group_nloc = band_group%nloc
  factor = tpiba*(0._DP,-1._DP)
  !
  jkb = 0
  DO nt = 1,ntyp
     nh_nt = nh(nt)
     DO ia = 1,nat
        IF(ityp(ia) == nt) THEN
           IF(ia == i_at) EXIT
           jkb = jkb+nh_nt
        ENDIF
     ENDDO
     IF(ia == i_at) EXIT
  ENDDO
  !
  IF(nh_nt < 1) RETURN
  !
  ALLOCATE(work(npwx,nh_nt))
  ALLOCATE(bec1(nh_nt,band_group%nlocx))
  ALLOCATE(bec2(nh_nt,band_group%nlocx))
  !
  DO ic = 1,3
     !
     ! first term: sum_l sum_G' [ i V_l(G) V^*_l(G') (G'*u) psi(G')
     !
     !$acc parallel loop collapse(2) present(work,vkb,g)
     DO ih = 1,nh_nt
        DO ig = 1,npw
           work(ig,ih) = vkb(ig,jkb+ih)*g(ic,ig)*factor
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(work,dvg_tmp,bec1)
     CALL DGEMM('C', 'N', nh_nt, band_group%nloc, 2*npw, 2._DP, work, 2*npwx, dvg_tmp, 2*npwx, &
     & 0._DP, bec1, nh_nt)
     !$acc end host_data
     !
     IF(gstart == 2) THEN
        !$acc parallel loop collapse(2) present(bec1,work,dvg_tmp)
        DO ib = 1,band_group_nloc
           DO ih = 1,nh_nt
              bec1(ih,ib) = bec1(ih,ib) - work(1,ih)*dvg_tmp(1,ib)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc host_data use_device(bec1)
     CALL mp_sum(bec1,intra_bgrp_comm)
     !$acc end host_data
     !
     !$acc parallel loop collapse(2) present(bec1,dvan)
     DO ib = 1,band_group_nloc
        DO ih = 1,nh_nt
           bec1(ih,ib) = dvan(ih,ih,nt)*bec1(ih,ib)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(vkb,bec1,dvpsi)
     CALL DGEMM('N', 'N', 2*npw, band_group%nloc, nh_nt, 1._DP, vkb(1,jkb+1), 2*npwx, bec1, nh_nt, &
     & 1._DP, dvpsi(1,1,ic), 2*npwx)
     !$acc end host_data
     !
     ! second term: sum_l sum_G' [-i (G*u) V_l(G) V^*_l(G') psi(G')
     !
     !$acc host_data use_device(vkb,dvg_tmp,bec2)
     CALL DGEMM('C', 'N', nh_nt, band_group%nloc, 2*npw, 2._DP, vkb(1,jkb+1), 2*npwx, dvg_tmp, &
     & 2*npwx, 0._DP, bec2, nh_nt)
     !$acc end host_data
     !
     IF(gstart == 2) THEN
        !$acc parallel loop collapse(2) present(bec2,vkb,dvg_tmp)
        DO ib = 1,band_group_nloc
           DO ih = 1,nh_nt
              bec2(ih,ib) = bec2(ih,ib) - vkb(1,jkb+ih)*dvg_tmp(1,ib)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc host_data use_device(bec2)
     CALL mp_sum(bec2,intra_bgrp_comm)
     !$acc end host_data
     !
     !$acc parallel loop collapse(2) present(bec2,dvan)
     DO ib = 1,band_group_nloc
        DO ih = 1,nh_nt
           bec2(ih,ib) = dvan(ih,ih,nt)*bec2(ih,ib)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(work,bec2,dvpsi)
     CALL DGEMM('N', 'N', 2*npw, band_group%nloc, nh_nt, 1._DP, work, 2*npwx, bec2, nh_nt, 1._DP, &
     & dvpsi(1,1,ic), 2*npwx)
     !$acc end host_data
     !
  ENDDO
  !
  DEALLOCATE(work)
  DEALLOCATE(bec1)
  DEALLOCATE(bec2)
  !
END SUBROUTINE
