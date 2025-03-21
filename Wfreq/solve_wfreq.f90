!
! Copyright (C) 2015-2025 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE solve_wfreq(l_read_restart,l_generate_plot,l_QDET)
  !-----------------------------------------------------------------------
  !
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart,l_generate_plot,l_QDET
  !
  IF( gamma_only ) THEN
     CALL solve_wfreq_gamma( l_read_restart,l_generate_plot,l_QDET )
  ELSE
     CALL solve_wfreq_k( l_read_restart,l_generate_plot )
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_wfreq_gamma(l_read_restart,l_generate_plot,l_QDET)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_lanczos,npwq,npwqx,l_macropol,l_frac_occ,&
                                 & occupation,nbnd_occ,nbnd_occ_full,n_bands,l_enable_lanczos,&
                                 & iuwfc,lrwfc,wfreq_eta,imfreq_list,refreq_list,wstat_save_dir,&
                                 & fftdriver,d_epsm1_ifr,z_epsm1_rfr,z_head_rfr,d_head_ifr,&
                                 & d_epsm1_ifr_a,d_head_ifr_a,z_epsm1_rfr_a,z_head_rfr_a,l_dc2025,&
                                 & d_epsm1_ifr_dc,d_head_ifr_dc,z_epsm1_rfr_dc,z_head_rfr_dc
  USE mp_global,            ONLY : inter_image_comm,my_image_id,nimage,inter_pool_comm,npool,&
                                 & inter_bgrp_comm,nbgrp,intra_bgrp_comm,me_bgrp
  USE mp,                   ONLY : mp_bcast,mp_sum
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE becmod,               ONLY : becp,allocate_bec_type_acc,deallocate_bec_type_acc
  USE uspp_init,            ONLY : init_us_2
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert,kpt_pool,band_group,macropert,ifr,rfr,occband
  USE class_idistribute,    ONLY : idistribute
  USE wfreq_io,             ONLY : write_wfreq,read_wfreq
  USE types_bz_grid,        ONLY : k_grid
  USE chi_invert,           ONLY : chi_invert_real,chi_invert_complex
  USE types_coulomb,        ONLY : pot3D
  USE uspp,                 ONLY : vkb,nkb
  USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu,reallocate_ps_gpu,allocate_gw_gpu,&
                                 & deallocate_gw_gpu,allocate_lanczos_gpu,deallocate_lanczos_gpu,&
                                 & allocate_chi_gpu,deallocate_chi_gpu,allocate_macropol_gpu,&
                                 & deallocate_macropol_gpu,ps_r
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart,l_generate_plot,l_QDET
  !
  ! Workspace
  !
  INTEGER :: ip,glob_ip,ig,ir,iv,ivloc,ivloc2,ifloc,iks,iks_g,ipol
  CHARACTER(LEN=25) :: filepot
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  INTEGER :: nbndval,nbndval_full
  INTEGER :: dffts_nnr,mypara_nloc,mypara_nglob,ifr_nloc,rfr_nloc
  REAL(DP),ALLOCATABLE :: subdiago(:,:),bnorm(:)
  REAL(DP),ALLOCATABLE :: diago(:,:),braket(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: braket
#endif
  COMPLEX(DP),ALLOCATABLE :: q_s(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dvpsi
#endif
  COMPLEX(DP),ALLOCATABLE :: phi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phis(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: phis
#endif
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: pertg_all
#endif
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
#if !defined(__CUDA)
  REAL(DP),ALLOCATABLE :: ps_r(:,:)
#endif
  TYPE(idistribute) :: mypara
  REAL(DP),ALLOCATABLE :: overlap(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: overlap
#endif
  REAL(DP) :: mwo,ecv,dfactor,frequency,dhead
  COMPLEX(DP) :: zmwo,zfactor,zm,zp,zhead
  INTEGER :: glob_jp,ic,ifreq,il
  INTEGER :: who
  REAL(DP),ALLOCATABLE :: dmatilda(:,:),dlambda(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dmatilda,dlambda
#endif
  COMPLEX(DP),ALLOCATABLE :: zmatilda(:,:),zlambda(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: zmatilda,zlambda
#endif
  REAL(DP),ALLOCATABLE :: dmati(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr(:,:,:)
  REAL(DP),ALLOCATABLE :: dmati_a(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr_a(:,:,:)
  REAL(DP),ALLOCATABLE :: dmati_freq0(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr_freq0(:,:)
  REAL(DP),ALLOCATABLE :: dmati_a_freq0(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr_a_freq0(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dmati,zmatr,dmati_a,zmatr_a
#endif
  REAL(DP) :: this_et,this_occ,docc
  INTEGER,ALLOCATABLE :: l2g(:)
  COMPLEX(DP),ALLOCATABLE :: evc_copy(:,:)
  !
  CALL io_push_title('(W)-Lanczos')
  !
  ! DISTRIBUTING...
  !
  mypara = idistribute()
  IF(l_macropol) THEN
     CALL macropert%copyto(mypara)
  ELSE
     CALL pert%copyto(mypara)
  ENDIF
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type_acc( becp )
  IF(l_macropol) THEN
     CALL allocate_bec_type_acc( nkb, MAX(mypara%nloc,3), becp ) ! I just need 2 becp at a time
  ELSE
     CALL allocate_bec_type_acc( nkb, mypara%nloc, becp ) ! I just need 2 becp at a time
  ENDIF
  !
  ! ALLOCATE dmati, zmatr, where chi0 is stored
  !
  ALLOCATE( dmati( mypara%nglob, mypara%nloc, ifr%nloc) )
  ALLOCATE( zmatr( mypara%nglob, mypara%nloc, rfr%nloc) )
  IF(.NOT. l_QDET) THEN
     !$acc enter data create(dmati,zmatr)
  ENDIF
  !
  IF(l_QDET) THEN
     ALLOCATE( evc_copy(npwx*npol, nbnd ) )
     ALLOCATE( dmati_a(mypara%nglob, mypara%nloc, ifr%nloc) )
     ALLOCATE( zmatr_a(mypara%nglob, mypara%nloc, rfr%nloc) )
     !$acc enter data create(evc_copy,dmati_a,zmatr_a)
     !
     !$acc kernels present(dmati_a)
     dmati_a(:,:,:) = 0._DP
     !$acc end kernels
     !
     !$acc kernels present(zmatr_a)
     zmatr_a(:,:,:) = 0._DP
     !$acc end kernels
  ENDIF
  !
  IF(l_read_restart) THEN
     CALL read_wfreq(dmati,zmatr,mypara%nglob,mypara%nloc)
     IF(.NOT. l_QDET) THEN
        !$acc update device(dmati,zmatr)
     ENDIF
  ELSE
     !$acc kernels present(dmati)
     dmati(:,:,:) = 0._DP
     !$acc end kernels
     !
     !$acc kernels present(zmatr)
     zmatr(:,:,:) = 0._DP
     !$acc end kernels
  ENDIF
  !
  barra_load = 0
  DO iks = 1,kpt_pool%nloc
     CALL band_group%init(nbnd_occ(iks),'b','band_group',.FALSE.)
     DO ivloc = 1,band_group%nloc
        barra_load = barra_load+1
     ENDDO
  ENDDO
  IF(l_read_restart .AND. .NOT. l_QDET) barra_load = 0
  !
  IF(barra_load == 0) THEN
     CALL start_bar_type ( barra, 'wlanczos', 1 )
     CALL update_bar_type( barra, 'wlanczos', 1 )
  ELSE
     CALL start_bar_type ( barra, 'wlanczos', barra_load )
  ENDIF
  !
  CALL pot3D%init('Wave',.FALSE.,'default')
  !
  !$acc enter data copyin(pot3D)
  !$acc enter data copyin(pot3D%sqvc)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
  CALL allocate_gw_gpu(mypara%nlocx,mypara%nloc)
#endif
  !
  dffts_nnr = dffts%nnr
  mypara_nloc = mypara%nloc
  mypara_nglob = mypara%nglob
  ifr_nloc = ifr%nloc
  rfr_nloc = rfr%nloc
  IF(l_frac_occ) THEN
     nbndval_full = MINVAL(nbnd_occ_full)
  ELSE
     nbndval_full = MINVAL(nbnd_occ)
  ENDIF
  !
  !$acc enter data copyin(imfreq_list,refreq_list)
  !
  ALLOCATE(dvpsi(npwx*npol,mypara%nlocx))
  ALLOCATE(overlap(mypara%nglob,nbnd-nbndval_full))
  ALLOCATE(pertr(dffts%nnr))
  ALLOCATE(pertg(npwqx))
  ALLOCATE(l2g(mypara%nloc))
  !$acc enter data create(dvpsi,overlap,pertr,pertg,l2g)
  !
  !$acc parallel loop present(l2g)
  DO ip = 1,mypara_nloc
     !
     ! l2g(ip) = mypara%l2g(ip)
     !
     l2g(ip) = nimage*(ip-1)+my_image_id+1
  ENDDO
  !$acc end parallel
  !
  ! Read PDEP
  !
  ALLOCATE(pertg_all(npwqx,mypara%nloc))
  pertg_all(:,:) = 0._DP
  !
  DO ip = 1,mypara%nloc
     glob_ip = mypara%l2g(ip)
     IF(glob_ip <= n_pdep_eigen_to_use) THEN
        CALL generate_pdep_fname(filepot,glob_ip)
        fname = TRIM(wstat_save_dir)//'/'//filepot
        CALL pdep_read_G_and_distribute(fname,pertg_all(:,ip))
     ENDIF
  ENDDO
  !
  ! LOOP
  !
  DO iks = 1,kpt_pool%nloc ! KPOINT-SPIN
     !
     ! Exit loop if no work to do
     !
     IF(barra_load == 0) EXIT
     !
     iks_g = kpt_pool%l2g(iks)
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
#if defined(__CUDA)
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb, .TRUE. )
#else
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb, .FALSE. )
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
     IF(l_QDET) THEN
        !
        !$acc kernels present(evc_copy,evc)
        evc_copy(:,:) = evc
        !$acc end kernels
        !
#if defined(__CUDA)
        CALL reallocate_ps_gpu(n_bands,nbnd)
#endif
        CALL apply_alpha_pa_to_m_wfcs(iks_g,nbnd,evc_copy,(1.0_DP,0.0_DP))
        !
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     IF(l_frac_occ) THEN
        nbndval_full = nbnd_occ_full(iks)
     ELSE
        nbndval_full = nbnd_occ(iks)
     ENDIF
     !
     mwo = -k_grid%weight(iks_g)/omega
     zmwo = CMPLX(mwo,KIND=DP)
     !
     ! Parallel macropol
     !
     IF(l_macropol) THEN
#if defined(__CUDA)
        CALL allocate_macropol_gpu(1)
        CALL reallocate_ps_gpu(nbndval_full,3)
#endif
        !
        ! PHI
        !
        CALL occband%init(nbndval,'i','occband',.FALSE.)
        !
        ALLOCATE(phis(npwx*npol,3,occband%nloc))
        !
        DO ivloc = 1,occband%nloc
           !
           iv = occband%l2g(ivloc)
           !
           IF(l_QDET) THEN
              CALL linsolve_commut_Hx(iks,nbndval_full,et(iv,iks),evc_copy(:,iv),phis(:,:,ivloc))
           ELSE
              CALL linsolve_commut_Hx(iks,nbndval_full,et(iv,iks),evc(:,iv),phis(:,:,ivloc))
           ENDIF
           !
        ENDDO
        !
#if defined(__CUDA)
        CALL deallocate_macropol_gpu()
#endif
        !
     ENDIF ! macropol
     !
     IF(l_enable_lanczos .AND. .NOT. l_QDET) THEN
#if defined(__CUDA)
        CALL allocate_lanczos_gpu(mypara%nloc)
#endif
        !
        ALLOCATE(bnorm(mypara%nloc))
        ALLOCATE(subdiago(n_lanczos-1,mypara%nloc))
        ALLOCATE(q_s(npwx*npol,mypara%nloc,n_lanczos))
        ALLOCATE(diago(n_lanczos,mypara%nloc))
        ALLOCATE(braket(mypara%nglob,n_lanczos,mypara%nloc))
        !$acc enter data create(q_s,diago,braket)
     ENDIF
     !
     CALL band_group%init(nbndval,'b','band_group',.FALSE.)
     !
     ! LOOP over band states
     !
     DO ivloc2 = 1,band_group%nloc
        !
        iv = band_group%l2g(ivloc2)
        !
        ! MACROPOL CASE
        !
        IF(l_macropol) THEN
           !
           ALLOCATE(phi(npwx*npol,3))
           phi(:,:) = 0._DP
           !
           DO ivloc = 1,occband%nloc
               IF(occband%l2g(ivloc) == iv) THEN
                  phi(:,:) = phis(:,:,ivloc)
               ENDIF
           ENDDO
           !
           CALL mp_sum(phi, inter_image_comm)
           !
        ENDIF
        !
        ! PSIC
        !
        IF(l_QDET) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc_copy(:,iv),psic,'Wave')
        ELSE
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,iv),psic,'Wave')
        ENDIF
        !
        !$acc kernels present(dvpsi)
        dvpsi(:,:) = 0._DP
        !$acc end kernels
        !
        DO ip = 1,mypara%nloc
           !
           glob_ip = mypara%l2g(ip)
           !
           ! Decide whether read dbs E or dhpi
           !
           IF(glob_ip <= n_pdep_eigen_to_use) THEN
              !
              pertg(:) = pertg_all(:,ip)
              !$acc update device(pertg)
              !
              ! Multiply by sqvc
              !
              !$acc parallel loop present(pertg,pot3D,pot3D%sqvc)
              DO ig = 1,npwq
                 pertg(ig) = pot3D%sqvc(ig)*pertg(ig)
              ENDDO
              !$acc end parallel
              !
              ! Bring it to R-space
              !
              CALL single_invfft_gamma(dffts,npwq,npwqx,pertg,pertr,TRIM(fftdriver))
              !
              !$acc parallel loop present(pertr)
              DO ir = 1,dffts_nnr
                 pertr(ir) = psic(ir)*pertr(ir)
              ENDDO
              !$acc end parallel
              !
              CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi(:,ip),'Wave')
              !
           ELSE
              !
              ipol = glob_ip-n_pdep_eigen_to_use
              !
              dvpsi(:,ip) = phi(:,ipol) * SQRT(fpi * e2)
              !$acc update device(dvpsi(:,ip))
              !
           ENDIF
           !
        ENDDO ! pert
        !
        IF(l_macropol) DEALLOCATE(phi)
        !
        IF(l_QDET) THEN
#if defined(__CUDA)
           CALL reallocate_ps_gpu(n_bands,mypara%nloc)
#endif
           CALL apply_alpha_pa_to_m_wfcs(iks_g,mypara%nloc,dvpsi,(1._DP,0._DP))
        ENDIF
        !
#if defined(__CUDA)
        CALL reallocate_ps_gpu(nbndval_full,mypara%nloc)
#endif
        CALL apply_alpha_pc_to_m_wfcs(nbndval_full,mypara%nloc,dvpsi,(1._DP,0._DP))
        !
        IF(nbnd > nbndval_full) THEN
           !
           ! OVERLAP( glob_ip, im=1:nbnd ) = < psi_im iks | dvpsi_glob_ip >
           !
#if defined(__CUDA)
           CALL reallocate_ps_gpu(nbnd-nbndval_full,mypara%nloc)
#else
           IF(ALLOCATED(ps_r)) DEALLOCATE(ps_r)
           ALLOCATE(ps_r(nbnd-nbndval_full,mypara%nloc))
#endif
           !
           IF(l_QDET) THEN
              CALL glbrak_gamma(evc_copy(:,nbndval_full+1:nbnd),dvpsi,ps_r,npw,npwx,&
              & nbnd-nbndval_full,mypara%nloc,nbnd-nbndval_full,npol)
           ELSE
              CALL glbrak_gamma(evc(:,nbndval_full+1:nbnd),dvpsi,ps_r,npw,npwx,&
              & nbnd-nbndval_full,mypara%nloc,nbnd-nbndval_full,npol)
           ENDIF
           !
           !$acc host_data use_device(ps_r)
           CALL mp_sum(ps_r,intra_bgrp_comm)
           !$acc end host_data
           !
           !$acc kernels present(overlap)
           overlap(:,1:nbnd-nbndval_full) = 0._DP
           !$acc end kernels
           !
           !$acc parallel loop collapse(2) present(overlap,l2g,ps_r)
           DO ic = 1,nbnd-nbndval_full
              DO ip = 1,mypara_nloc
                 overlap(l2g(ip),ic) = ps_r(ic,ip)
              ENDDO
           ENDDO
           !$acc end parallel
           !
#if !defined(__CUDA)
           DEALLOCATE(ps_r)
#endif
           !
           IF(nimage > 1) THEN
              !$acc update host(overlap)
              CALL mp_sum(overlap,inter_image_comm)
              !$acc update device(overlap)
           ENDIF
           !
           ! Update dmati with cond
           !
           DO ic = 1,nbnd-nbndval_full
              !
              ! Avoid double counting in frac occ case
              !
              IF(ic+nbndval_full <= iv) CYCLE
              !
              ecv = et(ic+nbndval_full,iks)-et(iv,iks)
              IF(l_frac_occ) docc = occupation(iv,iks)-occupation(ic+nbndval_full,iks)
              !
              IF(l_QDET) THEN
                 !$acc parallel loop collapse(3) present(imfreq_list,dmati_a,overlap,l2g)
                 DO ifreq = 1,ifr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = imfreq_list(ifreq)
                          dfactor = mwo*2._DP*ecv/(ecv**2+frequency**2)
                          IF(l_frac_occ) dfactor = dfactor*docc
                          dmati_a(glob_jp,ip,ifreq) = dmati_a(glob_jp,ip,ifreq) &
                          & +overlap(glob_jp,ic)*overlap(l2g(ip),ic)*dfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
              ELSE
                 !$acc parallel loop collapse(3) present(imfreq_list,dmati,overlap,l2g)
                 DO ifreq = 1,ifr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = imfreq_list(ifreq)
                          dfactor = mwo*2._DP*ecv/(ecv**2+frequency**2)
                          IF(l_frac_occ) dfactor = dfactor*docc
                          dmati(glob_jp,ip,ifreq) = dmati(glob_jp,ip,ifreq) &
                          & +overlap(glob_jp,ic)*overlap(l2g(ip),ic)*dfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
              ENDIF
              !
           ENDDO ! ic
           !
           ! Update zmatr with cond
           !
           DO ic = 1,nbnd-nbndval_full
              !
              ! Avoid double counting in frac occ case
              !
              IF(ic+nbndval_full <= iv) CYCLE
              !
              ecv = et(ic+nbndval_full,iks)-et(iv,iks)
              IF(l_frac_occ) docc = occupation(iv,iks)-occupation(ic+nbndval_full,iks)
              !
              IF(l_QDET) THEN
                 !$acc parallel loop collapse(3) present(refreq_list,zmatr_a,overlap,l2g)
                 DO ifreq = 1,rfr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = refreq_list(ifreq)
                          zp = CMPLX(ecv+frequency,-wfreq_eta,KIND=DP)
                          zm = CMPLX(ecv-frequency,-wfreq_eta,KIND=DP)
                          zfactor = zmwo/zp+zmwo/zm
                          IF(l_frac_occ) zfactor = zfactor*docc
                          zmatr_a(glob_jp,ip,ifreq) = zmatr_a(glob_jp,ip,ifreq) &
                          & +overlap(glob_jp,ic)*overlap(l2g(ip),ic)*zfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
              ELSE
                 !$acc parallel loop collapse(3) present(refreq_list,zmatr,overlap,l2g)
                 DO ifreq = 1,rfr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = refreq_list(ifreq)
                          zp = CMPLX(ecv+frequency,-wfreq_eta,KIND=DP)
                          zm = CMPLX(ecv-frequency,-wfreq_eta,KIND=DP)
                          zfactor = zmwo/zp+zmwo/zm
                          IF(l_frac_occ) zfactor = zfactor*docc
                          zmatr(glob_jp,ip,ifreq) = zmatr(glob_jp,ip,ifreq) &
                          & +overlap(glob_jp,ic)*overlap(l2g(ip),ic)*zfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
              ENDIF
              !
           ENDDO ! ic
           !
        ENDIF
        !
        IF(l_QDET) THEN
           CALL update_bar_type( barra, 'wlanczos', 1 )
           CYCLE
        ENDIF
        !
        ! Apply Pc, to be sure
        !
#if defined(__CUDA)
        CALL reallocate_ps_gpu(nbnd,mypara%nloc)
#endif
        CALL apply_alpha_pc_to_m_wfcs(nbnd,mypara%nloc,dvpsi,(1._DP,0._DP))
        !
        !$acc update host(dvpsi)
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF(l_enable_lanczos) THEN
           !
           CALL solve_deflated_lanczos_w_full_ortho(nbnd,mypara%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
           CALL get_brak_hyper_parallel(dvpsi,mypara%nloc,n_lanczos,q_s,braket,mypara)
           !
           DO ip = 1,mypara%nloc
              CALL diago_lanczos(bnorm(ip),diago(:,ip),subdiago(:,ip),braket(:,:,ip),mypara%nglob)
           ENDDO
           !
           !$acc update device(diago,braket)
           !
           this_et = et(iv,iks)
           this_occ = occupation(iv,iks)
           !
           ! Update dmati with lanczos
           !
           DO il = 1,n_lanczos
              !
              !$acc parallel loop collapse(3) present(imfreq_list,diago,dmati,braket)
              DO ifreq = 1,ifr_nloc
                 DO ip = 1,mypara_nloc
                    DO glob_jp = 1,mypara_nglob
                       frequency = imfreq_list(ifreq)
                       ecv = diago(il,ip)-this_et
                       dfactor = mwo*2._DP*ecv/(ecv**2+frequency**2)
                       IF(l_frac_occ) dfactor = dfactor*this_occ
                       dmati(glob_jp,ip,ifreq) = dmati(glob_jp,ip,ifreq)+braket(glob_jp,il,ip)*dfactor
                    ENDDO
                 ENDDO
              ENDDO ! ifreq
              !$acc end parallel
              !
           ENDDO ! il
           !
           ! Update zmatr with lanczos
           !
           DO il = 1,n_lanczos
              !
              !$acc parallel loop collapse(3) present(refreq_list,diago,zmatr,braket)
              DO ifreq = 1,rfr_nloc
                 DO ip = 1,mypara_nloc
                    DO glob_jp = 1,mypara_nglob
                       frequency = refreq_list(ifreq)
                       ecv = diago(il,ip)-this_et
                       zp = CMPLX(ecv+frequency,-wfreq_eta,KIND=DP)
                       zm = CMPLX(ecv-frequency,-wfreq_eta,KIND=DP)
                       zfactor = zmwo/zp+zmwo/zm
                       IF(l_frac_occ) zfactor = zfactor*this_occ
                       zmatr(glob_jp,ip,ifreq) = zmatr(glob_jp,ip,ifreq)+braket(glob_jp,il,ip)*zfactor
                    ENDDO
                 ENDDO
              ENDDO ! ifreq
              !$acc end parallel
              !
           ENDDO ! il
           !
        ENDIF ! l_enable_lanczos
        !
        CALL update_bar_type( barra, 'wlanczos', 1 )
        !
     ENDDO ! BANDS
     !
     IF(l_macropol) DEALLOCATE(phis)
     !
     IF(l_enable_lanczos .AND. .NOT. l_QDET) THEN
#if defined(__CUDA)
        CALL deallocate_lanczos_gpu()
#endif
        !
        DEALLOCATE(bnorm)
        DEALLOCATE(subdiago)
        !$acc exit data delete(q_s,diago,braket)
        DEALLOCATE(q_s)
        DEALLOCATE(diago)
        DEALLOCATE(braket)
     ENDIF
     !
  ENDDO ! KPOINT-SPIN
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
  CALL deallocate_gw_gpu()
#endif
  !
  IF(l_QDET) THEN
     !$acc exit data copyout(dmati_a,zmatr_a)
  ELSE
     !$acc exit data copyout(dmati,zmatr)
  ENDIF
  !$acc exit data delete(imfreq_list,refreq_list)
  !$acc exit data delete(dvpsi,overlap,pertr,pertg,l2g)
  DEALLOCATE(dvpsi)
  DEALLOCATE(overlap)
  DEALLOCATE(pertr)
  DEALLOCATE(pertg)
  DEALLOCATE(l2g)
  DEALLOCATE(pertg_all)
  !
  !$acc exit data delete(pot3D%sqvc)
  !$acc exit data delete(pot3D)
  !
  ! Synchronize and write data
  !
  IF(.NOT. l_read_restart .AND. .NOT. l_QDET) THEN
     IF(npool > 1) THEN
        CALL mp_sum(dmati,inter_pool_comm)
        CALL mp_sum(zmatr,inter_pool_comm)
     ENDIF
     IF(nbgrp > 1) THEN
        CALL mp_sum(dmati,inter_bgrp_comm)
        CALL mp_sum(zmatr,inter_bgrp_comm)
     ENDIF
     CALL write_wfreq(dmati,zmatr,mypara%nglob,mypara%nloc)
  ENDIF
  !
  IF(l_QDET) THEN
     IF(npool > 1) THEN
        CALL mp_sum(dmati_a,inter_pool_comm)
        CALL mp_sum(zmatr_a,inter_pool_comm)
     ENDIF
     IF(nbgrp > 1) THEN
        CALL mp_sum(dmati_a,inter_bgrp_comm)
        CALL mp_sum(zmatr_a,inter_bgrp_comm)
     ENDIF
  ENDIF
  !
  IF(l_QDET) THEN
     !$acc exit data delete(evc_copy)
     DEALLOCATE(evc_copy)
  ENDIF
  !
  CALL stop_bar_type( barra, 'wlanczos' )
  !
  CALL start_clock('chi_invert')
  !
  ! EPS-1 imfreq
  !
  ALLOCATE(dmatilda(mypara%nglob,mypara%nglob))
  ALLOCATE(dlambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
  !$acc enter data create(dmatilda,dlambda)
  IF(l_QDET) THEN
     ALLOCATE(d_epsm1_ifr_a(pert%nglob,pert%nloc,ifr%nloc))
     d_epsm1_ifr_a(:,:,:) = 0._DP
     IF(l_dc2025) THEN
        ALLOCATE(d_epsm1_ifr_dc(pert%nglob,pert%nloc,ifr%nloc))
        d_epsm1_ifr_dc(:,:,:) = 0._DP
        ALLOCATE(dmati_freq0(mypara%nglob,mypara%nloc))
        ALLOCATE(dmati_a_freq0(mypara%nglob,mypara%nloc))
     ENDIF
  ELSE
     ALLOCATE(d_epsm1_ifr(pert%nglob,pert%nloc,ifr%nloc))
     d_epsm1_ifr(:,:,:) = 0._DP
  ENDIF
  IF(l_macropol) THEN
     IF(l_QDET) THEN
        ALLOCATE(d_head_ifr_a(ifr%nloc))
        d_head_ifr_a(:) = 0._DP
        IF(l_dc2025) THEN
           ALLOCATE(d_head_ifr_dc(ifr%nloc))
           d_head_ifr_dc(:) = 0._DP
        ENDIF
     ELSE
        ALLOCATE(d_head_ifr(ifr%nloc))
        d_head_ifr(:) = 0._DP
     ENDIF
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_chi_gpu(.TRUE.)
#endif
  !
  CALL band_group%init(ifr%nloc,'b','band_group',.FALSE.)
  !
  IF(l_QDET .AND. l_dc2025) THEN
     !
     CALL ifr%g2l(1,ifloc,who)
     !
     IF(me_bgrp == who) THEN
        dmati_freq0(:,:) = dmati(:,:,1)
        dmati_a_freq0(:,:) = dmati_a(:,:,1)
     ENDIF
     !
     CALL mp_bcast(dmati_freq0,who,intra_bgrp_comm)
     CALL mp_bcast(dmati_a_freq0,who,intra_bgrp_comm)
     !
  ENDIF
  !
  DO ifloc = 1,band_group%nloc
     !
     ifreq = band_group%l2g(ifloc)
     !
     dmatilda(:,:) = 0._DP
     DO ip = 1,mypara%nloc
        glob_ip = mypara%l2g(ip)
        IF(l_QDET) THEN
           dmatilda(:,glob_ip) = dmati(:,ip,ifreq)-dmati_a(:,ip,ifreq)
        ELSE
           dmatilda(:,glob_ip) = dmati(:,ip,ifreq)
        ENDIF
     ENDDO
     !
     CALL mp_sum(dmatilda,inter_image_comm)
     CALL chi_invert_real(dmatilda,dhead,dlambda,mypara%nglob)
     !
     DO ip = 1,pert%nloc
        glob_ip = pert%l2g(ip)
        IF(l_QDET) THEN
           d_epsm1_ifr_a(1:n_pdep_eigen_to_use,ip,ifreq) = dlambda(1:n_pdep_eigen_to_use,glob_ip)
        ELSE
           d_epsm1_ifr(1:n_pdep_eigen_to_use,ip,ifreq) = dlambda(1:n_pdep_eigen_to_use,glob_ip)
        ENDIF
     ENDDO
     !
     IF(l_macropol) THEN
        IF(l_QDET) THEN
           d_head_ifr_a(ifreq) = dhead
        ELSE
           d_head_ifr(ifreq) = dhead
        ENDIF
     ENDIF
     !
     IF(l_QDET .AND. l_dc2025) THEN
        !
        ! Double counting
        !
        dmatilda(:,:) = 0._DP
        DO ip = 1,mypara%nloc
           glob_ip = mypara%l2g(ip)
           dmatilda(:,glob_ip) = dmati_freq0(:,ip)-dmati_a_freq0(:,ip)+dmati_a(:,ip,ifreq)
        ENDDO
        !
        CALL mp_sum(dmatilda,inter_image_comm)
        CALL chi_invert_real(dmatilda,dhead,dlambda,mypara%nglob)
        !
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           d_epsm1_ifr_dc(1:n_pdep_eigen_to_use,ip,ifreq) = dlambda(1:n_pdep_eigen_to_use,glob_ip)
        ENDDO
        !
        IF(l_macropol) d_head_ifr_dc(ifreq) = dhead
        !
     ENDIF
     !
  ENDDO
  !
#if defined(__CUDA)
  CALL deallocate_chi_gpu()
#endif
  !
  !$acc exit data delete(dmatilda,dlambda)
  DEALLOCATE(dlambda)
  DEALLOCATE(dmatilda)
  DEALLOCATE(dmati)
  IF(l_QDET) DEALLOCATE(dmati_a)
  IF(l_QDET .AND. l_dc2025) THEN
     DEALLOCATE(dmati_freq0)
     DEALLOCATE(dmati_a_freq0)
  ENDIF
  !
  IF(l_QDET) THEN
     CALL mp_sum(d_epsm1_ifr_a,inter_bgrp_comm)
     IF(l_macropol) CALL mp_sum(d_head_ifr_a,inter_bgrp_comm)
     IF(l_dc2025) THEN
        CALL mp_sum(d_epsm1_ifr_dc,inter_bgrp_comm)
        IF(l_macropol) CALL mp_sum(d_head_ifr_dc,inter_bgrp_comm)
     ENDIF
  ELSE
     CALL mp_sum(d_epsm1_ifr,inter_bgrp_comm)
     IF(l_macropol) CALL mp_sum(d_head_ifr,inter_bgrp_comm)
  ENDIF
  !
  ! EPS-1 refreq
  !
  ALLOCATE(zmatilda(mypara%nglob,mypara%nglob))
  ALLOCATE(zlambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
  !$acc enter data create(zmatilda,zlambda)
  IF(l_QDET) THEN
     ALLOCATE(z_epsm1_rfr_a(pert%nglob,pert%nloc,rfr%nloc))
     z_epsm1_rfr_a(:,:,:) = 0._DP
     IF(l_dc2025) THEN
        ALLOCATE(z_epsm1_rfr_dc(pert%nglob,pert%nloc,rfr%nloc))
        z_epsm1_rfr_dc(:,:,:) = 0._DP
        ALLOCATE(zmatr_freq0(mypara%nglob,mypara%nloc))
        ALLOCATE(zmatr_a_freq0(mypara%nglob,mypara%nloc))
     ENDIF
  ELSE
     ALLOCATE(z_epsm1_rfr(pert%nglob,pert%nloc,rfr%nloc))
     z_epsm1_rfr(:,:,:) = 0._DP
  ENDIF
  IF(l_macropol) THEN
     IF(l_QDET) THEN
        ALLOCATE(z_head_rfr_a(rfr%nloc))
        z_head_rfr_a(:) = 0._DP
        IF(l_dc2025) THEN
           ALLOCATE(z_head_rfr_dc(rfr%nloc))
           z_head_rfr_dc(:) = 0._DP
        ENDIF
     ELSE
        ALLOCATE(z_head_rfr(rfr%nloc))
        z_head_rfr(:) = 0._DP
     ENDIF
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_chi_gpu(.FALSE.)
#endif
  !
  CALL band_group%init(rfr%nloc,'b','band_group',.FALSE.)
  !
  IF(l_QDET .AND. l_dc2025) THEN
     !
     CALL rfr%g2l(1,ifloc,who)
     !
     IF(me_bgrp == who) THEN
        zmatr_freq0(:,:) = zmatr(:,:,1)
        zmatr_a_freq0(:,:) = zmatr_a(:,:,1)
     ENDIF
     !
     CALL mp_bcast(zmatr_freq0,who,intra_bgrp_comm)
     CALL mp_bcast(zmatr_a_freq0,who,intra_bgrp_comm)
     !
  ENDIF
  !
  DO ifloc = 1,band_group%nloc
     !
     ifreq = band_group%l2g(ifloc)
     !
     zmatilda(:,:) = 0._DP
     DO ip = 1,mypara%nloc
        glob_ip = mypara%l2g(ip)
        IF(l_QDET) THEN
           zmatilda(:,glob_ip) = zmatr(:,ip,ifreq)-zmatr_a(:,ip,ifreq)
        ELSE
           zmatilda(:,glob_ip) = zmatr(:,ip,ifreq)
        ENDIF
     ENDDO
     !
     CALL mp_sum(zmatilda,inter_image_comm)
     CALL chi_invert_complex(zmatilda,zhead,zlambda,mypara%nglob)
     !
     DO ip = 1,pert%nloc
        glob_ip = pert%l2g(ip)
        IF(l_QDET) THEN
           z_epsm1_rfr_a(1:n_pdep_eigen_to_use,ip,ifreq) = zlambda(1:n_pdep_eigen_to_use,glob_ip)
        ELSE
           z_epsm1_rfr(1:n_pdep_eigen_to_use,ip,ifreq) = zlambda(1:n_pdep_eigen_to_use,glob_ip)
        ENDIF
     ENDDO
     !
     IF(l_macropol) THEN
        IF(l_QDET) THEN
           z_head_rfr_a(ifreq) = zhead
        ELSE
           z_head_rfr(ifreq) = zhead
        ENDIF
     ENDIF
     !
     IF(l_QDET .AND. l_dc2025) THEN
        !
        ! Double counting
        !
        zmatilda(:,:) = 0._DP
        DO ip = 1,mypara%nloc
           glob_ip = mypara%l2g(ip)
           zmatilda(:,glob_ip) = zmatr_freq0(:,ip)-zmatr_a_freq0(:,ip)+zmatr_a(:,ip,ifreq)
        ENDDO
        !
        CALL mp_sum(zmatilda,inter_image_comm)
        CALL chi_invert_complex(zmatilda,zhead,zlambda,mypara%nglob)
        !
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           z_epsm1_rfr_dc(1:n_pdep_eigen_to_use,ip,ifreq) = zlambda(1:n_pdep_eigen_to_use,glob_ip)
        ENDDO
        !
        IF(l_macropol) z_head_rfr_dc(ifreq) = zhead
        !
     ENDIF
     !
  ENDDO
  !
#if defined(__CUDA)
  CALL deallocate_chi_gpu()
#endif
  !
  !$acc exit data delete(zmatilda,zlambda)
  DEALLOCATE(zlambda)
  DEALLOCATE(zmatilda)
  DEALLOCATE(zmatr)
  IF(l_QDET) DEALLOCATE(zmatr_a)
  IF(l_QDET .AND. l_dc2025) THEN
     DEALLOCATE(zmatr_freq0)
     DEALLOCATE(zmatr_a_freq0)
  ENDIF
  !
  IF(l_QDET) THEN
     CALL mp_sum(z_epsm1_rfr_a,inter_bgrp_comm)
     IF(l_macropol) CALL mp_sum(z_head_rfr_a,inter_bgrp_comm)
     IF(l_dc2025) THEN
        CALL mp_sum(z_epsm1_rfr_dc,inter_bgrp_comm)
        IF(l_macropol) CALL mp_sum(z_head_rfr_dc,inter_bgrp_comm)
     ENDIF
  ELSE
     CALL mp_sum(z_epsm1_rfr,inter_bgrp_comm)
     IF(l_macropol) CALL mp_sum(z_head_rfr,inter_bgrp_comm)
  ENDIF
  !
  CALL stop_clock('chi_invert')
  !
  IF(l_generate_plot) THEN
     CALL output_eps_head()
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_wfreq_k(l_read_restart,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_lanczos,npwq,npwqx,l_macropol,nbnd_occ,&
                                 & l_enable_lanczos,iuwfc,lrwfc,wfreq_eta,imfreq_list,refreq_list,&
                                 & wstat_save_dir,ngq,igq_q,z_epsm1_ifr_q,z_epsm1_rfr_q,z_head_rfr,&
                                 & z_head_ifr
  USE mp_global,            ONLY : my_image_id,inter_image_comm,nimage,inter_bgrp_comm,nbgrp,&
                                 & intra_bgrp_comm
  USE mp,                   ONLY : mp_bcast,mp_sum
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,current_spin,isk,nbnd,lsda,igk_k,current_k,ngk
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE becmod,               ONLY : becp,allocate_bec_type_acc,deallocate_bec_type_acc
  USE uspp_init,            ONLY : init_us_2
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert,macropert,ifr,rfr,occband,band_group
  USE class_idistribute,    ONLY : idistribute
  USE wfreq_io,             ONLY : write_wfreq,read_wfreq
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE chi_invert,           ONLY : chi_invert_complex
  USE types_coulomb,        ONLY : pot3D
  USE uspp,                 ONLY : vkb,nkb
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu,reallocate_ps_gpu,allocate_gw_gpu,&
                                 & deallocate_gw_gpu,allocate_lanczos_gpu,deallocate_lanczos_gpu,&
                                 & allocate_chi_gpu,deallocate_chi_gpu,allocate_macropol_gpu,&
                                 & deallocate_macropol_gpu,ps_c
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart,l_generate_plot
  !
  ! Workspace
  !
  INTEGER :: ip,glob_ip,ig,ir,iv,ivloc,ivloc2,ifloc,iks,ik,is,iq,ikqs,ikq,ipol
  CHARACTER(LEN=25) :: filepot
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  INTEGER :: nbndval
  INTEGER :: dffts_nnr,mypara_nloc,mypara_nglob,ifr_nloc,rfr_nloc
  REAL(DP),ALLOCATABLE :: subdiago(:,:),bnorm(:)
  REAL(DP),ALLOCATABLE :: diago(:,:)
  COMPLEX(DP),ALLOCATABLE :: braket(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: braket
#endif
  COMPLEX(DP),ALLOCATABLE :: q_s(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dvpsi
#endif
  COMPLEX(DP),ALLOCATABLE :: phi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phis(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: phis
#endif
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
  COMPLEX(DP),ALLOCATABLE :: evckpq(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: pertg_all,evckpq
#endif
  COMPLEX(DP),ALLOCATABLE :: phase(:)
  COMPLEX(DP),ALLOCATABLE :: psick(:),psick_nc(:,:)
  INTEGER :: npwkq
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
#if !defined(__CUDA)
  COMPLEX(DP),ALLOCATABLE :: ps_c(:,:)
#endif
  TYPE(idistribute) :: mypara
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: overlap
#endif
  REAL(DP) :: mwo,ecv,dfactor,frequency
  COMPLEX(DP) :: zmwo,zfactor,zm,zp,zhead
  INTEGER :: glob_jp,ic,ifreq,il
  COMPLEX(DP),ALLOCATABLE :: zmatilda(:,:),zlambda(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmati_q(:,:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr_q(:,:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: zmatilda,zlambda,zmati_q,zmatr_q
#endif
  LOGICAL :: l_gammaq
  REAL(DP) :: g0(3)
  REAL(DP) :: this_et
  INTEGER,ALLOCATABLE :: l2g(:)
  !
  CALL io_push_title('(W)-Lanczos')
  !
  ! DISTRIBUTING...
  !
  mypara = idistribute()
  IF(l_macropol) THEN
     CALL macropert%copyto(mypara)
  ELSE
     CALL pert%copyto(mypara)
  ENDIF
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type_acc( becp )
  IF(l_macropol) THEN
     CALL allocate_bec_type_acc( nkb, MAX(mypara%nloc,3), becp ) ! I just need 2 becp at a time
  ELSE
     CALL allocate_bec_type_acc( nkb, mypara%nloc, becp ) ! I just need 2 becp at a time
  ENDIF
  !
  ! ALLOCATE zmati_q, zmatr_q, where chi0 is stored
  !
  ALLOCATE( zmati_q( mypara%nglob, mypara%nloc, ifr%nloc, q_grid%np) )
  ALLOCATE( zmatr_q( mypara%nglob, mypara%nloc, rfr%nloc, q_grid%np) )
  !$acc enter data create(zmati_q,zmatr_q)
  !
  IF(noncolin) THEN
     ALLOCATE( psick_nc(dffts%nnr,npol) )
     !$acc enter data create(psick_nc)
  ELSE
     ALLOCATE( psick(dffts%nnr) )
     !$acc enter data create(psick)
  ENDIF
  ALLOCATE( phase(dffts%nnr) )
  ALLOCATE( evckpq(npwx*npol,nbnd) )
  !$acc enter data create(phase,evckpq)
  !
  IF(l_read_restart) THEN
     CALL read_wfreq(zmati_q,zmatr_q,mypara%nglob,mypara%nloc)
     !$acc update device(zmati_q,zmatr_q)
  ELSE
     !$acc kernels present(zmati_q)
     zmati_q(:,:,:,:) = 0._DP
     !$acc end kernels
     !
     !$acc kernels present(zmatr_q)
     zmatr_q(:,:,:,:) = 0._DP
     !$acc end kernels
  ENDIF
  !
  barra_load = 0
  DO iq = 1,q_grid%np
     DO iks = 1,k_grid%nps
        CALL band_group%init(nbnd_occ(iks),'b','band_group',.FALSE.)
        DO ivloc = 1,band_group%nloc
           barra_load = barra_load+1
        ENDDO
     ENDDO
  ENDDO
  IF(l_read_restart) barra_load = 0
  !
  IF(barra_load == 0) THEN
     CALL start_bar_type ( barra, 'wlanczos', 1 )
     CALL update_bar_type( barra, 'wlanczos', 1 )
  ELSE
     CALL start_bar_type ( barra, 'wlanczos', barra_load )
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_gpu()
  CALL allocate_gw_gpu(mypara%nlocx,mypara%nloc)
#endif
  !
  dffts_nnr = dffts%nnr
  mypara_nloc = mypara%nloc
  mypara_nglob = mypara%nglob
  ifr_nloc = ifr%nloc
  rfr_nloc = rfr%nloc
  nbndval = MINVAL(nbnd_occ)
  !
  !$acc enter data copyin(imfreq_list,refreq_list)
  !
  ALLOCATE(dvpsi(npwx*npol,mypara%nlocx))
  ALLOCATE(overlap(mypara%nglob,nbnd-nbndval))
  ALLOCATE(pertr(dffts%nnr))
  ALLOCATE(pertg(npwqx))
  ALLOCATE(l2g(mypara%nloc))
  !$acc enter data create(dvpsi,overlap,pertr,pertg,l2g)
  !
  !$acc parallel loop present(l2g)
  DO ip = 1,mypara_nloc
     !
     ! l2g(ip) = mypara%l2g(ip)
     !
     l2g(ip) = nimage*(ip-1)+my_image_id+1
  ENDDO
  !$acc end parallel
  !
  ALLOCATE(pertg_all(npwqx,mypara%nloc))
  !
  ! LOOP
  !
  DO iq = 1,q_grid%np   ! Q-POINT
     !
     ! Exit loop if no work to do
     !
     IF(barra_load == 0) EXIT
     !
     npwq = ngq(iq)
     l_gammaq = q_grid%l_pIsGamma(iq)
     !
     CALL pot3D%init('Wave',.TRUE.,'default',iq)
     !
     !$acc enter data copyin(pot3D)
     !$acc enter data copyin(pot3D%sqvc)
     !
     ! Read PDEP
     !
     pertg_all = 0._DP
     !
     DO ip = 1,mypara%nloc
        glob_ip = mypara%l2g(ip)
        IF(glob_ip <= n_pdep_eigen_to_use) THEN
           CALL generate_pdep_fname(filepot,glob_ip,iq)
           fname = TRIM(wstat_save_dir)//'/'//filepot
           CALL pdep_read_G_and_distribute(fname,pertg_all(:,ip),iq)
        ENDIF
     ENDDO
     !
     DO iks = 1,k_grid%nps   ! KPOINT-SPIN
        !
        ik = k_grid%ip(iks)
        is = k_grid%is(iks)
        !
        ! ... Set k-point, spin, kinetic energy, needed by Hpsi
        !
        current_k = iks
        IF ( lsda ) current_spin = isk(iks)
        call g2_kin( iks )
        !
        ! ... More stuff needed by the hamiltonian: nonlocal projectors
        !
#if defined(__CUDA)
        IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), k_grid%p_cart(1,ik), vkb, .TRUE. )
#else
        IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), k_grid%p_cart(1,ik), vkb, .FALSE. )
#endif
        npw = ngk(iks)
        !
        ! ... read in wavefunctions from the previous iteration
        !
        IF(k_grid%nps > 1) THEN
           IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc,0,inter_image_comm)
           !$acc update device(evc)
        ENDIF
        !
        CALL k_grid%find( k_grid%p_cart(:,ik) + q_grid%p_cart(:,iq), 'cart', ikq, g0 )
        ikqs = k_grid%ipis2ips(ikq,is)
        !
        npwkq = ngk(ikqs)
        !
        CALL compute_phase( g0, 'cart', phase )
        !
        ! Set wavefunctions at [k+q] in G space, for all bands,
        ! and store them in evckpq
        !
        IF(my_image_id == 0) CALL get_buffer( evckpq, lrwfc, iuwfc, ikqs )
        CALL mp_bcast( evckpq, 0, inter_image_comm )
        !
        !$acc update device(evckpq,phase)
        !
        nbndval = nbnd_occ(iks)
        !
        mwo = -k_grid%weight(iks)/omega
        zmwo = CMPLX(mwo,KIND=DP)
        !
        ! MACROPOL CASE
        !
        IF(l_macropol .AND. l_gammaq) THEN
#if defined(__CUDA)
           CALL allocate_macropol_gpu(1)
           CALL reallocate_ps_gpu(nbndval,3)
#endif
           !
           ! PHI
           !
           CALL occband%init(nbndval,'i','occband',.FALSE.)
           !
           ALLOCATE(phis(npwx*npol,3,occband%nloc))
           !
           DO ivloc = 1,occband%nloc
              !
              iv = occband%l2g(ivloc)
              !
              CALL linsolve_commut_Hx(iks,nbndval,et(iv,iks),evc(:,iv),phis(:,:,ivloc))
              !
           ENDDO
           !
#if defined(__CUDA)
           CALL deallocate_macropol_gpu()
#endif
           !
        ENDIF ! macropol
        !
        IF(l_enable_lanczos) THEN
#if defined(__CUDA)
           CALL allocate_lanczos_gpu(mypara%nloc)
#endif
           !
           ALLOCATE(bnorm(mypara%nloc))
           ALLOCATE(subdiago(n_lanczos-1,mypara%nloc))
           ALLOCATE(q_s(npwx*npol,mypara%nloc,n_lanczos))
           ALLOCATE(diago(n_lanczos,mypara%nloc))
           ALLOCATE(braket(mypara%nglob,n_lanczos,mypara%nloc))
           !$acc enter data create(q_s,diago,braket)
        ENDIF
        !
        CALL band_group%init(nbndval,'b','band_group',.FALSE.)
        !
        ! LOOP over band states
        !
        DO ivloc2 = 1,band_group%nloc
           !
           iv = band_group%l2g(ivloc2)
           !
           ! MACROPOL CASE
           !
           IF(l_macropol .AND. l_gammaq) THEN
              !
              ALLOCATE(phi(npwx*npol,3))
              phi = 0._DP
              !
              DO ivloc = 1,occband%nloc
                 IF(occband%l2g(ivloc) == iv) THEN
                    phi(:,:) = phis(:,:,ivloc)
                 ENDIF
              ENDDO
              !
              CALL mp_sum(phi, inter_image_comm)
              !
           ENDIF
           !
           ! PSIC
           !
           IF(noncolin) THEN
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(1:npwx,iv),psick_nc(:,1),'Wave',igk_k(:,ikqs))
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(npwx+1:npwx*2,iv),psick_nc(:,2),'Wave',igk_k(:,ikqs))
           ELSE
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(:,iv),psick,'Wave',igk_k(:,ikqs))
           ENDIF
           !
           !$acc kernels present(dvpsi)
           dvpsi(:,:) = 0._DP
           !$acc end kernels
           !
           DO ip = 1,mypara%nloc
              !
              glob_ip = mypara%l2g(ip)
              !
              ! Decide whether read dbs E or dhpi
              !
              IF(glob_ip <= n_pdep_eigen_to_use) THEN
                 !
                 pertg(:) = pertg_all(:,ip)
                 !$acc update device(pertg)
                 !
                 ! Multiply by sqvc
                 !
                 !$acc parallel loop present(pertg,pot3D,pot3D%sqvc)
                 DO ig = 1,npwq
                    pertg(ig) = pot3D%sqvc(ig)*pertg(ig)
                 ENDDO
                 !$acc end parallel
                 !
                 ! Bring it to R-space
                 !
                 IF(noncolin) THEN
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                    !$acc parallel loop present(pertr,phase,psick_nc)
                    DO ir = 1,dffts_nnr
                       pertr(ir) = phase(ir)*psick_nc(ir,1)*CONJG(pertr(ir))
                    ENDDO
                    !$acc end parallel
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1:npwx,ip),'Wave',igk_k(:,current_k))
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                    !$acc parallel loop present(pertr,phase,psick_nc)
                    DO ir = 1,dffts_nnr
                       pertr(ir) = phase(ir)*psick_nc(ir,2)*CONJG(pertr(ir))
                    ENDDO
                    !$acc end parallel
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(npwx+1:npwx*2,ip),'Wave',igk_k(:,current_k))
                 ELSE
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                    !$acc parallel loop present(pertr,phase,psick)
                    DO ir = 1,dffts_nnr
                       pertr(ir) = phase(ir)*psick(ir)*CONJG(pertr(ir))
                    ENDDO
                    !$acc end parallel
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(:,ip),'Wave',igk_k(:,current_k))
                 ENDIF
                 !
              ELSE
                 !
                 IF(l_gammaq) THEN
                    ipol = glob_ip-n_pdep_eigen_to_use
                    dvpsi(:,ip) = phi(:,ipol) * SQRT(fpi * e2)
                    !$acc update device(dvpsi(:,ip))
                 ENDIF
                 !
              ENDIF
              !
           ENDDO ! pert
           !
           IF(l_macropol .AND. l_gammaq) DEALLOCATE(phi)
           !
#if defined(__CUDA)
           CALL reallocate_ps_gpu(nbndval,mypara%nloc)
#endif
           CALL apply_alpha_pc_to_m_wfcs(nbndval,mypara%nloc,dvpsi,(1._DP,0._DP))
           !
           IF(nbnd > nbndval) THEN
              !
              ! OVERLAP( glob_ip, im=1:nbnd ) = < psi_im iks | dvpsi_glob_ip >
              !
#if defined(__CUDA)
              CALL reallocate_ps_gpu(nbnd-nbndval,mypara%nloc)
#else
              IF(ALLOCATED(ps_c)) DEALLOCATE(ps_c)
              ALLOCATE(ps_c(nbnd-nbndval,mypara%nloc))
#endif
              !
              CALL glbrak_k(evc(:,nbndval+1:nbnd),dvpsi,ps_c,npw,npwx,nbnd-nbndval,&
              & mypara%nloc,nbnd-nbndval,npol)
              !
              !$acc host_data use_device(ps_c)
              CALL mp_sum(ps_c,intra_bgrp_comm)
              !$acc end host_data
              !
              !$acc kernels present(overlap)
              overlap(:,1:nbnd-nbndval) = 0._DP
              !$acc end kernels
              !
              !$acc parallel loop collapse(2) present(overlap,l2g,ps_c)
              DO ic = 1,nbnd-nbndval
                 DO ip = 1,mypara_nloc
                    overlap(l2g(ip),ic) = ps_c(ic,ip)
                 ENDDO
              ENDDO
              !$acc end parallel
              !
#if !defined(__CUDA)
              DEALLOCATE(ps_c)
#endif
              !
              IF(nimage > 1) THEN
                 !$acc update host(overlap)
                 CALL mp_sum(overlap,inter_image_comm)
                 !$acc update device(overlap)
              ENDIF
              !
              ! Update zmati with cond
              !
              DO ic = 1,nbnd-nbndval
                 !
                 ecv = et(ic+nbndval,iks)-et(iv,ikqs)
                 !
                 !$acc parallel loop collapse(3) present(imfreq_list,zmati_q,overlap,l2g)
                 DO ifreq = 1,ifr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = imfreq_list(ifreq)
                          dfactor = mwo*2._DP*ecv/(ecv**2+frequency**2)
                          zmati_q(glob_jp,ip,ifreq,iq) = zmati_q(glob_jp,ip,ifreq,iq) &
                          & +CONJG(overlap(l2g(ip),ic))*overlap(glob_jp,ic)*dfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ENDDO ! ic
              !
              ! Update zmatr with cond
              !
              DO ic = 1,nbnd-nbndval
                 !
                 ecv = et(ic+nbndval,iks)-et(iv,ikqs)
                 !
                 !$acc parallel loop collapse(3) present(refreq_list,zmatr_q,overlap,l2g)
                 DO ifreq = 1,rfr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = refreq_list(ifreq)
                          zp = CMPLX(ecv+frequency,-wfreq_eta,KIND=DP)
                          zm = CMPLX(ecv-frequency,-wfreq_eta,KIND=DP)
                          zfactor = zmwo/zp+zmwo/zm
                          zmatr_q(glob_jp,ip,ifreq,iq) = zmatr_q(glob_jp,ip,ifreq,iq) &
                          & +CONJG(overlap(l2g(ip),ic))*overlap(glob_jp,ic)*zfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ENDDO ! ic
              !
           ENDIF
           !
           ! Apply Pc, to be sure
           !
#if defined(__CUDA)
           CALL reallocate_ps_gpu(nbnd,mypara%nloc)
#endif
           CALL apply_alpha_pc_to_m_wfcs(nbnd,mypara%nloc,dvpsi,(1._DP,0._DP))
           !
           !$acc update host(dvpsi)
           !
           ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
           ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
           !
           IF(l_enable_lanczos) THEN
              !
              CALL solve_deflated_lanczos_w_full_ortho(nbnd,mypara%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
              CALL get_brak_hyper_parallel_complex(dvpsi,mypara%nloc,n_lanczos,q_s,braket,mypara)
              !
              DO ip = 1,mypara%nloc
                 CALL diago_lanczos_complex(bnorm(ip),diago(:,ip),subdiago(:,ip),braket(:,:,ip),mypara%nglob)
              ENDDO
              !
              !$acc update device(diago,braket)
              !
              this_et = et(iv,ikqs)
              !
              ! Update zmati with lanczos
              !
              DO il = 1,n_lanczos
                 !
                 !$acc parallel loop collapse(3) present(imfreq_list,diago,zmati_q,braket)
                 DO ifreq = 1,ifr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = imfreq_list(ifreq)
                          ecv = diago(il,ip)-this_et
                          dfactor = mwo*2._DP*ecv/(ecv**2+frequency**2)
                          zmati_q(glob_jp,ip,ifreq,iq) = zmati_q(glob_jp,ip,ifreq,iq) &
                          & +CONJG(braket(glob_jp,il,ip))*dfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ENDDO ! il
              !
              ! Update zmatr with lanczos
              !
              DO il = 1,n_lanczos
                 !
                 !$acc parallel loop collapse(3) present(refreq_list,diago,zmatr_q,braket)
                 DO ifreq = 1,rfr_nloc
                    DO ip = 1,mypara_nloc
                       DO glob_jp = 1,mypara_nglob
                          frequency = refreq_list(ifreq)
                          ecv = diago(il,ip)-this_et
                          zp = CMPLX(ecv+frequency,-wfreq_eta,KIND=DP)
                          zm = CMPLX(ecv-frequency,-wfreq_eta,KIND=DP)
                          zfactor = zmwo/zp + zmwo/zm
                          zmatr_q(glob_jp,ip,ifreq,iq) = zmatr_q(glob_jp,ip,ifreq,iq) &
                          & +CONJG(braket(glob_jp,il,ip))*zfactor
                       ENDDO
                    ENDDO
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ENDDO ! il
              !
           ENDIF ! l_enable_lanczos
           !
           CALL update_bar_type( barra, 'wlanczos', 1 )
           !
        ENDDO ! BANDS
        !
        IF(l_macropol .AND. l_gammaq) DEALLOCATE(phis)
        !
        IF(l_enable_lanczos) THEN
#if defined(__CUDA)
           CALL deallocate_lanczos_gpu()
#endif
           !
           DEALLOCATE(bnorm)
           DEALLOCATE(subdiago)
           !$acc exit data delete(q_s,diago,braket)
           DEALLOCATE(q_s)
           DEALLOCATE(diago)
           DEALLOCATE(braket)
        ENDIF
        !
     ENDDO ! KPOINT-SPIN
     !
     !$acc exit data delete(pot3D%sqvc)
     !$acc exit data delete(pot3D)
     !
  ENDDO ! QPOINT
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
  CALL deallocate_gw_gpu()
#endif
  !
  IF(noncolin) THEN
     !$acc exit data delete(psick_nc)
     DEALLOCATE(psick_nc)
  ELSE
     !$acc exit data delete(psick)
     DEALLOCATE(psick)
  ENDIF
  !$acc exit data delete(phase,evckpq)
  DEALLOCATE(phase)
  DEALLOCATE(evckpq)
  !
  !$acc exit data copyout(zmati_q,zmatr_q)
  !$acc exit data delete(imfreq_list,refreq_list)
  !$acc exit data delete(dvpsi,overlap,pertr,pertg,l2g)
  DEALLOCATE(dvpsi)
  DEALLOCATE(overlap)
  DEALLOCATE(pertr)
  DEALLOCATE(pertg)
  DEALLOCATE(l2g)
  DEALLOCATE(pertg_all)
  !
  ! Synchronize and write data
  !
  IF(.NOT. l_read_restart) THEN
     IF(nbgrp > 1) THEN
        CALL mp_sum(zmati_q,inter_bgrp_comm)
        CALL mp_sum(zmatr_q,inter_bgrp_comm)
     ENDIF
     CALL write_wfreq(zmati_q,zmatr_q,mypara%nglob,mypara%nloc)
  ENDIF
  !
  CALL stop_bar_type( barra, 'wlanczos' )
  !
  CALL start_clock('chi_invert')
  !
  ! EPS-1 imfreq
  !
  ALLOCATE(zmatilda(mypara%nglob,mypara%nglob))
  ALLOCATE(zlambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
  !$acc enter data create(zmatilda,zlambda)
  ALLOCATE(z_epsm1_ifr_q(pert%nglob,pert%nloc,ifr%nloc,q_grid%np))
  z_epsm1_ifr_q = 0._DP
  IF(l_macropol) THEN
     ALLOCATE(z_head_ifr(ifr%nloc))
     z_head_ifr = 0._DP
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_chi_gpu(.FALSE.)
#endif
  !
  CALL band_group%init(ifr%nloc,'b','band_group',.FALSE.)
  !
  DO iq = 1,q_grid%np
     !
     l_gammaq = q_grid%l_pIsGamma(iq)
     !
     DO ifloc = 1,band_group%nloc
        !
        ifreq = band_group%l2g(ifloc)
        !
        zmatilda = 0._DP
        DO ip = 1,mypara%nloc
           glob_ip = mypara%l2g(ip)
           zmatilda(:,glob_ip) = zmati_q(:,ip,ifreq,iq)
        ENDDO
        !
        CALL mp_sum(zmatilda,inter_image_comm)
        CALL chi_invert_complex(zmatilda,zhead,zlambda,mypara%nglob,l_gammaq)
        !
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           z_epsm1_ifr_q(1:n_pdep_eigen_to_use,ip,ifreq,iq) = zlambda(1:n_pdep_eigen_to_use,glob_ip)
        ENDDO
        IF(l_macropol .AND. l_gammaq) z_head_ifr(ifreq) = zhead
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(zmati_q)
  !
  CALL mp_sum(z_epsm1_ifr_q,inter_bgrp_comm)
  IF(l_macropol) CALL mp_sum(z_head_ifr,inter_bgrp_comm)
  !
  ! EPS-1 refreq
  !
  ALLOCATE(z_epsm1_rfr_q(pert%nglob,pert%nloc,rfr%nloc,q_grid%np))
  z_epsm1_rfr_q = 0._DP
  IF(l_macropol) THEN
     ALLOCATE(z_head_rfr(rfr%nloc))
     z_head_rfr = 0._DP
  ENDIF
  !
  CALL band_group%init(rfr%nloc,'b','band_group',.FALSE.)
  !
  DO iq = 1,q_grid%np
     !
     l_gammaq = q_grid%l_pIsGamma(iq)
     !
     DO ifloc = 1,band_group%nloc
        !
        ifreq = band_group%l2g(ifloc)
        !
        zmatilda = 0._DP
        DO ip = 1,mypara%nloc
           glob_ip = mypara%l2g(ip)
           zmatilda(:,glob_ip) = zmatr_q(:,ip,ifreq,iq)
        ENDDO
        !
        CALL mp_sum(zmatilda,inter_image_comm)
        CALL chi_invert_complex(zmatilda,zhead,zlambda,mypara%nglob,l_gammaq)
        !
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           z_epsm1_rfr_q(1:n_pdep_eigen_to_use,ip,ifreq,iq) = zlambda(1:n_pdep_eigen_to_use,glob_ip)
        ENDDO
        IF(l_macropol .AND. l_gammaq) z_head_rfr(ifreq) = zhead
        !
     ENDDO
     !
  ENDDO
  !
#if defined(__CUDA)
  CALL deallocate_chi_gpu()
#endif
  !
  !$acc exit data delete(zmatilda,zlambda)
  DEALLOCATE(zlambda)
  DEALLOCATE(zmatilda)
  DEALLOCATE(zmatr_q)
  !
  CALL mp_sum(z_epsm1_rfr_q,inter_bgrp_comm)
  IF(l_macropol) CALL mp_sum(z_head_rfr,inter_bgrp_comm)
  !
  CALL stop_clock('chi_invert')
  !
  IF(l_generate_plot) THEN
     CALL output_eps_head()
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE output_eps_head( )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : z_head_rfr,refreq_list,l_macropol,wfreq_save_dir
  USE constants,            ONLY : rytoev,fpi
  USE mp_world,             ONLY : mpime,root
  USE distribution_center,  ONLY : rfr
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : intra_bgrp_comm
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE json_module,          ONLY : json_file
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: out_tabella(:,:)
  INTEGER :: ifreq,glob_ifreq
  REAL(DP) :: lf, rf, ep1, ep2, nn, kk, rr
  TYPE(bar_type) :: barra
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  CHARACTER(20),EXTERNAL :: human_readable_time
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  IF(l_macropol) THEN
     !
     CALL io_push_title('(O)ptics')
     !
     CALL start_bar_type ( barra, 'optics', rfr%nloc )
     !
     time_spent(1) = get_clock( 'optics' )
     !
     ! head_rfr
     !
     ALLOCATE( out_tabella(rfr%nglob,8) )
     !
     out_tabella = 0._DP
     DO ifreq = 1,rfr%nloc
        glob_ifreq = rfr%l2g(ifreq)
        rf =  REAL( z_head_rfr( ifreq ), KIND=DP ) + 1._DP
        lf = -AIMAG( z_head_rfr( ifreq ) )
        ep1 = rf / (rf**2+lf**2)
        ep2 = lf / (rf**2+lf**2)
        nn = SQRT( 0.5_DP * ( ep1 + SQRT( ep1**2 + ep2**2 ) ) )
        kk = SQRT( 0.5_DP * (-ep1 + SQRT( ep1**2 + ep2**2 ) ) )
        rr = ((1._DP-nn)**2+kk**2)/((1._DP+nn)**2+kk**2)
        !
        out_tabella( glob_ifreq, 1 ) = refreq_list( ifreq )*rytoev
        out_tabella( glob_ifreq, 2 ) = ep1
        out_tabella( glob_ifreq, 3 ) = ep2
        out_tabella( glob_ifreq, 4 ) = lf
        out_tabella( glob_ifreq, 5 ) = nn
        out_tabella( glob_ifreq, 6 ) = kk
        out_tabella( glob_ifreq, 7 ) = rr
        out_tabella( glob_ifreq, 8 ) = (ep1-1._DP)*3._DP*omega/fpi/(ep1+2._DP)
        !
        CALL update_bar_type( barra, 'optics', 1 )
        !
     ENDDO
     !
     CALL stop_bar_type( barra, 'optics' )
     !
     CALL mp_sum( out_tabella, intra_bgrp_comm )
     !
     IF(mpime == root) THEN
        !
        CALL json%initialize()
        !
        CALL json%add('e',out_tabella(:,1))
        CALL json%add('eps1',out_tabella(:,2))
        CALL json%add('eps2',out_tabella(:,3))
        CALL json%add('EELF',out_tabella(:,4))
        CALL json%add('n',out_tabella(:,5))
        CALL json%add('k',out_tabella(:,6))
        CALL json%add('refl',out_tabella(:,7))
        CALL json%add('pol',out_tabella(:,8))
        !
        OPEN( NEWUNIT=iunit, FILE=TRIM(wfreq_save_dir)//'/optics.json' )
        CALL json%print( iunit )
        CLOSE( iunit )
        !
        CALL json%destroy()
        !
     ENDIF
     !
     time_spent(2) = get_clock( 'optics' )
     !
     WRITE(stdout,*)
     CALL io_push_bar()
     WRITE(stdout,'(5x, "File ",a," written in ",a)') TRIM(wfreq_save_dir)//'/optics.json',&
     & TRIM(human_readable_time(time_spent(2)-time_spent(1)))
     CALL io_push_bar()
     !
     DEALLOCATE( out_tabella )
     !
  ENDIF
  !
END SUBROUTINE
