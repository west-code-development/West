!
! Copyright (C) 2015-2021 M. Govoni
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
SUBROUTINE solve_gfreq(l_read_restart)
  !-----------------------------------------------------------------------
  !
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart
  !
  IF( gamma_only ) THEN
     CALL solve_gfreq_gamma( l_read_restart )
  ELSE
     CALL solve_gfreq_k( l_read_restart )
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_gfreq_gamma(l_read_restart)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_lanczos,npwq,qp_bandrange,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,&
                                 & o_restart_time,npwqx,fftdriver,wstat_save_dir
  USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,npool,intra_bgrp_comm,nproc_bgrp,nbgrp
  USE mp,                   ONLY : mp_bcast,mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE becmod,               ONLY : becp,allocate_bec_type,deallocate_bec_type
  USE uspp,                 ONLY : vkb,nkb
  USE uspp_init,            ONLY : init_us_2
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert,band_group,kpt_pool
  USE wfreq_restart,        ONLY : solvegfreq_restart_write,solvegfreq_restart_read,bks_type
  USE wfreq_io,             ONLY : writeout_overlap,writeout_solvegfreq
  USE types_coulomb,        ONLY : pot3D
  USE types_bz_grid,        ONLY : k_grid
#if defined(__CUDA)
  USE wavefunctions,        ONLY : evc_host=>evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE west_gpu,             ONLY : ps_r,allocate_gpu,deallocate_gpu,allocate_gw_gpu,deallocate_gw_gpu,&
                                 & reallocate_ps_gpu,memcpy_H2D
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart
  !
  ! Workspace
  !
  LOGICAL :: l_write_restart
  INTEGER :: ip,ig,glob_ip,ir,ib,ibloc,iks,im,iks_g
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  CHARACTER(LEN=25) :: filepot
  INTEGER :: nbndval
  INTEGER :: dffts_nnr,pert_nloc
  REAL(DP),ALLOCATABLE :: diago(:,:),subdiago(:,:),bnorm(:)
  REAL(DP),ALLOCATABLE :: braket(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: braket
#endif
  COMPLEX(DP),ALLOCATABLE :: q_s(:,:,:)
  !$acc declare device_resident(q_s)
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dvpsi
#endif
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  !$acc declare device_resident(pertg,pertr)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: pertg_all
#endif
#if !defined(__CUDA)
  REAL(DP),ALLOCATABLE :: ps_r(:,:)
#endif
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  REAL(DP),ALLOCATABLE :: overlap(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: overlap
#endif
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bks_type) :: bks
  INTEGER,ALLOCATABLE :: l2g(:)
  REAL(DP),ALLOCATABLE :: sqvc(:)
  !$acc declare device_resident(l2g,sqvc)
  !
  CALL io_push_title('(G)-Lanczos')
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type( becp )
  CALL allocate_bec_type ( nkb, pert%nloc, becp ) ! I just need 2 becp at a time
  !
  CALL pot3D%init('Wave',.FALSE.,'default')
  CALL band_group%init(qp_bandrange(2)-qp_bandrange(1)+1,'b','band_group',.FALSE.)
  !
  IF(l_read_restart) THEN
     CALL solvegfreq_restart_read( bks )
  ELSE
     bks%lastdone_ks   = 0
     bks%lastdone_band = 0
     bks%old_ks        = 0
     bks%old_band      = 0
     bks%max_ks        = k_grid%nps
     bks%min_ks        = 1
  ENDIF
  !
  barra_load = 0
  DO iks = 1,kpt_pool%nloc
     IF(iks < bks%lastdone_ks) CYCLE
     !
     DO ibloc = 1,band_group%nloc
        ib = band_group%l2g(ibloc)+qp_bandrange(1)-1
        !
        IF(iks == bks%lastdone_ks .AND. ib <= bks%lastdone_band) CYCLE
        !
        barra_load = barra_load+1
     ENDDO
  ENDDO
  !
  IF( barra_load == 0 ) THEN
     CALL start_bar_type( barra, 'glanczos', 1 )
     CALL update_bar_type( barra, 'glanczos', 1 )
  ELSE
     CALL start_bar_type( barra, 'glanczos', barra_load )
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_gpu()
  CALL allocate_gw_gpu(pert%nlocx,pert%nloc)
  CALL reallocate_ps_gpu(nbnd,pert%nloc)
#endif
  !
  dffts_nnr = dffts%nnr
  pert_nloc = pert%nloc
  !
#if !defined(__CUDA)
  ALLOCATE(ps_r(nbnd,pert%nloc))
#endif
  ALLOCATE(dvpsi(npwx*npol,pert%nlocx))
  !$acc enter data create(dvpsi)
  ALLOCATE(sqvc(npwqx))
#if defined(__CUDA)
  CALL memcpy_H2D(sqvc,pot3D%sqvc,npwqx)
#endif
  ALLOCATE(overlap(pert%nglob,nbnd))
  !$acc enter data create(overlap)
  ALLOCATE(pertr(dffts%nnr))
  ALLOCATE(pertg(npwqx))
  IF(l_enable_lanczos) THEN
     ALLOCATE(bnorm(pert%nloc))
     ALLOCATE(diago(n_lanczos,pert%nloc))
     ALLOCATE(subdiago(n_lanczos-1,pert%nloc))
     ALLOCATE(q_s(npwx*npol,pert%nloc,n_lanczos))
     ALLOCATE(braket(pert%nglob,n_lanczos,pert%nloc))
     !$acc enter data create(braket)
  ENDIF
  ALLOCATE(l2g(pert%nloc))
  !
  !$acc parallel loop present(l2g)
  DO ip = 1,pert_nloc
     !
     ! l2g(ip) = pert%l2g(ip)
     !
     l2g(ip) = nimage*(ip-1)+my_image_id+1
  ENDDO
  !$acc end parallel
  !
  ! Read PDEP
  !
  ALLOCATE(pertg_all(npwqx,pert%nloc))
  pertg_all = 0._DP
  !
  DO ip = 1,pert%nloc
     glob_ip = pert%l2g(ip)
     CALL generate_pdep_fname(filepot,glob_ip)
     fname = TRIM(wstat_save_dir)//'/'//filepot
     CALL pdep_read_G_and_distribute(fname,pertg_all(:,ip))
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
     IF(iks < bks%lastdone_ks) CYCLE
     !
     iks_g = kpt_pool%l2g(iks)
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
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
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_host,0,inter_image_comm)
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
     CALL using_evc(2)
     CALL using_evc_d(0)
#endif
     !
     nbndval = nbnd_occ(iks)
     !
     bks%max_band = nbndval
     bks%min_band = 1
     !
     time_spent(1) = get_clock( 'glanczos' )
     !
     ! LOOP over band states
     !
     DO ibloc = 1,band_group%nloc
        ib = band_group%l2g(ibloc)+qp_bandrange(1)-1
        !
        IF(iks == bks%lastdone_ks .AND. ib <= bks%lastdone_band) CYCLE
        !
        ! PSIC
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc_work(:,ib),psic,'Wave')
        !
        !$acc kernels present(dvpsi)
        dvpsi(:,:) = 0._DP
        !$acc end kernels
        !
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           !
#if defined(__CUDA)
           CALL memcpy_H2D(pertg,pertg_all(:,ip),npwqx)
#else
           pertg = pertg_all(:,ip)
#endif
           !
           ! Multiply by sqvc
           !
           !$acc parallel loop present(pertg,sqvc)
           DO ig = 1,npwq
#if defined(__CUDA)
              pertg(ig) = sqvc(ig)*pertg(ig)
#else
              pertg(ig) = pot3D%sqvc(ig)*pertg(ig)
#endif
           ENDDO
           !$acc end parallel
           !
           ! Bring it to R-space
           !
           !$acc host_data use_device(pertg,pertr)
           CALL single_invfft_gamma(dffts,npwq,npwqx,pertg,pertr,TRIM(fftdriver))
           !$acc end host_data
           !
           !$acc parallel loop present(pertr)
           DO ir = 1,dffts_nnr
              pertr(ir) = psic(ir)*pertr(ir)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(pertr,dvpsi)
           CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi(:,ip),'Wave')
           !$acc end host_data
           !
        ENDDO ! pert
        !
        ! OVERLAP( glob_ip, im=1:n_hstates ) = < psi_im iks | dvpsi_glob_ip >
        !
#if defined(__CUDA)
        !$acc host_data use_device(dvpsi,ps_r)
        CALL glbrak_gamma_gpu(evc_work,dvpsi,ps_r,npw,npwx,nbnd,pert%nloc,nbnd,npol)
        !$acc end host_data
#else
        CALL glbrak_gamma(evc_work,dvpsi,ps_r,npw,npwx,nbnd,pert%nloc,nbnd,npol)
#endif
        IF(nproc_bgrp > 1) THEN
           !$acc host_data use_device(ps_r)
           CALL mp_sum(ps_r,intra_bgrp_comm)
           !$acc end host_data
        ENDIF
        !
        !$acc kernels present(overlap)
        overlap(:,:) = 0._DP
        !$acc end kernels
        !
        !$acc parallel loop collapse(2) present(overlap,ps_r)
        DO im = 1,nbnd
           DO ip = 1,pert_nloc
              overlap(l2g(ip),im) = ps_r(im,ip)
           ENDDO
        ENDDO
        !$acc end parallel
        !
        !$acc update host(overlap)
        CALL mp_sum(overlap,inter_image_comm)
        CALL writeout_overlap( 'g', kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
        !
        !$acc update host(dvpsi)
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF( l_enable_lanczos ) THEN
           !
#if defined(__CUDA)
           CALL solve_deflated_lanczos_w_full_ortho_gpu(nbnd,pert%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
#else
           CALL solve_deflated_lanczos_w_full_ortho(nbnd,pert%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
#endif
           CALL get_brak_hyper_parallel(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert)
           !
           DO ip = 1, pert%nloc
              CALL diago_lanczos(bnorm(ip),diago(:,ip),subdiago(:,ip),braket(:,:,ip),pert%nglob)
           ENDDO
           !
           ! MPI-IO
           !
           CALL writeout_solvegfreq( kpt_pool%l2g(iks), ib, diago, braket, pert%nloc, pert%nglob, pert%myoffset )
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'glanczos' )
        l_write_restart = .FALSE.
        !
        IF( o_restart_time >= 0._DP ) THEN
           IF( time_spent(2)-time_spent(1) > o_restart_time*60._DP ) l_write_restart = .TRUE.
           IF( ib == qp_bandrange(2) ) l_write_restart = .TRUE.
        ENDIF
        !
        ! Write final restart file
        !
        IF( iks == k_grid%nps .AND. ib == qp_bandrange(2) ) l_write_restart = .TRUE.
        !
        ! But do not write here when using pool or band group
        !
        IF( npool*nbgrp > 1 ) l_write_restart = .FALSE.
        !
        IF( l_write_restart ) THEN
           bks%lastdone_ks = iks
           bks%lastdone_band = ib
           CALL solvegfreq_restart_write( bks )
           bks%old_ks = iks
           bks%old_band = ib
           time_spent(1) = get_clock( 'glanczos' )
        ENDIF
        !
        CALL update_bar_type( barra, 'glanczos', 1 )
        !
     ENDDO ! BANDS
     !
  ENDDO ! KPOINT-SPIN
  !
#if defined(__CUDA)
  CALL deallocate_gw_gpu()
  CALL deallocate_gpu()
#endif
  !
#if !defined(__CUDA)
  DEALLOCATE(ps_r)
#endif
  !$acc exit data delete(dvpsi)
  DEALLOCATE(dvpsi)
  DEALLOCATE(sqvc)
  !$acc exit data delete(overlap)
  DEALLOCATE(overlap)
  DEALLOCATE(pertr)
  DEALLOCATE(pertg)
  IF(l_enable_lanczos) THEN
     DEALLOCATE(bnorm)
     DEALLOCATE(diago)
     DEALLOCATE(subdiago)
     DEALLOCATE(q_s)
     !$acc exit data delete(braket)
     DEALLOCATE(braket)
  ENDIF
  DEALLOCATE(l2g)
  DEALLOCATE(pertg_all)
  !
  ! Write final restart file when using pool or band group
  !
  IF( npool*nbgrp > 1 ) THEN
     bks%lastdone_ks = k_grid%nps
     bks%lastdone_band = qp_bandrange(2)
     CALL solvegfreq_restart_write( bks )
  ENDIF
  !
  CALL stop_bar_type( barra, 'glanczos' )
  !
  CALL mp_barrier( world_comm )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_gfreq_k(l_read_restart)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_lanczos,npwq,qp_bandrange,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,&
                                 & o_restart_time,npwqx,wstat_save_dir,ngq,igq_q
  USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm,nproc_bgrp,nbgrp
  USE mp,                   ONLY : mp_bcast,mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE becmod,               ONLY : becp,allocate_bec_type,deallocate_bec_type
  USE uspp,                 ONLY : vkb,nkb
  USE uspp_init,            ONLY : init_us_2
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert,band_group,kpt_pool
  USE wfreq_restart,        ONLY : solvegfreq_restart_write_q,solvegfreq_restart_read_q,bksks_type
  USE wfreq_io,             ONLY : writeout_overlap,writeout_solvegfreq
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE types_coulomb,        ONLY : pot3D
#if defined(__CUDA)
  USE wavefunctions,        ONLY : evc_host=>evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE west_gpu,             ONLY : ps_c,allocate_gpu,deallocate_gpu,allocate_gw_gpu,deallocate_gw_gpu,&
                                 & reallocate_ps_gpu,memcpy_H2D
#else
  USE wavefunctions,        ONLY : evc_work=>evc
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart
  !
  ! Workspace
  !
  LOGICAL :: l_write_restart
  INTEGER :: ip,ig,glob_ip,ir,ib,ibloc,iks,ik,im,ikks,ikk,iq
  INTEGER :: npwk
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  CHARACTER(LEN=25) :: filepot
  INTEGER :: nbndval
  INTEGER :: dffts_nnr,pert_nloc
  REAL(DP) :: g0(3)
  REAL(DP),ALLOCATABLE :: diago(:,:),subdiago(:,:),bnorm(:)
  COMPLEX(DP),ALLOCATABLE :: braket(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: braket
#endif
  COMPLEX(DP),ALLOCATABLE :: q_s(:,:,:)
  !$acc declare device_resident(q_s)
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dvpsi
#endif
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  !$acc declare device_resident(pertg,pertr)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
  COMPLEX(DP),ALLOCATABLE :: evck(:,:),phase(:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: pertg_all,evck,phase
#endif
  COMPLEX(DP),ALLOCATABLE :: psick(:),psick_nc(:,:)
  !$acc declare device_resident(psick,psick_nc)
#if !defined(__CUDA)
  COMPLEX(DP),ALLOCATABLE :: ps_c(:,:)
#endif
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: overlap
#endif
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bksks_type) :: bksks
  INTEGER,ALLOCATABLE :: l2g(:)
  REAL(DP),ALLOCATABLE :: sqvc(:)
  !$acc declare device_resident(l2g,sqvc)
  !
  CALL io_push_title('(G)-Lanczos')
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type( becp )
  CALL allocate_bec_type ( nkb, pert%nloc, becp ) ! I just need 2 becp at a time
  !
  CALL band_group%init(qp_bandrange(2)-qp_bandrange(1)+1,'b','band_group',.FALSE.)
  !
  IF(l_read_restart) THEN
     CALL solvegfreq_restart_read_q( bksks )
  ELSE
     bksks%lastdone_ks     = 0
     bksks%lastdone_kks    = 0
     bksks%lastdone_band   = 0
     bksks%old_ks          = 0
     bksks%old_kks         = 0
     bksks%old_band        = 0
     bksks%max_ks          = k_grid%nps
     bksks%min_ks          = 1
     bksks%max_kks         = k_grid%nps
     bksks%min_kks         = 1
  ENDIF
  !
  barra_load = 0
  DO ikks = 1,k_grid%nps
     IF(ikks < bksks%lastdone_ks) CYCLE
     !
     DO ibloc = 1,band_group%nloc
        ib = band_group%l2g(ibloc)+qp_bandrange(1)-1
        !
        IF(ikks == bksks%lastdone_ks .AND. ib < bksks%lastdone_band) CYCLE
        !
        DO iks = 1,k_grid%nps
           IF(ikks == bksks%lastdone_ks .AND. ib == bksks%lastdone_band .AND. iks <= bksks%lastdone_kks) CYCLE
           barra_load = barra_load+1
        ENDDO
     ENDDO
  ENDDO
  !
  IF( barra_load == 0 ) THEN
     CALL start_bar_type( barra, 'glanczos', 1 )
     CALL update_bar_type( barra, 'glanczos', 1 )
  ELSE
     CALL start_bar_type( barra, 'glanczos', barra_load )
  ENDIF
  !
  IF(noncolin) THEN
     ALLOCATE(psick_nc(dffts%nnr,npol))
  ELSE
     ALLOCATE(psick(dffts%nnr))
  ENDIF
  ALLOCATE(phase(dffts%nnr))
  ALLOCATE(evck(npwx*npol,nbnd))
  !$acc enter data create(phase,evck)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
  CALL allocate_gw_gpu(pert%nlocx,pert%nloc)
  CALL reallocate_ps_gpu(nbnd,pert%nloc)
#endif
  !
  dffts_nnr = dffts%nnr
  pert_nloc = pert%nloc
  !
#if !defined(__CUDA)
  ALLOCATE(ps_c(nbnd,pert%nloc))
#endif
  ALLOCATE(dvpsi(npwx*npol,pert%nlocx))
  !$acc enter data create(dvpsi)
  ALLOCATE(sqvc(npwqx))
  ALLOCATE(overlap(pert%nglob,nbnd))
  !$acc enter data create(overlap)
  ALLOCATE(pertr(dffts%nnr))
  ALLOCATE(pertg(npwqx))
  IF(l_enable_lanczos) THEN
     ALLOCATE(bnorm(pert%nloc))
     ALLOCATE(diago(n_lanczos,pert%nloc))
     ALLOCATE(subdiago(n_lanczos-1,pert%nloc))
     ALLOCATE(q_s(npwx*npol,pert%nloc,n_lanczos))
     ALLOCATE(braket(pert%nglob,n_lanczos,pert%nloc))
     !$acc enter data create(braket)
  ENDIF
  ALLOCATE(l2g(pert%nloc))
  !
  !$acc parallel loop present(l2g)
  DO ip = 1,pert_nloc
     !
     ! l2g(ip) = pert%l2g(ip)
     !
     l2g(ip) = nimage*(ip-1)+my_image_id+1
  ENDDO
  !$acc end parallel
  !
  ALLOCATE(pertg_all(npwqx,pert%nloc))
  !
  ! LOOP
  !
  ! ... Outer k-point loop (wfc matrix element): ikks, npwk, evck, psick
  ! ... Inner k-point loop (wfc summed over k'): iks, npw, evc (passed to h_psi: current_k = iks)
  !
  DO ikks = 1, k_grid%nps ! KPOINT-SPIN (MATRIX ELEMENT)
     !
     ! Exit loop if no work to do
     !
     IF(barra_load == 0) EXIT
     !
     IF(ikks < bksks%lastdone_ks) CYCLE
     !
     ikk = k_grid%ip(ikks)
     !
     npwk = ngk(ikks)
     !
     IF(my_image_id == 0) CALL get_buffer(evck,lrwfc,iuwfc,ikks)
     CALL mp_bcast(evck,0,inter_image_comm)
     !
     !$acc update device(evck)
     !
     nbndval = nbnd_occ(ikks)
     !
     bksks%max_band = nbndval
     bksks%min_band = 1
     !
     DO iks = 1, k_grid%nps ! KPOINT-SPIN (INTEGRAL OVER K')
        IF(ikks == bksks%lastdone_ks .AND. iks < bksks%lastdone_kks) CYCLE
        !
        ik = k_grid%ip(iks)
        !
        time_spent(1) = get_clock( 'glanczos' )
        !
        CALL q_grid%find( k_grid%p_cart(:,ikk) - k_grid%p_cart(:,ik), 'cart', iq, g0 )
        !
        npwq = ngq(iq)
        !
        ! compute Coulomb potential
        !
        CALL pot3D%init('Wave',.TRUE.,'default',iq)
        !
#if defined(__CUDA)
        CALL memcpy_H2D(sqvc,pot3D%sqvc,npwqx)
#endif
        !
        ! The Hamiltonian is evaluated at k'
        !
        npw = ngk(iks)
        !
        ! ... Set k-point, spin, kinetic energy, needed by Hpsi
        !
        current_k = iks
        IF ( lsda ) current_spin = isk(iks)
#if defined(__CUDA)
        call g2_kin_gpu(iks)
        !
        ! ... More stuff needed by the hamiltonian: nonlocal projectors
        !
        IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
        call g2_kin(iks)
        !
        ! ... More stuff needed by the hamiltonian: nonlocal projectors
        !
        IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
        !
!       ! ... Needed for LDA+U
!       !
!       IF ( kpt_pool%nloc > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
!            CALL get_buffer ( wfcU, nwordwfcU, iunhub, iks )
!       !
!       current_k = iks
!       current_spin = isk(iks)
!       !
!       CALL gk_sort(xk(1,iks),ngm,g,gcutw,npw,igk,g2kin)
!       g2kin=g2kin*tpiba2
!       !
!       ! reads unperturbed wavefuctions psi_k in G_space, for all bands
!       !
!       CALL init_us_2 (npw, igk, xk (1, iks), vkb)
!       !
        CALL compute_phase( g0, 'cart', phase )
        !
        !$acc update device(phase)
        !
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_host,0,inter_image_comm)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
        !
#if defined(__CUDA)
        !
        ! ... Sync GPU
        !
        CALL using_becp_auto(2)
        CALL using_becp_d_auto(0)
        CALL using_evc(2)
        CALL using_evc_d(0)
#endif
        !
        ! Read PDEP
        !
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           CALL generate_pdep_fname(filepot,glob_ip,iq)
           fname = TRIM(wstat_save_dir)//'/'//filepot
           CALL pdep_read_G_and_distribute(fname,pertg_all(:,ip),iq)
        ENDDO
        !
        ! LOOP over band states
        !
        DO ibloc = 1,band_group%nloc
           ib = band_group%l2g(ibloc)+qp_bandrange(1)-1
           !
           IF(ikks == bksks%lastdone_ks .AND. iks == bksks%lastdone_kks .AND. ib <= bksks%lastdone_band) CYCLE
           !
           ! PSIC
           !
           IF(noncolin) THEN
              !$acc host_data use_device(evck,psick_nc)
              CALL single_invfft_k(dffts,npwk,npwx,evck(1:npwx,ib),psick_nc(:,1),'Wave',igk_k(:,ikks))
              CALL single_invfft_k(dffts,npwk,npwx,evck(1+npwx:npwx*2,ib),psick_nc(:,2),'Wave',igk_k(:,ikks))
              !$acc end host_data
           ELSE
              !$acc host_data use_device(evck,psick)
              CALL single_invfft_k(dffts,npwk,npwx,evck(:,ib),psick,'Wave',igk_k(:,ikks))
              !$acc end host_data
           ENDIF
           !
           !$acc kernels present(dvpsi)
           dvpsi(:,:) = 0._DP
           !$acc end kernels
           !
           DO ip = 1,pert%nloc
              glob_ip = pert%l2g(ip)
              !
#if defined(__CUDA)
              CALL memcpy_H2D(pertg,pertg_all(:,ip),npwqx)
#else
              pertg = pertg_all(:,ip)
#endif
              !
              ! Multiply by sqvc
              !
              !$acc parallel loop present(pertg,sqvc)
              DO ig = 1,npwq
#if defined(__CUDA)
                 pertg(ig) = sqvc(ig)*pertg(ig)
#else
                 pertg(ig) = pot3D%sqvc(ig)*pertg(ig)
#endif
              ENDDO
              !$acc end parallel
              !
              ! Bring it to R-space
              !
              IF(noncolin) THEN
                 !$acc host_data use_device(pertg,pertr)
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                 !$acc end host_data
                 !$acc parallel loop present(pertr,phase,psick_nc)
                 DO ir = 1,dffts_nnr
                    pertr(ir) = CONJG(phase(ir))*psick_nc(ir,1)*CONJG(pertr(ir))
                 ENDDO
                 !$acc end parallel
                 !$acc host_data use_device(pertr,dvpsi,pertg)
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1:npwx,ip),'Wave',igk_k(:,current_k))
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                 !$acc end host_data
                 !$acc parallel loop present(pertr,phase,psick_nc)
                 DO ir = 1,dffts_nnr
                    pertr(ir) = CONJG(phase(ir))*psick_nc(ir,2)*CONJG(pertr(ir))
                 ENDDO
                 !$acc end parallel
                 !$acc host_data use_device(pertr,dvpsi)
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1+npwx:npwx*2,ip),'Wave',igk_k(:,current_k))
                 !$acc end host_data
              ELSE
                 !$acc host_data use_device(pertg,pertr)
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                 !$acc end host_data
                 !$acc parallel loop present(pertr,phase,psick)
                 DO ir = 1,dffts_nnr
                    pertr(ir) = CONJG(phase(ir))*psick(ir)*CONJG(pertr(ir))
                 ENDDO
                 !$acc end parallel
                 !$acc host_data use_device(pertr,dvpsi)
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(:,ip),'Wave',igk_k(:,current_k))
                 !$acc end host_data
              ENDIF
              !
           ENDDO ! pert
           !
           ! OVERLAP( glob_ip, im=1:n_hstates ) = < psi_im iks | dvpsi_glob_ip >
           !
#if defined(__CUDA)
           !$acc host_data use_device(dvpsi,ps_c)
           CALL glbrak_k_gpu(evc_work,dvpsi,ps_c,npw,npwx,nbnd,pert%nloc,nbnd,npol)
           !$acc end host_data
#else
           CALL glbrak_k(evc_work,dvpsi,ps_c,npw,npwx,nbnd,pert%nloc,nbnd,npol)
#endif
           IF(nproc_bgrp > 1) THEN
              !$acc host_data use_device(ps_c)
              CALL mp_sum(ps_c,intra_bgrp_comm)
              !$acc end host_data
           ENDIF
           !
           !$acc kernels present(overlap)
           overlap(:,:) = 0._DP
           !$acc end kernels
           !
           !$acc parallel loop collapse(2) present(overlap,ps_c)
           DO im = 1,nbnd
              DO ip = 1,pert_nloc
                 overlap(l2g(ip),im) = ps_c(im,ip)
              ENDDO
           ENDDO
           !$acc end parallel
           !
           !$acc update host(overlap)
           CALL mp_sum(overlap,inter_image_comm)
           CALL writeout_overlap( 'g', kpt_pool%l2g(ikks), kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
           !
           CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
           !
           !$acc update host(dvpsi)
           !
           ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
           ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
           !
           IF( l_enable_lanczos ) THEN
              !
#if defined(__CUDA)
              CALL solve_deflated_lanczos_w_full_ortho_gpu(nbnd,pert%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
#else
              CALL solve_deflated_lanczos_w_full_ortho(nbnd,pert%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
#endif
              CALL get_brak_hyper_parallel_complex(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert)
              !
              DO ip = 1, pert%nloc
                 CALL diago_lanczos_complex(bnorm(ip),diago(:,ip),subdiago(:,ip),braket(:,:,ip),pert%nglob)
              ENDDO
              !
              ! MPI-IO
              !
              CALL writeout_solvegfreq( kpt_pool%l2g(ikks), kpt_pool%l2g(iks), ib, diago, braket, pert%nloc, &
              & pert%nglob, pert%myoffset )
              !
           ENDIF ! l_enable_lanczos
           !
           time_spent(2) = get_clock( 'glanczos' )
           l_write_restart = .FALSE.
           !
           IF( o_restart_time >= 0._DP ) THEN
              IF( time_spent(2)-time_spent(1) > o_restart_time*60._DP ) l_write_restart = .TRUE.
              IF( ib == qp_bandrange(2) ) l_write_restart = .TRUE.
           ENDIF
           !
           ! Write final restart file
           !
           IF( ikks == k_grid%nps .AND. iks == k_grid%nps .AND. ib == qp_bandrange(2) ) l_write_restart = .TRUE.
           !
           ! But do not write here when using band group
           !
           IF( nbgrp > 1 ) l_write_restart = .FALSE.
           !
           IF( l_write_restart ) THEN
              bksks%lastdone_ks = ikks
              bksks%lastdone_kks = iks
              bksks%lastdone_band = ib
              CALL solvegfreq_restart_write_q( bksks )
              bksks%old_ks = ikks
              bksks%old_kks = iks
              bksks%old_band = ib
              time_spent(1) = get_clock( 'glanczos' )
           ENDIF
           !
           CALL update_bar_type( barra, 'glanczos', 1 )
           !
        ENDDO ! BANDS
        !
     ENDDO ! KPOINT-SPIN (INTEGRAL OVER K')
     !
  ENDDO ! KPOINT-SPIN (MATRIX ELEMENT)
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
  CALL deallocate_gw_gpu()
#endif
  !
  IF(noncolin) THEN
     DEALLOCATE(psick_nc)
  ELSE
     DEALLOCATE(psick)
  ENDIF
  !$acc exit data delete(phase,evck)
  DEALLOCATE(phase)
  DEALLOCATE(evck)
  !
#if !defined(__CUDA)
  DEALLOCATE(ps_c)
#endif
  !$acc exit data delete(dvpsi)
  DEALLOCATE(dvpsi)
  DEALLOCATE(sqvc)
  !$acc exit data delete(overlap)
  DEALLOCATE(overlap)
  DEALLOCATE(pertr)
  DEALLOCATE(pertg)
  IF(l_enable_lanczos) THEN
     DEALLOCATE(bnorm)
     DEALLOCATE(diago)
     DEALLOCATE(subdiago)
     DEALLOCATE(q_s)
     !$acc exit data delete(braket)
     DEALLOCATE(braket)
  ENDIF
  DEALLOCATE(l2g)
  DEALLOCATE(pertg_all)
  !
  ! Write final restart file when using band group
  !
  IF( nbgrp > 1 ) THEN
     bksks%lastdone_ks = k_grid%nps
     bksks%lastdone_kks = k_grid%nps
     bksks%lastdone_band = qp_bandrange(2)
     CALL solvegfreq_restart_write_q( bksks )
  ENDIF
  !
  CALL stop_bar_type( barra, 'glanczos' )
  !
  CALL mp_barrier( world_comm )
  !
END SUBROUTINE
