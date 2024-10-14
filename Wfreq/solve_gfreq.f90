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
  USE westcom,              ONLY : n_lanczos,npwq,npwqx,qp_bands,n_bands,nbnd_occ,l_enable_lanczos,&
                                 & l_enable_off_diagonal,iuwfc,lrwfc,wstat_save_dir,fftdriver,ijpmap,&
                                 & d_body2_ifr,d_diago,d_body2_ifr_full,d_diago_full,d_epsm1_ifr
  USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_pool_comm,npool,&
                                 & intra_bgrp_comm
  USE mp,                   ONLY : mp_bcast,mp_sum
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE becmod,               ONLY : becp,allocate_bec_type_acc,deallocate_bec_type_acc
  USE uspp,                 ONLY : vkb,nkb
  USE uspp_init,            ONLY : init_us_2
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert,kpt_pool,band_group,ifr
  USE wfreq_io,             ONLY : write_overlap,write_gfreq,read_gfreq
  USE types_coulomb,        ONLY : pot3D
  USE types_bz_grid,        ONLY : k_grid
  USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
  USE west_gpu,             ONLY : ps_r,allocate_gpu,deallocate_gpu,allocate_gw_gpu,deallocate_gw_gpu,&
                                 & allocate_lanczos_gpu,deallocate_lanczos_gpu,reallocate_ps_gpu
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
  INTEGER :: ip,glob_ip,glob_jp,il,ig,ir,ib,ibloc,ib_index,jb,jb_index,im,iks,iks_g,is,ifreq
  INTEGER :: ipair,iloc_pair,nloc_pairs
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  CHARACTER(LEN=25) :: filepot
  INTEGER :: nbndval
  INTEGER :: dffts_nnr,pert_nloc,pert_nglob,ifr_nloc
  REAL(DP) :: reduce
  REAL(DP),ALLOCATABLE :: diago(:,:),subdiago(:,:),bnorm(:)
  REAL(DP),ALLOCATABLE :: braket(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: braket
#endif
  REAL(DP),ALLOCATABLE :: diago1(:,:),subdiago1(:,:)
  COMPLEX(DP),ALLOCATABLE :: q_s(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dvpsi
#endif
  COMPLEX(DP),ALLOCATABLE :: psic1(:),dvpsi1(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
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
  INTEGER,ALLOCATABLE :: l2g(:)
  !
  CALL io_push_title('(G)-Lanczos')
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type_acc( becp )
  CALL allocate_bec_type_acc( nkb, pert%nloc, becp ) ! I just need 2 becp at a time
  !
  CALL pot3D%init('Wave',.FALSE.,'default')
  !
  !$acc enter data copyin(pot3D)
  !$acc enter data copyin(pot3D%sqvc)
  !
  CALL band_group%init(n_bands,'b','band_group',.FALSE.)
  !
  barra_load = 0
  DO iks = 1,kpt_pool%nloc
     iks_g = kpt_pool%l2g(iks)
     is = k_grid%is(iks_g)
     DO ibloc = 1,band_group%nloc
        ib_index = band_group%l2g(ibloc)
        IF(l_enable_off_diagonal) THEN
           barra_load = barra_load+ib_index
        ELSE
           barra_load = barra_load+1
        ENDIF
     ENDDO
  ENDDO
  IF(l_read_restart) barra_load = 0
  !
  IF(barra_load == 0) THEN
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
  IF(l_enable_lanczos) CALL allocate_lanczos_gpu(pert%nloc)
#endif
  !
  !$acc enter data copyin(d_epsm1_ifr)
  !
  dffts_nnr = dffts%nnr
  pert_nloc = pert%nloc
  pert_nglob = pert%nglob
  ifr_nloc = ifr%nloc
  !
#if !defined(__CUDA)
  ALLOCATE(ps_r(nbnd,pert%nloc))
#endif
  ALLOCATE(dvpsi(npwx*npol,pert%nlocx))
  ALLOCATE(overlap(pert%nglob,nbnd))
  ALLOCATE(pertr(dffts%nnr))
  ALLOCATE(pertg(npwqx))
  !$acc enter data create(dvpsi,overlap,pertr,pertg)
  IF(l_enable_lanczos) THEN
     ALLOCATE(bnorm(pert%nloc))
     ALLOCATE(diago(n_lanczos,pert%nloc))
     ALLOCATE(subdiago(n_lanczos-1,pert%nloc))
     ALLOCATE(q_s(npwx*npol,pert%nloc,n_lanczos))
     ALLOCATE(braket(pert%nglob,n_lanczos,pert%nloc))
     !$acc enter data create(q_s,braket)
     IF(l_enable_off_diagonal) THEN
        ALLOCATE(dvpsi1(npwx*npol,pert%nlocx))
        ALLOCATE(psic1(dffts%nnr))
        !$acc enter data create(dvpsi1,psic1)
        ALLOCATE(diago1(n_lanczos,pert%nloc))
        ALLOCATE(subdiago1(n_lanczos-1,pert%nloc))
        !
        nloc_pairs = 0
        DO iks = 1,k_grid%nps
           im = 0
           DO ibloc = 1,band_group%nloc
              ib_index = band_group%l2g(ibloc)
              ib = qp_bands(ib_index,iks)
              DO jb_index = 1,n_bands
                 jb = qp_bands(jb_index,iks)
                 IF(jb <= ib) im = im+1
              ENDDO
           ENDDO
           nloc_pairs = MAX(nloc_pairs,im)
        ENDDO
        !
        ALLOCATE(d_body2_ifr_full(n_lanczos,pert%nloc,ifr%nloc,nloc_pairs,k_grid%nps))
        ALLOCATE(d_diago_full(n_lanczos,pert%nloc,nloc_pairs,k_grid%nps))
        d_body2_ifr_full(:,:,:,:,:) = 0._DP
        d_diago_full(:,:,:,:) = 0._DP
     ELSE
        ALLOCATE(d_body2_ifr(n_lanczos,pert%nloc,ifr%nloc,band_group%nloc,k_grid%nps))
        ALLOCATE(d_diago(n_lanczos,pert%nloc,band_group%nloc,k_grid%nps))
        d_body2_ifr(:,:,:,:,:) = 0._DP
        d_diago(:,:,:,:) = 0._DP
     ENDIF
  ENDIF
  ALLOCATE(l2g(pert%nloc))
  !$acc enter data create(l2g)
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
  IF(l_read_restart .AND. l_enable_lanczos) THEN
     IF(l_enable_off_diagonal) THEN
        CALL read_gfreq(d_diago_full,d_body2_ifr_full,pert%nloc,nloc_pairs)
     ELSE
        CALL read_gfreq(d_diago,d_body2_ifr,pert%nloc,band_group%nloc)
     ENDIF
  ENDIF
  !
  ! Read PDEP
  !
  ALLOCATE(pertg_all(npwqx,pert%nloc))
  pertg_all(:,:) = 0._DP
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
     iks_g = kpt_pool%l2g(iks)
     is = k_grid%is(iks_g)
     iloc_pair = 0
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     CALL g2_kin( iks )
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
     nbndval = nbnd_occ(iks)
     !
     ! LOOP over band states
     !
     DO ibloc = 1,band_group%nloc
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index,is)
        !
        ! PSIC
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ib),psic,'Wave')
        !
        !$acc kernels present(dvpsi)
        dvpsi(:,:) = 0._DP
        !$acc end kernels
        !
        DO ip = 1,pert%nloc
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
        ENDDO ! pert
        !
        ! OVERLAP( glob_ip, im=1:n_hstates ) = < psi_im iks | dvpsi_glob_ip >
        !
        CALL glbrak_gamma(evc,dvpsi,ps_r,npw,npwx,nbnd,pert%nloc,nbnd,npol)
        !
        !$acc host_data use_device(ps_r)
        CALL mp_sum(ps_r,intra_bgrp_comm)
        !$acc end host_data
        !
        !$acc kernels present(overlap)
        overlap(:,:) = 0._DP
        !$acc end kernels
        !
        !$acc parallel loop collapse(2) present(overlap,l2g,ps_r)
        DO im = 1,nbnd
           DO ip = 1,pert_nloc
              overlap(l2g(ip),im) = ps_r(im,ip)
           ENDDO
        ENDDO
        !$acc end parallel
        !
        !$acc update host(overlap)
        CALL mp_sum(overlap,inter_image_comm)
        CALL write_overlap( 'g', kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
        !
        !$acc update host(dvpsi)
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF(l_enable_lanczos) THEN
           !
           CALL solve_deflated_lanczos_w_full_ortho(nbnd,pert%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
           !
           DO jb_index = 1,n_bands
              !
              jb = qp_bands(jb_index,is)
              !
              IF(l_enable_off_diagonal .AND. jb <= ib) THEN
                 !
                 ipair = ijpmap(jb_index,ib_index)
                 iloc_pair = iloc_pair+1
                 !
                 ! PSIC
                 !
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jb),psic1,'Wave')
                 !
                 !$acc kernels present(dvpsi1)
                 dvpsi1(:,:) = 0._DP
                 !$acc end kernels
                 !
                 DO ip = 1,pert%nloc
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
                    !$acc parallel loop present(pertr,psic1)
                    DO ir = 1,dffts_nnr
                       pertr(ir) = psic1(ir)*pertr(ir)
                    ENDDO
                    !$acc end parallel
                    !
                    CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi1(:,ip),'Wave')
                    !
                 ENDDO ! pert
                 !
                 CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi1,(1._DP,0._DP))
                 !
                 !$acc update host(dvpsi1)
                 !
                 CALL get_brak_hyper_parallel(dvpsi1,pert%nloc,n_lanczos,q_s,braket,pert)
                 !
                 diago1(:,:) = diago
                 subdiago1(:,:) = subdiago
                 !
                 DO ip = 1,pert%nloc
                    CALL diago_lanczos(bnorm(ip),diago1(:,ip),subdiago1(:,ip),braket(:,:,ip),pert%nglob)
                 ENDDO
                 !
                 d_diago_full(:,:,iloc_pair,iks_g) = diago1
                 !
                 !$acc enter data create(d_body2_ifr_full(:,:,:,iloc_pair,iks_g))
                 !$acc update device(braket)
                 !
                 !$acc parallel present(braket,d_epsm1_ifr,d_body2_ifr_full(:,:,:,iloc_pair,iks_g))
                 !$acc loop collapse(3)
                 DO ifreq = 1,ifr_nloc
                    DO ip = 1,pert_nloc
                       DO il = 1,n_lanczos
                          reduce = 0._DP
                          !$acc loop reduction(+:reduce)
                          DO glob_jp = 1,pert_nglob
                             reduce = reduce+braket(glob_jp,il,ip)*d_epsm1_ifr(glob_jp,ip,ifreq)
                          ENDDO
                          d_body2_ifr_full(il,ip,ifreq,iloc_pair,iks_g) = reduce
                       ENDDO
                    ENDDO
                 ENDDO
                 !$acc end parallel
                 !
                 !$acc exit data copyout(d_body2_ifr_full(:,:,:,iloc_pair,iks_g))
                 !
                 CALL update_bar_type( barra, 'glanczos', 1 )
                 !
              ELSEIF(.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                 !
                 CALL get_brak_hyper_parallel(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert)
                 !
                 DO ip = 1,pert%nloc
                    CALL diago_lanczos(bnorm(ip),diago(:,ip),subdiago(:,ip),braket(:,:,ip),pert%nglob)
                 ENDDO
                 !
                 d_diago(:,:,ibloc,iks_g) = diago
                 !
                 !$acc enter data create(d_body2_ifr(:,:,:,ibloc,iks_g))
                 !$acc update device(braket)
                 !
                 !$acc parallel present(braket,d_epsm1_ifr,d_body2_ifr(:,:,:,ibloc,iks_g))
                 !$acc loop collapse(3)
                 DO ifreq = 1,ifr_nloc
                    DO ip = 1,pert_nloc
                       DO il = 1,n_lanczos
                          reduce = 0._DP
                          !$acc loop reduction(+:reduce)
                          DO glob_jp = 1,pert_nglob
                             reduce = reduce+braket(glob_jp,il,ip)*d_epsm1_ifr(glob_jp,ip,ifreq)
                          ENDDO
                          d_body2_ifr(il,ip,ifreq,ibloc,iks_g) = reduce
                       ENDDO
                    ENDDO
                 ENDDO
                 !$acc end parallel
                 !
                 !$acc exit data copyout(d_body2_ifr(:,:,:,ibloc,iks_g))
                 !
                 CALL update_bar_type( barra, 'glanczos', 1 )
                 !
              ENDIF
              !
           ENDDO
           !
        ENDIF ! l_enable_lanczos
        !
     ENDDO ! BANDS
     !
  ENDDO ! KPOINT-SPIN
  !
#if defined(__CUDA)
  CALL deallocate_gw_gpu()
  CALL deallocate_gpu()
  IF(l_enable_lanczos) CALL deallocate_lanczos_gpu()
#endif
  !
  !$acc exit data delete(d_epsm1_ifr)
  !
#if !defined(__CUDA)
  DEALLOCATE(ps_r)
#endif
  !$acc exit data delete(dvpsi,overlap,pertr,pertg)
  DEALLOCATE(dvpsi)
  DEALLOCATE(overlap)
  DEALLOCATE(pertr)
  DEALLOCATE(pertg)
  IF(l_enable_lanczos) THEN
     DEALLOCATE(bnorm)
     DEALLOCATE(diago)
     DEALLOCATE(subdiago)
     !$acc exit data delete(q_s,braket)
     DEALLOCATE(q_s)
     DEALLOCATE(braket)
     IF(l_enable_off_diagonal) THEN
        !$acc exit data delete(dvpsi1,psic1)
        DEALLOCATE(dvpsi1)
        DEALLOCATE(psic1)
        DEALLOCATE(subdiago1)
        DEALLOCATE(diago1)
     ENDIF
  ENDIF
  !$acc exit data delete(l2g)
  DEALLOCATE(l2g)
  DEALLOCATE(pertg_all)
  !
  !$acc exit data delete(pot3D%sqvc)
  !$acc exit data delete(pot3D)
  !
  ! Synchronize and write data
  !
  IF(.NOT. l_read_restart .AND. l_enable_lanczos) THEN
     IF(l_enable_off_diagonal) THEN
        IF(npool > 1) THEN
           CALL mp_sum(d_body2_ifr_full,inter_pool_comm)
           CALL mp_sum(d_diago_full,inter_pool_comm)
        ENDIF
        CALL write_gfreq(d_diago_full,d_body2_ifr_full,pert%nloc,nloc_pairs)
     ELSE
        IF(npool > 1) THEN
           CALL mp_sum(d_body2_ifr,inter_pool_comm)
           CALL mp_sum(d_diago,inter_pool_comm)
        ENDIF
        CALL write_gfreq(d_diago,d_body2_ifr,pert%nloc,band_group%nloc)
     ENDIF
  ENDIF
  !
  CALL stop_bar_type( barra, 'glanczos' )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_gfreq_k(l_read_restart)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_lanczos,npwq,npwqx,qp_bands,n_bands,nbnd_occ,l_enable_lanczos,&
                                 & iuwfc,lrwfc,wstat_save_dir,ngq,igq_q,z_body2_ifr_q,d_diago_q,&
                                 & z_epsm1_ifr_q
  USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm
  USE mp,                   ONLY : mp_bcast,mp_sum
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE becmod,               ONLY : becp,allocate_bec_type_acc,deallocate_bec_type_acc
  USE uspp,                 ONLY : vkb,nkb
  USE uspp_init,            ONLY : init_us_2
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert,kpt_pool,band_group,ifr
  USE wfreq_io,             ONLY : write_overlap,write_gfreq,read_gfreq
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE types_coulomb,        ONLY : pot3D
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE west_gpu,             ONLY : ps_c,allocate_gpu,deallocate_gpu,allocate_gw_gpu,deallocate_gw_gpu,&
                                 & allocate_lanczos_gpu,deallocate_lanczos_gpu,reallocate_ps_gpu
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
  INTEGER :: ip,glob_ip,glob_jp,il,ig,ir,ib,ibloc,ib_index,im,iks,ik,is,ikks,ikk,iq,ifreq
  INTEGER :: npwk
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  CHARACTER(LEN=25) :: filepot
  INTEGER :: nbndval
  INTEGER :: dffts_nnr,pert_nloc,pert_nglob,ifr_nloc
  COMPLEX(DP) :: reduce
  REAL(DP) :: g0(3)
  REAL(DP),ALLOCATABLE :: diago(:,:),subdiago(:,:),bnorm(:)
  COMPLEX(DP),ALLOCATABLE :: braket(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: braket
#endif
  COMPLEX(DP),ALLOCATABLE :: q_s(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dvpsi
#endif
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
  COMPLEX(DP),ALLOCATABLE :: evck(:,:),phase(:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: pertg_all,evck
#endif
  COMPLEX(DP),ALLOCATABLE :: psick(:),psick_nc(:,:)
#if !defined(__CUDA)
  COMPLEX(DP),ALLOCATABLE :: ps_c(:,:)
#endif
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: overlap
#endif
  INTEGER,ALLOCATABLE :: l2g(:)
  !
  CALL io_push_title('(G)-Lanczos')
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type_acc( becp )
  CALL allocate_bec_type_acc( nkb, pert%nloc, becp ) ! I just need 2 becp at a time
  !
  CALL band_group%init(n_bands,'b','band_group',.FALSE.)
  !
  barra_load = k_grid%nps*band_group%nloc*k_grid%nps
  IF(l_read_restart) barra_load = 0
  !
  IF(barra_load == 0) THEN
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
  IF(l_enable_lanczos) CALL allocate_lanczos_gpu(pert%nloc)
#endif
  !
  dffts_nnr = dffts%nnr
  pert_nloc = pert%nloc
  pert_nglob = pert%nglob
  ifr_nloc = ifr%nloc
  !
  IF(noncolin) THEN
     ALLOCATE(psick_nc(dffts%nnr,npol))
     !$acc enter data create(psick_nc)
  ELSE
     ALLOCATE(psick(dffts%nnr))
     !$acc enter data create(psick)
  ENDIF
  ALLOCATE(phase(dffts%nnr))
  ALLOCATE(evck(npwx*npol,nbnd))
  !$acc enter data create(phase,evck)
  !
#if !defined(__CUDA)
  ALLOCATE(ps_c(nbnd,pert%nloc))
#endif
  ALLOCATE(dvpsi(npwx*npol,pert%nlocx))
  ALLOCATE(overlap(pert%nglob,nbnd))
  ALLOCATE(pertr(dffts%nnr))
  ALLOCATE(pertg(npwqx))
  !$acc enter data create(dvpsi,overlap,pertr,pertg)
  IF(l_enable_lanczos) THEN
     ALLOCATE(bnorm(pert%nloc))
     ALLOCATE(diago(n_lanczos,pert%nloc))
     ALLOCATE(subdiago(n_lanczos-1,pert%nloc))
     ALLOCATE(q_s(npwx*npol,pert%nloc,n_lanczos))
     ALLOCATE(braket(pert%nglob,n_lanczos,pert%nloc))
     !$acc enter data create(q_s,braket)
     ALLOCATE(z_body2_ifr_q(n_lanczos,pert%nloc,ifr%nloc,band_group%nloc,k_grid%nps,q_grid%nps))
     ALLOCATE(d_diago_q(n_lanczos,pert%nloc,band_group%nloc,k_grid%nps,q_grid%nps))
  ENDIF
  ALLOCATE(l2g(pert%nloc))
  !$acc enter data create(l2g)
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
  IF(l_read_restart .AND. l_enable_lanczos) &
  & CALL read_gfreq(d_diago_q,z_body2_ifr_q,pert%nloc,band_group%nloc)
  !
  ALLOCATE(pertg_all(npwqx,pert%nloc))
  !
  ! LOOP
  !
  ! ... Outer k-point loop (wfc matrix element): ikks, npwk, evck, psick
  ! ... Inner k-point loop (wfc summed over k'): iks, npw, evc (passed to h_psi: current_k = iks)
  !
  DO ikks = 1,k_grid%nps ! KPOINT-SPIN (MATRIX ELEMENT)
     !
     ! Exit loop if no work to do
     !
     IF(barra_load == 0) EXIT
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
     DO iks = 1,k_grid%nps ! KPOINT-SPIN (INTEGRAL OVER K')
        !
        ik = k_grid%ip(iks)
        is = k_grid%is(iks)
        !
        CALL q_grid%find( k_grid%p_cart(:,ikk) - k_grid%p_cart(:,ik), 'cart', iq, g0 )
        !
        npwq = ngq(iq)
        !
        ! compute Coulomb potential
        !
        CALL pot3D%init('Wave',.TRUE.,'default',iq)
        !
        !$acc enter data copyin(pot3D)
        !$acc enter data copyin(pot3D%sqvc)
        !
        ! The Hamiltonian is evaluated at k'
        !
        npw = ngk(iks)
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
!       ! ... Needed for LDA+U
!       !
!       IF ( kpt_pool%nloc > 1 .AND. lda_plus_u .AND. (U_projection /= 'pseudo') ) &
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
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
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
           !
           ib_index = band_group%l2g(ibloc)
           ib = qp_bands(ib_index,is)
           !
           ! PSIC
           !
           IF(noncolin) THEN
              CALL single_invfft_k(dffts,npwk,npwx,evck(1:npwx,ib),psick_nc(:,1),'Wave',igk_k(:,ikks))
              CALL single_invfft_k(dffts,npwk,npwx,evck(1+npwx:npwx*2,ib),psick_nc(:,2),'Wave',igk_k(:,ikks))
           ELSE
              CALL single_invfft_k(dffts,npwk,npwx,evck(:,ib),psick,'Wave',igk_k(:,ikks))
           ENDIF
           !
           !$acc kernels present(dvpsi)
           dvpsi(:,:) = 0._DP
           !$acc end kernels
           !
           DO ip = 1,pert%nloc
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
                    pertr(ir) = CONJG(phase(ir))*psick_nc(ir,1)*CONJG(pertr(ir))
                 ENDDO
                 !$acc end parallel
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1:npwx,ip),'Wave',igk_k(:,current_k))
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                 !$acc parallel loop present(pertr,phase,psick_nc)
                 DO ir = 1,dffts_nnr
                    pertr(ir) = CONJG(phase(ir))*psick_nc(ir,2)*CONJG(pertr(ir))
                 ENDDO
                 !$acc end parallel
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1+npwx:npwx*2,ip),'Wave',igk_k(:,current_k))
              ELSE
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg,pertr,'Wave',igq_q(:,iq))
                 !$acc parallel loop present(pertr,phase,psick)
                 DO ir = 1,dffts_nnr
                    pertr(ir) = CONJG(phase(ir))*psick(ir)*CONJG(pertr(ir))
                 ENDDO
                 !$acc end parallel
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(:,ip),'Wave',igk_k(:,current_k))
              ENDIF
              !
           ENDDO ! pert
           !
           ! OVERLAP( glob_ip, im=1:n_hstates ) = < psi_im iks | dvpsi_glob_ip >
           !
           CALL glbrak_k(evc,dvpsi,ps_c,npw,npwx,nbnd,pert%nloc,nbnd,npol)
           !
           !$acc host_data use_device(ps_c)
           CALL mp_sum(ps_c,intra_bgrp_comm)
           !$acc end host_data
           !
           !$acc kernels present(overlap)
           overlap(:,:) = 0._DP
           !$acc end kernels
           !
           !$acc parallel loop collapse(2) present(overlap,l2g,ps_c)
           DO im = 1,nbnd
              DO ip = 1,pert_nloc
                 overlap(l2g(ip),im) = ps_c(im,ip)
              ENDDO
           ENDDO
           !$acc end parallel
           !
           !$acc update host(overlap)
           CALL mp_sum(overlap,inter_image_comm)
           CALL write_overlap( 'g', kpt_pool%l2g(ikks), kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
           !
           CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
           !
           !$acc update host(dvpsi)
           !
           ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
           ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
           !
           IF(l_enable_lanczos) THEN
              !
              CALL solve_deflated_lanczos_w_full_ortho(nbnd,pert%nloc,n_lanczos,dvpsi,diago,subdiago,q_s,bnorm)
              CALL get_brak_hyper_parallel_complex(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert)
              !
              DO ip = 1,pert%nloc
                 CALL diago_lanczos_complex(bnorm(ip),diago(:,ip),subdiago(:,ip),braket(:,:,ip),pert%nglob)
              ENDDO
              !
              d_diago_q(:,:,ibloc,ikks,iq) = diago
              !
              !$acc enter data create(z_body2_ifr_q(:,:,:,ibloc,ikks,iq)) copyin(z_epsm1_ifr_q(:,:,:,iq))
              !$acc update device(braket)
              !
              !$acc parallel present(braket,z_epsm1_ifr_q(:,:,:,iq),z_body2_ifr_q(:,:,:,ibloc,ikks,iq))
              !$acc loop collapse(3)
              DO ifreq = 1,ifr_nloc
                 DO ip = 1,pert_nloc
                    DO il = 1,n_lanczos
                       reduce = 0._DP
                       !$acc loop reduction(+:reduce)
                       DO glob_jp = 1,pert_nglob
                          reduce = reduce+braket(glob_jp,il,ip)*z_epsm1_ifr_q(glob_jp,ip,ifreq,iq)
                       ENDDO
                       z_body2_ifr_q(il,ip,ifreq,ibloc,ikks,iq) = reduce
                    ENDDO
                 ENDDO
              ENDDO
              !$acc end parallel
              !
              !$acc exit data delete(z_epsm1_ifr_q(:,:,:,iq)) copyout(z_body2_ifr_q(:,:,:,ibloc,ikks,iq))
              !
           ENDIF ! l_enable_lanczos
           !
           CALL update_bar_type( barra, 'glanczos', 1 )
           !
        ENDDO ! BANDS
        !
        !$acc exit data delete(pot3D%sqvc)
        !$acc exit data delete(pot3D)
        !
     ENDDO ! KPOINT-SPIN (INTEGRAL OVER K')
     !
  ENDDO ! KPOINT-SPIN (MATRIX ELEMENT)
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
  CALL deallocate_gw_gpu()
  IF(l_enable_lanczos) CALL deallocate_lanczos_gpu()
#endif
  !
  IF(noncolin) THEN
     !$acc exit data delete(psick_nc)
     DEALLOCATE(psick_nc)
  ELSE
     !$acc exit data delete(psick)
     DEALLOCATE(psick)
  ENDIF
  !$acc exit data delete(phase,evck)
  DEALLOCATE(phase)
  DEALLOCATE(evck)
  !
#if !defined(__CUDA)
  DEALLOCATE(ps_c)
#endif
  !$acc exit data delete(dvpsi,overlap,pertr,pertg)
  DEALLOCATE(dvpsi)
  DEALLOCATE(overlap)
  DEALLOCATE(pertr)
  DEALLOCATE(pertg)
  IF(l_enable_lanczos) THEN
     DEALLOCATE(bnorm)
     DEALLOCATE(diago)
     DEALLOCATE(subdiago)
     !$acc exit data delete(q_s,braket)
     DEALLOCATE(q_s)
     DEALLOCATE(braket)
  ENDIF
  !$acc exit data delete(l2g)
  DEALLOCATE(l2g)
  DEALLOCATE(pertg_all)
  !
  ! Write data
  !
  IF(.NOT. l_read_restart .AND. l_enable_lanczos) &
  & CALL write_gfreq(d_diago_q,z_body2_ifr_q,pert%nloc,band_group%nloc)
  !
  CALL stop_bar_type( barra, 'glanczos' )
  !
END SUBROUTINE
