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
  USE westcom,              ONLY : n_lanczos,npwq,qp_bands,n_bands,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,&
                                 & o_restart_time,npwqx,fftdriver,wstat_save_dir,l_enable_off_diagonal,&
                                 & ijpmap
  USE mp_global,            ONLY : my_image_id,inter_image_comm,npool,intra_bgrp_comm,nbgrp
  USE mp,                   ONLY : mp_bcast,mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : fpi,e2
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
  USE wavefunctions,        ONLY : evc,psic
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
  INTEGER :: ip,ig,glob_ip,ir,ib,ibloc,iks,im,iks_g,jb,ib_index,jb_index,ipair
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  CHARACTER(LEN=25) :: filepot
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:), braket(:, :, :)
  REAL(DP),ALLOCATABLE :: diago1( :, : ), subdiago1( :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:), psic1(:), dvpsi1(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
  REAL(DP),ALLOCATABLE :: ps_r(:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  REAL(DP),ALLOCATABLE :: overlap(:,:)
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bks_type) :: bks
  !
  CALL io_push_title("(G)-Lanczos")
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type( becp )
  CALL allocate_bec_type ( nkb, pert%nloc, becp ) ! I just need 2 becp at a time
  !
  CALL pot3D%init('Wave',.FALSE.,'default')
  CALL band_group%init(n_bands,'b','band_group',.FALSE.)
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
        ib = qp_bands(band_group%l2g(ibloc))
        !
        IF(iks == bks%lastdone_ks .AND. ib <= bks%lastdone_band) CYCLE
        !
        IF (l_enable_off_diagonal) THEN
           barra_load = barra_load + n_bands - ib_index + 1
        ELSE
           barra_load = barra_load + 1
        ENDIF
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
  ! Read PDEP
  !
  ALLOCATE(pertg_all(npwqx,pert%nloc))
  pertg_all = 0._DP
  !
  DO ip = 1,pert%nloc
     glob_ip = pert%l2g(ip)
     CALL generate_pdep_fname(filepot,glob_ip)
     fname = TRIM(wstat_save_dir)//"/"//filepot
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
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
     !npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     !
     bks%max_band = nbndval
     bks%min_band = 1
     !
     ALLOCATE(dvpsi(npwx*npol,pert%nlocx))
     IF (l_enable_off_diagonal) ALLOCATE(dvpsi1(npwx*npol,pert%nlocx))
     !
     time_spent(1) = get_clock( 'glanczos' )
     !
     ! LOOP over band states
     !
     DO ibloc = 1,band_group%nloc
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index)
        !
        IF(iks == bks%lastdone_ks .AND. ib <= bks%lastdone_band) CYCLE
        !
        ! PSIC
        !
        CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ib),psic,'Wave')
        !
        ! ZEROS
        !
        dvpsi = 0._DP
        !
        ALLOCATE( pertg( npwqx ) )
        ALLOCATE( pertr( dffts%nnr ) )
        !
        DO ip=1,pert%nloc
           glob_ip = pert%l2g(ip)
           pertg = pertg_all(:,ip)
           !
           ! Multiply by sqvc
           !pertg(:) = sqvc(:) * pertg(:) ! / SQRT(fpi*e2)     ! CONTROLLARE QUESTO
           DO ig = 1, npwq
              pertg(ig) = pot3D%sqvc(ig) * pertg(ig)
           ENDDO
           !
           ! Bring it to R-space
           CALL single_invfft_gamma(dffts,npwq,npwqx,pertg(1),pertr,TRIM(fftdriver))
           DO ir=1,dffts%nnr
              pertr(ir)=psic(ir)*pertr(ir)
           ENDDO
           CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave')
           !
        ENDDO ! pert
        !
        DEALLOCATE(pertr)
        DEALLOCATE(pertg)
        !
        ! OVERLAP( glob_ip, im=1:n_hstates ) = < psi_im iks | dvpsi_glob_ip >
        !
        IF(ALLOCATED(ps_r)) DEALLOCATE(ps_r)
        ALLOCATE(ps_r(nbnd,pert%nloc))
        CALL glbrak_gamma(evc,dvpsi,ps_r,npw,npwx,nbnd,pert%nloc,nbnd,npol)
        CALL mp_sum(ps_r,intra_bgrp_comm)
        !
        IF(ALLOCATED(overlap)) DEALLOCATE(overlap)
        ALLOCATE(overlap(pert%nglob, nbnd ) )
        overlap = 0._DP
        DO im = 1, nbnd
           DO ip = 1, pert%nloc
              glob_ip = pert%l2g(ip)
              overlap(glob_ip,im) = ps_r(im,ip)
           ENDDO
        ENDDO
        !
        DEALLOCATE(ps_r)
        CALL mp_sum(overlap,inter_image_comm)
        CALL writeout_overlap( 'g', kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
        DEALLOCATE(overlap)
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF( l_enable_lanczos ) THEN
           !
           ALLOCATE( bnorm    (                         pert%nloc ) )
           ALLOCATE( diago    (            n_lanczos  , pert%nloc ) )
           ALLOCATE( subdiago (            n_lanczos-1, pert%nloc ) )
           ALLOCATE( q_s      ( npwx*npol, pert%nloc  , n_lanczos ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM
           !
           CALL solve_deflated_lanczos_w_full_ortho ( nbnd, pert%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
           !
           DO jb_index = 1, n_bands
              !
              jb = qp_bands(jb_index)
              !
              IF (l_enable_off_diagonal .and. jb <= ib) THEN     
                 !
                 ipair = ijpmap(jb_index,ib_index)
                 !
                 ! PSIC
                 !
                 ALLOCATE( psic1(dffts%nnr) )
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(1,jb),psic1,'Wave') 
                 !
                 ! ZEROS
                 !
                 dvpsi1 = 0._DP   
                 !
                 ! Read PDEP
                 !
                 ALLOCATE( pertg( npwqx ) )
                 ALLOCATE( pertr( dffts%nnr ) )
                 !
                 DO ip=1,pert%nloc
                    glob_ip = pert%l2g(ip)
                    pertg = pertg_all(:,ip)
                    !
                    ! Multiply by sqvc
                    DO ig = 1, npwq
                       pertg(ig) = pot3D%sqvc(ig) * pertg(ig)
                    ENDDO
                    !
                    ! Bring it to R-space
                    CALL single_invfft_gamma(dffts,npwq,npwqx,pertg(1),pertr,TRIM(fftdriver))
                    DO ir=1,dffts%nnr
                       pertr(ir)=psic1(ir)*pertr(ir)
                    ENDDO
                    CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi1(1,ip),'Wave')
                    !
                 ENDDO ! pert
                 ! 
                 DEALLOCATE(pertr)
                 DEALLOCATE(pertg)
                 DEALLOCATE(psic1)
                 !
                 CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi1,(1._DP,0._DP))
                 !
                 ALLOCATE( diago1    (               n_lanczos   , pert%nloc ) )
                 ALLOCATE( subdiago1 (               n_lanczos-1 , pert%nloc ) )               
                 ALLOCATE( braket   ( pert%nglob, n_lanczos   , pert%nloc ) )
                 !
                 diago1 = diago
                 subdiago1 = subdiago
                 !
                 CALL get_brak_hyper_parallel(dvpsi1,pert%nloc,n_lanczos,q_s,braket,pert)
                 !
                 DO ip = 1, pert%nloc
                    CALL diago_lanczos( bnorm(ip), diago1( :, ip), subdiago1( :, ip), braket(:,:,ip), pert%nglob )
                 ENDDO
                 !
                 DEALLOCATE( subdiago1 )
                 !
                 CALL writeout_solvegfreq( kpt_pool%l2g(iks), ipair, diago1, braket, pert%nloc, pert%nglob, pert%myoffset )
                 ! 
                 DEALLOCATE( diago1, braket )
                 !
                 CALL update_bar_type( barra, 'glanczos', 1 )
                 !
              ELSEIF ( .NOT. l_enable_off_diagonal .AND. jb == ib ) THEN
                 !
                 ALLOCATE( braket   ( pert%nglob, n_lanczos   , pert%nloc ) )
                 CALL get_brak_hyper_parallel(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert)
                 DEALLOCATE( q_s )
                 !
                 DO ip = 1, pert%nloc
                    CALL diago_lanczos( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), pert%nglob )
                 ENDDO
                 !
                 DEALLOCATE( bnorm )
                 DEALLOCATE( subdiago )
                 !
                 ! MPI-IO
                 !
                 CALL writeout_solvegfreq( kpt_pool%l2g(iks), ib, diago, braket, pert%nloc, pert%nglob, pert%myoffset )
                 !
                 DEALLOCATE( diago )
                 DEALLOCATE( braket )
                 !
                 CALL update_bar_type( barra, 'glanczos', 1 )
                 !
              ENDIF
              !
           ENDDO
           !
           IF (l_enable_off_diagonal) DEALLOCATE(bnorm, diago, subdiago, q_s)
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'glanczos' )
        l_write_restart = .FALSE.
        !
        IF( o_restart_time >= 0._DP ) THEN
           IF( time_spent(2)-time_spent(1) > o_restart_time*60._DP ) l_write_restart = .TRUE.
           IF( ib == qp_bands(n_bands) ) l_write_restart = .TRUE.
        ENDIF
        !
        ! Write final restart file
        !
        IF( iks == k_grid%nps .AND. ib == qp_bands(n_bands) ) l_write_restart = .TRUE.
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
     ENDDO ! BANDS
     !
     DEALLOCATE(dvpsi)
     IF (l_enable_off_diagonal) DEALLOCATE(dvpsi1)
     !
  ENDDO ! KPOINT-SPIN
  !
  DEALLOCATE(pertg_all)
  !
  ! Write final restart file when using pool or band group
  !
  IF( npool*nbgrp > 1 ) THEN
     bks%lastdone_ks = k_grid%nps
     bks%lastdone_band = qp_bands(n_bands)
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
  USE westcom,              ONLY : n_lanczos,npwq,qp_bands,n_bands,l_enable_lanczos,&
                                 & nbnd_occ,iuwfc,lrwfc,o_restart_time,npwqx,&
                                 & wstat_save_dir,ngq,igq_q
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nbgrp
  USE mp,                   ONLY : mp_bcast,mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : fpi,e2
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
  USE wavefunctions,        ONLY : evc
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
  INTEGER :: ip,ig,glob_ip,ir,ib,ibloc,iks,ik,im,ikks,ikk,iq,ib_index
  INTEGER :: npwk
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  CHARACTER(LEN=25)     :: filepot
  INTEGER :: nbndval
  REAL(DP) :: g0(3)
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:)
  COMPLEX(DP),ALLOCATABLE :: braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:), pertr(:)
  COMPLEX(DP), ALLOCATABLE :: evck(:,:), phase(:)
  COMPLEX(DP), ALLOCATABLE :: psick(:), psick_nc(:,:)
  COMPLEX(DP),ALLOCATABLE :: ps_c(:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bksks_type) :: bksks
  !
  CALL io_push_title("(G)-Lanczos")
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type( becp )
  CALL allocate_bec_type ( nkb, pert%nloc, becp ) ! I just need 2 becp at a time
  !
  CALL band_group%init(n_bands,'b','band_group',.FALSE.)
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
  ALLOCATE( evck(npwx*npol,nbnd) )
  IF ( noncolin ) THEN
     ALLOCATE( psick_nc( dffts%nnr, npol ) )
  ELSE
     ALLOCATE( psick(dffts%nnr) )
  ENDIF
  ALLOCATE( phase(dffts%nnr) )
  !
  barra_load = 0
  DO ikks = 1,k_grid%nps
     IF(ikks < bksks%lastdone_ks) CYCLE
     !
     DO ibloc = 1,band_group%nloc
        ib = qp_bands(band_group%l2g(ibloc))
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
  ! LOOP
  !
  ! ... Outer k-point loop (wfc matrix element): ikks, npwk, evck, psick
  ! ... Inner k-point loop (wfc summed over k'): iks, npw, evc (passed to h_psi: current_k = iks)
  !
  DO ikks = 1, k_grid%nps   ! KPOINT-SPIN (MATRIX ELEMENT)
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
     IF(my_image_id==0) CALL get_buffer( evck, lrwfc, iuwfc, ikks )
     CALL mp_bcast(evck,0,inter_image_comm)
     !
     nbndval = nbnd_occ(ikks)
     !
     bksks%max_band=nbndval
     bksks%min_band=1
     !
     ALLOCATE(dvpsi(npwx*npol,pert%nlocx))
     !
     ! LOOP over band states
     !
     DO ibloc = 1,band_group%nloc
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index)
        !
        IF(ikks == bksks%lastdone_ks .AND. ib < bksks%lastdone_band) CYCLE
        !
        ! PSIC
        !
        IF(noncolin) THEN
           CALL single_invfft_k(dffts,npwk,npwx,evck(1     ,ib),psick_nc(1,1),'Wave',igk_k(1,ikks))
           CALL single_invfft_k(dffts,npwk,npwx,evck(1+npwx,ib),psick_nc(1,2),'Wave',igk_k(1,ikks))
        ELSE
           CALL single_invfft_k(dffts,npwk,npwx,evck(1,ib),psick,'Wave',igk_k(1,ikks))
        ENDIF
        !
        DO iks = 1, k_grid%nps ! KPOINT-SPIN (INTEGRAL OVER K')
           IF (ikks==bksks%lastdone_ks .AND. ib==bksks%lastdone_band .AND. iks <= bksks%lastdone_kks) CYCLE
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
           IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
           !
!          !
!          ! ... Needed for LDA+U
!          !
!          IF ( kpt_pool%nloc > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
!               CALL get_buffer ( wfcU, nwordwfcU, iunhub, iks )
!          !
!          current_k = iks
!          current_spin = isk(iks)
!          !
!          CALL gk_sort(xk(1,iks),ngm,g,gcutw,npw,igk,g2kin)
!          g2kin=g2kin*tpiba2
!          !
!          ! reads unperturbed wavefuctions psi_k in G_space, for all bands
!          !
!          !
!          CALL init_us_2 (npw, igk, xk (1, iks), vkb)
           !
           CALL compute_phase( g0, 'cart', phase )
           !
           IF ( my_image_id == 0 ) CALL get_buffer( evc, lrwfc, iuwfc, iks )
           CALL mp_bcast( evc, 0, inter_image_comm )
           !
           ! ZEROS
           !
           dvpsi = 0._DP
           !
           ! Read PDEP
           !
           ALLOCATE( pertg( npwqx ) )
           ALLOCATE( pertr( dffts%nnr ) )
           !
           DO ip=1,pert%nloc
              glob_ip = pert%l2g(ip)
              !
              ! Exhume dbs eigenvalue
              !
              CALL generate_pdep_fname( filepot, glob_ip, iq )
              fname = TRIM( wstat_save_dir ) // "/"// filepot
              CALL pdep_read_G_and_distribute(fname,pertg,iq)
              !
              ! Multiply by sqvc
              !pertg(:) = sqvc(:) * pertg(:) ! / SQRT(fpi*e2)     ! CONTROLLARE QUESTO
              DO ig = 1, npwq
                 pertg(ig) = pot3D%sqvc(ig) * pertg(ig)
              ENDDO
              !
              ! Bring it to R-space
              IF(noncolin) THEN
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                 DO ir=1,dffts%nnr
                    pertr(ir)=CONJG(phase(ir))*psick_nc(ir,1)*CONJG(pertr(ir))
                 ENDDO
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                 DO ir=1,dffts%nnr
                    pertr(ir)=CONJG(phase(ir))*psick_nc(ir,2)*CONJG(pertr(ir))
                 ENDDO
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1+npwx,ip),'Wave',igk_k(1,current_k))
              ELSE
                 CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                 DO ir=1,dffts%nnr
                    pertr(ir)=CONJG(phase(ir))*psick(ir)*CONJG(pertr(ir))
                 ENDDO
                 CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
              ENDIF
              !
           ENDDO ! pert
           !
           DEALLOCATE(pertr)
           DEALLOCATE(pertg)
           !
           ! OVERLAP( glob_ip, im=1:n_hstates ) = < psi_im iks | dvpsi_glob_ip >
           !
           IF(ALLOCATED(ps_c)) DEALLOCATE(ps_c)
           ALLOCATE(ps_c(nbnd,pert%nloc))
           CALL glbrak_k(evc,dvpsi,ps_c,npw,npwx,nbnd,pert%nloc,nbnd,npol)
           CALL mp_sum(ps_c,intra_bgrp_comm)
           !
           IF(ALLOCATED(overlap)) DEALLOCATE(overlap)
           ALLOCATE(overlap(pert%nglob, nbnd ) )
           overlap = 0._DP
           DO im = 1, nbnd
              DO ip = 1, pert%nloc
                 glob_ip = pert%l2g(ip)
                 overlap(glob_ip,im) = ps_c(im,ip)
              ENDDO
           ENDDO
           !
           DEALLOCATE(ps_c)
           CALL mp_sum(overlap,inter_image_comm)
           CALL writeout_overlap( 'g', kpt_pool%l2g(ikks), kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
           DEALLOCATE(overlap)
           !
           CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
           !
           ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
           ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
           !
           IF( l_enable_lanczos ) THEN
              !
              ALLOCATE( bnorm    (                         pert%nloc ) )
              ALLOCATE( diago    (            n_lanczos  , pert%nloc ) )
              ALLOCATE( subdiago (            n_lanczos-1, pert%nloc ) )
              ALLOCATE( q_s      ( npwx*npol, pert%nloc  , n_lanczos ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM
              !
              CALL solve_deflated_lanczos_w_full_ortho(nbnd, pert%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
              ALLOCATE( braket( pert%nglob, n_lanczos, pert%nloc) )
              CALL get_brak_hyper_parallel_complex(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert)
              !
              DO ip = 1, pert%nloc
                 CALL diago_lanczos_complex( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), pert%nglob )
              ENDDO
              !
              DEALLOCATE( q_s )
              DEALLOCATE( bnorm )
              DEALLOCATE( subdiago )
              !
              ! MPI-IO
              !
              CALL writeout_solvegfreq( kpt_pool%l2g(ikks), kpt_pool%l2g(iks), ib, diago, braket, pert%nloc, &
              & pert%nglob, pert%myoffset )
              !
              DEALLOCATE( diago )
              DEALLOCATE( braket )
              !
           ENDIF ! l_enable_lanczos
           !
           time_spent(2) = get_clock( 'glanczos' )
           l_write_restart = .FALSE.
           !
           IF( o_restart_time >= 0._DP ) THEN
              IF( time_spent(2)-time_spent(1) > o_restart_time*60._DP ) l_write_restart = .TRUE.
              IF( ib == qp_bands(n_bands) ) l_write_restart = .TRUE.
           ENDIF
           !
           ! Write final restart file
           !
           IF( ikks == k_grid%nps .AND. iks == k_grid%nps .AND. ib == qp_bands(n_bands) ) l_write_restart = .TRUE.
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
        ENDDO ! KPOINT-SPIN (INTEGRAL OVER K')
        !
     ENDDO ! BANDS
     !
     DEALLOCATE(dvpsi)
     !
  ENDDO ! KPOINT-SPIN (MATRIX ELEMENT)
  !
  DEALLOCATE( phase )
  IF ( noncolin ) THEN
     DEALLOCATE( psick_nc )
  ELSE
     DEALLOCATE( psick )
  ENDIF
  DEALLOCATE( evck )
  !
  ! Write final restart file when using band group
  !
  IF( nbgrp > 1 ) THEN
     bksks%lastdone_ks = k_grid%nps
     bksks%lastdone_kks = k_grid%nps
     bksks%lastdone_band = qp_bands(n_bands)
     CALL solvegfreq_restart_write_q( bksks )
  ENDIF
  !
  CALL stop_bar_type( barra, 'glanczos' )
  !
  CALL mp_barrier( world_comm )
  !
END SUBROUTINE
