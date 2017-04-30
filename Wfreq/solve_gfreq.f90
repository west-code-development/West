!
! Copyright (C) 2015-2016 M. Govoni 
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
  USE westcom,              ONLY : sqvc,west_prefix,n_pdep_eigen_to_use,n_lanczos,npwq0,qp_bandrange,iks_l2g,&
                                 & wfreq_dirname,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,o_restart_time,npwq0x,fftdriver
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm,inter_pool_comm
  USE mp,                   ONLY : mp_bcast,mp_barrier,mp_sum
  USE io_global,            ONLY : stdout, ionode
  USE gvect,                ONLY : g,ngm,gstart
  USE gvecw,                ONLY : gcutw
  USE cell_base,            ONLY : tpiba2,bg
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : tpi,fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,nkstot,current_k,ngk
  USE wavefunctions_module, ONLY : evc,psic,psic_nc
  USE io_files,             ONLY : tmp_dir,nwordwfc,iunwfc
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
!  USE fft_at_k,             ONLY : SINGLEBAND_INVFFT_k,SINGLEBAND_FWFFT_k
  USE becmod,               ONLY : becp,allocate_bec_type,deallocate_bec_type
  USE uspp,                 ONLY : vkb,nkb
  USE pdep_io,              ONLY : pdep_read_G_and_distribute 
  USE io_push,              ONLY : io_push_title
!  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol 
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert
  USE wfreq_restart,        ONLY : solvegfreq_restart_write,solvegfreq_restart_read,bks_type
  USE wfreq_io,             ONLY : writeout_overlap,writeout_solvegfreq,preallocate_solvegfreq
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i3,ip,ig,glob_ip,ir,ib,iks,m,im
  CHARACTER(LEN=256)    :: wstat_dirname, fname
  CHARACTER(LEN=6)      :: my_label_b
  COMPLEX(DP),ALLOCATABLE :: auxr(:)
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:), braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  REAL(DP),ALLOCATABLE :: ps_r(:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  REAL(DP),ALLOCATABLE :: overlap(:,:)
  LOGICAL :: l_iks_skip, l_ib_skip
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
  ! Remember the directory names
  !
  wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
  !
  IF(l_read_restart) THEN
     CALL solvegfreq_restart_read( bks )
  ELSE
     bks%lastdone_ks   = 0 
     bks%lastdone_band = 0 
     bks%old_ks        = 0 
     bks%old_band      = 0 
     bks%max_ks        = nks 
     bks%min_ks        = 1 
  ENDIF
  !
  barra_load = 0
  DO iks = 1, nks
     IF(iks<bks%lastdone_ks) CYCLE
     DO ib = qp_bandrange(1), qp_bandrange(2)
        IF(iks==bks%lastdone_ks .AND. ib <= bks%lastdone_band ) CYCLE
        barra_load = barra_load + 1
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
  DO iks = 1, nks   ! KPOINT-SPIN
     IF(iks<bks%lastdone_ks) CYCLE
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
     IF(nks>1) THEN
        !iuwfc = 20
        !lrwfc = nbnd * npwx * npol 
        !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        !CALL mp_bcast(evc,0,inter_image_comm)
        !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
!     !
!     ! ... Needed for LDA+U
!     !
!     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
!          CALL get_buffer ( wfcU, nwordwfcU, iunhub, iks )
!     !
!     current_k = iks
!     current_spin = isk(iks)
!     !
!     CALL gk_sort(xk(1,iks),ngm,g,gcutw,npw,igk,g2kin)
!     g2kin=g2kin*tpiba2
!     !
!     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
!     !
!     !
!     CALL init_us_2 (npw, igk, xk (1, iks), vkb)
     !
     nbndval = nbnd_occ(iks)
     !
     bks%max_band=nbndval
     bks%min_band=1
     !
     ALLOCATE(dvpsi(npwx*npol,pert%nlocx))   
     CALL preallocate_solvegfreq( iks_l2g(iks), qp_bandrange(1), qp_bandrange(2), pert )
     !
     time_spent(1) = get_clock( 'glanczos' ) 
     !
     ! LOOP over band states 
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        IF(iks==bks%lastdone_ks .AND. ib <= bks%lastdone_band ) CYCLE
        !
        ! PSIC
        !
!        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ib),psic,'Wave') 
!        ELSEIF(noncolin) THEN
!           CALL SINGLEBAND_invfft_k(npw,evc(1     ,ib),npwx,psic_nc(1,1),dffts%nnr,.TRUE.)
!           CALL SINGLEBAND_invfft_k(npw,evc(1+npwx,ib),npwx,psic_nc(1,2),dffts%nnr,.TRUE.)
!        ELSE
!           CALL SINGLEBAND_invfft_k(npw,evc(1,ib),npwx,psic,dffts%nnr,.TRUE.)
!        ENDIF
        !
        ! ZEROS
        !
        dvpsi = 0._DP   
        !
        ! Read PDEP
        !
        ALLOCATE( pertg(npwq0x) )
        ALLOCATE( pertr( dffts%nnr ) )
        !
        DO ip=1,pert%nloc
           glob_ip = pert%l2g(ip)
           !
           ! Exhume dbs eigenvalue
           !
           WRITE(my_label_b,'(i6.6)') glob_ip
           fname = TRIM( wstat_dirname ) // "/E"//TRIM(ADJUSTL(my_label_b))//".dat"
           CALL pdep_read_G_and_distribute(fname,pertg)
           !
           ! Multiply by sqvc
           !pertg(:) = sqvc(:) * pertg(:) ! / SQRT(fpi*e2)     ! CONTROLLARE QUESTO
           DO ig = 1, npwq0
              pertg(ig) = sqvc(ig) * pertg(ig)
           ENDDO
           !
           ! Bring it to R-space
!           IF(gamma_only) THEN
              CALL single_invfft_gamma(dffts,npwq0,npwq0x,pertg(1),pertr,TRIM(fftdriver))
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic(ir)*pertr(ir)
              ENDDO
              CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave')
!           ELSEIF(noncolin) THEN
!              CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
!              DO ir=1,dffts%nnr 
!                 pertr(ir)=psic_nc(ir,1)*DCONJG(pertr(ir))
!              ENDDO
!              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
!              CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
!              DO ir=1,dffts%nnr 
!                 pertr(ir)=psic_nc(ir,2)*DCONJG(pertr(ir))
!              ENDDO
!              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1+npwx,ip),npwx,.TRUE.)
!           ELSE
!              CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
!              DO ir=1,dffts%nnr 
!                 pertr(ir)=psic(ir)*DCONJG(pertr(ir))
!              ENDDO
!              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
!           ENDIF 
           !
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
        CALL writeout_overlap( 'g', iks_l2g(iks), ib, overlap, pert%nglob, nbnd )
        DEALLOCATE(overlap)
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF( l_enable_lanczos ) THEN  
           !
           ALLOCATE( bnorm    (                             pert%nloc ) )
           ALLOCATE( diago    (               n_lanczos   , pert%nloc ) )
           ALLOCATE( subdiago (               n_lanczos-1 , pert%nloc ) )
           ALLOCATE( q_s      ( npwx*npol   , pert%nloc   , n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM 
           !
           CALL solve_deflated_lanczos_w_full_ortho ( nbnd, pert%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
           !
           ALLOCATE( braket   ( pert%nglob, n_lanczos   , pert%nloc ) )
           CALL get_brak_hyper_parallel(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert%nloc,pert%nlocx,pert%nglob)
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
           CALL writeout_solvegfreq( iks_l2g(iks), ib, diago, braket, pert%nloc, pert%nglob, pert%myoffset )
           !
           DEALLOCATE( diago ) 
           DEALLOCATE( braket )
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'glanczos' ) 
        !
        IF( o_restart_time >= 0._DP ) THEN 
           IF( (time_spent(2)-time_spent(1)) > o_restart_time*60._DP .OR. ib == qp_bandrange(2) ) THEN 
              bks%lastdone_ks=iks
              bks%lastdone_band=ib 
              CALL solvegfreq_restart_write( bks )
              bks%old_ks=iks
              bks%old_band=ib
              time_spent(1) = get_clock( 'glanczos' )
           ENDIF
        ENDIF
        !
        CALL update_bar_type( barra, 'glanczos', 1 )
        !
     ENDDO ! BANDS
     !
     DEALLOCATE(dvpsi)
     !
  ENDDO ! KPOINT-SPIN 
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
  USE westcom,              ONLY : sqvc,west_prefix,n_pdep_eigen_to_use,n_lanczos,npwq0,qp_bandrange,iks_l2g,&
                                 & wfreq_dirname,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,o_restart_time,npwq0x,fftdriver
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm,inter_pool_comm
  USE mp,                   ONLY : mp_bcast,mp_barrier,mp_sum
  USE io_global,            ONLY : stdout, ionode
  USE gvect,                ONLY : g,ngm,gstart
  USE gvecw,                ONLY : gcutw
  USE cell_base,            ONLY : tpiba2,bg
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : tpi,fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,nkstot,current_k,ngk
  USE wavefunctions_module, ONLY : evc,psic,psic_nc
  USE io_files,             ONLY : tmp_dir,nwordwfc,iunwfc
!  USE fft_at_gamma,         ONLY : DOUBLEBAND_INVFFT,SINGLEBAND_INVFFT,DOUBLEBAND_FWFFT,SINGLEBAND_FWFFT
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE becmod,               ONLY : becp,allocate_bec_type,deallocate_bec_type
  USE uspp,                 ONLY : vkb,nkb
  USE pdep_io,              ONLY : pdep_read_G_and_distribute 
  USE io_push,              ONLY : io_push_title
!  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol 
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert
  USE wfreq_restart,        ONLY : solvegfreq_restart_write,solvegfreq_restart_read,bks_type
  USE wfreq_io,             ONLY : writeout_overlap,writeout_solvegfreq,preallocate_solvegfreq
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i3,ip,ig,glob_ip,ir,ib,iks,m,im
  CHARACTER(LEN=256)    :: wstat_dirname, fname
  CHARACTER(LEN=6)      :: my_label_b
  COMPLEX(DP),ALLOCATABLE :: auxr(:)
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:)
  COMPLEX(DP),ALLOCATABLE :: braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP),ALLOCATABLE :: ps_c(:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
  LOGICAL :: l_iks_skip, l_ib_skip
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
  ! Remember the directory names
  !
  wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
  !
  IF(l_read_restart) THEN
     CALL solvegfreq_restart_read( bks )
  ELSE
     bks%lastdone_ks   = 0 
     bks%lastdone_band = 0 
     bks%old_ks        = 0 
     bks%old_band      = 0 
     bks%max_ks        = nks 
     bks%min_ks        = 1 
  ENDIF
  !
  barra_load = 0
  DO iks = 1, nks
     IF(iks<bks%lastdone_ks) CYCLE
     DO ib = qp_bandrange(1), qp_bandrange(2)
        IF(iks==bks%lastdone_ks .AND. ib <= bks%lastdone_band ) CYCLE
        barra_load = barra_load + 1
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
  DO iks = 1, nks   ! KPOINT-SPIN
     IF(iks<bks%lastdone_ks) CYCLE
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
     IF(nks>1) THEN
        !iuwfc = 20
        !lrwfc = nbnd * npwx * npol 
        !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        !CALL mp_bcast(evc,0,inter_image_comm)
        !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
!     !
!     ! ... Needed for LDA+U
!     !
!     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
!          CALL get_buffer ( wfcU, nwordwfcU, iunhub, iks )
!     !
!     current_k = iks
!     current_spin = isk(iks)
!     !
!     CALL gk_sort(xk(1,iks),ngm,g,gcutw,npw,igk,g2kin)
!     g2kin=g2kin*tpiba2
!     !
!     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
!     !
!     !
!     CALL init_us_2 (npw, igk, xk (1, iks), vkb)
     !
     nbndval = nbnd_occ(iks)
     !
     bks%max_band=nbndval
     bks%min_band=1
     !
     ALLOCATE(dvpsi(npwx*npol,pert%nlocx))   
     CALL preallocate_solvegfreq( iks_l2g(iks), qp_bandrange(1), qp_bandrange(2), pert )
     !
     time_spent(1) = get_clock( 'glanczos' ) 
     !
     ! LOOP over band states 
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        IF(iks==bks%lastdone_ks .AND. ib <= bks%lastdone_band ) CYCLE
        !
        ! PSIC
        !
!        IF(gamma_only) THEN
!           CALL SINGLEBAND_invfft(npw,evc(1,ib),npwx,psic,dffts%nnr) 
!        ELSE
        IF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1     ,ib),psic_nc(1,1),'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx,ib),psic_nc(1,2),'Wave',igk_k(1,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic,'Wave',igk_k(1,current_k))
        ENDIF
        !
        ! ZEROS
        !
        dvpsi = 0._DP   
        !
        ! Read PDEP
        !
        ALLOCATE( pertg(npwq0x) )
        ALLOCATE( pertr( dffts%nnr ) )
        !
        DO ip=1,pert%nloc
           glob_ip = pert%l2g(ip)
           !
           ! Exhume dbs eigenvalue
           !
           WRITE(my_label_b,'(i6.6)') glob_ip
           fname = TRIM( wstat_dirname ) // "/E"//TRIM(ADJUSTL(my_label_b))//".dat"
           CALL pdep_read_G_and_distribute(fname,pertg)
           !
           ! Multiply by sqvc
           !pertg(:) = sqvc(:) * pertg(:) ! / SQRT(fpi*e2)     ! CONTROLLARE QUESTO
           DO ig = 1, npwq0
              pertg(ig) = sqvc(ig) * pertg(ig)
           ENDDO
           !
           ! Bring it to R-space
!           IF(gamma_only) THEN
!              CALL SINGLEBAND_invfft(npwq0,pertg(1),npwx,pertr,dffts%nnr)
!              DO ir=1,dffts%nnr 
!                 pertr(ir)=psic(ir)*pertr(ir)
!              ENDDO
!              CALL SINGLEBAND_fwfft(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx)
!           ELSE
            IF(noncolin) THEN
              CALL single_invfft_k(dffts,npwq0,npwq0x,pertg(1),pertr,TRIM(fftdriver)) ! no igk
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic_nc(ir,1)*DCONJG(pertr(ir))
              ENDDO
              CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
              CALL single_invfft_k(dffts,npwq0,npwq0x,pertg(1),pertr,TRIM(fftdriver),igk_k(1,current_k))
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic_nc(ir,2)*DCONJG(pertr(ir))
              ENDDO
              CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1+npwx,ip),'Wave',igk_k(1,current_k))
           ELSE
              CALL single_invfft_k(dffts,npwq0,npwq0x,pertg(1),pertr,TRIM(fftdriver)) ! no igk
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic(ir)*DCONJG(pertr(ir))
              ENDDO
              CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
           ENDIF 
           !
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
        CALL writeout_overlap( 'g', iks_l2g(iks), ib, overlap, pert%nglob, nbnd )
        DEALLOCATE(overlap)
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF( l_enable_lanczos ) THEN  
           !
           ALLOCATE( bnorm    (                             pert%nloc ) )
           ALLOCATE( diago    (               n_lanczos   , pert%nloc ) )
           ALLOCATE( subdiago (               n_lanczos-1 , pert%nloc ) )
           ALLOCATE( q_s      ( npwx*npol   , pert%nloc   , n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM 
           !
           CALL solve_deflated_lanczos_w_full_ortho ( nbnd, pert%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
           !
           ALLOCATE( braket   ( pert%nglob, n_lanczos   , pert%nloc ) )
           CALL get_brak_hyper_parallel_complex(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert%nloc,pert%nlocx,pert%nglob)
           DEALLOCATE( q_s )
           !
           DO ip = 1, pert%nloc
              CALL diago_lanczos_complex( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), pert%nglob )
           ENDDO
           !
           DEALLOCATE( bnorm )
           DEALLOCATE( subdiago )
           !
           ! MPI-IO 
           !
           CALL writeout_solvegfreq( iks_l2g(iks), ib, diago, braket, pert%nloc, pert%nglob, pert%myoffset )
           !
           DEALLOCATE( diago ) 
           DEALLOCATE( braket )
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'glanczos' ) 
        !
        IF( o_restart_time >= 0._DP ) THEN 
           IF( (time_spent(2)-time_spent(1)) > o_restart_time*60._DP .OR. ib == qp_bandrange(2) ) THEN 
              bks%lastdone_ks=iks
              bks%lastdone_band=ib 
              CALL solvegfreq_restart_write( bks )
              bks%old_ks=iks
              bks%old_band=ib
              time_spent(1) = get_clock( 'glanczos' )
           ENDIF
        ENDIF
        !
        CALL update_bar_type( barra, 'glanczos', 1 )
        !
     ENDDO ! BANDS
     !
     DEALLOCATE(dvpsi)
     !
  ENDDO ! KPOINT-SPIN 
  ! 
  CALL stop_bar_type( barra, 'glanczos' )
  !
END SUBROUTINE 
