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
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_lanczos,npwq,l_macropol,d_epsm1_ifr,z_epsm1_rfr,&
                                 & l_enable_lanczos,iuwfc,lrwfc,wfreq_eta,imfreq_list,refreq_list,tr2_dfpt,&
                                 & z_head_rfr,d_head_ifr,o_restart_time,l_skip_nl_part_of_hcomr,npwqx,&
                                 & fftdriver,wstat_save_dir,l_frac_occ,occupation,nbnd_occ,nbnd_occ_full,&
                                 & qp_bands,d_epsm1_ifr_a,z_epsm1_rfr_a,z_head_rfr_a,d_head_ifr_a
  USE mp_global,            ONLY : my_image_id,inter_image_comm,inter_pool_comm,npool,intra_bgrp_comm,&
                                 & inter_bgrp_comm,nbgrp
  USE mp,                   ONLY : mp_bcast,mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : bg,omega
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,current_spin,isk,xk,nbnd,lsda,igk_k,current_k,ngk
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
  USE distribution_center,  ONLY : pert,macropert,ifr,rfr,occband,band_group,kpt_pool
  USE class_idistribute,    ONLY : idistribute
  USE wfreq_restart,        ONLY : solvewfreq_restart_write,solvewfreq_restart_read,bks_type
  USE types_bz_grid,        ONLY : k_grid
  USE chi_invert,           ONLY : chi_invert_real,chi_invert_complex
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart,l_generate_plot,l_QDET
  !
  ! Workspace
  !
  LOGICAL :: l_write_restart
  INTEGER :: i1,i2,ip,ig,glob_ip,ir,iv,ivloc,ivloc2,ifloc,iks,ipol,iks_g
  CHARACTER(LEN=25) :: filepot
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  INTEGER :: nbndval, nbndval_full
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:), braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi(:,:), phis(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: phi_tmp(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
  COMPLEX(DP) :: zkonstant
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  REAL(DP),ALLOCATABLE :: ps_r(:,:)
  TYPE(idistribute) :: mypara
  REAL(DP),ALLOCATABLE :: overlap(:,:)
  REAL(DP) :: mwo, ecv, dfactor, frequency, dhead
  COMPLEX(DP) :: zmwo, zfactor,zm,zp, zhead, zhead_a
  INTEGER :: glob_jp,ic,ifreq,il
  REAL(DP),ALLOCATABLE :: dmatilda(:,:), dlambda(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatilda(:,:), zlambda(:,:), zlambda_a(:,:)
  REAL(DP),ALLOCATABLE :: dmati(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr(:,:,:)
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bks_type) :: bks
  REAL(DP),ALLOCATABLE :: eprec(:)
  INTEGER :: ierr
  REAL(DP),ALLOCATABLE :: e(:)
  REAL(DP) :: docc
  !
  COMPLEX(DP),ALLOCATABLE :: evc_a(:,:)
  REAL(DP),ALLOCATABLE :: dmati_a(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr_a(:,:,:)
  REAL(DP),ALLOCATABLE :: dmatilda_a(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatilda_a(:,:), zmatilda_diff(:,:)
  !
  CALL io_push_title("(W)-Lanczos")
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
  CALL deallocate_bec_type( becp )
  IF(l_macropol) THEN
     CALL allocate_bec_type ( nkb, MAX(mypara%nloc,3), becp ) ! I just need 2 becp at a time
  ELSE
     CALL allocate_bec_type ( nkb, mypara%nloc, becp ) ! I just need 2 becp at a time
  ENDIF
  !
  ! ALLOCATE dmati, zmatr, where chi0 is stored
  !
  ALLOCATE( dmati( mypara%nglob, mypara%nloc, ifr%nloc) )
  ALLOCATE( zmatr( mypara%nglob, mypara%nloc, rfr%nloc) )
  dmati = 0._DP
  zmatr = 0._DP
  !
  ALLOCATE( evc_a(npwx*npol, nbnd) )
  IF(l_QDET) THEN
     ALLOCATE( dmati_a( mypara%nglob, mypara%nloc, ifr%nloc) )
     ALLOCATE( zmatr_a( mypara%nglob, mypara%nloc, rfr%nloc) )
     dmati_a = 0._DP
     zmatr_a = 0._DP
  ENDIF
  !
  IF(l_read_restart) THEN
     CALL solvewfreq_restart_read( bks, dmati, zmatr, mypara%nglob, mypara%nloc )
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
     CALL band_group%init(nbnd_occ(iks),'b','band_group',.FALSE.)
     !
     DO ivloc = 1,band_group%nloc
        iv = band_group%l2g(ivloc)
        !
        IF(iks == bks%lastdone_ks .AND. iv <= bks%lastdone_band) CYCLE
        !
        barra_load = barra_load+1
     ENDDO
  ENDDO
  !
  IF( barra_load == 0 ) THEN
     CALL start_bar_type ( barra, 'wlanczos', 1 )
     CALL update_bar_type( barra, 'wlanczos', 1 )
  ELSE
     CALL start_bar_type ( barra, 'wlanczos', barra_load )
  ENDIF
  !
  CALL pot3D%init('Wave',.FALSE.,'default')
  !
  ! Read PDEP
  !
  ALLOCATE(pertg_all(npwqx,mypara%nloc))
  pertg_all = 0._DP
  !
  DO ip = 1,mypara%nloc
     glob_ip = mypara%l2g(ip)
     IF(glob_ip <= n_pdep_eigen_to_use) THEN
        CALL generate_pdep_fname(filepot,glob_ip)
        fname = TRIM(wstat_save_dir)//"/"//filepot
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
     IF(.NOT. l_QDET .AND. barra_load == 0) EXIT
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
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     evc_a = evc
     IF(l_QDET) THEN
        CALL apply_alpha_pa_to_m_wfcs(nbnd,evc_a,(1.0_DP,0.0_DP))
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     IF (l_frac_occ) THEN
        nbndval_full = nbnd_occ_full(iks)
     ELSE
        nbndval_full = nbnd_occ(iks)
     ENDIF
     !
     mwo = -k_grid%weight(iks_g)/omega
     zmwo = CMPLX(-k_grid%weight(iks_g)/omega,0._DP,KIND=DP)
     !
     bks%max_band = nbndval
     bks%min_band = 1
     !
     ! Parallel macropol
     IF(l_macropol) THEN
        !
        ! PHI
        !
        CALL occband%init(nbndval,'i','occband',.FALSE.)
        !
        ALLOCATE(phis(npwx*npol,3,occband%nloc))
        !
        phis = 0._DP
        !
        ALLOCATE(phi(npwx*npol,3))
        ALLOCATE(phi_tmp(npwx*npol,3))
        !
        DO ivloc = 1, occband%nloc
           !
           iv = occband%l2g(ivloc)
           !
           CALL commut_Hx_psi (iks, 1, 1, evc_a(1,iv), phi_tmp(1,1), l_skip_nl_part_of_hcomr)
           CALL commut_Hx_psi (iks, 1, 2, evc_a(1,iv), phi_tmp(1,2), l_skip_nl_part_of_hcomr)
           CALL commut_Hx_psi (iks, 1, 3, evc_a(1,iv), phi_tmp(1,3), l_skip_nl_part_of_hcomr)
           phi = 0._DP
           DO i1 = 1, 3
              DO i2 = 1, 3
                 zkonstant = CMPLX( bg(i1,i2) , 0._DP, KIND=DP )
                 CALL ZAXPY(npwx*npol,zkonstant,phi_tmp(1,i2),1,phi(1,i1),1)
              ENDDO
           ENDDO
           !
           phi_tmp = phi
           CALL apply_alpha_pc_to_m_wfcs(nbndval_full,3,phi_tmp,(1._DP,0._DP))
           !
           ALLOCATE( eprec(3) )
           ALLOCATE( e(3) )
           !
           CALL set_eprec( 1, evc_a(1,iv), eprec(1))
           eprec(2) = eprec(1)
           eprec(3) = eprec(1)
           e(1) = et(iv,iks)
           e(2) = et(iv,iks)
           e(3) = et(iv,iks)
           !
           CALL precondition_m_wfcts( 3, phi_tmp, phi, eprec )
           !
           CALL linsolve_sternheimer_m_wfcts (nbndval_full, 3, phi_tmp, phi, e, eprec, tr2_dfpt, ierr )
           !
           IF(ierr/=0) THEN
              WRITE(stdout, '(7X,"** WARNING : MACROPOL not converged, ierr = ",i8)') ierr
           ENDIF
           !
           phis(:,:,ivloc) = phi(:,:)
           !
           DEALLOCATE( eprec )
           DEALLOCATE( e )
           !
        END DO
        !
        DEALLOCATE( phi, phi_tmp )
        !
     ENDIF ! macropol
     !
     ALLOCATE(dvpsi(npwx*npol,mypara%nlocx))
     !
     time_spent(1) = get_clock( 'wlanczos' )
     !
     CALL band_group%init(nbndval,'b','band_group',.FALSE.)
     !
     ! LOOP over band states
     !
     DO ivloc2 = 1,band_group%nloc
        !
        iv = band_group%l2g(ivloc2)
        !
        IF(iks == bks%lastdone_ks .AND. iv <= bks%lastdone_band .AND. .NOT. l_QDET) CYCLE
        !
        ! MACROPOL CASE
        !
        IF(l_macropol) THEN
           !
           ALLOCATE(phi(npwx*npol,3))
           phi = 0._DP
           !
           DO ivloc = 1, occband%nloc
               IF( occband%l2g(ivloc) == iv ) THEN
                   phi(:,:) = phis(:,:,ivloc)
               ENDIF
           END DO
           !
           CALL mp_sum(phi, inter_image_comm)
           !
        ENDIF
        !
        ! PSIC
        !
        CALL single_invfft_gamma(dffts,npwq,npwqx,evc_a(1,iv),psic,TRIM(fftdriver))
        !
        ! ZEROS
        !
        dvpsi = 0._DP
        !
        ALLOCATE( pertg( npwqx ) )
        ALLOCATE( pertr( dffts%nnr ) )
        !
        DO ip=1,mypara%nloc
           !
           glob_ip = mypara%l2g(ip)
           !
           ! Decide whether read dbs E or dhpi
           !
           IF(glob_ip<=n_pdep_eigen_to_use) THEN
              !
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
                 pertr(ir)=psic(ir)*pertr(ir)
              ENDDO
              CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave')
              !
           ELSE
              !
              ipol = glob_ip-n_pdep_eigen_to_use
              !
              dvpsi(:,ip) = phi(:,ipol) * SQRT(fpi * e2)
              !
           ENDIF
           !
        ENDDO ! pert
        !
        DEALLOCATE(pertr)
        DEALLOCATE(pertg)
        IF(l_macropol) THEN
           DEALLOCATE(phi)
        ENDIF
        !
        IF(l_QDET) CALL apply_alpha_pa_to_m_wfcs(mypara%nlocx,dvpsi,(1._DP,0._DP))
        CALL apply_alpha_pc_to_m_wfcs(nbndval_full,mypara%nloc,dvpsi,(1._DP,0._DP))
        !
        ! OVERLAP( glob_ip, im=1:nbnd ) = < psi_im iks | dvpsi_glob_ip >
        !
        IF(ALLOCATED(ps_r)) DEALLOCATE(ps_r)
        ALLOCATE(ps_r(nbnd-nbndval_full,mypara%nloc))
        CALL glbrak_gamma(evc_a(1,nbndval_full+1),dvpsi(1,1),ps_r,npw,npwx,nbnd-nbndval_full,mypara%nloc,nbnd-nbndval_full,npol)
        CALL mp_sum(ps_r,intra_bgrp_comm)
        !
        IF(ALLOCATED(overlap)) DEALLOCATE(overlap)
        ALLOCATE(overlap(mypara%nglob, nbndval_full+1:nbnd ) )
        overlap = 0._DP
        DO ic = nbndval_full+1, nbnd
           DO ip = 1, mypara%nloc
              glob_ip = mypara%l2g(ip)
              overlap(glob_ip,ic) = ps_r(ic-nbndval_full,ip)
           ENDDO
        ENDDO
        !
        DEALLOCATE(ps_r)
        CALL mp_sum(overlap,inter_image_comm)
        !
        ! Update dmati with cond
        !
        DO ifreq = 1, ifr%nloc
           !
           frequency = imfreq_list( ifreq )
           !
           DO ic = nbndval_full+1, nbnd
              !
              ecv = et(ic,iks)-et(iv,iks)
              dfactor = mwo * 2._DP * ecv / ( ecv**2 + frequency**2 )
              !
              IF (ic <= iv) CYCLE  ! avoid double counting in frac occ case
              !
              IF (l_frac_occ) THEN
                 !
                 docc = occupation(iv,iks) - occupation(ic,iks)
                 dfactor = dfactor * docc
                 !
              ENDIF
              !
              DO ip = 1, mypara%nloc
                 glob_ip = mypara%l2g(ip)
                 DO glob_jp = 1, mypara%nglob
                    !
                    IF(l_QDET) THEN
                       dmati_a(glob_jp,ip,ifreq) = dmati_a(glob_jp,ip,ifreq) &
                       & + overlap( glob_jp, ic) * overlap( glob_ip, ic) * dfactor
                    ELSE
                       dmati(glob_jp,ip,ifreq) = dmati(glob_jp,ip,ifreq) &
                       & + overlap( glob_jp, ic) * overlap( glob_ip, ic) * dfactor
                    ENDIF
                    !
                 ENDDO
              ENDDO
              !
           ENDDO ! ic
        ENDDO ! ifreq
        !
        ! Update zmatr with cond
        !
        DO ifreq = 1, rfr%nloc
           !
           frequency = refreq_list( ifreq )
           !
           DO ic = nbndval_full+1, nbnd
              !
              ecv = et(ic,iks)-et(iv,iks)
              zp = CMPLX( ecv + frequency, - wfreq_eta, KIND=DP )
              zm = CMPLX( ecv - frequency, - wfreq_eta, KIND=DP )
              zfactor = zmwo / zp + zmwo / zm
              !
              IF (ic <= iv) CYCLE  ! avoid double counting in frac occ case
              !
              IF (l_frac_occ) THEN
                 docc = occupation(iv,iks) - occupation(ic,iks)
                 zfactor = zfactor * docc
                 !
              ENDIF
              !
              DO ip = 1, mypara%nloc
                 glob_ip = mypara%l2g(ip)
                 DO glob_jp = 1, mypara%nglob
                    !
                    IF(l_QDET) THEN
                       zmatr_a(glob_jp,ip,ifreq) = zmatr_a(glob_jp,ip,ifreq) &
                       + CMPLX( overlap( glob_jp, ic) * overlap( glob_ip, ic), 0._DP, KIND=DP ) * zfactor
                    ELSE
                       zmatr(glob_jp,ip,ifreq) = zmatr(glob_jp,ip,ifreq) &
                       + CMPLX( overlap( glob_jp, ic) * overlap( glob_ip, ic), 0._DP, KIND=DP ) * zfactor
                    ENDIF
                    !
                 ENDDO
              ENDDO
              !
           ENDDO ! ic
        ENDDO ! ifreq
        !
        DEALLOCATE(overlap)
        !
        IF(l_QDET) CYCLE
        !   
        ! Apply Pc, to be sure
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,mypara%nloc,dvpsi,(1._DP,0._DP))
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF( l_enable_lanczos ) THEN
           !
           ALLOCATE( bnorm    (                         mypara%nloc ) )
           ALLOCATE( diago    (            n_lanczos  , mypara%nloc ) )
           ALLOCATE( subdiago (            n_lanczos-1, mypara%nloc ) )
           ALLOCATE( q_s      ( npwx*npol, mypara%nloc, n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM
           !
           CALL solve_deflated_lanczos_w_full_ortho(nbnd, mypara%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
           !
           ALLOCATE( braket( mypara%nglob, n_lanczos, mypara%nloc ) )
           CALL get_brak_hyper_parallel(dvpsi,mypara%nloc,n_lanczos,q_s,braket,mypara)
           !
           DO ip = 1, mypara%nloc
              CALL diago_lanczos( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), mypara%nglob )
           ENDDO
           !
           DEALLOCATE( q_s )
           DEALLOCATE( bnorm )
           DEALLOCATE( subdiago )
           !
           ! Update dmati with lanczos
           !
           DO ifreq = 1, ifr%nloc
              !
              frequency = imfreq_list( ifreq )
              !
              DO il = 1, n_lanczos
                 !
                 DO ip = 1, mypara%nloc
                    glob_ip = mypara%l2g(ip)
                    ecv = diago( il, ip ) - et(iv,iks)
                    dfactor = mwo * 2._DP * ecv / ( ecv**2 + frequency**2 )
                    IF (l_frac_occ) dfactor = dfactor * occupation(iv,iks)
                    !
                    DO glob_jp = 1, mypara%nglob
                       !
                       dmati(glob_jp,ip,ifreq) = dmati(glob_jp,ip,ifreq) + braket( glob_jp, il, ip) * dfactor
                       !
                    ENDDO
                 ENDDO
                 !
              ENDDO ! il
           ENDDO ! ifreq
           !
           ! Update zmatr with lanczos
           !
           DO ifreq = 1, rfr%nloc
              !
              frequency = refreq_list( ifreq )
              !
              DO il = 1, n_lanczos
                 !
                 DO ip = 1, mypara%nloc
                    glob_ip = mypara%l2g(ip)
                    ecv = diago( il, ip ) - et(iv,iks)
                    zp = CMPLX( ecv + frequency, - wfreq_eta, KIND=DP )
                    zm = CMPLX( ecv - frequency, - wfreq_eta, KIND=DP )
                    zfactor = zmwo / zp + zmwo / zm
                    IF (l_frac_occ) zfactor = zfactor * occupation(iv,iks)
                    !
                    DO glob_jp = 1, mypara%nglob
                       !
                       zmatr(glob_jp,ip,ifreq) = zmatr(glob_jp,ip,ifreq) + &
                       & CMPLX( braket( glob_jp, il, ip), 0._DP, KIND=DP ) * zfactor
                       !
                    ENDDO
                 ENDDO
                 !
              ENDDO ! il
           ENDDO ! ifreq
           !
           DEALLOCATE( diago )
           DEALLOCATE( braket )
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'wlanczos' )
        l_write_restart = .FALSE.
        !
        IF( o_restart_time >= 0._DP ) THEN
           IF( time_spent(2)-time_spent(1) > o_restart_time*60._DP ) l_write_restart = .TRUE.
           IF( iv == nbndval ) l_write_restart = .TRUE.
        ENDIF
        !
        ! Write final restart file
        !
        IF( iks == k_grid%nps .AND. iv == nbndval .AND. .NOT. l_QDET ) l_write_restart = .TRUE.
        !
        ! But do not write here when using pool or band group
        !
        IF( npool*nbgrp > 1 ) l_write_restart = .FALSE.
        !
        IF( l_write_restart ) THEN
           bks%lastdone_ks = iks
           bks%lastdone_band = iv
           CALL solvewfreq_restart_write(bks,dmati,zmatr,mypara%nglob,mypara%nloc)
           bks%old_ks = iks
           bks%old_band = iv
           time_spent(1) = get_clock( 'wlanczos' )
        ENDIF
        !
        CALL update_bar_type( barra, 'wlanczos', 1 )
        !
     ENDDO ! BANDS
     !
     IF(l_macropol) DEALLOCATE(phis)
     !
     DEALLOCATE(dvpsi)
     !
  ENDDO ! KPOINT-SPIN
  !
  DEALLOCATE(pertg_all)
  DEALLOCATE(evc_a)
  !
  ! Synchronize and write final restart file when using pool or band group
  !
  IF(npool*nbgrp > 1 .AND. .NOT. l_read_restart .AND. .NOT. l_QDET) THEN
     bks%lastdone_ks = k_grid%nps
     bks%lastdone_band = nbndval
     CALL mp_sum(dmati,inter_bgrp_comm)
     CALL mp_sum(dmati,inter_pool_comm)
     CALL mp_sum(zmatr,inter_bgrp_comm)
     CALL mp_sum(zmatr,inter_pool_comm)
     CALL solvewfreq_restart_write(bks,dmati,zmatr,mypara%nglob,mypara%nloc)
  ENDIF
  !
  IF(l_QDET) THEN
     CALL mp_sum(dmati_a,inter_bgrp_comm)
     CALL mp_sum(dmati_a,inter_pool_comm)
     CALL mp_sum(zmatr_a,inter_bgrp_comm)
     CALL mp_sum(zmatr_a,inter_pool_comm)
  ENDIF
  !
  CALL stop_bar_type( barra, 'wlanczos' )
  !
  CALL start_clock('chi_invert')
  !
  ! EPS-1 imfreq
  !
  IF (ALLOCATED(d_epsm1_ifr)) DEALLOCATE(d_epsm1_ifr)
  IF (ALLOCATED(d_head_ifr)) DEALLOCATE(d_head_ifr)
  ALLOCATE(dmatilda(mypara%nglob,mypara%nglob))
  ALLOCATE(dlambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
  ALLOCATE(d_epsm1_ifr(pert%nglob,pert%nloc,ifr%nloc))
  d_epsm1_ifr = 0._DP
  IF (l_QDET) THEN
     ALLOCATE(dmatilda_a(mypara%nglob,mypara%nglob))
     ALLOCATE(d_epsm1_ifr_a(pert%nglob,pert%nloc,ifr%nloc))
     d_epsm1_ifr_a = 0._DP
  ENDIF
  IF(l_macropol) THEN
     ALLOCATE(d_head_ifr(ifr%nloc))
     d_head_ifr = 0._DP
     IF(l_QDET) THEN
        ALLOCATE(d_head_ifr_a(ifr%nloc))
        d_head_ifr_a = 0._DP  
     ENDIF     
  ENDIF
  !
  CALL band_group%init(ifr%nloc,'b','band_group',.FALSE.)
  !
  DO ifloc = 1,band_group%nloc
     !
     ifreq = band_group%l2g(ifloc)
     !
     dmatilda = 0._DP
     IF (l_QDET) dmatilda_a = 0._DP
     DO ip = 1,mypara%nloc
        glob_ip = mypara%l2g(ip)
        dmatilda(:,glob_ip) = dmati(:,ip,ifreq)
        IF (l_QDET) dmatilda_a(:,glob_ip) = dmati_a(:,ip,ifreq)
     ENDDO
     !
     CALL mp_sum(dmatilda,inter_image_comm)
     IF(l_QDET) THEN
        CALL mp_sum(dmatilda_a,inter_image_comm)
        CALL chi_invert_real(dmatilda-dmatilda_a,dhead,dlambda,mypara%nglob)
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           d_epsm1_ifr_a(1:n_pdep_eigen_to_use,ip,ifreq) = dlambda(1:n_pdep_eigen_to_use,glob_ip)
        ENDDO
        IF(l_macropol) d_head_ifr_a(ifreq) = dhead
     ENDIF 
     !
     CALL chi_invert_real(dmatilda,dhead,dlambda,mypara%nglob)
     !
     DO ip = 1,pert%nloc
        glob_ip = pert%l2g(ip)
        d_epsm1_ifr(1:n_pdep_eigen_to_use,ip,ifreq) = dlambda(1:n_pdep_eigen_to_use,glob_ip)
     ENDDO
     IF(l_macropol) d_head_ifr(ifreq) = dhead
     !
  ENDDO
  !
  DEALLOCATE(dlambda)
  DEALLOCATE(dmatilda)
  DEALLOCATE(dmati)
  IF(l_QDET) THEN
     DEALLOCATE(dmatilda_a)
     DEALLOCATE(dmati_a)
  ENDIF     
  !
  CALL mp_sum(d_epsm1_ifr,inter_bgrp_comm)
  IF(l_macropol) CALL mp_sum(d_head_ifr,inter_bgrp_comm)
  IF(l_QDET) THEN
     CALL mp_sum(d_epsm1_ifr_a,inter_bgrp_comm)
     IF(l_macropol) CALL mp_sum(d_head_ifr_a,inter_bgrp_comm)
  ENDIF
  !
  ! EPS-1 refreq
  !
  IF (ALLOCATED(z_epsm1_rfr)) DEALLOCATE(z_epsm1_rfr)
  IF (ALLOCATED(z_head_rfr)) DEALLOCATE(z_head_rfr)
  ALLOCATE(zmatilda(mypara%nglob,mypara%nglob))
  ALLOCATE(zlambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
  ALLOCATE(zlambda_a(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
  ALLOCATE(z_epsm1_rfr(pert%nglob,pert%nloc,rfr%nloc))
  z_epsm1_rfr = 0._DP
  IF (l_QDET) THEN
     ALLOCATE(zmatilda_a(mypara%nglob,mypara%nglob))
     ALLOCATE(zmatilda_diff(mypara%nglob,mypara%nglob))
     ALLOCATE(z_epsm1_rfr_a(pert%nglob,pert%nloc,rfr%nloc))
     z_epsm1_rfr_a = 0._DP
  ENDIF
  IF(l_macropol) THEN
     ALLOCATE(z_head_rfr(rfr%nloc))
     z_head_rfr = 0._DP
     IF(l_QDET) THEN
        ALLOCATE(z_head_rfr_a(rfr%nloc))
        z_head_rfr_a = 0._DP  
     ENDIF 
  ENDIF
  !
  CALL band_group%init(rfr%nloc,'b','band_group',.FALSE.)
  !
  DO ifloc = 1,band_group%nloc
     !
     ifreq = band_group%l2g(ifloc)
     !
     zmatilda = 0._DP
     IF (l_QDET) zmatilda_a = 0._DP
     DO ip = 1,mypara%nloc
        glob_ip = mypara%l2g(ip)
        zmatilda(:,glob_ip) = zmatr(:,ip,ifreq)
        IF (l_QDET) zmatilda_a(:,glob_ip) = zmatr_a(:,ip,ifreq)
     ENDDO
     !
     CALL mp_sum(zmatilda,inter_image_comm)
     CALL chi_invert_complex(zmatilda,zhead,zlambda,mypara%nglob)
     !
     DO ip = 1,pert%nloc
        glob_ip = pert%l2g(ip)
        z_epsm1_rfr(1:n_pdep_eigen_to_use,ip,ifreq) = zlambda(1:n_pdep_eigen_to_use,glob_ip)
     ENDDO
     IF(l_macropol) z_head_rfr(ifreq) = zhead
     !
     IF(l_QDET) THEN
        CALL mp_sum(zmatilda_a,inter_image_comm)
        zmatilda_diff = zmatilda - zmatilda_a
        CALL chi_invert_complex(zmatilda_diff,zhead_a,zlambda_a,mypara%nglob)
        DO ip = 1,pert%nloc
           glob_ip = pert%l2g(ip)
           z_epsm1_rfr_a(1:n_pdep_eigen_to_use,ip,ifreq) = zlambda_a(1:n_pdep_eigen_to_use,glob_ip)
        ENDDO
        IF(l_macropol) z_head_rfr_a(ifreq) = zhead_a
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE(zlambda)
  DEALLOCATE(zmatilda)
  DEALLOCATE(zmatr)
  IF(l_QDET) THEN
     DEALLOCATE(zlambda_a)
     DEALLOCATE(zmatilda_a, zmatilda_diff)
     DEALLOCATE(zmatr_a)
  ENDIF  
  !
  CALL mp_sum(z_epsm1_rfr,inter_bgrp_comm)
  IF(l_macropol) CALL mp_sum(z_head_rfr,inter_bgrp_comm)
  IF(l_QDET) THEN
     CALL mp_sum(z_epsm1_rfr_a,inter_bgrp_comm)
     IF(l_macropol) CALL mp_sum(z_head_rfr_a,inter_bgrp_comm)
  ENDIF
  !
  CALL stop_clock('chi_invert')
  !
  IF(l_generate_plot) THEN
     CALL output_eps_head()
  ENDIF
  !
  CALL mp_barrier(world_comm)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_wfreq_k(l_read_restart,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_lanczos,npwq,l_macropol,z_epsm1_ifr_q,&
                                 & z_epsm1_rfr_q,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,wfreq_eta,&
                                 & imfreq_list,refreq_list,tr2_dfpt,z_head_rfr,z_head_ifr,&
                                 & o_restart_time,l_skip_nl_part_of_hcomr,npwqx,wstat_save_dir,&
                                 & ngq,igq_q
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,inter_bgrp_comm,nbgrp
  USE mp,                   ONLY : mp_bcast,mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : bg,omega
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,current_spin,isk,nbnd,lsda,igk_k,current_k,ngk
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
  USE distribution_center,  ONLY : pert,macropert,ifr,rfr,occband,band_group
  USE class_idistribute,    ONLY : idistribute
  USE wfreq_restart,        ONLY : solvewfreq_restart_write,solvewfreq_restart_read,bksq_type
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE chi_invert,           ONLY : chi_invert_complex
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart,l_generate_plot
  !
  ! Workspace
  !
  LOGICAL :: l_write_restart
  INTEGER :: i1,i2,ip,ig,glob_ip,ir,iv,ivloc,ivloc2,ifloc,iks,ik,is,iq,ikqs,ikq,ipol
  CHARACTER(LEN=25) :: filepot
  CHARACTER(LEN=:),ALLOCATABLE    :: fname
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:)
  COMPLEX(DP),ALLOCATABLE :: braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi(:,:), phis(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: phi_tmp(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP),ALLOCATABLE :: pertg_all(:,:)
  COMPLEX(DP),ALLOCATABLE :: evckpq(:,:)
  COMPLEX(DP),ALLOCATABLE :: psick(:),psick_nc(:,:)
  COMPLEX(DP),ALLOCATABLE :: phase(:)
  INTEGER :: npwkq
  COMPLEX(DP) :: zkonstant
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP),ALLOCATABLE :: ps_c(:,:)
  TYPE(idistribute) :: mypara
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
  REAL(DP) :: mwo, ecv, dfactor, frequency
  COMPLEX(DP) :: zmwo, zfactor,zm,zp, zhead
  INTEGER :: glob_jp,ic,ifreq,il
  COMPLEX(DP),ALLOCATABLE :: zmatilda(:,:), zlambda(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmati_q(:,:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr_q(:,:,:,:)
  LOGICAL :: l_gammaq
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bksq_type) :: bksq
  REAL(DP),ALLOCATABLE :: eprec(:)
  INTEGER :: ierr
  REAL(DP),ALLOCATABLE :: e(:)
  REAL(DP) :: g0(3)
  !
  CALL io_push_title("(W)-Lanczos")
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
  CALL deallocate_bec_type( becp )
  IF(l_macropol) THEN
     CALL allocate_bec_type ( nkb, MAX(mypara%nloc,3), becp ) ! I just need 2 becp at a time
  ELSE
     CALL allocate_bec_type ( nkb, mypara%nloc, becp ) ! I just need 2 becp at a time
  ENDIF
  !
  ! ALLOCATE zmati_q, zmatr_q, where chi0 is stored
  !
  ALLOCATE( zmati_q( mypara%nglob, mypara%nloc, ifr%nloc, q_grid%np) )
  ALLOCATE( zmatr_q( mypara%nglob, mypara%nloc, rfr%nloc, q_grid%np) )
  zmati_q = 0._DP
  zmatr_q = 0._DP
  !
  ALLOCATE( evckpq(npwx*npol,nbnd) )
  IF (noncolin) THEN
     ALLOCATE( psick_nc(dffts%nnr,npol) )
  ELSE
     ALLOCATE( psick(dffts%nnr) )
  ENDIF
  ALLOCATE( phase(dffts%nnr) )
  !
  IF(l_read_restart) THEN
     CALL solvewfreq_restart_read( bksq, zmati_q, zmatr_q, mypara%nglob, mypara%nloc )
  ELSE
     bksq%lastdone_ks   = 0
     bksq%lastdone_q    = 0
     bksq%lastdone_band = 0
     bksq%old_ks        = 0
     bksq%old_q         = 0
     bksq%old_band      = 0
     bksq%max_q         = q_grid%np
     bksq%max_ks        = k_grid%nps
     bksq%min_q         = 1
     bksq%min_ks        = 1
  ENDIF
  !
  barra_load = 0
  DO iq = 1,q_grid%np
     IF(iq < bksq%lastdone_q) CYCLE
     !
     DO iks = 1,k_grid%nps
        IF(iq == bksq%lastdone_q .AND. iks < bksq%lastdone_ks) CYCLE
        !
        CALL band_group%init(nbnd_occ(iks),'b','band_group',.FALSE.)
        !
        DO ivloc = 1,band_group%nloc
           iv = band_group%l2g(ivloc)
           !
           IF(iq == bksq%lastdone_q .AND. iks == bksq%lastdone_ks .AND. iv <= bksq%lastdone_band) CYCLE
           !
           barra_load = barra_load+1
        ENDDO
     ENDDO
  ENDDO
  !
  IF( barra_load == 0 ) THEN
     CALL start_bar_type ( barra, 'wlanczos', 1 )
     CALL update_bar_type( barra, 'wlanczos', 1 )
  ELSE
     CALL start_bar_type ( barra, 'wlanczos', barra_load )
  ENDIF
  !
  ! LOOP
  !
  DO iq = 1, q_grid%np   ! Q-POINT
     !
     ! Exit loop if no work to do
     !
     IF(barra_load == 0) EXIT
     !
     IF(iq < bksq%lastdone_q) CYCLE
     !
     npwq = ngq(iq)
     l_gammaq = q_grid%l_pIsGamma(iq)
     !
     CALL pot3D%init('Wave',.TRUE.,'default',iq)
     !
     ! Read PDEP
     !
     ALLOCATE(pertg_all(npwqx,mypara%nloc))
     pertg_all = 0._DP
     !
     DO ip = 1,mypara%nloc
        glob_ip = mypara%l2g(ip)
        IF(glob_ip <= n_pdep_eigen_to_use) THEN
           CALL generate_pdep_fname(filepot,glob_ip,iq)
           fname = TRIM(wstat_save_dir)//"/"//filepot
           CALL pdep_read_G_and_distribute(fname,pertg_all(:,ip),iq)
        ENDIF
     ENDDO
     !
     DO iks = 1, k_grid%nps   ! KPOINT-SPIN
        !
        ik = k_grid%ip(iks)
        is = k_grid%is(iks)
        !
        IF(iq==bksq%lastdone_q .AND. iks<bksq%lastdone_ks) CYCLE
        !
        ! ... Set k-point, spin, kinetic energy, needed by Hpsi
        !
        current_k = iks
        IF ( lsda ) current_spin = isk(iks)
        call g2_kin( iks )
        !
        ! ... More stuff needed by the hamiltonian: nonlocal projectors
        !
        IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), k_grid%p_cart(1,ik), vkb )
        npw = ngk(iks)
        !
        ! ... read in wavefunctions from the previous iteration
        !
        IF(k_grid%nps>1) THEN
           IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
           CALL mp_bcast(evc,0,inter_image_comm)
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
        IF ( my_image_id == 0 ) CALL get_buffer( evckpq, lrwfc, iuwfc, ikqs )
        CALL mp_bcast( evckpq, 0, inter_image_comm )
        !
        nbndval = nbnd_occ(iks)
        !
        mwo = - k_grid%weight(iks) / omega
        zmwo = CMPLX( - k_grid%weight(iks) / omega, 0._DP, KIND=DP)
        !
        bksq%max_band = nbndval
        bksq%min_band = 1
        !
        ! MACROPOL CASE
        !
        IF(l_macropol .AND. l_gammaq) THEN
           !
           ! PHI
           !
           CALL occband%init(nbndval,'i','occband',.FALSE.)
           !
           ALLOCATE(phis(npwx*npol,3,occband%nloc))
           !
           phis = 0._DP
           !
           ALLOCATE(phi(npwx*npol,3))
           ALLOCATE(phi_tmp(npwx*npol,3))
           !
           DO ivloc = 1, occband%nloc
              !
              iv = occband%l2g(ivloc)
              !
              CALL commut_Hx_psi (iks, 1, 1, evc(1,iv), phi_tmp(1,1), l_skip_nl_part_of_hcomr)
              CALL commut_Hx_psi (iks, 1, 2, evc(1,iv), phi_tmp(1,2), l_skip_nl_part_of_hcomr)
              CALL commut_Hx_psi (iks, 1, 3, evc(1,iv), phi_tmp(1,3), l_skip_nl_part_of_hcomr)
              phi = 0._DP
              DO i1 = 1, 3
                 DO i2 = 1, 3
                    zkonstant = CMPLX( bg(i1,i2) , 0._DP, KIND=DP )
                    CALL ZAXPY(npwx*npol,zkonstant,phi_tmp(1,i2),1,phi(1,i1),1)
                 ENDDO
              ENDDO
              !
              phi_tmp = phi
              CALL apply_alpha_pc_to_m_wfcs(nbndval,3,phi_tmp,(1._DP,0._DP))
              !
              ALLOCATE( eprec(3) )
              ALLOCATE( e(3) )
              !
              CALL set_eprec( 1, evc(1,iv), eprec(1))
              eprec(2) = eprec(1)
              eprec(3) = eprec(1)
              e(1) = et(iv,iks)
              e(2) = et(iv,iks)
              e(3) = et(iv,iks)
              !
              CALL precondition_m_wfcts( 3, phi_tmp, phi, eprec )
              !
              CALL linsolve_sternheimer_m_wfcts (nbndval, 3, phi_tmp, phi, e, eprec, tr2_dfpt, ierr )
              !
              IF(ierr/=0) THEN
                 WRITE(stdout, '(7X,"** WARNING : MACROPOL not converged, ierr = ",i8)') ierr
              ENDIF
              !
              phis(:,:,ivloc) = phi(:,:)
              !
              DEALLOCATE( eprec )
              DEALLOCATE( e )
              !
           END DO
           !
           DEALLOCATE( phi, phi_tmp )
           !
        ENDIF ! macropol
        !
        ALLOCATE(dvpsi(npwx*npol,mypara%nlocx))
        !
        time_spent(1) = get_clock( 'wlanczos' )
        !
        CALL band_group%init(nbndval,'b','band_group',.FALSE.)
        !
        ! LOOP over band states
        !
        DO ivloc2 = 1,band_group%nloc
           !
           iv = band_group%l2g(ivloc2)
           !
           IF(iq == bksq%lastdone_q .AND. iks == bksq%lastdone_ks .AND. iv <= bksq%lastdone_band) CYCLE
           !
           ! MACROPOL CASE
           !
           IF(l_macropol .AND. l_gammaq) THEN
              !
              ALLOCATE(phi(npwx*npol,3))
              phi = 0._DP
              !
              DO ivloc = 1, occband%nloc
                 IF( occband%l2g(ivloc) == iv ) THEN
                    phi(:,:) = phis(:,:,ivloc)
                 ENDIF
              END DO
              !
              CALL mp_sum(phi, inter_image_comm)
              !
           ENDIF
           !
           ! PSIC
           !
           IF(noncolin) THEN
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(1     ,iv),psick_nc(1,1),'Wave',igk_k(1,ikqs))
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(npwx+1,iv),psick_nc(1,2),'Wave',igk_k(1,ikqs))
           ELSE
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(1,iv),psick,'Wave',igk_k(1,ikqs))
           ENDIF
           !
           ! ZEROS
           !
           dvpsi = 0._DP
           !
           ALLOCATE( pertg( npwqx ) )
           ALLOCATE( pertr( dffts%nnr ) )
           !
           DO ip=1,mypara%nloc
              !
              glob_ip = mypara%l2g(ip)
              !
              ! Decide whether read dbs E or dhpi
              !
              IF(glob_ip<=n_pdep_eigen_to_use) THEN
                 !
                 pertg = pertg_all(:,ip)
                 !
                 ! Multiply by sqvc
                 DO ig = 1, npwq
                    pertg(ig) = pot3D%sqvc(ig) * pertg(ig)
                 ENDDO
                 !
                 ! Bring it to R-space
                 IF(noncolin) THEN
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                    DO ir=1,dffts%nnr
                       pertr(ir)=phase(ir)*psick_nc(ir,1)*CONJG(pertr(ir))
                    ENDDO
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                    DO ir=1,dffts%nnr
                       pertr(ir)=phase(ir)*psick_nc(ir,2)*CONJG(pertr(ir))
                    ENDDO
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(npwx+1,ip),'Wave',igk_k(1,current_k))
                 ELSE
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                    DO ir=1,dffts%nnr
                       pertr(ir)=phase(ir)*psick(ir)*CONJG(pertr(ir))
                    ENDDO
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
                 ENDIF
                 !
              ELSE
                 !
                 IF (l_gammaq) THEN
                    ipol = glob_ip-n_pdep_eigen_to_use
                    dvpsi(:,ip) = phi(:,ipol) * SQRT(fpi * e2)
                 ENDIF
                 !
              ENDIF
              !
           ENDDO ! pert
           !
           DEALLOCATE(pertr)
           DEALLOCATE(pertg)
           IF(l_macropol.AND.l_gammaq) THEN
              DEALLOCATE(phi)
           ENDIF
           !
           CALL apply_alpha_pc_to_m_wfcs(nbndval,mypara%nloc,dvpsi,(1._DP,0._DP))
           !
           ! OVERLAP( glob_ip, im=1:nbnd ) = < psi_im iks | dvpsi_glob_ip >
           !
           IF(ALLOCATED(ps_c)) DEALLOCATE(ps_c)
           ALLOCATE(ps_c(nbnd-nbndval,mypara%nloc))
           CALL glbrak_k(evc(1,nbndval+1),dvpsi(1,1),ps_c,npw,npwx,nbnd-nbndval,mypara%nloc,nbnd-nbndval,npol)
           CALL mp_sum(ps_c,intra_bgrp_comm)
           !
           IF(ALLOCATED(overlap)) DEALLOCATE(overlap)
           ALLOCATE(overlap(mypara%nglob, nbndval+1:nbnd ) )
           overlap = 0._DP
           DO ic = nbndval+1, nbnd
              DO ip = 1, mypara%nloc
                 glob_ip = mypara%l2g(ip)
                 overlap(glob_ip,ic) = ps_c(ic-nbndval,ip)
              ENDDO
           ENDDO
           !
           DEALLOCATE(ps_c)
           CALL mp_sum(overlap,inter_image_comm)
           !
           ! Update zmati with cond
           !
           DO ifreq = 1, ifr%nloc
              !
              frequency = imfreq_list( ifreq )
              !
              DO ic = nbndval+1, nbnd
                 !
                 ecv = et(ic,iks)-et(iv,ikqs)
                 dfactor = mwo * 2._DP * ecv / ( ecv**2 + frequency**2 )
                 zfactor = CMPLX( dfactor, 0._DP, KIND=DP )
                 !
                 DO ip = 1, mypara%nloc
                    glob_ip = mypara%l2g(ip)
                    DO glob_jp = 1, mypara%nglob
                       !
                       zmati_q(glob_jp,ip,ifreq,iq) = zmati_q(glob_jp,ip,ifreq,iq) + CONJG(overlap( glob_ip, ic)) * &
                       & overlap( glob_jp, ic) * zfactor
                       !
                    ENDDO
                 ENDDO
                 !
              ENDDO ! ic
           ENDDO ! ifreq
           !
           ! Update zmatr with cond
           !
           DO ifreq = 1, rfr%nloc
              !
              frequency = refreq_list( ifreq )
              !
              DO ic = nbndval+1, nbnd
                 !
                 ecv = et(ic,iks)-et(iv,ikqs)
                 zp = CMPLX( ecv + frequency, - wfreq_eta, KIND=DP )
                 zm = CMPLX( ecv - frequency, - wfreq_eta, KIND=DP )
                 zfactor = zmwo / zp + zmwo / zm
                 !
                 DO ip = 1, mypara%nloc
                    glob_ip = mypara%l2g(ip)
                    DO glob_jp = 1, mypara%nglob
                       !
                       zmatr_q(glob_jp,ip,ifreq,iq) = zmatr_q(glob_jp,ip,ifreq,iq) + CONJG( overlap( glob_ip, ic) ) * &
                       & overlap( glob_jp, ic) * zfactor
                       !
                    ENDDO
                 ENDDO
                 !
              ENDDO ! ic
           ENDDO ! ifreq
           !
           DEALLOCATE(overlap)
           !
           ! Apply Pc, to be sure
           !
           CALL apply_alpha_pc_to_m_wfcs(nbnd,mypara%nloc,dvpsi,(1._DP,0._DP))
           !
           ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
           ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
           !
           IF( l_enable_lanczos ) THEN
              !
              ALLOCATE( bnorm    (                         mypara%nloc ) )
              ALLOCATE( diago    (            n_lanczos  , mypara%nloc ) )
              ALLOCATE( subdiago (            n_lanczos-1, mypara%nloc ) )
              ALLOCATE( q_s      ( npwx*npol, mypara%nloc, n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM
              !
              CALL solve_deflated_lanczos_w_full_ortho(nbnd, mypara%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
              !
              ALLOCATE( braket( mypara%nglob, n_lanczos  , mypara%nloc ) )
              CALL get_brak_hyper_parallel_complex(dvpsi,mypara%nloc,n_lanczos,q_s,braket,mypara)
              !
              DO ip = 1, mypara%nloc
                 CALL diago_lanczos_complex( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), mypara%nglob )
              ENDDO
              !
              DEALLOCATE( q_s )
              DEALLOCATE( bnorm )
              DEALLOCATE( subdiago )
              !
              ! Update zmati with lanczos
              !
              DO ifreq = 1, ifr%nloc
                 !
                 frequency = imfreq_list( ifreq )
                 !
                 DO il = 1, n_lanczos
                    !
                    DO ip = 1, mypara%nloc
                       glob_ip = mypara%l2g(ip)
                       ecv = diago( il, ip ) - et(iv,ikqs)
                       dfactor = mwo * 2._DP * ecv / ( ecv**2 + frequency**2 )
                       zfactor = CMPLX( dfactor, 0._DP, KIND=DP )
                       DO glob_jp = 1, mypara%nglob
                          !
                          zmati_q(glob_jp,ip,ifreq,iq) = zmati_q(glob_jp,ip,ifreq,iq) + CONJG(braket( glob_jp, il, ip)) * zfactor
                          !
                       ENDDO
                    ENDDO
                    !
                 ENDDO ! il
              ENDDO ! ifreq
              !
              ! Update zmatr with lanczos
              !
              DO ifreq = 1, rfr%nloc
                 !
                 frequency = refreq_list( ifreq )
                 !
                 DO il = 1, n_lanczos
                    !
                    DO ip = 1, mypara%nloc
                       glob_ip = mypara%l2g(ip)
                       ecv = diago( il, ip ) - et(iv,ikqs)
                       zp = CMPLX( ecv + frequency, - wfreq_eta, KIND=DP )
                       zm = CMPLX( ecv - frequency, - wfreq_eta, KIND=DP )
                       zfactor = zmwo / zp + zmwo / zm
                       DO glob_jp = 1, mypara%nglob
                          !
                          zmatr_q(glob_jp,ip,ifreq,iq) = zmatr_q(glob_jp,ip,ifreq,iq) + CONJG(braket( glob_jp, il, ip)) * zfactor
                          !
                       ENDDO
                    ENDDO
                    !
                 ENDDO ! il
              ENDDO ! ifreq
              !
              DEALLOCATE( diago )
              DEALLOCATE( braket )
              !
           ENDIF ! l_enable_lanczos
           !
           time_spent(2) = get_clock( 'wlanczos' )
           l_write_restart = .FALSE.
           !
           IF( o_restart_time >= 0._DP ) THEN
              IF( time_spent(2)-time_spent(1) > o_restart_time*60._DP ) l_write_restart = .TRUE.
              IF( iv == nbndval ) l_write_restart = .TRUE.
           ENDIF
           !
           ! Write final restart file
           !
           IF( iq == q_grid%np .AND. iks == k_grid%nps .AND. iv == nbndval ) l_write_restart = .TRUE.
           !
           ! But do not write here when using band group
           !
           IF( nbgrp > 1 ) l_write_restart = .FALSE.
           !
           IF( l_write_restart ) THEN
              bksq%lastdone_q = iq
              bksq%lastdone_ks = iks
              bksq%lastdone_band = iv
              CALL solvewfreq_restart_write(bksq,zmati_q,zmatr_q,mypara%nglob,mypara%nloc)
              bksq%old_q = iq
              bksq%old_ks = iks
              bksq%old_band = iv
              time_spent(1) = get_clock( 'wlanczos' )
           ENDIF
           !
           CALL update_bar_type( barra, 'wlanczos', 1 )
           !
        ENDDO ! BANDS
        !
        IF(l_macropol .AND. l_gammaq) DEALLOCATE(phis)
        !
        DEALLOCATE(dvpsi)
        !
     ENDDO ! KPOINT-SPIN
     !
     DEALLOCATE(pertg_all)
     !
  ENDDO ! QPOINT
  !
  DEALLOCATE( evckpq )
  IF (noncolin) THEN
     DEALLOCATE( psick_nc )
  ELSE
     DEALLOCATE( psick )
  ENDIF
  DEALLOCATE( phase )
  !
  ! Synchronize and write final restart file when using band group
  !
  IF(nbgrp > 1 .AND. .NOT. l_read_restart) THEN
     bksq%lastdone_q = q_grid%np
     bksq%lastdone_ks = k_grid%nps
     bksq%lastdone_band = nbndval
     CALL mp_sum(zmati_q,inter_bgrp_comm)
     CALL mp_sum(zmatr_q,inter_bgrp_comm)
     CALL solvewfreq_restart_write(bksq,zmati_q,zmatr_q,mypara%nglob,mypara%nloc)
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
  ALLOCATE(z_epsm1_ifr_q(pert%nglob,pert%nloc,ifr%nloc,q_grid%np))
  z_epsm1_ifr_q = 0._DP
  IF(l_macropol) THEN
     ALLOCATE(z_head_ifr(ifr%nloc))
     z_head_ifr = 0._DP
  ENDIF
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
        !
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
  DEALLOCATE(zlambda)
  DEALLOCATE(zmatilda)
  DEALLOCATE(zmati_q)
  !
  CALL mp_sum(z_epsm1_ifr_q,inter_bgrp_comm)
  IF(l_macropol) CALL mp_sum(z_head_ifr,inter_bgrp_comm)
  !
  ! EPS-1 refreq
  !
  ALLOCATE(zmatilda(mypara%nglob,mypara%nglob))
  ALLOCATE(zlambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use))
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
  CALL mp_barrier(world_comm)
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
     CALL io_push_title("(O)ptics")
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
     DO ifreq = 1, rfr%nloc
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
     !CALL serial_table_output(mpime==root,4000,'optics',out_tabella,&
     !& rfr%nglob,8,&
     !& (/'     E[eV]','      eps1','      eps2','      EELF','         n','         k','      Refl',' pol[au^3]'/))
     !
     IF( mpime == root ) THEN
        !
        CALL json%initialize()
        !
        CALL json%add("e",out_tabella(:,1))
        CALL json%add("eps1",out_tabella(:,2))
        CALL json%add("eps2",out_tabella(:,3))
        CALL json%add("EELF",out_tabella(:,4))
        CALL json%add("n",out_tabella(:,5))
        CALL json%add("k",out_tabella(:,6))
        CALL json%add("refl",out_tabella(:,7))
        CALL json%add("pol",out_tabella(:,8))
        !
        OPEN( NEWUNIT=iunit, FILE=TRIM(wfreq_save_dir)//"/optics.json" )
        CALL json%print( iunit )
        CLOSE( iunit )
        !
        CALL json%destroy()
        !
     ENDIF
     !
     time_spent(2) = get_clock( 'optics' )
     !
     WRITE(stdout,'(  5x," ")')
     CALL io_push_bar()
     WRITE(stdout, "(5x, 'File ',a,' written in ',a20)") TRIM(wfreq_save_dir)//"/optics.json",&
     & human_readable_time(time_spent(2)-time_spent(1))
     CALL io_push_bar()
     !
     DEALLOCATE( out_tabella )
     !
  ENDIF
  !
END SUBROUTINE
