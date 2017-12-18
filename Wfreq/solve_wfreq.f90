!
! Copyright (C) 2015-2017 M. Govoni 
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
SUBROUTINE solve_wfreq(l_read_restart,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  USE control_flags,        ONLY : gamma_only 
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  LOGICAL,INTENT(IN) :: l_read_restart,l_generate_plot
  !
  IF( gamma_only ) THEN 
    CALL solve_wfreq_gamma( l_read_restart,l_generate_plot )
  ELSE
    CALL solve_wfreq_k( l_read_restart,l_generate_plot )
  ENDIF
  !
END SUBROUTINE 
!
!-----------------------------------------------------------------------
SUBROUTINE solve_wfreq_gamma(l_read_restart,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP 
  USE westcom,              ONLY : west_prefix,n_pdep_eigen_to_use,n_lanczos,npwq,l_macropol,iks_l2g,d_epsm1_ifr,z_epsm1_rfr,&
                                 & l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,wfreq_eta,imfreq_list,refreq_list,tr2_dfpt,&
                                 & z_head_rfr,d_head_ifr,o_restart_time,l_skip_nl_part_of_hcomr,npwqx,fftdriver, wstat_save_dir
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE mp_world,             ONLY : mpime
  USE mp,                   ONLY : mp_bcast,mp_barrier,mp_sum
  USE io_global,            ONLY : stdout,ionode
  USE gvect,                ONLY : g,ngm,gstart
  USE gvecw,                ONLY : gcutw
  USE cell_base,            ONLY : tpiba2,bg,omega
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : tpi,fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,current_k,wk,ngk
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
  USE distribution_center,  ONLY : pert,macropert,ifr,rfr
  USE class_idistribute,    ONLY : idistribute 
  USE wfreq_restart,        ONLY : solvewfreq_restart_write,solvewfreq_restart_read,bks_type
  USE class_bz_grid,        ONLY : bz_grid
  USE types_bz_grid,        ONLY : k_grid, q_grid
  USE chi_invert,           ONLY : chi_invert_real, chi_invert_complex
  USE class_coulomb,        ONLY : coulomb
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
  INTEGER :: i1,i2,i3,im,ip,ig,glob_ip,ir,iv,iks,ipol,m
  CHARACTER(LEN=512)    :: fname
  CHARACTER(LEN=6)      :: my_label_b
  COMPLEX(DP),ALLOCATABLE :: auxr(:)
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:), braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi_tmp(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP) :: zkonstant
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  REAL(DP),ALLOCATABLE :: ps_r(:,:)
  TYPE(idistribute) :: mypara
  REAL(DP),ALLOCATABLE :: overlap(:,:)
  REAL(DP) :: mwo, ecv, dfactor, frequency, dhead
  COMPLEX(DP) :: zmwo, zfactor,zm,zp, zhead
  INTEGER :: glob_jp,ic,ifreq,il
  REAL(DP),ALLOCATABLE :: dmatilda(:,:), dlambda(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatilda(:,:), zlambda(:,:)
  REAL(DP),ALLOCATABLE :: dmati(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatr(:,:,:)
  LOGICAL :: l_iks_skip, l_iv_skip
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bks_type) :: bks
  REAL(DP),ALLOCATABLE :: eprec(:)
  INTEGER :: ierr
  REAL(DP),ALLOCATABLE :: e(:)
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
!  ALLOCATE( sqvc(npwqx) )
  !
  IF(l_read_restart) THEN
     CALL solvewfreq_restart_read( bks, dmati, zmatr, mypara%nglob, mypara%nloc )
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
     DO iv = 1, nbnd_occ(iks)
        IF(iks==bks%lastdone_ks .AND. iv <= bks%lastdone_band ) CYCLE
        barra_load = barra_load + 1
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
  CALL pot3D%init('Wave','default')
!  CALL store_sqvc( sqvc, npwq, 'spherical', 1, .FALSE., isz )
  !CALL store_sqvc(sqvc,npwq,1,isz,.FALSE.)
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
     mwo = - wk(iks) / omega
     zmwo = CMPLX( - wk(iks) / omega, 0._DP, KIND=DP)
     !
     bks%max_band = nbndval
     bks%min_band = 1
     !
     ALLOCATE(dvpsi(npwx*npol,mypara%nlocx)) 
     !
     time_spent(1) = get_clock( 'wlanczos' ) 
     !
     ! LOOP over band states 
     !
     DO iv = 1, nbndval
        IF(iks==bks%lastdone_ks .AND. iv <= bks%lastdone_band ) CYCLE
        !
        ! MACROPOL CASE
        !
        IF(l_macropol) THEN
           !
           ! PHI 
           !
           ALLOCATE(phi(npwx*npol,3))
           ALLOCATE(phi_tmp(npwx*npol,3))
           CALL commutator_Hx_psi (iks, 1, 1, evc(1,iv), phi_tmp(1,1), l_skip_nl_part_of_hcomr)
           CALL commutator_Hx_psi (iks, 1, 2, evc(1,iv), phi_tmp(1,2), l_skip_nl_part_of_hcomr)
           CALL commutator_Hx_psi (iks, 1, 3, evc(1,iv), phi_tmp(1,3), l_skip_nl_part_of_hcomr)
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
           DEALLOCATE( eprec )
           DEALLOCATE( e )
           DEALLOCATE( phi_tmp )
           !
        ENDIF
        !
        ! PSIC
        !
!        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,iv),psic,'Wave') 
!        ELSEIF(noncolin) THEN
!           CALL SINGLEBAND_invfft_k(npw,evc(1     ,iv),npwx,psic_nc(1,1),dffts%nnr,.TRUE.)
!           CALL SINGLEBAND_invfft_k(npw,evc(1+npwx,iv),npwx,psic_nc(1,2),dffts%nnr,.TRUE.)
!        ELSE
!           CALL SINGLEBAND_invfft_k(npw,evc(1,iv),npwx,psic,dffts%nnr,.TRUE.)
!        ENDIF
        !
        ! ZEROS
        !
        dvpsi = 0._DP   
        !
        ! Read PDEP
        !
        ALLOCATE( pertg(npwqx) )
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
              ! Exhume dbs eigenvalue
              !
              WRITE(my_label_b,'(i6.6)') glob_ip
              fname = TRIM( wstat_save_dir ) // "/E"//TRIM(ADJUSTL(my_label_b))//".json"
              CALL pdep_read_G_and_distribute(fname,pertg)
              !
              ! Multiply by sqvc
              DO ig = 1, npwq
                 pertg(ig) = pot3D%sqvc(ig) * pertg(ig)
              ENDDO
              !
              ! Bring it to R-space
!              IF(gamma_only) THEN
                 CALL single_invfft_gamma(dffts,npwq,npwqx,pertg(1),pertr,TRIM(fftdriver))
                 DO ir=1,dffts%nnr 
                    pertr(ir)=psic(ir)*pertr(ir)
                 ENDDO
                 CALL single_fwfft_gamma(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave')
!              ELSEIF(noncolin) THEN
!                 CALL SINGLEBAND_invfft_k(npwq,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
!                 DO ir=1,dffts%nnr 
!                    pertr(ir)=psic_nc(ir,1)*pertr(ir)
!                 ENDDO
!                 CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
!                 CALL SINGLEBAND_invfft_k(npwq,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
!                 DO ir=1,dffts%nnr 
!                    pertr(ir)=psic_nc(ir,2)*pertr(ir)
!                 ENDDO
!                 CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1+npwx,ip),npwx,.TRUE.)
!              ELSE
!                 CALL SINGLEBAND_invfft_k(npwq,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
!                 DO ir=1,dffts%nnr 
!                    pertr(ir)=psic(ir)*pertr(ir)
!                 ENDDO
!                 CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
!              ENDIF 
              !
           ELSE
              !
              ipol = glob_ip-n_pdep_eigen_to_use
              !
              dvpsi(:,ip) = phi(:,ipol) * DSQRT(fpi * e2)
              !
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
        CALL apply_alpha_pc_to_m_wfcs(nbndval,mypara%nloc,dvpsi,(1._DP,0._DP))
        !
        ! OVERLAP( glob_ip, im=1:nbnd ) = < psi_im iks | dvpsi_glob_ip >
        !
        IF(ALLOCATED(ps_r)) DEALLOCATE(ps_r)
        ALLOCATE(ps_r(nbnd-nbndval,mypara%nloc))
        CALL glbrak_gamma(evc(1,nbndval+1),dvpsi(1,1),ps_r,npw,npwx,nbnd-nbndval,mypara%nloc,nbnd-nbndval,npol)
        CALL mp_sum(ps_r,intra_bgrp_comm) 
        !
        IF(ALLOCATED(overlap)) DEALLOCATE(overlap)
        ALLOCATE(overlap(mypara%nglob, nbndval+1:nbnd ) )
        overlap = 0._DP
        DO ic = nbndval+1, nbnd
           DO ip = 1, mypara%nloc
              glob_ip = mypara%l2g(ip)
              overlap(glob_ip,ic) = ps_r(ic-nbndval,ip)
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
           DO ic = nbndval+1, nbnd 
              !
              ecv = et(ic,iks)-et(iv,iks)
              dfactor = mwo * 2._DP * ecv / ( ecv**2 + frequency**2 )
              !
              DO ip = 1, mypara%nloc
                 glob_ip = mypara%l2g(ip)
                 DO glob_jp = 1, mypara%nglob
                    !
                    dmati(glob_jp,ip,ifreq) = dmati(glob_jp,ip,ifreq) + overlap( glob_jp, ic) * overlap( glob_ip, ic) * dfactor
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
              ecv = et(ic,iks)-et(iv,iks)
              zp = CMPLX( ecv + frequency, - wfreq_eta, KIND=DP )
              zm = CMPLX( ecv - frequency, - wfreq_eta, KIND=DP )
              zfactor = zmwo / zp + zmwo / zm 
              !
              DO ip = 1, mypara%nloc
                 glob_ip = mypara%l2g(ip)
                 DO glob_jp = 1, mypara%nglob
                    !
                    zmatr(glob_jp,ip,ifreq) = zmatr(glob_jp,ip,ifreq) + CMPLX( overlap( glob_jp, ic) * &
                    & overlap( glob_ip, ic), 0._DP, KIND=DP ) * zfactor
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
           ALLOCATE( bnorm    (                             mypara%nloc ) )
           ALLOCATE( diago    (               n_lanczos   , mypara%nloc ) )
           ALLOCATE( subdiago (               n_lanczos-1 , mypara%nloc ) )
           ALLOCATE( q_s      ( npwx*npol   , mypara%nloc , n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM 
           !
           CALL solve_deflated_lanczos_w_full_ortho ( nbnd, mypara%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
           !
           ALLOCATE( braket   ( mypara%nglob, n_lanczos   , mypara%nloc ) )
           CALL get_brak_hyper_parallel(dvpsi,mypara%nloc,n_lanczos,q_s,braket,mypara%nloc,mypara%nlocx,mypara%nglob)
           DEALLOCATE( q_s )
           !
           DO ip = 1, mypara%nloc
              CALL diago_lanczos( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), mypara%nglob )
           ENDDO
           !
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
           ! MPI-IO 
           !
           !CALL writeout_solvewfreq( iks_l2g(iks), iv, diago, braket, io_comm, mypara%nloc, mypara%nglob, mypara%myoffset )
           !
           DEALLOCATE( diago ) 
           DEALLOCATE( braket )
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'wlanczos' ) 
        !
        IF( o_restart_time >= 0._DP ) THEN 
           IF( (time_spent(2)-time_spent(1)) > o_restart_time*60._DP .OR. iv == nbndval ) THEN 
              bks%lastdone_ks=iks
              bks%lastdone_band=iv
              CALL solvewfreq_restart_write(bks,dmati,zmatr,mypara%nglob,mypara%nloc)
              bks%old_ks = iks 
              bks%old_band = iv
              time_spent(1) = get_clock( 'wlanczos' )
           ENDIF
        ENDIF
        !
        CALL update_bar_type( barra, 'wlanczos', 1 )
        !
     ENDDO ! BANDS
     !
     DEALLOCATE(dvpsi)
     !
  ENDDO ! KPOINT-SPIN 
  !
  CALL stop_bar_type( barra, 'wlanczos' )
  !
  ! EPS-1 imfreq
  !
  ALLOCATE( dmatilda( mypara%nglob, mypara%nglob ) )
  ALLOCATE( dlambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use) ) 
  ALLOCATE( d_epsm1_ifr( pert%nglob, pert%nloc, ifr%nloc) )
  d_epsm1_ifr = 0._DP
  IF(l_macropol) ALLOCATE( d_head_ifr( ifr%nloc) )
  !
  CALL start_clock('chi0_diag')
  !
  DO ifreq = 1, ifr%nloc
     !
     dmatilda = 0._DP
     DO ip = 1, mypara%nloc
        glob_ip = mypara%l2g(ip) 
        dmatilda( :, glob_ip) = dmati( :, ip, ifreq )
     ENDDO
     !
     CALL mp_sum( dmatilda, inter_image_comm )
     ! 
     CALL chi_invert_real( dmatilda, dhead, dlambda, mypara%nglob)
     !
     DO ip = 1, pert%nloc
        glob_ip = pert%l2g(ip)
        d_epsm1_ifr(1:n_pdep_eigen_to_use,ip,ifreq) = dlambda( 1:n_pdep_eigen_to_use, glob_ip)
     ENDDO 
     IF( l_macropol ) d_head_ifr( ifreq) = dhead
     !
  ENDDO
  !
  CALL stop_clock('chi0_diag')
  !
  DEALLOCATE( dlambda )
  DEALLOCATE( dmatilda )
  DEALLOCATE(dmati)
  !
  ! EPS-1 refreq
  !
  ALLOCATE( zmatilda( mypara%nglob, mypara%nglob ) )
  ALLOCATE( zlambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use ) ) 
  ALLOCATE( z_epsm1_rfr( pert%nglob, pert%nloc, rfr%nloc) )
  z_epsm1_rfr = 0._DP
  IF(l_macropol) ALLOCATE( z_head_rfr( rfr%nloc) )
  !
  DO ifreq = 1, rfr%nloc
     !
     zmatilda = 0._DP
     DO ip = 1, mypara%nloc
        glob_ip = mypara%l2g(ip) 
        zmatilda( :, glob_ip) = zmatr( :, ip, ifreq )
     ENDDO
     !
     CALL mp_sum( zmatilda, inter_image_comm ) 
     CALL chi_invert_complex( zmatilda, zhead, zlambda, mypara%nglob)
     !
     DO ip = 1, pert%nloc
        glob_ip = pert%l2g(ip)
        z_epsm1_rfr(1:n_pdep_eigen_to_use,ip,ifreq) = zlambda( 1:n_pdep_eigen_to_use, glob_ip)
     ENDDO 
     IF( l_macropol ) z_head_rfr( ifreq) = zhead
     !
  ENDDO
  !
  DEALLOCATE( zlambda )
  DEALLOCATE( zmatilda )
  DEALLOCATE( zmatr ) 
  !
  IF(l_generate_plot) THEN
     CALL output_eps_head( )
  ENDIF
  !
END SUBROUTINE 
!
!-----------------------------------------------------------------------
SUBROUTINE solve_wfreq_k(l_read_restart,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP 
  USE westcom,              ONLY : west_prefix,n_pdep_eigen_to_use,n_lanczos,npwq,l_macropol,iks_l2g,z_epsm1_ifr_q,&
                                 &  z_epsm1_rfr_q,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,wfreq_eta,imfreq_list,refreq_list,tr2_dfpt,&
                                 & z_head_rfr,z_head_ifr,o_restart_time,l_skip_nl_part_of_hcomr,npwqx,fftdriver, wstat_save_dir,&
                                 & ngq, igq_q
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE mp_world,             ONLY : mpime
  USE mp,                   ONLY : mp_bcast,mp_barrier,mp_sum
  USE io_global,            ONLY : stdout,ionode
  USE gvect,                ONLY : g,ngm,gstart
  USE gvecw,                ONLY : gcutw
  USE cell_base,            ONLY : tpiba2,bg,omega
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : tpi,fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,current_k,wk,ngk
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
  USE distribution_center,  ONLY : pert,macropert,ifr,rfr
  USE class_idistribute,    ONLY : idistribute 
  USE wfreq_restart,        ONLY : solvewfreq_restart_write,solvewfreq_restart_read,bks_type
  USE class_bz_grid,        ONLY : bz_grid
  USE types_bz_grid,        ONLY : k_grid, q_grid, compute_phase
  USE chi_invert,           ONLY : chi_invert_complex
  USE class_coulomb,        ONLY : coulomb
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
  INTEGER :: i1,i2,i3,im,ip,ig,glob_ip,ir,iv,iks,ik,is,iq,ikqs,ipol,m
  CHARACTER(LEN=512)    :: fname
  CHARACTER(LEN=6)      :: my_label_b
  CHARACTER(LEN=5)      :: my_label_q
  COMPLEX(DP),ALLOCATABLE :: auxr(:)
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:)
  COMPLEX(DP),ALLOCATABLE :: braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi_tmp(:,:)
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
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
  LOGICAL :: l_iks_skip, l_iv_skip
  LOGICAL :: l_gammaq
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bks_type) :: bks
  REAL(DP),ALLOCATABLE :: eprec(:)
  INTEGER :: ierr
  REAL(DP),ALLOCATABLE :: e(:)
  REAL(DP) :: kpq(3), g0(3) 
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
!  ALLOCATE( sqvc(npwqx) )
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
     CALL solvewfreq_restart_read( bks, zmati_q, zmatr_q, mypara%nglob, mypara%nloc )
  ELSE
     bks%lastdone_ks   = 0 
     bks%lastdone_q    = 0 
     bks%lastdone_band = 0 
     bks%old_ks        = 0 
     bks%old_q         = 0 
     bks%old_band      = 0 
     bks%max_ks        = k_grid%np
     bks%min_ks        = 1 
  ENDIF
  !
  barra_load = 0
  DO iq = 1, q_grid%np
     IF (iq<bks%lastdone_q) CYCLE
     DO iks = 1, k_grid%nps
        IF(iks<bks%lastdone_ks) CYCLE
        DO iv = 1, nbnd_occ(iks)
           IF(iks==bks%lastdone_ks .AND. iv <= bks%lastdone_band ) CYCLE
           barra_load = barra_load + 1
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
     IF (iq<bks%lastdone_q) CYCLE
     !
     npwq = ngq(iq)
     l_gammaq = q_grid%l_pIsGamma(iq)
     !
     CALL pot3D%init('Wave','default',iq)
!     CALL store_sqvc( sqvc, npwq, 'spherical', iq, .TRUE., isz )
     !IF (l_gammaq) THEN
     !   CALL store_sqvc(sqvc,npwq,1,isz,.FALSE.)
     !ELSE
     !   CALL store_sqvc_q(sqvc,npwq,1,iq,.TRUE.)
     !ENDIF
     !
     DO iks = 1, k_grid%nps   ! KPOINT-SPIN
        !
        ik = k_grid%ip(iks) 
        is = k_grid%is(iks) 
        !
        IF(iq==bks%lastdone_q .AND. iks<bks%lastdone_ks) CYCLE
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
           !iuwfc = 20
           !lrwfc = nbnd * npwx * npol 
           !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
           IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
           !CALL mp_bcast(evc,0,inter_image_comm)
           !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
           CALL mp_bcast(evc,0,inter_image_comm)
        ENDIF
!       !
!       ! ... Needed for LDA+U
!       !
!       IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
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
!       !
!       CALL init_us_2 (npw, igk, xk (1, iks), vkb)
        !
        !ikqs = kpq_grid%index_kq(iks,iq)
        !npwkq = ngk(ikqs)
        !
        CALL k_grid%find( k_grid%p_cart(:,ik) + q_grid%p_cart(:,iq), is, 'cart', ikqs, g0 )
        !CALL k_grid%add( k_grid%p_cart(:,ik), q_grid%p_cart(:,iq), kpq, g0, 'cart' ) 
        !ikqs = k_grid%find( kpq, 'cart' )
        !
        npwkq = ngk(ikqs)
        !
        CALL compute_phase( g0, 'cart', phase )
        !!
        !! computes the phase needed to bring the wavefunction at k+q 
        !! to the equivalent [k+q] point in the first BZ
        !!
        !CALL kpq_grid%get_phase(iks,iq)
        !phase = kpq_grid%phase
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
        bks%max_band = nbndval
        bks%min_band = 1
        !
        ALLOCATE(dvpsi(npwx*npol,mypara%nlocx)) 
        !
        time_spent(1) = get_clock( 'wlanczos' ) 
        !
        ! LOOP over band states 
        !
        DO iv = 1, nbndval
           IF(iq==bks%lastdone_q .AND. iks==bks%lastdone_ks .AND. iv <= bks%lastdone_band ) CYCLE
           !
           ! MACROPOL CASE
           !
           IF(l_macropol .AND. l_gammaq) THEN
              !
              ! PHI 
              !
              ALLOCATE(phi(npwx*npol,3))
              ALLOCATE(phi_tmp(npwx*npol,3))
              CALL commutator_Hx_psi (iks, 1, 1, evc(1,iv), phi_tmp(1,1), l_skip_nl_part_of_hcomr)
              CALL commutator_Hx_psi (iks, 1, 2, evc(1,iv), phi_tmp(1,2), l_skip_nl_part_of_hcomr)
              CALL commutator_Hx_psi (iks, 1, 3, evc(1,iv), phi_tmp(1,3), l_skip_nl_part_of_hcomr)
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
              DEALLOCATE( eprec )
              DEALLOCATE( e )
              DEALLOCATE( phi_tmp )
              !
              !
           ENDIF
           !
           ! PSIC
           !
!        IF(gamma_only) THEN
!           CALL SINGLEBAND_invfft(npw,evc(1,iv),npwx,psic,dffts%nnr) 
!        ELSE
           IF(noncolin) THEN
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(1     ,iv),psick_nc(1,1),'Wave',igk_k(1,ikqs))
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(1+npwx,iv),psick_nc(1,2),'Wave',igk_k(1,ikqs))
           ELSE
              CALL single_invfft_k(dffts,npwkq,npwx,evckpq(1,iv),psick,'Wave',igk_k(1,ikqs))
           ENDIF
           !
           ! ZEROS
           !
           dvpsi = 0._DP   
           !
           ! Read PDEP
           !
           ALLOCATE( pertg(npwqx) )
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
                 ! Exhume dbs eigenvalue
                 !
                 WRITE(my_label_b,'(i6.6)') glob_ip
                 WRITE(my_label_q,'(i5.5)') iq
                 fname = TRIM( wstat_save_dir ) // "/EQ"//TRIM(ADJUSTL(my_label_q))//"_"//TRIM(ADJUSTL(my_label_b))//".json"
                 CALL pdep_read_G_and_distribute(fname,pertg,iq)
                 !
                 ! Multiply by sqvc
                 DO ig = 1, npwq
                    pertg(ig) = pot3D%sqvc(ig) * pertg(ig)
                 ENDDO
                 !
                 ! Bring it to R-space
!                 IF(gamma_only) THEN
!                    CALL SINGLEBAND_invfft(npwq,pertg(1),npwx,pertr,dffts%nnr)
!                    DO ir=1,dffts%nnr 
!                       pertr(ir)=psic(ir)*pertr(ir)
!                    ENDDO
!                    CALL SINGLEBAND_fwfft(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx)
!                 ELSE
                 IF(noncolin) THEN
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                    DO ir=1,dffts%nnr 
                       pertr(ir)=phase(ir)*psick_nc(ir,1)*DCONJG(pertr(ir))
                    ENDDO
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,TRIM(fftdriver))
                    DO ir=1,dffts%nnr 
                       pertr(ir)=phase(ir)*psick_nc(ir,2)*DCONJG(pertr(ir))
                    ENDDO
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1+npwx,ip),'Wave',igk_k(1,current_k))
                 ELSE
                    CALL single_invfft_k(dffts,npwq,npwqx,pertg(1),pertr,'Wave',igq_q(1,iq))
                    DO ir=1,dffts%nnr 
                       pertr(ir)=phase(ir)*psick(ir)*DCONJG(pertr(ir))
                    ENDDO
                    CALL single_fwfft_k(dffts,npw,npwx,pertr,dvpsi(1,ip),'Wave',igk_k(1,current_k))
                 ENDIF 
                 !
              ELSE
                 !
                 IF (l_gammaq) THEN 
                    ipol = glob_ip-n_pdep_eigen_to_use
                    dvpsi(:,ip) = phi(:,ipol) * DSQRT(fpi * e2)
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
                       zmati_q(glob_jp,ip,ifreq,iq) = zmati_q(glob_jp,ip,ifreq,iq) + DCONJG(overlap( glob_ip, ic)) * &
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
                       zmatr_q(glob_jp,ip,ifreq,iq) = zmatr_q(glob_jp,ip,ifreq,iq) + DCONJG( overlap( glob_ip, ic) ) * &
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
              ALLOCATE( bnorm    (                             mypara%nloc ) )
              ALLOCATE( diago    (               n_lanczos   , mypara%nloc ) )
              ALLOCATE( subdiago (               n_lanczos-1 , mypara%nloc ) )
              ALLOCATE( q_s      ( npwx*npol   , mypara%nloc , n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM 
              !
              CALL solve_deflated_lanczos_w_full_ortho ( nbnd, mypara%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
              !
              ALLOCATE( braket   ( mypara%nglob, n_lanczos   , mypara%nloc ) )
              CALL get_brak_hyper_parallel_complex(dvpsi,mypara%nloc,n_lanczos,q_s,braket,mypara%nloc,mypara%nlocx,mypara%nglob)
              DEALLOCATE( q_s )
              !
              DO ip = 1, mypara%nloc
                 CALL diago_lanczos_complex( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), mypara%nglob )
              ENDDO
              !
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
                          zmati_q(glob_jp,ip,ifreq,iq) = zmati_q(glob_jp,ip,ifreq,iq) + DCONJG(braket( glob_jp, il, ip)) * zfactor
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
                          zmatr_q(glob_jp,ip,ifreq,iq) = zmatr_q(glob_jp,ip,ifreq,iq) + DCONJG(braket( glob_jp, il, ip)) * zfactor
                          !
                       ENDDO
                    ENDDO
                    !
                 ENDDO ! il
              ENDDO ! ifreq
              !
              ! MPI-IO 
              !
              !CALL writeout_solvewfreq( iks_l2g(iks), iv, diago, braket, io_comm, mypara%nloc, mypara%nglob, mypara%myoffset )
              !
              DEALLOCATE( diago ) 
              DEALLOCATE( braket )
              !
           ENDIF ! l_enable_lanczos
           !
           time_spent(2) = get_clock( 'wlanczos' ) 
           !
           IF( o_restart_time >= 0._DP ) THEN 
              IF( (time_spent(2)-time_spent(1)) > o_restart_time*60._DP .OR. iv == nbndval) THEN
                 bks%lastdone_q=iq
                 bks%lastdone_ks=iks
                 bks%lastdone_band=iv
                 CALL solvewfreq_restart_write(bks,zmati_q(:,:,:,iq),zmatr_q(:,:,:,iq),mypara%nglob,mypara%nloc)
                 bks%old_q = iq
                 bks%old_ks = iks 
                 bks%old_band = iv
                 time_spent(1) = get_clock( 'wlanczos' )
              ENDIF
           ENDIF
           !
           CALL update_bar_type( barra, 'wlanczos', 1 )
           !
        ENDDO ! BANDS
        !
        DEALLOCATE(dvpsi)
        !
     ENDDO ! KPOINT-SPIN 
     !
  ENDDO ! QPOINT
  !
  CALL stop_bar_type( barra, 'wlanczos' )
  !
  DEALLOCATE( evckpq )
  DEALLOCATE( phase )
  !
  ! EPS-1 imfreq
  !
  ALLOCATE( zmatilda( mypara%nglob, mypara%nglob ) )
  ALLOCATE( zlambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use) ) 
  ALLOCATE( z_epsm1_ifr_q( pert%nglob, pert%nloc, ifr%nloc, q_grid%np) )
  z_epsm1_ifr_q = 0._DP
  IF(l_macropol) ALLOCATE( z_head_ifr( ifr%nloc) )
  !
  CALL start_clock('chi0_diag')
  !
  DO iq = 1, q_grid%np
     !
     l_gammaq = q_grid%l_pIsGamma(iq)
     !
     DO ifreq = 1, ifr%nloc
        !
        zmatilda = 0._DP
        DO ip = 1, mypara%nloc
           glob_ip = mypara%l2g(ip) 
           zmatilda( :, glob_ip) = zmati_q( :, ip, ifreq, iq )
        ENDDO
        !
        CALL mp_sum( zmatilda, inter_image_comm )
        ! 
        CALL chi_invert_complex( zmatilda, zhead, zlambda, mypara%nglob, l_gammaq)
        !
        DO ip = 1, pert%nloc
           glob_ip = pert%l2g(ip)
           z_epsm1_ifr_q(1:n_pdep_eigen_to_use,ip,ifreq,iq) = zlambda( 1:n_pdep_eigen_to_use, glob_ip)
        ENDDO 
        IF( l_macropol .AND. l_gammaq ) z_head_ifr( ifreq) = zhead
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_clock('chi0_diag')
  !
  DEALLOCATE( zlambda )
  DEALLOCATE( zmatilda )
  DEALLOCATE( zmati_q )
  !
  ! EPS-1 refreq
  !
  ALLOCATE( zmatilda( mypara%nglob, mypara%nglob ) )
  ALLOCATE( zlambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use ) ) 
  ALLOCATE( z_epsm1_rfr_q( pert%nglob, pert%nloc, rfr%nloc, q_grid%np) )
  z_epsm1_rfr_q = 0._DP
  IF(l_macropol) ALLOCATE( z_head_rfr( rfr%nloc) )
  !
  DO iq = 1, q_grid%np
     !
     l_gammaq = q_grid%l_pIsGamma(iq)
     !
     DO ifreq = 1, rfr%nloc
        !
        zmatilda = 0._DP
        DO ip = 1, mypara%nloc
           glob_ip = mypara%l2g(ip) 
           zmatilda( :, glob_ip) = zmatr_q( :, ip, ifreq, iq )
        ENDDO
        !
        CALL mp_sum( zmatilda, inter_image_comm ) 
        CALL chi_invert_complex( zmatilda, zhead, zlambda, mypara%nglob, l_gammaq)
        !
        DO ip = 1, pert%nloc
           glob_ip = pert%l2g(ip)
           z_epsm1_rfr_q(1:n_pdep_eigen_to_use,ip,ifreq,iq) = zlambda( 1:n_pdep_eigen_to_use, glob_ip)
        ENDDO 
        IF( l_macropol .AND. l_gammaq ) z_head_rfr( ifreq ) = zhead
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE( zlambda )
  DEALLOCATE( zmatilda )
  DEALLOCATE( zmatr_q ) 
  !
  IF(l_generate_plot) THEN 
     CALL output_eps_head( )
  ENDIF
  !
END SUBROUTINE 
!
!
SUBROUTINE output_eps_head( )
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : d_head_ifr,z_head_ifr,z_head_rfr,refreq_list,l_macropol,imfreq_list,wfreq_save_dir
  USE constants,            ONLY : rytoev,fpi
  !USE west_io,              ONLY : serial_table_output
  USE mp_world,             ONLY : mpime,root
  USE distribution_center,  ONLY : ifr,rfr
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
  CHARACTER(LEN=9) :: prefisso
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
        nn = DSQRT( 0.5_DP * ( ep1 + DSQRT( ep1**2 + ep2**2 ) ) ) 
        kk = DSQRT( 0.5_DP * (-ep1 + DSQRT( ep1**2 + ep2**2 ) ) )
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
        CALL json%print_file( iunit )
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
