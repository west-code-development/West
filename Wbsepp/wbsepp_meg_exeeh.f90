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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
SUBROUTINE wbsepp_meg ( )
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : e2, fpi
  USE cell_base,            ONLY : alat, tpiba2, omega
  USE distribution_center,  ONLY : pert
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE westcom,              ONLY : nbnd_occ,ev,npwq0x,npwq0,dvg,sqvc,fftdriver,&
                                 & dvg_exc,d0psi,n_plep_read_from_file
  USE pwcom,                ONLY : wk,nks,nelup,neldw,isk,g,igk_k,ngm,tpiba2,xk,omega,npw,npwx,lsda,nkstot,&
                                 & current_k,ngk,nbnd,wg
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_at_gamma,         ONLY : double_fwfft_gamma,double_invfft_gamma
  USE plep_db,              ONLY : plep_db_read
  USE pdep_db,              ONLY : pdep_db_read
  USE mp_world,             ONLY : mpime
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_image_comm, intra_bgrp_comm
  USE wavefunctions_module, ONLY : evc, psic
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart
  USE bse_module,           ONLY : bseparal
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER    :: nbndx_occ,num_triplepair,ibnd, jbnd, alnd, iks
  INTEGER    :: il1, ig1, il, jl, kl, index_i, ig
  INTEGER    :: n_pdep, n_plep, n_pdep_read_from_file
  INTEGER    :: exc1, exc2, exc3
  INTEGER(DP), ALLOCATABLE :: M_ijk(:,:)
  REAL(DP)   :: ev_tmp, delta_e, qg2
  REAL(DP)   :: exchange_eeh_tmp
  REAL(DP)   :: alpha_ija_vx, alpha_ija_vc, alpha_aij_vx, alpha_aij_vc
  REAL(DP)   :: xq(3)
  REAL(DP), ALLOCATABLE :: coeff_aij_vx(:,:,:),coeff_aij_vc(:,:,:)
  REAL(DP), ALLOCATABLE :: coeff_ija_vx(:,:,:),coeff_ija_vc(:,:,:)
  REAL(DP), ALLOCATABLE :: ev_exc(:), exchange_eeh(:)
  COMPLEX(DP), ALLOCATABLE :: psic_i(:), psic_j(:)
  COMPLEX(DP), ALLOCATABLE :: phi_c(:), phi_x(:)
  REAL(DP), EXTERNAL       :: DDOT
  !
  CALL start_clock( 'meg_calculation')
  !
  ! ... INITIALIZATION
  !
  n_pdep_read_from_file = 10
  n_plep_read_from_file = 10
  !
  n_pdep  = n_pdep_read_from_file
  n_plep  = n_plep_read_from_file
  nbndx_occ = MAXVAL(nbnd_occ)
  delta_e = 0.05
  xq(:)   = 0.0_DP
  iks     = 1
  !  
  pert = idistribute()
  CALL pert%init(num_triplepair,'b','num_triple_pairs',.TRUE.)
  !
  ALLOCATE( dvg_exc( npwx, nbndx_occ, n_plep, nks))
  ALLOCATE( ev(n_plep))
  ALLOCATE( ev_exc(n_plep))
  !
  CALL plep_db_read( n_plep )
  ev_exc(:) = ev(:)  
  !
  DEALLOCATE(ev)
  !
  ALLOCATE( dvg( npwq0x, n_pdep))
  ALLOCATE( ev(  n_pdep ))
  !
  CALL pdep_db_read( n_pdep )
  !
  ! Main program
  !
  ALLOCATE (M_ijk(n_plep*n_plep*n_plep,3))
  !
  M_ijk(:,:) = 0
  !
  index_i = 0
  !
  DO il = 1, n_plep
     ! 
     DO jl = 1, n_plep
        !
        DO kl = 1, n_plep
           !
           ev_tmp = ev_exc( il ) - (ev_exc( jl ) + ev_exc( kl )) 
           !
           IF ( abs(ev_tmp) <= delta_e) THEN
              !
              index_i = index_i + 1
              !
              M_ijk(index_i,1) = il
              M_ijk(index_i,2) = jl
              M_ijk(index_i,3) = kl
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO 
  !
  num_triplepair = index_i
  !
  ALLOCATE (exchange_eeh(num_triplepair))
  !
  exchange_eeh(:) = 0.0_DP
  !
  ! ... DISTRIBUTE num_triplepair 
  !
  bseparal = idistribute()
  CALL bseparal%init(num_triplepair,'i','num_triple_pairs',.TRUE.)
  CALL wbse_memory_report() ! Before allocating I report the memory required. 
  ! 
  DO il1 = 1, bseparal%nloc
     !
     ig1 = bseparal%l2g(il1)
     ! 
     IF( ig1 < 1 .OR. ig1 > num_triplepair ) CYCLE
     !
     exc1 = M_ijk(ig1, 1)
     exc2 = M_ijk(ig1, 2)
     exc3 = M_ijk(ig1, 3)
     !
     ! EXCHANGE EEH CHANNEL
     !
     ALLOCATE(coeff_aij_vx(nbndx_occ,nbndx_occ,alnd))
     ALLOCATE(coeff_aij_vc(nbndx_occ,nbndx_occ,alnd))
     ALLOCATE(coeff_ija_vx(nbndx_occ,nbndx_occ,alnd))
     ALLOCATE(coeff_ija_vc(nbndx_occ,nbndx_occ,alnd))
     !
     DO ibnd = 1, nbndx_occ
        !
        DO jbnd = 1, nbndx_occ
           !
           ALLOCATE (psic_i(dfftp%nnr), psic_j(dfftp%nnr))
           !
           psic_i (:) = (0.0d0,0.0d0)
           psic_j (:) = (0.0d0,0.0d0)
           !
           IF (gamma_only) THEN
              !           
              CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc(1,ibnd,exc1,iks),dvg_exc(1,jbnd,exc3,iks), psic_i,'Wave')
              !
              CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc(1,ibnd,exc2,iks), evc(1,jbnd), psic_j,'Wave')
              !
           ELSE
              !
              STOP
              !
           ENDIF
           !
           DO alnd = 1, n_pdep
              !
              ALLOCATE(phi_c(npwq0x),phi_x(npwq0x),psic(dffts%nnr))
              !
              psic(:) = (0.0d0,0.0d0)
              ! 
              IF (gamma_only) THEN
                 !
                 ! 1) compute phi_c = phi_i * sqrt(v_c*)
                 !    G->0 : v_c(G) = 0.0
                 ! 
                 phi_c(:) = (0.d0,0.d0)
                 DO ig = 1, npwq0x
                    ! 
                    qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
                    !  
                    IF (qg2 > 1.d-8) THEN
                       !
                       phi_c(ig) = dvg(ig,alnd) * DSQRT(e2*fpi/(tpiba2*qg2))
                       !
                    ENDIF
                    !
                 ENDDO
                 !
                 ! 2) compute phi_x = phi_i * sqrt(v_c)
                 ! using F-B correction
                 !
                 phi_x(:) = (0.d0,0.d0)
                 phi_x(:) = dvg(:,alnd) * sqvc(:)
                 !
                 CALL double_invfft_gamma(dffts,npwq0,npwq0x,phi_x,phi_c,psic,TRIM(fftdriver))
                 !
                 alpha_ija_vx  = 0.0_DP
                 alpha_ija_vc  = 0.0_DP
                 alpha_ija_vx  = SUM(DBLE(psic_i(:)) * AIMAG(psic_i(:)) *  DBLE(psic(:)))
                 alpha_ija_vc  = SUM(DBLE(psic_i(:)) * AIMAG(psic_i(:)) * AIMAG(psic(:))) &
                                                     * ev(alnd)/(1.0_DP-ev(alnd))
                 !
                 alpha_ija_vx  = alpha_ija_vx * wg(ibnd,iks)/(dffts%nr1*dffts%nr2*dffts%nr3)
                 alpha_ija_vc  = alpha_ija_vc * wg(ibnd,iks)/(dffts%nr1*dffts%nr2*dffts%nr3)
                 !
                 CALL mp_sum(alpha_ija_vx,intra_bgrp_comm)
                 CALL mp_sum(alpha_ija_vc,intra_bgrp_comm)
                 !
                 coeff_ija_vx(ibnd,jbnd,alnd) =  alpha_ija_vx
                 coeff_ija_vc(ibnd,jbnd,alnd) =  alpha_ija_vc
                 !
                 alpha_aij_vx  = 0.0_DP
                 alpha_aij_vc  = 0.0_DP
                 alpha_aij_vx  = SUM(DBLE(psic_j(:)) * AIMAG(psic_j(:)) *  DBLE(psic(:)))
                 alpha_aij_vc  = SUM(DBLE(psic_j(:)) * AIMAG(psic_j(:)) * AIMAG(psic(:))) &
                                                     * ev(alnd)/(1.0_DP-ev(alnd))
                 !
                 alpha_aij_vx  = alpha_aij_vx * wg(ibnd,iks)/(dffts%nr1*dffts%nr2*dffts%nr3)
                 alpha_aij_vc  = alpha_aij_vc * wg(ibnd,iks)/(dffts%nr1*dffts%nr2*dffts%nr3)
                 !
                 CALL mp_sum(alpha_aij_vx,intra_bgrp_comm)
                 CALL mp_sum(alpha_aij_vc,intra_bgrp_comm)
                 !
                 coeff_aij_vx(ibnd,jbnd,alnd) =  alpha_aij_vx
                 coeff_aij_vc(ibnd,jbnd,alnd) =  alpha_aij_vc
                 !
              ELSE
                 STOP
              ENDIF
              !
              DEALLOCATE (phi_c,phi_x,psic)
              !
           ENDDO
           !
           DEALLOCATE (psic_i, psic_j)
           !
        ENDDO
        !
     ENDDO
     !  
     exchange_eeh_tmp = 0.0d0
     !
     DO ibnd = 1, nbndx_occ 
        !
        DO jbnd = 1, nbndx_occ
           !
           DO alnd = 1, n_pdep
              !
              exchange_eeh_tmp = exchange_eeh_tmp + coeff_ija_vc(ibnd,jbnd,alnd) * coeff_aij_vc(ibnd,jbnd,alnd) &
                                                  + coeff_ija_vx(ibnd,jbnd,alnd) * coeff_aij_vx(ibnd,jbnd,alnd)
              !
           ENDDO 
           ! 
        ENDDO  
        ! 
     ENDDO 
     !
     exchange_eeh(ig1) = exchange_eeh_tmp 
     !
     DEALLOCATE(coeff_aij_vx,coeff_aij_vc)
     DEALLOCATE(coeff_ija_vx,coeff_ija_vc)
     !
  ENDDO
  !
  CALL mp_sum(exchange_eeh, inter_image_comm)
  !
  DEALLOCATE( exchange_eeh )
  DEALLOCATE( M_ijk )
  DEALLOCATE( ev )
  DEALLOCATE( dvg )
  DEALLOCATE( ev_exc )
  DEALLOCATE( dvg_exc )
  !
  CALL stop_clock( 'meg_calculation' )
  !
  RETURN
  !
END SUBROUTINE
