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
SUBROUTINE solve_qp(l_secant,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_secant,l_generate_plot
  !
  IF( gamma_only ) THEN
     CALL solve_qp_gamma( l_secant, l_generate_plot )
  ELSE
     CALL solve_qp_k( l_secant, l_generate_plot )
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_qp_gamma(l_secant,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies.
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_lanczos,qp_bandrange,imfreq_list_integrate,&
                                 & n_secant_maxiter,trev_secant,l_enable_lanczos,imfreq_list,&
                                 & n_imfreq,d_epsm1_ifr,z_epsm1_rfr,n_spectralf,ecut_spectralf,&
                                 & d_body1_ifr,d_body2_ifr,d_diago,z_body_rfr,sigma_z,sigma_eqplin,&
                                 & sigma_eqpsec,sigma_sc_eks,sigma_sc_eqplin,sigma_sc_eqpsec,&
                                 & sigma_diff,sigma_spectralf,sigma_freq,l_enable_off_diagonal,&
                                 & ijpmap,npair,d_body1_ifr_full,d_body2_ifr_full,d_diago_full,&
                                 & z_body_rfr_full
  USE mp_global,            ONLY : inter_image_comm,inter_pool_comm,my_pool_id,intra_bgrp_comm,&
                                 & inter_bgrp_comm
  USE mp_world,             ONLY : mpime,root
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE pwcom,                ONLY : et,nbnd
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE constants,            ONLY : rytoev,pi
  USE west_io,              ONLY : serial_table_output
  USE distribution_center,  ONLY : pert,ifr,rfr,aband,band_group,kpt_pool
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wfreq_io,             ONLY : readin_overlap,readin_solvegfreq,readin_solvehf
  USE wfreq_db,             ONLY : wfreq_db_write
  USE types_bz_grid,        ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_secant,l_generate_plot
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: sigma_hf(:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_cor_in (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_cor_out(:,:)
  REAL(DP),ALLOCATABLE :: z_in(:,:)
  REAL(DP),ALLOCATABLE :: qp_energy(:,:)
  COMPLEX(DP),ALLOCATABLE :: sc(:,:,:)
  REAL(DP),ALLOCATABLE :: en(:,:,:)
  LOGICAL,ALLOCATABLE :: l_conv(:,:)
  REAL(DP),PARAMETER :: eshift = 0.007349862_DP ! = 0.1 eV
  INTEGER :: ib,ibloc,iks,ifixed,ip,glob_ip,ifreq,il,im,glob_im,glob_jp,glob_ifreq,iks_g,jb,index
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  CHARACTER(LEN=5) :: myglobk
  INTEGER :: notconv
  REAL(DP),ALLOCATABLE :: dtemp(:)
  REAL(DP),ALLOCATABLE :: dtemp2(:,:)
  COMPLEX(DP),ALLOCATABLE :: ztemp2(:,:)
  REAL(DP),ALLOCATABLE :: diago(:,:)
  REAL(DP),ALLOCATABLE :: braket(:,:,:)
  REAL(DP),ALLOCATABLE :: overlap(:,:), overlap1(:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  !
  CALL start_clock('solve_qp')
  !
  ALLOCATE( imfreq_list_integrate( 2, ifr%nloc ) )
  ALLOCATE( dtemp( n_imfreq ) )
  !
  dtemp = 0._DP
  DO ifreq = 1, ifr%nloc
     glob_ifreq = ifr%l2g(ifreq)
     dtemp( glob_ifreq ) = imfreq_list( ifreq )
  ENDDO
  CALL mp_sum( dtemp, intra_bgrp_comm )
  !
  DO ifreq = 1, ifr%nloc
     glob_ifreq = ifr%l2g(ifreq)
     IF( glob_ifreq == 1 ) THEN
        imfreq_list_integrate(1,ifreq) = 0._DP
     ELSE
        imfreq_list_integrate(1,ifreq) = ( dtemp(glob_ifreq) + dtemp(glob_ifreq-1) ) * 0.5_DP
     ENDIF
     !
     IF( glob_ifreq == n_imfreq ) THEN
        imfreq_list_integrate(2,ifreq) = dtemp(n_imfreq)
     ELSE
        imfreq_list_integrate(2,ifreq) = ( dtemp(glob_ifreq) + dtemp(glob_ifreq+1) ) * 0.5_DP
     ENDIF
  ENDDO
  !
  DEALLOCATE( dtemp )
  !
  ! TEMP
  ! 
  ALLOCATE( overlap( pert%nglob, nbnd ) )
  ALLOCATE( dtemp2( nbnd, ifr%nloc ) )
  ALLOCATE( ztemp2( nbnd, rfr%nloc ) )
  IF (l_enable_lanczos) THEN
     ALLOCATE( braket( pert%nglob, n_lanczos, pert%nloc ) )
     ALLOCATE( diago( n_lanczos, pert%nloc ) )
  ENDIF
  IF (l_enable_off_diagonal) THEN
     ALLOCATE( overlap1( pert%nglob, nbnd ) )
     ALLOCATE( d_body1_ifr_full( aband%nloc, ifr%nloc, 1:npair, k_grid%nps ) )
     ALLOCATE( z_body_rfr_full( aband%nloc, rfr%nloc, 1:npair, k_grid%nps ) )
     IF ( l_enable_lanczos ) THEN
        ALLOCATE( d_body2_ifr_full( n_lanczos, pert%nloc, ifr%nloc, 1:npair, k_grid%nps ) )
        ALLOCATE( d_diago_full( n_lanczos, pert%nloc, 1:npair, k_grid%nps ) )
     ENDIF
     !
     d_body1_ifr_full = 0._DP
     z_body_rfr_full = 0._DP
     IF ( l_enable_lanczos ) THEN
        d_body2_ifr_full = 0._DP    
        d_diago_full = 0._DP
     ENDIF
  ELSE
     !
     ALLOCATE( d_body1_ifr( aband%nloc, ifr%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps ) )
     ALLOCATE( z_body_rfr( aband%nloc, rfr%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps ) )
     IF( l_enable_lanczos ) THEN
        ALLOCATE( d_body2_ifr( n_lanczos, pert%nloc, ifr%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps ) )
        ALLOCATE( d_diago( n_lanczos, pert%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps ) )
     ENDIF
     !
     d_body1_ifr = 0._DP
     z_body_rfr = 0._DP
     IF( l_enable_lanczos ) THEN
        d_body2_ifr = 0._DP
        d_diago = 0._DP
     ENDIF
     !
  ENDIF
  !
  ! d_body1_ifr, d_body2_ifr, z_diago_rfr, d_diago
  !
  CALL io_push_title('Collecting results from W and G')
  !
  CALL band_group%init(qp_bandrange(2)-qp_bandrange(1)+1,'b','band_group',.FALSE.)
  !
  IF (l_enable_off_diagonal) THEN
     barra_load = kpt_pool%nloc * band_group%nloc * band_group%nglob
  ELSE
     barra_load = kpt_pool%nloc * band_group%nloc
  ENDIF
  CALL start_bar_type( barra, 'coll_gw', barra_load )
  !
  ! LOOP
  !
  DO iks = 1,kpt_pool%nloc ! KPOINT-SPIN
     !
     iks_g = kpt_pool%l2g(iks)
     !
     DO ibloc = 1, band_group%nloc
        ib = band_group%l2g(ibloc)+qp_bandrange(1)-1
        !
        CALL readin_overlap( 'g', kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
        !
        DO jb = qp_bandrange(1), qp_bandrange(2)
           !
           IF ( l_enable_off_diagonal .AND. jb <= ib ) THEN
              !
              index = ijpmap(jb,ib)
              !
              CALL readin_overlap( 'g', kpt_pool%l2g(iks), jb, overlap1, pert%nglob, nbnd )
           ELSEIF ( .NOT. l_enable_off_diagonal .AND. jb == ib ) THEN
              CONTINUE
           ELSE
              CYCLE
           ENDIF
           !
           ! ------
           ! d_body1_ifr
           ! ------
           !
           dtemp2 = 0._DP
           !
           DO ifreq = 1, ifr%nloc 
              !
              DO im = 1, nbnd 
                 !
                 DO glob_jp = 1, n_pdep_eigen_to_use
                    DO ip = 1, pert%nloc
                       glob_ip = pert%l2g(ip)
                       IF (l_enable_off_diagonal .AND. jb <= ib) THEN
                          dtemp2(im,ifreq) = dtemp2(im,ifreq) + overlap1(glob_jp,im)&
                          &*overlap(glob_ip,im)*d_epsm1_ifr(glob_jp,ip,ifreq)
                       ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                          dtemp2(im,ifreq) = dtemp2(im,ifreq) + overlap(glob_jp,im)*overlap(glob_ip,im)*d_epsm1_ifr(glob_jp,ip,ifreq)
                       ENDIF
                    ENDDO
                 ENDDO 
                 !
              ENDDO ! im
              !
              CALL mp_sum( dtemp2, inter_image_comm ) 
              !
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 IF (l_enable_off_diagonal .AND. jb <= ib) THEN
                    d_body1_ifr_full(im,ifreq,index,iks) = dtemp2(glob_im,ifreq)
                 ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                    d_body1_ifr(im,ifreq,ib,iks) = dtemp2(glob_im,ifreq)
                 ENDIF
              ENDDO
              !
           ENDDO ! ifreq
           !
           ! -----
           ! z_body_rfr
           ! -----
           !
           ztemp2 = 0._DP
           !
           DO ifreq = 1, rfr%nloc 
              !
              DO im = 1, nbnd 
                 !
                 DO glob_jp = 1, n_pdep_eigen_to_use
                    DO ip = 1, pert%nloc
                       glob_ip = pert%l2g(ip)
                       IF (l_enable_off_diagonal .AND. jb <= ib) THEN
                          ztemp2(im,ifreq) = ztemp2(im,ifreq) + overlap1(glob_jp,im)&
                          &*overlap(glob_ip,im)*z_epsm1_rfr(glob_jp,ip,ifreq)
                       ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                          ztemp2(im,ifreq) = ztemp2(im,ifreq) + overlap(glob_jp,im)*overlap(glob_ip,im)*z_epsm1_rfr(glob_jp,ip,ifreq) 
                       ENDIF
                    ENDDO
                 ENDDO 
                 !
              ENDDO ! im
              !
              CALL mp_sum( ztemp2, inter_image_comm ) 
              !
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 IF (l_enable_off_diagonal .AND. jb <= ib) THEN
                    z_body_rfr_full(im,ifreq,index,iks) = ztemp2(glob_im,ifreq)
                 ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                    z_body_rfr(im,ifreq,ib,iks) = ztemp2(glob_im,ifreq)
                 ENDIF
              ENDDO
              !
           ENDDO ! ifreq
           !
           ! -----------------------------
           ! LANCZOS part : d_diago, d_body2_ifr
           ! -----------------------------
           !
           IF( l_enable_lanczos ) THEN
              !
              IF (l_enable_off_diagonal .AND. jb <= ib) THEN
                 CALL readin_solvegfreq( kpt_pool%l2g(iks), index, diago, braket, pert%nloc, pert%nglob, pert%myoffset )
              ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                 CALL readin_solvegfreq( kpt_pool%l2g(iks), ib, diago, braket, pert%nloc, pert%nglob, pert%myoffset )
              ENDIF
              !
              DO ip = 1, pert%nloc
                 DO il = 1, n_lanczos 
                    IF (l_enable_off_diagonal .AND. jb <= ib) THEN
                       d_diago_full(il,ip,index,iks_g) = diago(il,ip) 
                    ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                       d_diago(il,ip,ib,iks_g) = diago(il,ip)
                    ENDIF
                 ENDDO
              ENDDO
              !
              DO ifreq = 1,ifr%nloc
                 DO ip = 1, pert%nloc
                    DO il = 1, n_lanczos
                       DO glob_jp = 1, pert%nglob
                          IF (l_enable_off_diagonal .AND. jb <= ib) THEN
                             d_body2_ifr_full(il,ip,ifreq,index,iks_g) = &
                             & d_body2_ifr_full(il,ip,ifreq,index,iks_g) + braket(glob_jp,il,ip) &
                             & *d_epsm1_ifr(glob_jp,ip,ifreq)
                          ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                             d_body2_ifr(il,ip,ifreq,ib,iks_g) = d_body2_ifr(il,ip,ifreq,ib,iks_g) + &
                             & braket(glob_jp,il,ip)*d_epsm1_ifr(glob_jp,ip,ifreq)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              !
           ENDIF
           !
        ENDDO ! jbnd
        !
        CALL update_bar_type( barra, 'coll_gw', qp_bandrange(2)-qp_bandrange(1)+1 )
        !
     ENDDO ! ibnd
     !
  ENDDO ! iks
  !
  CALL mp_sum(d_body1_ifr,inter_bgrp_comm)
  CALL mp_sum(d_body1_ifr,inter_pool_comm)
  CALL mp_sum(z_body_rfr,inter_bgrp_comm)
  CALL mp_sum(z_body_rfr,inter_pool_comm)
  !
  IF(l_enable_lanczos) THEN
     CALL mp_sum(d_body2_ifr,inter_bgrp_comm)
     CALL mp_sum(d_body2_ifr,inter_pool_comm)
     CALL mp_sum(d_diago,inter_bgrp_comm)
     CALL mp_sum(d_diago,inter_pool_comm)
  ENDIF
  !
  CALL stop_bar_type( barra, 'coll_gw' )
  !
  DEALLOCATE( overlap )
  DEALLOCATE( d_epsm1_ifr )
  DEALLOCATE( z_epsm1_rfr )
  DEALLOCATE( dtemp2 )
  DEALLOCATE( ztemp2 )
  IF( l_enable_lanczos ) THEN
     DEALLOCATE( braket )
     DEALLOCATE( diago )
  ENDIF
  !
  ! Get Sigma_X
  !
  ALLOCATE( sigma_hf (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  CALL readin_solvehf( sigma_hf(qp_bandrange(1),1), qp_bandrange(2)-qp_bandrange(1)+1, k_grid%nps )
  !
  ! For CORR
  !
  IF( l_secant ) THEN
     !
     CALL io_push_title('(Q)uasiparticle energies')
     !
     ALLOCATE( en(qp_bandrange(1):qp_bandrange(2),k_grid%nps,2) )
     ALLOCATE( sc(qp_bandrange(1):qp_bandrange(2),k_grid%nps,2) )
     ALLOCATE( l_conv(qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( sigma_cor_in (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( sigma_cor_out(qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( z_in (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( qp_energy (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     !
     en = 0._DP
     !
     ! 1st step of secant solver : E_KS - 0.5 * eshift
     ! 1st step of secant solver : E_KS + 0.5 * eshift
     !
     DO iks = 1, kpt_pool%nloc
        !
        iks_g = kpt_pool%l2g(iks)
        !
        DO ib = qp_bandrange(1), qp_bandrange(2)
           en(ib,iks_g,1) = et(ib,iks) - eshift*0.5_DP
           en(ib,iks_g,2) = et(ib,iks) + eshift*0.5_DP
        ENDDO
     ENDDO
     CALL mp_sum(en,inter_pool_comm)
     CALL calc_corr_gamma( sc(:,:,1), en(:,:,1), .TRUE., .FALSE.)
     CALL calc_corr_gamma( sc(:,:,2), en(:,:,2), .TRUE., .FALSE.)
     !
     ! Stage sigma_corr_in
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           sigma_cor_in(ib,iks) = ( sc(ib,iks,2) + sc(ib,iks,1) ) * 0.5_DP
        ENDDO
     ENDDO
     !
     ! Stage z_in
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           z_in(ib,iks) = 1._DP / ( 1._DP-REAL( sc(ib,iks,2)-sc(ib,iks,1), KIND=DP ) / eshift  )
        ENDDO
     ENDDO
     !
     ! en 1 = EKS, sc 1 = sigma_corr_in
     ! en 2 = EKS + Z * ( sigmax - vxc + sigmacorrin)
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           en(ib,iks,1) = et(ib,iks)
           sc(ib,iks,1) = sigma_cor_in(ib,iks)
           en(ib,iks,2) = et(ib,iks) + z_in(ib,iks) * ( sigma_hf(ib,iks) + REAL( sigma_cor_in(ib,iks),KIND=DP) )
        ENDDO
     ENDDO
     sigma_z           = z_in
     sigma_eqplin(:,:) = en(:,:,2)
     sigma_sc_eks      = sigma_cor_in
     CALL output_eqp_report(0,en(:,:,1),en(:,:,2),sigma_cor_in(:,:))
     !
     ! nth step of the secant solver
     !
     l_conv = .FALSE.
     notconv = k_grid%nps * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
     DO ifixed = 1, n_secant_maxiter
        !
        CALL calc_corr_gamma( sc(:,:,2), en(:,:,2), .TRUE., .FALSE.)
        !
        IF( my_pool_id == 0 ) THEN
           !
           ! Resulting new energy:
           !
           DO iks = 1, k_grid%nps
              DO ib = qp_bandrange(1), qp_bandrange(2)
                  IF( .NOT. l_conv(ib,iks) ) THEN
                     qp_energy(ib,iks) = en(ib,iks,2) &
                     & +(et(ib,iks)+sigma_hf(ib,iks)+REAL(sc(ib,iks,2),KIND=DP)-en(ib,iks,2)) &
                     & /(1._DP-REAL(sc(ib,iks,2)-sc(ib,iks,1),KIND=DP)/(en(ib,iks,2)-en(ib,iks,1)))
                  ENDIF
              ENDDO
           ENDDO
           !
           ! Estimate l_conv
           !
           DO iks = 1, k_grid%nps
              DO ib = qp_bandrange(1), qp_bandrange(2)
                 l_conv(ib,iks) = ( ABS( qp_energy(ib,iks) - en(ib,iks,2) ) < trev_secant )
              ENDDO
           ENDDO
           !
           ! Count the number of notconverged QP energies
           !
           notconv = 0
           DO iks = 1, k_grid%nps
              DO ib = qp_bandrange(1), qp_bandrange(2)
                 IF( .NOT.l_conv(ib,iks) ) notconv = notconv + 1
              ENDDO
           ENDDO
           !
           IF( ifixed == 1 ) THEN
              sigma_sc_eqplin(:,:) = sc(:,:,2)
           ENDIF
           !
        ENDIF
        !
        CALL mp_bcast(en(:,:,2),0,inter_pool_comm)
        CALL mp_bcast(sc(:,:,2),0,inter_pool_comm)
        CALL mp_bcast(qp_energy,0,inter_pool_comm)
        CALL mp_bcast(notconv,0,inter_pool_comm)
        !
        sigma_eqpsec = qp_energy
        sigma_sc_eqpsec = sc(:,:,2)
        sigma_diff(:,:) = qp_energy(:,:) - en(:,:,2)
        !
        IF( notconv == 0 ) THEN
           CALL io_push_title('CONVERGENCE ACHIEVED !!!')
           CALL output_eqp_report(-1,en(:,:,2),qp_energy,sc(:,:,2))
           EXIT
        ELSE
           CALL output_eqp_report(ifixed,en(:,:,2),qp_energy,sc(:,:,2))
        ENDIF
        !
        ! Iterate
        !
        en(:,:,1) = en(:,:,2)
        sc(:,:,1) = sc(:,:,2)
        en(:,:,2) = qp_energy(:,:)
        !
     ENDDO
     !   
     sigma_cor_out(:,:) = sc(:,:,2)
     !
     IF (l_enable_off_diagonal) CALL calc_corr_gamma( sc(:,:,2), en(:,:,2), .TRUE., .TRUE.) 
     !
     DEALLOCATE( en, sc, l_conv )
     !
     IF( notconv /= 0 ) THEN
        CALL io_push_title('CONVERGENCE **NOT** ACHIEVED !!!')
     ENDIF
     !
     ! Output it per k-point
     !
     ALLOCATE(out_tab(qp_bandrange(2)-qp_bandrange(1)+1,7))
     !
     DO iks=1,k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           out_tab( ib - qp_bandrange(1) + 1, 1) = REAL( ib, KIND=DP)
           out_tab( ib - qp_bandrange(1) + 1, 2) = et(ib,iks) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 3) = (et(ib,iks)+sigma_hf(ib,iks)) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 4) = qp_energy(ib,iks) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 5) = (qp_energy(ib,iks) - et(ib,iks) ) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 6) = REAL(  sigma_cor_out(ib,iks), KIND=DP ) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 7) = AIMAG( sigma_cor_out(ib,iks) ) * rytoev
        ENDDO
        WRITE(myglobk,'(i5.5)') iks
        !
        CALL serial_table_output(mpime==root,'eqp_K'//myglobk,out_tab,&
           & qp_bandrange(2)-qp_bandrange(1)+1,7,&
           & (/'      band','    E0[eV]','   EHF[eV]','   Eqp[eV]','Eqp-E0[eV]','Sc_Eqp[eV]',' Width[eV]'/))
     ENDDO
     !
     DEALLOCATE( out_tab )
     DEALLOCATE( sigma_cor_in )
     DEALLOCATE( sigma_cor_out )
     DEALLOCATE( z_in )
     DEALLOCATE( qp_energy )
     !
  ENDIF
  !
  IF( l_generate_plot ) THEN
     !
     CALL io_push_title('(P)lotting the QP corrections...')
     !
     CALL start_bar_type( barra, 'qplot', n_spectralf )
     !
     ALLOCATE( en(qp_bandrange(1):qp_bandrange(2),k_grid%nps,1) )
     ALLOCATE( sc(qp_bandrange(1):qp_bandrange(2),k_grid%nps,1) )
     !
     DO glob_ifreq = 1, n_spectralf
        sigma_freq(glob_ifreq) = (ecut_spectralf(2)-ecut_spectralf(1))/REAL(n_spectralf-1,KIND=DP)*REAL(glob_ifreq-1,KIND=DP) &
        & +ecut_spectralf(1)
     ENDDO
     !
     DO glob_ifreq = 1, n_spectralf
        en = (ecut_spectralf(2)-ecut_spectralf(1))/REAL(n_spectralf-1,KIND=DP)*REAL(glob_ifreq-1,KIND=DP) &
        & +ecut_spectralf(1)
        CALL calc_corr_gamma( sc(:,:,1), en(:,:,1), .FALSE.)
        DO iks=1,k_grid%nps
           DO ib = qp_bandrange(1), qp_bandrange(2)
              sigma_spectralf(glob_ifreq,ib,iks) = sc(ib,iks,1)
           ENDDO
        ENDDO
        CALL update_bar_type(barra,'qplot',1)
     ENDDO
     !
     CALL stop_bar_type(barra,'qplot')
     !
     DEALLOCATE( en )
     DEALLOCATE( sc )
     !
  ENDIF
  !
  DEALLOCATE( sigma_hf )
  !
  CALL wfreq_db_write( )
  !
  CALL stop_clock( 'solve_qp' )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_qp_k(l_secant,l_generate_plot)
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies.
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_lanczos,qp_bandrange,imfreq_list_integrate,&
                                 & n_secant_maxiter,trev_secant,l_enable_lanczos,imfreq_list,n_imfreq,&
                                 & z_epsm1_ifr_q,z_epsm1_rfr_q,n_spectralf,ecut_spectralf,&
                                 & z_body1_ifr_q,z_body2_ifr_q,d_diago_q,z_body_rfr_q,sigma_z,&
                                 & sigma_eqplin,sigma_eqpsec,sigma_sc_eks,sigma_sc_eqplin,&
                                 & sigma_sc_eqpsec,sigma_diff,sigma_spectralf,sigma_freq
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm,inter_bgrp_comm
  USE mp_world,             ONLY : mpime,root
  USE mp,                   ONLY : mp_sum
  USE pwcom,                ONLY : et,nbnd
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE constants,            ONLY : rytoev,pi
  USE west_io,              ONLY : serial_table_output
  USE distribution_center,  ONLY : pert,ifr,rfr,aband,band_group,kpt_pool
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wfreq_io,             ONLY : readin_overlap,readin_solvegfreq,readin_solvehf
  USE wfreq_db,             ONLY : wfreq_db_write
  USE types_bz_grid,        ONLY : k_grid,q_grid
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_secant,l_generate_plot
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: sigma_hf(:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_cor_in (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_cor_out(:,:)
  REAL(DP),ALLOCATABLE :: z_in(:,:)
  REAL(DP),ALLOCATABLE :: qp_energy(:,:)
  COMPLEX(DP),ALLOCATABLE :: sc(:,:,:)
  REAL(DP),ALLOCATABLE :: en(:,:,:)
  LOGICAL,ALLOCATABLE :: l_conv(:,:)
  REAL(DP),PARAMETER :: eshift = 0.007349862_DP ! = 0.1 eV
  INTEGER :: ib,ibloc,iks,ik,ikks,ikk,iq,ifixed,ip,glob_ip,ifreq,il,im,glob_im,glob_jp,glob_ifreq
  INTEGER :: is,iss
  REAL(DP) :: g0(3)
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  CHARACTER(LEN=5) :: myglobk
  INTEGER :: notconv
  REAL(DP),ALLOCATABLE :: dtemp(:)
  COMPLEX(DP),ALLOCATABLE :: ztemp2(:,:)
  REAL(DP),ALLOCATABLE :: diago(:,:)
  COMPLEX(DP),ALLOCATABLE :: braket(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  !
  CALL start_clock('solve_qp')
  !
  ALLOCATE( imfreq_list_integrate( 2, ifr%nloc ) )
  ALLOCATE( dtemp( n_imfreq ) )
  !
  dtemp = 0._DP
  DO ifreq = 1, ifr%nloc
     glob_ifreq = ifr%l2g(ifreq)
     dtemp( glob_ifreq ) = imfreq_list( ifreq )
  ENDDO
  CALL mp_sum( dtemp, intra_bgrp_comm )
  !
  DO ifreq = 1, ifr%nloc
     glob_ifreq = ifr%l2g(ifreq)
     IF( glob_ifreq == 1 ) THEN
        imfreq_list_integrate(1,ifreq) = 0._DP
     ELSE
        imfreq_list_integrate(1,ifreq) = ( dtemp(glob_ifreq) + dtemp(glob_ifreq-1) ) * 0.5_DP
     ENDIF
     !
     IF( glob_ifreq == n_imfreq ) THEN
        imfreq_list_integrate(2,ifreq) = dtemp(n_imfreq)
     ELSE
        imfreq_list_integrate(2,ifreq) = ( dtemp(glob_ifreq) + dtemp(glob_ifreq+1) ) * 0.5_DP
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE( dtemp )
  !
  ! TEMP
  !
  ALLOCATE( overlap( pert%nglob, nbnd ) )
  ALLOCATE( z_body1_ifr_q( aband%nloc, ifr%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps, q_grid%nps ) )
  ALLOCATE( z_body_rfr_q( aband%nloc, rfr%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps, q_grid%nps ) )
  ALLOCATE( ztemp2( nbnd, MAX(ifr%nloc,rfr%nloc) ) )
  IF( l_enable_lanczos ) THEN
     ALLOCATE( z_body2_ifr_q( n_lanczos, pert%nloc, ifr%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps, q_grid%nps ) )
     ALLOCATE( d_diago_q( n_lanczos, pert%nloc, qp_bandrange(1):qp_bandrange(2), k_grid%nps, q_grid%nps ) )
     ALLOCATE( braket( pert%nglob, n_lanczos, pert%nloc ) )
     ALLOCATE( diago( n_lanczos, pert%nloc ) )
  ENDIF
  !
  z_body1_ifr_q = 0._DP
  z_body_rfr_q = 0._DP
  IF( l_enable_lanczos ) THEN
     z_body2_ifr_q = 0._DP
     d_diago_q = 0._DP
  ENDIF
  !
  ! d_body1_ifr_q, d_body2_ifr_q, z_diago_rfr_q, d_diago_q
  !
  CALL io_push_title('Collecting results from W and G')
  !
  CALL band_group%init(qp_bandrange(2)-qp_bandrange(1)+1,'b','band_group',.FALSE.)
  !
  barra_load = k_grid%nps * band_group%nloc * q_grid%nps
  CALL start_bar_type( barra, 'coll_gw', barra_load )
  !
  ! LOOP
  !
  ! ... Outer k-point loop (wfc matrix element): iks
  ! ... Inner k-point loop (wfc summed over k'): ikks
  ! ... BEWARE: iks and ikks are switched w.r.t. solve_gfreq_k
  !
  DO iks = 1, k_grid%nps   ! KPOINT-SPIN (MATRIX ELEMENT)
     !
     ik = k_grid%ip(iks)
     is = k_grid%is(iks)
     !
     DO ibloc = 1, band_group%nloc
        ib = band_group%l2g(ibloc)+qp_bandrange(1)-1
        !
        DO ikks = 1, k_grid%nps   ! KPOINT-SPIN (INTEGRAL OVER K')
           !
           ikk = k_grid%ip(ikks)
           iss = k_grid%is(ikks)
           IF( is /= iss ) CYCLE
           !
           CALL q_grid%find( k_grid%p_cart(:,ik) - k_grid%p_cart(:,ikk), 'cart', iq, g0 )
           !
           CALL readin_overlap( 'g', kpt_pool%l2g(iks), kpt_pool%l2g(ikks), ib, overlap, pert%nglob, nbnd )
           !
           ! ------
           ! z_body1_ifr_q
           ! ------
           !
           ztemp2 = 0._DP
           !
           DO ifreq = 1, ifr%nloc
              !
              DO im = 1, nbnd
                 !
                 DO glob_jp = 1, n_pdep_eigen_to_use
                    DO ip = 1, pert%nloc
                       glob_ip = pert%l2g(ip)
                       ztemp2(im,ifreq) = ztemp2(im,ifreq)+CONJG(overlap(glob_jp,im))*overlap(glob_ip,im) &
                       & *z_epsm1_ifr_q(glob_jp,ip,ifreq,iq)
                    ENDDO
                 ENDDO
                 !
              ENDDO ! im
              !
           ENDDO ! ifreq
           !
           CALL mp_sum(ztemp2,inter_image_comm)
           !
           DO ifreq = 1, ifr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 z_body1_ifr_q(im,ifreq,ib,iks,iq) = ztemp2(glob_im,ifreq)
              ENDDO
           ENDDO
           !
           ! -----
           ! z_body_rfr_q
           ! -----
           !
           ztemp2 = 0._DP
           !
           DO ifreq = 1, rfr%nloc
              !
              DO im = 1, nbnd
                 !
                 DO glob_jp = 1, n_pdep_eigen_to_use
                    DO ip = 1, pert%nloc
                       glob_ip = pert%l2g(ip)
                       ztemp2(im,ifreq) = ztemp2(im,ifreq)+CONJG(overlap(glob_jp,im))*overlap(glob_ip,im) &
                       & *z_epsm1_rfr_q(glob_jp,ip,ifreq,iq)
                    ENDDO
                 ENDDO
                 !
              ENDDO ! im
              !
           ENDDO ! ifreq
           !
           CALL mp_sum(ztemp2,inter_image_comm)
           !
           DO ifreq = 1, rfr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 z_body_rfr_q(im,ifreq,ib,iks,iq) = ztemp2(glob_im,ifreq)
              ENDDO
           ENDDO
           !
           ! -----------------------------
           ! LANCZOS part : d_diago_q, z_body2_ifr_q
           ! -----------------------------
           !
           IF( l_enable_lanczos ) THEN
              !
              CALL readin_solvegfreq( kpt_pool%l2g(iks), kpt_pool%l2g(ikks), ib, diago, braket, &
                 & pert%nloc, pert%nglob, pert%myoffset )
              !
              DO ip = 1, pert%nloc
                 DO il = 1, n_lanczos
                    d_diago_q(il,ip,ib,iks,iq) = diago(il,ip)
                 ENDDO
              ENDDO
              !
              DO ifreq = 1, ifr%nloc
                 DO ip = 1, pert%nloc
                    DO il = 1, n_lanczos
                       DO glob_jp = 1, pert%nglob
                          z_body2_ifr_q(il,ip,ifreq,ib,iks,iq) = z_body2_ifr_q(il,ip,ifreq,ib,iks,iq) + &
                          & braket(glob_jp,il,ip)*z_epsm1_ifr_q(glob_jp,ip,ifreq,iq)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              !
           ENDIF
           !
           CALL update_bar_type( barra, 'coll_gw', 1 )
           !
        ENDDO ! ibnd
        !
     ENDDO ! ikks
     !
  ENDDO ! iks
  !
  CALL mp_sum(z_body1_ifr_q,inter_bgrp_comm)
  CALL mp_sum(z_body_rfr_q,inter_bgrp_comm)
  !
  IF(l_enable_lanczos) THEN
     CALL mp_sum(z_body2_ifr_q,inter_bgrp_comm)
     CALL mp_sum(d_diago_q,inter_bgrp_comm)
  ENDIF
  !
  CALL stop_bar_type( barra, 'coll_gw' )
  !
  DEALLOCATE( overlap )
  DEALLOCATE( z_epsm1_ifr_q )
  DEALLOCATE( z_epsm1_rfr_q )
  DEALLOCATE( ztemp2 )
  IF( l_enable_lanczos ) THEN
     DEALLOCATE( braket )
     DEALLOCATE( diago )
  ENDIF
  !
  ! Get Sigma_X
  !
  ALLOCATE( sigma_hf (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  CALL readin_solvehf( sigma_hf(qp_bandrange(1),1), qp_bandrange(2)-qp_bandrange(1)+1, k_grid%nps )
  !
  ! For CORR
  !
  IF( l_secant ) THEN
     !
     CALL io_push_title('(Q)uasiparticle energies')
     !
     ALLOCATE( en(qp_bandrange(1):qp_bandrange(2),k_grid%nps,2) )
     ALLOCATE( sc(qp_bandrange(1):qp_bandrange(2),k_grid%nps,2) )
     ALLOCATE( l_conv(qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( sigma_cor_in (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( sigma_cor_out(qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( z_in (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( qp_energy (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     !
     ! 1st step of secant solver : E_KS - 0.5 * eshift
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           en(ib,iks,1) = et(ib,iks) - eshift*0.5_DP
        ENDDO
     ENDDO
     CALL calc_corr_k( sc(:,:,1), en(:,:,1), .TRUE.)
     !
     ! 1st step of secant solver : E_KS + 0.5 * eshift
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           en(ib,iks,2) = et(ib,iks) + eshift*0.5_DP
        ENDDO
     ENDDO
     CALL calc_corr_k( sc(:,:,2), en(:,:,2), .TRUE.)
     !
     ! Stage sigma_corr_in
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           sigma_cor_in(ib,iks) = ( sc(ib,iks,2) + sc(ib,iks,1) ) * 0.5_DP
        ENDDO
     ENDDO
     !
     ! Stage z_in
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           z_in(ib,iks) = 1._DP / ( 1._DP-REAL( sc(ib,iks,2)-sc(ib,iks,1), KIND=DP ) / eshift  )
        ENDDO
     ENDDO
     !
     ! en 1 = EKS, sc 1 = sigma_corr_in
     ! en 2 = EKS + Z * ( sigmax - vxc + sigmacorrin)
     !
     DO iks = 1, k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           en(ib,iks,1) = et(ib,iks)
           sc(ib,iks,1) = sigma_cor_in(ib,iks)
           en(ib,iks,2) = et(ib,iks) + z_in(ib,iks) * ( sigma_hf(ib,iks) + REAL( sigma_cor_in(ib,iks),KIND=DP) )
        ENDDO
     ENDDO
     sigma_z           = z_in
     sigma_eqplin(:,:) = en(:,:,2)
     sigma_sc_eks      = sigma_cor_in
     CALL output_eqp_report(0,en(:,:,1),en(:,:,2),sigma_cor_in(:,:))
     !
     ! nth step of the secant solver
     !
     l_conv = .FALSE.
     notconv =  k_grid%nps * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
     DO ifixed = 1, n_secant_maxiter
        !
        CALL calc_corr_k( sc(:,:,2), en(:,:,2), .TRUE.)
        !
        ! Resulting new energy:
        !
        DO iks = 1, k_grid%nps
           DO ib = qp_bandrange(1), qp_bandrange(2)
               IF( .NOT. l_conv(ib,iks) ) THEN
                  qp_energy(ib,iks) = en(ib,iks,2) + &
                         & ( et(ib,iks) + sigma_hf(ib,iks) + REAL(sc(ib,iks,2),KIND=DP) - en(ib,iks,2) ) / &
                         & ( 1._DP - REAL( sc(ib,iks,2) - sc(ib,iks,1), KIND=DP ) / ( en(ib,iks,2) - en(ib,iks,1) ) )
               ENDIF
           ENDDO
        ENDDO
        !
        ! Estimate l_conv
        !
        DO iks = 1, k_grid%nps
           DO ib = qp_bandrange(1), qp_bandrange(2)
              l_conv(ib,iks) = ( ABS( qp_energy(ib,iks) - en(ib,iks,2) ) < trev_secant )
           ENDDO
        ENDDO
        !
        ! Count the number of notconverged QP energies
        !
        notconv = 0
        DO iks = 1, k_grid%nps
           DO ib = qp_bandrange(1), qp_bandrange(2)
              IF( .NOT.l_conv(ib,iks) ) notconv = notconv + 1
           ENDDO
        ENDDO
        !
        IF( ifixed == 1 ) THEN
           sigma_sc_eqplin(:,:) = sc(:,:,2)
        ENDIF
        !
        sigma_eqpsec = qp_energy
        sigma_sc_eqpsec = sc(:,:,2)
        sigma_diff(:,:) = qp_energy(:,:) - en(:,:,2)
        !
        IF( notconv == 0 ) THEN
           !
           CALL io_push_title('CONVERGENCE ACHIEVED !!!')
           CALL output_eqp_report(-1,en(:,:,2),qp_energy,sc(:,:,2))
           EXIT
        ELSE
           CALL output_eqp_report(ifixed,en(:,:,2),qp_energy,sc(:,:,2))
        ENDIF
        !
        ! Iterate
        !
        en(:,:,1) = en(:,:,2)
        sc(:,:,1) = sc(:,:,2)
        en(:,:,2) = qp_energy(:,:)
        !
     ENDDO
     !
     sigma_cor_out(:,:) = sc(:,:,2)
     DEALLOCATE( en, sc, l_conv )
     !
     IF( notconv /= 0 ) THEN
        CALL io_push_title('CONVERGENCE **NOT** ACHIEVED !!!')
     ENDIF
     !
     ! Output it per k-point
     !
     ALLOCATE(out_tab(qp_bandrange(2)-qp_bandrange(1)+1,7))
     !
     DO iks=1,k_grid%nps
        DO ib = qp_bandrange(1), qp_bandrange(2)
           out_tab( ib - qp_bandrange(1) + 1, 1) = REAL( ib, KIND=DP)
           out_tab( ib - qp_bandrange(1) + 1, 2) = et(ib,iks) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 3) = (et(ib,iks)+sigma_hf(ib,iks)) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 4) = qp_energy(ib,iks) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 5) = (qp_energy(ib,iks) - et(ib,iks) ) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 6) = REAL(  sigma_cor_out(ib,iks), KIND=DP ) * rytoev
           out_tab( ib - qp_bandrange(1) + 1, 7) = AIMAG( sigma_cor_out(ib,iks) ) * rytoev
        ENDDO
        WRITE(myglobk,'(i5.5)') iks
        !
        CALL serial_table_output(mpime==root,'eqp_K'//myglobk,out_tab,&
           & qp_bandrange(2)-qp_bandrange(1)+1,7,&
           & (/'      band','    E0[eV]','   EHF[eV]','   Eqp[eV]','Eqp-E0[eV]','Sc_Eqp[eV]',' Width[eV]'/))
     ENDDO
     !
     DEALLOCATE( out_tab )
     DEALLOCATE( sigma_cor_in )
     DEALLOCATE( sigma_cor_out )
     DEALLOCATE( z_in )
     DEALLOCATE( qp_energy )
     !
  ENDIF
  !
  IF( l_generate_plot ) THEN
     !
     CALL io_push_title('(P)lotting the QP corrections...')
     !
     CALL start_bar_type( barra, 'qplot', n_spectralf )
     !
     ALLOCATE( en(qp_bandrange(1):qp_bandrange(2),k_grid%nps,1) )
     ALLOCATE( sc(qp_bandrange(1):qp_bandrange(2),k_grid%nps,1) )
     !
     DO glob_ifreq = 1, n_spectralf
        sigma_freq(glob_ifreq) = (ecut_spectralf(2)-ecut_spectralf(1)) / REAL(n_spectralf-1,KIND=DP) * REAL(glob_ifreq-1,KIND=DP) &
                               & + ecut_spectralf(1)
     ENDDO
     !
     DO glob_ifreq = 1, n_spectralf
        en = (ecut_spectralf(2)-ecut_spectralf(1)) / REAL(n_spectralf-1,KIND=DP) * REAL(glob_ifreq-1,KIND=DP) + ecut_spectralf(1)
        CALL calc_corr_k( sc(:,:,1), en(:,:,1), .FALSE.)
        DO iks=1, k_grid%nps
           DO ib = qp_bandrange(1), qp_bandrange(2)
              sigma_spectralf(glob_ifreq,ib,iks) = sc(ib,iks,1)
           ENDDO
        ENDDO
        CALL update_bar_type(barra,'qplot',1)
     ENDDO
     !
     CALL stop_bar_type(barra,'qplot')
     !
     DEALLOCATE( en )
     DEALLOCATE( sc )
     !
  ENDIF
  !
  DEALLOCATE( sigma_hf )
  !
  CALL wfreq_db_write( )
  !
  CALL stop_clock( 'solve_qp' )
  !
END SUBROUTINE
!
SUBROUTINE output_eqp_report(iteration,en1,en2,sc1)
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : qp_bandrange,trev_secant,logfile
  USE constants,            ONLY : rytoev
  USE mp_world,             ONLY : mpime,root
  USE io_global,            ONLY : stdout
  USE io_push,              ONLY : io_push_bar
  USE json_module,          ONLY : json_file
  USE types_bz_grid,        ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: iteration
  REAL(DP),INTENT(IN) :: en1(qp_bandrange(1):qp_bandrange(2),k_grid%nps)
  REAL(DP),INTENT(IN) :: en2(qp_bandrange(1):qp_bandrange(2),k_grid%nps)
  COMPLEX(DP),INTENT(IN) :: sc1(qp_bandrange(1):qp_bandrange(2),k_grid%nps)
  !
  ! Workspace
  !
  INTEGER :: ib, iks, ik, is
  CHARACTER(LEN=4) :: symb(2)
  TYPE(json_file) :: json
  INTEGER :: secitr, iunit
  LOGICAL :: found
  LOGICAL :: lnospin
  CHARACTER(LEN=10) :: ccounter, csecitr
  INTEGER :: counter
  !
  ! STDOUT
  !
  lnospin = ( k_grid%nps == k_grid%np )
  WRITE(stdout,'(5X)')
  CALL io_push_bar()
  IF( iteration >= 0 ) WRITE(stdout,"(5X,'Iter: ',i6.6)") iteration
  DO ik = 1, k_grid%np
     WRITE( stdout, '(5x,"k(",i6.6,") = (",3f12.7,") cryst. coord.")') ik, k_grid%p_cryst(1:3,ik)
     IF( lnospin ) THEN
        WRITE(stdout,'(5X,a,1X,a,1X,a)') 'band  ', '   QP en. [eV]', 'conv'
        iks = k_grid%ipis2ips(ik,1)
        DO ib = qp_bandrange(1), qp_bandrange(2)
           symb(1)='  no'
           IF( (iteration /= 0) .AND. (ABS(en2(ib,iks)-en1(ib,iks)) < trev_secant) ) symb(1)=' yes'
           WRITE(stdout,'(5X,i6.6,1X,1f14.6,1X,a)') ib, en2(ib,iks)*rytoev, symb(1)
        ENDDO
     ELSE
        WRITE(stdout,'(5X,a,1X,a,1X,a,1X,a,1X,a)') 'band  ', '   QP en. [eV]', 'conv', '   QP en. [eV]', 'conv'
        DO ib = qp_bandrange(1), qp_bandrange(2)
           symb(1:2)='  no'
           IF( (iteration /= 0) .AND. (ABS(en2(ib,k_grid%ipis2ips(ik,1))-en1(ib,k_grid%ipis2ips(ik,1))) &
             & < trev_secant) ) symb(1)=' yes'
           IF( (iteration /= 0) .AND. (ABS(en2(ib,k_grid%ipis2ips(ik,2))-en1(ib,k_grid%ipis2ips(ik,2))) &
             & < trev_secant) ) symb(2)=' yes'
           WRITE(stdout,'(5X,i6.6,1X,1f14.6,1X,a,1f14.6,1X,a)') ib, en2(ib,k_grid%ipis2ips(ik,1))*rytoev, symb(1), &
           & en2(ib,k_grid%ipis2ips(ik,2))*rytoev, symb(2)
        ENDDO
     ENDIF
     IF (k_grid%np>1.AND.ik<k_grid%np)  WRITE(stdout,'(5X, 33(a))') '---------------------------------'
  ENDDO
  CALL io_push_bar()
  !
  ! LOGFILE
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     !
     CALL json%load(filename=TRIM(logfile))
     !
     CALL json%get('exec.Q.secitr',secitr, found )
     !
     IF( found ) THEN
        secitr = secitr+1
     ELSE
        secitr = 1
     ENDIF
     !
     WRITE(csecitr,'(i10)') secitr
     !
     CALL json%update('exec.Q.secitr', secitr, found )
     !
     counter = 0
     DO iks = 1, k_grid%nps
        ik = k_grid%ip(iks)
        is = k_grid%is(iks)
        DO ib = qp_bandrange(1), qp_bandrange(2)
           counter = counter + 1
           WRITE( ccounter, '(i10)') counter
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').ksb',(/ik,is,ib/))
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').ein('//TRIM(ADJUSTL(csecitr))//')',en1(ib,iks)*rytoev)
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').eout('//TRIM(ADJUSTL(csecitr))//')',en2(ib,iks)*rytoev)
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').sc_ein('//TRIM(ADJUSTL(csecitr))//')',&
              &(/REAL(sc1(ib,iks)*rytoev,KIND=DP),AIMAG(sc1(ib,iks)*rytoev)/))
        ENDDO
     ENDDO
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
END SUBROUTINE
