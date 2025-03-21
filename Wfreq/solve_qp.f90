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
SUBROUTINE solve_qp(l_secant,l_generate_plot,l_QDET)
  !-----------------------------------------------------------------------
  !
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_secant,l_generate_plot
  LOGICAL,INTENT(IN) :: l_QDET ! True if QDET double-counting term is calculated.
  !
  IF( gamma_only ) THEN
     CALL solve_qp_gamma( l_secant, l_generate_plot, l_QDET )
  ELSE
     CALL solve_qp_k( l_secant, l_generate_plot )
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_qp_gamma(l_secant,l_generate_plot,l_QDET)
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies.
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_pdep_eigen_to_use,qp_bands,n_bands,imfreq_list_integrate,&
                                 & n_secant_maxiter,trev_secant,imfreq_list,n_imfreq,d_epsm1_ifr,&
                                 & z_epsm1_rfr,n_spectralf,ecut_spectralf,d_body1_ifr,z_body_rfr,&
                                 & sigma_z,sigma_eqplin,sigma_eqpsec,sigma_sc_eks,sigma_sc_eqplin,&
                                 & sigma_sc_eqpsec,sigma_diff,sigma_spectralf,sigma_freq,&
                                 & l_enable_off_diagonal,ijpmap,d_body1_ifr_full,z_body_rfr_full,&
                                 & sigma_sc_eks_full,sigma_sc_eqplin_full,sigma_corr_full,&
                                 & l_dc2025,d_epsm1_ifr_dc,z_epsm1_rfr_dc
  USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_pool_comm,my_pool_id,&
                                 & intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE pwcom,                ONLY : et,nbnd,nspin
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE distribution_center,  ONLY : pert,kpt_pool,band_group,ifr,rfr,aband
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wfreq_io,             ONLY : read_overlap,read_hf
  USE wfreq_db,             ONLY : wfreq_db_write
  USE types_bz_grid,        ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_secant,l_generate_plot,l_QDET
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: sigma_hf(:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_cor_out(:,:)
  REAL(DP),ALLOCATABLE :: qp_energy(:,:)
  COMPLEX(DP),ALLOCATABLE :: sc(:,:,:)
  REAL(DP),ALLOCATABLE :: en(:,:,:)
  LOGICAL,ALLOCATABLE :: l_conv(:,:)
  REAL(DP),PARAMETER :: eshift = 0.007349862_DP ! = 0.1 eV
  INTEGER :: ib,ibloc,ib_index,jb,jb_index,ipair,iloc_pair,nloc_pairs,iks,iks_g,is
  INTEGER :: ifixed,ip,ifreq,im,im_index,glob_im,glob_jp,glob_ifreq
  INTEGER :: notconv
  REAL(DP),ALLOCATABLE :: dtemp(:)
  REAL(DP),ALLOCATABLE :: dtemp2(:,:)
  COMPLEX(DP),ALLOCATABLE :: ztemp2(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dtemp2,ztemp2
#endif
  REAL(DP),ALLOCATABLE :: overlap(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: overlap
#endif
  REAL(DP),ALLOCATABLE :: overlap_loc(:,:)
  REAL(DP),ALLOCATABLE :: d_epsm1_ifr_trans(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: z_epsm1_rfr_trans(:,:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  REAL(DP) :: reduce_r
  COMPLEX(DP) :: reduce_c
  INTEGER :: pert_nloc,ifr_nloc,rfr_nloc
  INTEGER,ALLOCATABLE :: l2g(:)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('solve_qp')
#else
  CALL start_clock('solve_qp')
#endif
  !
  CALL io_push_title('Collecting results from W and G')
  !
  CALL band_group%init(n_bands,'b','band_group',.FALSE.)
  !
  IF( .NOT. l_QDET ) THEN
     ALLOCATE( imfreq_list_integrate( 2, ifr%nloc ) )
  ENDIF
  ALLOCATE( dtemp( n_imfreq ) )
  !
  dtemp(:) = 0._DP
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
  ALLOCATE( overlap_loc(pert%nloc,nbnd ) )
  ALLOCATE( dtemp2( nbnd, ifr%nloc ) )
  ALLOCATE( ztemp2( nbnd, rfr%nloc ) )
  !$acc enter data create(overlap,overlap_loc,dtemp2,ztemp2) copyin(qp_bands)
  ALLOCATE( d_epsm1_ifr_trans( pert%nloc, pert%nglob, ifr%nloc ) )
  ALLOCATE( z_epsm1_rfr_trans( pert%nloc, pert%nglob, rfr%nloc ) )
  IF (l_QDET .AND. l_dc2025) THEN
     d_epsm1_ifr_trans(:,:,:) = RESHAPE( d_epsm1_ifr_dc, [pert%nloc,pert%nglob,ifr%nloc], ORDER=[2,1,3] )
     z_epsm1_rfr_trans(:,:,:) = RESHAPE( z_epsm1_rfr_dc, [pert%nloc,pert%nglob,rfr%nloc], ORDER=[2,1,3] )
  ELSE
     d_epsm1_ifr_trans(:,:,:) = RESHAPE( d_epsm1_ifr, [pert%nloc,pert%nglob,ifr%nloc], ORDER=[2,1,3] )
     z_epsm1_rfr_trans(:,:,:) = RESHAPE( z_epsm1_rfr, [pert%nloc,pert%nglob,rfr%nloc], ORDER=[2,1,3] )
  ENDIF
  !$acc enter data copyin(d_epsm1_ifr_trans,z_epsm1_rfr_trans)
  IF( l_enable_off_diagonal ) THEN
     !
     nloc_pairs = 0
     DO is = 1, nspin
        im = 0
        DO ibloc = 1, band_group%nloc
           ib_index = band_group%l2g(ibloc)
           ib = qp_bands(ib_index,is)
           DO jb_index = 1, n_bands
              jb = qp_bands(jb_index,is)
              IF(jb <= ib) im = im+1
           ENDDO
        ENDDO
        nloc_pairs = MAX(nloc_pairs,im)
     ENDDO
     !
     IF( .NOT. l_QDET ) THEN
        ALLOCATE( d_body1_ifr_full( aband%nloc, ifr%nloc, nloc_pairs, k_grid%nps ) )
        ALLOCATE( z_body_rfr_full( aband%nloc, rfr%nloc, nloc_pairs, k_grid%nps ) )
     ENDIF
     !
     d_body1_ifr_full(:,:,:,:) = 0._DP
     z_body_rfr_full(:,:,:,:) = 0._DP
     !
  ELSE
     !
     IF( .NOT. l_QDET ) THEN
        ALLOCATE( d_body1_ifr( aband%nloc, ifr%nloc, band_group%nloc, k_grid%nps ) )
        ALLOCATE( z_body_rfr( aband%nloc, rfr%nloc, band_group%nloc, k_grid%nps ) )
     ENDIF
     !
     d_body1_ifr(:,:,:,:) = 0._DP
     z_body_rfr(:,:,:,:) = 0._DP
     !
  ENDIF
  !
  pert_nloc = pert%nloc
  ifr_nloc = ifr%nloc
  rfr_nloc = rfr%nloc
  !
  ALLOCATE( l2g( pert%nloc ) )
  !$acc enter data create(l2g)
  !
  !$acc parallel loop present(l2g)
  DO ip = 1, pert_nloc
     !
     ! l2g(ip) = pert%l2g(ip)
     !
     l2g(ip) = nimage*(ip-1)+my_image_id+1
  ENDDO
  !$acc end parallel
  !
  ! d_body1_ifr, d_body2_ifr, z_diago_rfr, d_diago
  !
  barra_load = 0
  DO ibloc = 1,band_group%nloc
     ib_index = band_group%l2g(ibloc)
     !
     IF(l_enable_off_diagonal) THEN
        barra_load = barra_load+ib_index
     ELSE
        barra_load = barra_load+1
     ENDIF
  ENDDO
  barra_load = barra_load*kpt_pool%nloc
  !
  CALL start_bar_type( barra, 'coll_gw', barra_load )
  !
  ! LOOP
  !
  DO iks = 1,kpt_pool%nloc ! KPOINT-SPIN
     !
     iks_g = kpt_pool%l2g(iks)
     is = k_grid%is(iks_g)
     iloc_pair = 0
     !
     DO ibloc = 1, band_group%nloc
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index,is)
        !
        CALL read_overlap( 'g', kpt_pool%l2g(iks), ib, overlap, pert%nglob, nbnd )
        !
        !$acc update device(overlap)
        !
        !$acc parallel loop collapse(2) present(overlap_loc,overlap,l2g)
        DO im = 1, nbnd
           DO ip = 1, pert_nloc
              overlap_loc(ip,im) = overlap(l2g(ip),im)
           ENDDO
        ENDDO
        !$acc end parallel
        !
        DO jb_index = 1, n_bands
           !
           jb = qp_bands(jb_index,is)
           !
           IF(l_enable_off_diagonal .AND. jb <= ib) THEN
              !
              ipair = ijpmap(jb_index,ib_index)
              iloc_pair = iloc_pair+1
              !
              CALL read_overlap( 'g', kpt_pool%l2g(iks), jb, overlap, pert%nglob, nbnd )
              !
              !$acc update device(overlap)
              !
           ELSEIF(l_enable_off_diagonal .OR. jb /= ib) THEN
              CYCLE
           ENDIF
           !
           ! ------
           ! d_body1_ifr
           ! ------
           !
           IF((l_enable_off_diagonal .AND. jb <= ib) &
           & .OR. (.NOT. l_enable_off_diagonal .AND. jb == ib)) THEN
              IF(l_QDET) THEN
                 !
                 ! For the QDET double-counting term, all states need to be within qp_bands
                 !
                 !$acc parallel present(qp_bands,overlap,overlap_loc,d_epsm1_ifr_trans,dtemp2)
                 !$acc loop collapse(2)
                 DO ifreq = 1, ifr_nloc
                    DO im_index = 1, n_bands
                       im = qp_bands(im_index,is)
                       reduce_r = 0._DP
                       !$acc loop collapse(2) reduction(+:reduce_r)
                       DO glob_jp = 1, n_pdep_eigen_to_use
                          DO ip = 1, pert_nloc
                             reduce_r = reduce_r+overlap(glob_jp,im)*overlap_loc(ip,im) &
                             & *d_epsm1_ifr_trans(ip,glob_jp,ifreq)
                          ENDDO
                       ENDDO
                       dtemp2(im,ifreq) = reduce_r
                    ENDDO ! im_index
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ELSE
                 !
                 !$acc parallel present(overlap,overlap_loc,d_epsm1_ifr_trans,dtemp2)
                 !$acc loop collapse(2)
                 DO ifreq = 1, ifr_nloc
                    DO im = 1, nbnd
                       reduce_r = 0._DP
                       !$acc loop collapse(2) reduction(+:reduce_r)
                       DO glob_jp = 1, n_pdep_eigen_to_use
                          DO ip = 1, pert_nloc
                             reduce_r = reduce_r+overlap(glob_jp,im)*overlap_loc(ip,im) &
                             & *d_epsm1_ifr_trans(ip,glob_jp,ifreq)
                          ENDDO
                       ENDDO
                       dtemp2(im,ifreq) = reduce_r
                    ENDDO ! im
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ENDIF
           ENDIF
           !
           !$acc update host(dtemp2)
           !
           CALL mp_sum(dtemp2,inter_image_comm)
           !
           DO ifreq = 1, ifr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 IF(l_enable_off_diagonal .AND. jb <= ib) THEN
                    d_body1_ifr_full(im,ifreq,iloc_pair,iks_g) = dtemp2(glob_im,ifreq)
                 ELSEIF(.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                    d_body1_ifr(im,ifreq,ibloc,iks_g) = dtemp2(glob_im,ifreq)
                 ENDIF
              ENDDO
           ENDDO
           !
           ! -----
           ! z_body_rfr
           ! -----
           !
           IF((l_enable_off_diagonal .AND. jb <= ib) &
           & .OR. (.NOT. l_enable_off_diagonal .AND. jb == ib)) THEN
              IF(l_QDET) THEN
                 !
                 ! For the QDET double-counting term, all states need to be within qp_bands
                 !
                 !$acc parallel present(qp_bands,overlap,overlap_loc,z_epsm1_rfr_trans,ztemp2)
                 !$acc loop collapse(2)
                 DO ifreq = 1, rfr_nloc
                    DO im_index = 1, n_bands
                       im = qp_bands(im_index,is)
                       reduce_c = 0._DP
                       !$acc loop collapse(2) reduction(+:reduce_c)
                       DO glob_jp = 1, n_pdep_eigen_to_use
                          DO ip = 1, pert_nloc
                             reduce_c = reduce_c+overlap(glob_jp,im)*overlap_loc(ip,im) &
                             & *z_epsm1_rfr_trans(ip,glob_jp,ifreq)
                          ENDDO
                       ENDDO
                       ztemp2(im,ifreq) = reduce_c
                    ENDDO ! im_index
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ELSE
                 !
                 !$acc parallel present(overlap,overlap_loc,z_epsm1_rfr_trans,ztemp2)
                 !$acc loop collapse(2)
                 DO ifreq = 1, rfr_nloc
                    DO im = 1, nbnd
                       reduce_c = 0._DP
                       !$acc loop collapse(2) reduction(+:reduce_c)
                       DO glob_jp = 1, n_pdep_eigen_to_use
                          DO ip = 1, pert_nloc
                             reduce_c = reduce_c+overlap(glob_jp,im)*overlap_loc(ip,im) &
                             & *z_epsm1_rfr_trans(ip,glob_jp,ifreq)
                          ENDDO
                       ENDDO
                       ztemp2(im,ifreq) = reduce_c
                    ENDDO ! im
                 ENDDO ! ifreq
                 !$acc end parallel
                 !
              ENDIF
           ENDIF
           !
           !$acc update host(ztemp2)
           !
           CALL mp_sum(ztemp2,inter_image_comm)
           !
           DO ifreq = 1, rfr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 IF(l_enable_off_diagonal .AND. jb <= ib) THEN
                    z_body_rfr_full(im,ifreq,iloc_pair,iks_g) = ztemp2(glob_im,ifreq)
                 ELSEIF(.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                    z_body_rfr(im,ifreq,ibloc,iks_g) = ztemp2(glob_im,ifreq)
                 ENDIF
              ENDDO
           ENDDO
           !
           ! -----------------------------
           ! LANCZOS part : d_diago, d_body2_ifr (computed in solve_gfreq)
           ! -----------------------------
           !
           CALL update_bar_type( barra, 'coll_gw', 1 )
           !
        ENDDO ! jbnd
        !
     ENDDO ! ibnd
     !
  ENDDO ! iks
  !
  IF(l_enable_off_diagonal) THEN
     CALL mp_sum(d_body1_ifr_full,inter_pool_comm)
     CALL mp_sum(z_body_rfr_full,inter_pool_comm)
  ELSE
     CALL mp_sum(d_body1_ifr,inter_pool_comm)
     CALL mp_sum(z_body_rfr,inter_pool_comm)
  ENDIF
  !
  CALL stop_bar_type( barra, 'coll_gw' )
  !
  !$acc exit data delete(l2g,overlap,overlap_loc,dtemp2,ztemp2,qp_bands,d_epsm1_ifr_trans,z_epsm1_rfr_trans)
  DEALLOCATE( l2g )
  DEALLOCATE( overlap )
  DEALLOCATE( overlap_loc )
  DEALLOCATE( dtemp2 )
  DEALLOCATE( ztemp2 )
  DEALLOCATE( d_epsm1_ifr_trans )
  DEALLOCATE( z_epsm1_rfr_trans )
  !
  ! Get Sigma_X
  !
  ALLOCATE( sigma_hf (n_bands,k_grid%nps) )
  CALL read_hf( sigma_hf, n_bands )
  !
  ! For CORR
  !
  IF( l_secant ) THEN
     !
     CALL io_push_title('(Q)uasiparticle energies')
     !
     ALLOCATE( en(n_bands,k_grid%nps,2) )
     ALLOCATE( sc(n_bands,k_grid%nps,2) )
     ALLOCATE( l_conv(n_bands,k_grid%nps) )
     ALLOCATE( qp_energy (n_bands,k_grid%nps) )
     !
     ! 1st step of secant solver
     ! en 1 = E_KS - 0.5 * eshift
     ! en 2 = E_KS + 0.5 * eshift
     !
     IF( my_pool_id == 0 ) THEN
        DO iks = 1, k_grid%nps
           is = k_grid%is(iks)
           DO ib = 1, n_bands
              en(ib,iks,1) = et(qp_bands(ib,is),iks) - eshift*0.5_DP
              en(ib,iks,2) = et(qp_bands(ib,is),iks) + eshift*0.5_DP
           ENDDO
           !
        ENDDO
     ENDIF
     !
     CALL mp_bcast(en,0,inter_pool_comm)
     !
     IF(l_enable_off_diagonal) THEN
        CALL calc_corr_gamma( sc(:,:,1), en(:,:,1), .TRUE., .TRUE., .FALSE.)
        sigma_sc_eks_full(:,:) = sigma_corr_full
        CALL calc_corr_gamma( sc(:,:,2), en(:,:,2), .TRUE., .TRUE., .FALSE.)
        sigma_sc_eks_full(:,:) = ( sigma_sc_eks_full + sigma_corr_full ) * 0.5_DP
     ELSE
        CALL calc_corr_gamma( sc(:,:,1), en(:,:,1), .TRUE., .FALSE., .FALSE.)
        CALL calc_corr_gamma( sc(:,:,2), en(:,:,2), .TRUE., .FALSE., .FALSE.)
     ENDIF
     !
     ! Stage sigma_sc_eks
     !
     DO iks = 1, k_grid%nps
        DO ib = 1, n_bands
           sigma_sc_eks(ib,iks) = ( sc(ib,iks,2) + sc(ib,iks,1) ) * 0.5_DP
        ENDDO
     ENDDO
     !
     ! Stage sigma_z
     !
     DO iks = 1, k_grid%nps
        DO ib = 1, n_bands
           sigma_z(ib,iks) = 1._DP / ( 1._DP - REAL(sc(ib,iks,2)-sc(ib,iks,1),KIND=DP ) / eshift )
        ENDDO
     ENDDO
     !
     ! en 1 = EKS
     ! sc 1 = sigma_sc_eks
     ! en 2 = EKS + Z * ( sigmax - vxc + sigma_sc_eks)
     !
     IF( my_pool_id == 0 ) THEN
        DO iks = 1, k_grid%nps
           is = k_grid%is(iks)
           DO ib = 1, n_bands
              en(ib,iks,1) = et(qp_bands(ib,is),iks)
              sc(ib,iks,1) = sigma_sc_eks(ib,iks)
              en(ib,iks,2) = et(qp_bands(ib,is),iks) + sigma_z(ib,iks) &
              & * ( sigma_hf(ib,iks) + REAL(sigma_sc_eks(ib,iks),KIND=DP) )
           ENDDO
        ENDDO
     ENDIF
     !
     CALL mp_bcast(en,0,inter_pool_comm)
     CALL mp_bcast(sc(:,:,1),0,inter_pool_comm)
     !
     sigma_eqplin(:,:) = en(:,:,2)
     CALL output_eqp_report(0,en(:,:,1),en(:,:,2),sigma_sc_eks(:,:))
     !
     ! nth step of the secant solver
     !
     l_conv(:,:) = .FALSE.
     notconv = k_grid%nps * n_bands
     DO ifixed = 1, n_secant_maxiter
        !
        CALL calc_corr_gamma( sc(:,:,2), en(:,:,2), .TRUE., l_enable_off_diagonal, .FALSE. )
        !
        IF( ifixed == 1 ) THEN
           IF(l_enable_off_diagonal) THEN
              sigma_sc_eqplin_full(:,:) = sigma_corr_full
           ELSE
              sigma_sc_eqplin(:,:) = sc(:,:,2)
           ENDIF
        ENDIF
        !
        IF( my_pool_id == 0 ) THEN
           !
           ! Resulting new energy:
           !
           DO iks = 1, k_grid%nps
              is = k_grid%is(iks)
              DO ib = 1, n_bands
                  IF( .NOT. l_conv(ib,iks) ) THEN
                     qp_energy(ib,iks) = en(ib,iks,2) &
                     & + (et(qp_bands(ib,is),iks) + sigma_hf(ib,iks) + REAL(sc(ib,iks,2),KIND=DP) - en(ib,iks,2)) &
                     & / (1._DP - REAL(sc(ib,iks,2) - sc(ib,iks,1),KIND=DP) / (en(ib,iks,2) - en(ib,iks,1)))
                  ENDIF
              ENDDO
           ENDDO
           !
           ! Estimate l_conv
           !
           DO iks = 1, k_grid%nps
              DO ib = 1, n_bands
                 l_conv(ib,iks) = ( ABS( qp_energy(ib,iks) - en(ib,iks,2) ) < trev_secant )
              ENDDO
           ENDDO
           !
           ! Count the number of not converged QP energies
           !
           notconv = 0
           DO iks = 1, k_grid%nps
              DO ib = 1, n_bands
                 IF( .NOT. l_conv(ib,iks) ) notconv = notconv + 1
              ENDDO
           ENDDO
           !
        ENDIF
        !
        CALL mp_bcast(qp_energy,0,inter_pool_comm)
        CALL mp_bcast(notconv,0,inter_pool_comm)
        !
        sigma_eqpsec(:,:) = qp_energy
        IF(.NOT. l_enable_off_diagonal) sigma_sc_eqpsec(:,:) = sc(:,:,2)
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
     DEALLOCATE( en, sc, l_conv )
     !
     IF( notconv /= 0 ) THEN
        CALL io_push_title('CONVERGENCE **NOT** ACHIEVED !!!')
     ENDIF
     !
     DEALLOCATE( qp_energy )
     !
  ENDIF
  !
  IF(l_QDET) THEN
     !
     IF(l_enable_off_diagonal) THEN
        !
        ALLOCATE( sigma_cor_out(n_bands,k_grid%nps) )
        !
        ! sigma_cor_out is not required, the important output is contained in the
        ! global variable sigma_corr_full
        !
        CALL calc_corr_gamma( sigma_cor_out, sigma_eqpsec - sigma_diff, .TRUE., .TRUE., .TRUE.)
        !
        DEALLOCATE( sigma_cor_out )
        !
     ELSE
        !
        CALL calc_corr_gamma( sigma_sc_eqpsec, sigma_eqpsec - sigma_diff, .TRUE., .FALSE., .TRUE.)
        !
     ENDIF
     !
     CALL stop_clock( 'solve_qp' )
     !
     RETURN
     !
  ENDIF
  !
  IF( l_generate_plot ) THEN
     !
     CALL io_push_title('(P)lotting the QP corrections...')
     !
     CALL start_bar_type( barra, 'qplot', n_spectralf )
     !
     ALLOCATE( en(n_bands,k_grid%nps,1) )
     ALLOCATE( sc(n_bands,k_grid%nps,1) )
     !
     DO glob_ifreq = 1, n_spectralf
        sigma_freq(glob_ifreq) = &
        & (ecut_spectralf(2)-ecut_spectralf(1))/REAL(n_spectralf-1,KIND=DP)*REAL(glob_ifreq-1,KIND=DP) &
        & +ecut_spectralf(1)
     ENDDO
     !
     DO glob_ifreq = 1, n_spectralf
        en(:,:,:) = (ecut_spectralf(2)-ecut_spectralf(1))/REAL(n_spectralf-1,KIND=DP)*REAL(glob_ifreq-1,KIND=DP) &
        & +ecut_spectralf(1)
        CALL calc_corr_gamma( sc(:,:,1), en(:,:,1), .FALSE., .FALSE., .FALSE.)
        DO iks=1,k_grid%nps
           DO ib = 1, n_bands
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
#if defined(__CUDA)
  CALL stop_clock_gpu('solve_qp')
#else
  CALL stop_clock( 'solve_qp' )
#endif
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
  USE westcom,              ONLY : n_pdep_eigen_to_use,qp_bands,n_bands,imfreq_list_integrate,&
                                 & n_secant_maxiter,trev_secant,imfreq_list,n_imfreq,z_epsm1_ifr_q,&
                                 & z_epsm1_rfr_q,n_spectralf,ecut_spectralf,z_body1_ifr_q,z_body_rfr_q,&
                                 & sigma_z,sigma_eqplin,sigma_eqpsec,sigma_sc_eks,sigma_sc_eqplin,&
                                 & sigma_sc_eqpsec,sigma_diff,sigma_spectralf,sigma_freq
  USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE pwcom,                ONLY : et,nbnd
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE distribution_center,  ONLY : pert,kpt_pool,band_group,ifr,rfr,aband
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wfreq_io,             ONLY : read_overlap,read_hf
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
  REAL(DP),ALLOCATABLE :: qp_energy(:,:)
  COMPLEX(DP),ALLOCATABLE :: sc(:,:,:)
  REAL(DP),ALLOCATABLE :: en(:,:,:)
  LOGICAL,ALLOCATABLE :: l_conv(:,:)
  REAL(DP),PARAMETER :: eshift = 0.007349862_DP ! = 0.1 eV
  INTEGER :: ib,ibloc,ib_index,iks,ik,ikks,ikk,iq,is,iss
  INTEGER :: ifixed,ip,ifreq,im,glob_im,glob_jp,glob_ifreq
  REAL(DP) :: g0(3)
  INTEGER :: notconv
  REAL(DP),ALLOCATABLE :: dtemp(:)
  COMPLEX(DP),ALLOCATABLE :: ztemp2(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: ztemp2
#endif
  COMPLEX(DP),ALLOCATABLE :: overlap(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: overlap
#endif
  COMPLEX(DP),ALLOCATABLE :: overlap_loc(:,:)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_ifr_trans_q(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_rfr_trans_q(:,:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP) :: reduce
  INTEGER :: pert_nloc,ifr_nloc,rfr_nloc
  INTEGER,ALLOCATABLE :: l2g(:)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('solve_qp')
#else
  CALL start_clock('solve_qp')
#endif
  !
  CALL io_push_title('Collecting results from W and G')
  !
  CALL band_group%init(n_bands,'b','band_group',.FALSE.)
  !
  ALLOCATE( imfreq_list_integrate( 2, ifr%nloc ) )
  ALLOCATE( dtemp( n_imfreq ) )
  !
  dtemp(:) = 0._DP
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
  ALLOCATE( overlap_loc( pert%nloc,nbnd ) )
  ALLOCATE( z_body1_ifr_q( aband%nloc, ifr%nloc, band_group%nloc, k_grid%nps, q_grid%nps ) )
  ALLOCATE( z_body_rfr_q( aband%nloc, rfr%nloc, band_group%nloc, k_grid%nps, q_grid%nps ) )
  ALLOCATE( ztemp2( nbnd, MAX(ifr%nloc,rfr%nloc) ) )
  ALLOCATE( z_epsm1_ifr_trans_q( pert%nloc, pert%nglob, ifr%nloc ) )
  ALLOCATE( z_epsm1_rfr_trans_q( pert%nloc, pert%nglob, rfr%nloc ) )
  !$acc enter data create(overlap,overlap_loc,ztemp2,z_epsm1_ifr_trans_q,z_epsm1_rfr_trans_q)
  !
  pert_nloc = pert%nloc
  ifr_nloc = ifr%nloc
  rfr_nloc = rfr%nloc
  !
  ALLOCATE( l2g( pert%nloc ) )
  !$acc enter data create(l2g)
  !
  !$acc parallel loop present(l2g)
  DO ip = 1, pert_nloc
     !
     ! l2g(ip) = pert%l2g(ip)
     !
     l2g(ip) = nimage*(ip-1)+my_image_id+1
  ENDDO
  !$acc end parallel
  !
  z_body1_ifr_q(:,:,:,:,:) = 0._DP
  z_body_rfr_q(:,:,:,:,:) = 0._DP
  !
  ! d_body1_ifr_q, d_body2_ifr_q, z_diago_rfr_q, d_diago_q
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
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index,is)
        !
        DO ikks = 1, k_grid%nps   ! KPOINT-SPIN (INTEGRAL OVER K')
           !
           ikk = k_grid%ip(ikks)
           iss = k_grid%is(ikks)
           IF( is /= iss ) CYCLE
           !
           CALL q_grid%find( k_grid%p_cart(:,ik) - k_grid%p_cart(:,ikk), 'cart', iq, g0 )
           !
           z_epsm1_ifr_trans_q(:,:,:) = RESHAPE( z_epsm1_ifr_q(:,:,:,iq), &
                                               & [pert%nloc,pert%nglob,ifr%nloc], ORDER=[2,1,3] )
           z_epsm1_rfr_trans_q(:,:,:) = RESHAPE( z_epsm1_rfr_q(:,:,:,iq), &
                                               & [pert%nloc,pert%nglob,rfr%nloc], ORDER=[2,1,3] )
           !
           !$acc update device(z_epsm1_ifr_trans_q,z_epsm1_rfr_trans_q)
           !
           CALL read_overlap( 'g', kpt_pool%l2g(iks), kpt_pool%l2g(ikks), ib, overlap, pert%nglob, nbnd )
           !
           !$acc update device(overlap)
           !
           !$acc parallel loop collapse(2) present(overlap_loc,overlap,l2g)
           DO im = 1, nbnd
              DO ip = 1, pert_nloc
                 overlap_loc(ip,im) = overlap(l2g(ip),im)
              ENDDO
           ENDDO
           !$acc end parallel
           !
           ! ------
           ! z_body1_ifr_q
           ! ------
           !
           !$acc parallel present(overlap,overlap_loc,z_epsm1_ifr_trans_q,ztemp2)
           !$acc loop collapse(2)
           DO ifreq = 1, ifr_nloc
              DO im = 1, nbnd
                 reduce = 0._DP
                 !$acc loop collapse(2) reduction(+:reduce)
                 DO glob_jp = 1, n_pdep_eigen_to_use
                    DO ip = 1, pert_nloc
                       reduce = reduce+CONJG(overlap(glob_jp,im))*overlap_loc(ip,im) &
                       & *z_epsm1_ifr_trans_q(ip,glob_jp,ifreq)
                    ENDDO
                 ENDDO
                 ztemp2(im,ifreq) = reduce
              ENDDO ! im
           ENDDO ! ifreq
           !$acc end parallel
           !
           !$acc update host(ztemp2(:,1:ifr_nloc))
           !
           CALL mp_sum(ztemp2(:,1:ifr%nloc),inter_image_comm)
           !
           DO ifreq = 1, ifr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 z_body1_ifr_q(im,ifreq,ibloc,iks,iq) = ztemp2(glob_im,ifreq)
              ENDDO
           ENDDO
           !
           ! -----
           ! z_body_rfr_q
           ! -----
           !
           !$acc parallel present(overlap,overlap_loc,z_epsm1_rfr_trans_q,ztemp2)
           !$acc loop collapse(2)
           DO ifreq = 1, rfr_nloc
              DO im = 1, nbnd
                 reduce = 0._DP
                 !$acc loop collapse(2) reduction(+:reduce)
                 DO glob_jp = 1, n_pdep_eigen_to_use
                    DO ip = 1, pert_nloc
                       reduce = reduce+CONJG(overlap(glob_jp,im))*overlap_loc(ip,im) &
                       & *z_epsm1_rfr_trans_q(ip,glob_jp,ifreq)
                    ENDDO
                 ENDDO
                 ztemp2(im,ifreq) = reduce
              ENDDO ! im
           ENDDO ! ifreq
           !$acc end parallel
           !
           !$acc update host(ztemp2(:,1:rfr_nloc))
           !
           CALL mp_sum(ztemp2(:,1:rfr%nloc),inter_image_comm)
           !
           DO ifreq = 1, rfr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 z_body_rfr_q(im,ifreq,ibloc,iks,iq) = ztemp2(glob_im,ifreq)
              ENDDO
           ENDDO
           !
           ! -----------------------------
           ! LANCZOS part : d_diago_q, z_body2_ifr_q (computed in solve_gfreq)
           ! -----------------------------
           !
           CALL update_bar_type( barra, 'coll_gw', 1 )
           !
        ENDDO ! ibnd
        !
     ENDDO ! ikks
     !
  ENDDO ! iks
  !
  CALL stop_bar_type( barra, 'coll_gw' )
  !
  !$acc exit data delete(l2g,overlap,overlap_loc,ztemp2,z_epsm1_ifr_trans_q,z_epsm1_rfr_trans_q)
  DEALLOCATE( l2g )
  DEALLOCATE( overlap )
  DEALLOCATE( overlap_loc )
  DEALLOCATE( z_epsm1_ifr_q )
  DEALLOCATE( z_epsm1_rfr_q )
  DEALLOCATE( ztemp2 )
  DEALLOCATE( z_epsm1_ifr_trans_q )
  DEALLOCATE( z_epsm1_rfr_trans_q )
  !
  ! Get Sigma_X
  !
  ALLOCATE( sigma_hf (n_bands,k_grid%nps) )
  CALL read_hf( sigma_hf, n_bands )
  !
  ! For CORR
  !
  IF( l_secant ) THEN
     !
     CALL io_push_title('(Q)uasiparticle energies')
     !
     ALLOCATE( en(n_bands,k_grid%nps,2) )
     ALLOCATE( sc(n_bands,k_grid%nps,2) )
     ALLOCATE( l_conv(n_bands,k_grid%nps) )
     ALLOCATE( qp_energy (n_bands,k_grid%nps) )
     !
     ! 1st step of secant solver
     ! en 1 = E_KS - 0.5 * eshift
     ! en 2 = E_KS + 0.5 * eshift
     !
     DO iks = 1, k_grid%nps
        is = k_grid%is(iks)
        DO ib = 1, n_bands
           en(ib,iks,1) = et(qp_bands(ib,is),iks) - eshift*0.5_DP
           en(ib,iks,2) = et(qp_bands(ib,is),iks) + eshift*0.5_DP
        ENDDO
     ENDDO
     !
     CALL calc_corr_k( sc(:,:,1), en(:,:,1), .TRUE.)
     CALL calc_corr_k( sc(:,:,2), en(:,:,2), .TRUE.)
     !
     ! Stage sigma_sc_eks
     !
     DO iks = 1, k_grid%nps
        DO ib = 1, n_bands
           sigma_sc_eks(ib,iks) = ( sc(ib,iks,2) + sc(ib,iks,1) ) * 0.5_DP
        ENDDO
     ENDDO
     !
     ! Stage sigma_z
     !
     DO iks = 1, k_grid%nps
        DO ib = 1, n_bands
           sigma_z(ib,iks) = 1._DP / ( 1._DP - REAL(sc(ib,iks,2)-sc(ib,iks,1),KIND=DP) / eshift )
        ENDDO
     ENDDO
     !
     ! en 1 = EKS
     ! sc 1 = sigma_sc_eks
     ! en 2 = EKS + Z * ( sigmax - vxc + sigma_sc_eks)
     !
     DO iks = 1, k_grid%nps
        is = k_grid%is(iks)
        DO ib = 1, n_bands
           en(ib,iks,1) = et(qp_bands(ib,is),iks)
           sc(ib,iks,1) = sigma_sc_eks(ib,iks)
           en(ib,iks,2) = et(qp_bands(ib,is),iks) + sigma_z(ib,iks) &
           & * ( sigma_hf(ib,iks) + REAL(sigma_sc_eks(ib,iks),KIND=DP) )
        ENDDO
     ENDDO
     sigma_eqplin(:,:) = en(:,:,2)
     CALL output_eqp_report(0,en(:,:,1),en(:,:,2),sigma_sc_eks(:,:))
     !
     ! nth step of the secant solver
     !
     l_conv(:,:) = .FALSE.
     notconv = k_grid%nps * n_bands
     DO ifixed = 1, n_secant_maxiter
        !
        CALL calc_corr_k( sc(:,:,2), en(:,:,2), .TRUE.)
        !
        IF( ifixed == 1 ) THEN
           sigma_sc_eqplin(:,:) = sc(:,:,2)
        ENDIF
        !
        ! Resulting new energy:
        !
        DO iks = 1, k_grid%nps
           is = k_grid%is(iks)
           DO ib = 1, n_bands
              IF( .NOT. l_conv(ib,iks) ) THEN
                 qp_energy(ib,iks) = en(ib,iks,2) &
                 & + (et(qp_bands(ib,is),iks) + sigma_hf(ib,iks) + REAL(sc(ib,iks,2),KIND=DP) - en(ib,iks,2)) &
                 & / (1._DP - REAL(sc(ib,iks,2) - sc(ib,iks,1), KIND=DP) / (en(ib,iks,2) - en(ib,iks,1)))
              ENDIF
           ENDDO
        ENDDO
        !
        ! Estimate l_conv
        !
        DO iks = 1, k_grid%nps
           DO ib = 1, n_bands
              l_conv(ib,iks) = ( ABS( qp_energy(ib,iks) - en(ib,iks,2) ) < trev_secant )
           ENDDO
        ENDDO
        !
        ! Count the number of not converged QP energies
        !
        notconv = 0
        DO iks = 1, k_grid%nps
           DO ib = 1, n_bands
              IF( .NOT. l_conv(ib,iks) ) notconv = notconv + 1
           ENDDO
        ENDDO
        !
        sigma_eqpsec(:,:) = qp_energy
        sigma_sc_eqpsec(:,:) = sc(:,:,2)
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
     DEALLOCATE( en, sc, l_conv )
     !
     IF( notconv /= 0 ) THEN
        CALL io_push_title('CONVERGENCE **NOT** ACHIEVED !!!')
     ENDIF
     !
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
     ALLOCATE( en(n_bands,k_grid%nps,1) )
     ALLOCATE( sc(n_bands,k_grid%nps,1) )
     !
     DO glob_ifreq = 1, n_spectralf
        sigma_freq(glob_ifreq) = &
        & (ecut_spectralf(2)-ecut_spectralf(1)) / REAL(n_spectralf-1,KIND=DP) * REAL(glob_ifreq-1,KIND=DP) &
        & + ecut_spectralf(1)
     ENDDO
     !
     DO glob_ifreq = 1, n_spectralf
        en(:,:,:) = (ecut_spectralf(2)-ecut_spectralf(1))/REAL(n_spectralf-1,KIND=DP)*REAL(glob_ifreq-1,KIND=DP) &
        & +ecut_spectralf(1)
        CALL calc_corr_k( sc(:,:,1), en(:,:,1), .FALSE.)
        DO iks=1, k_grid%nps
           DO ib = 1, n_bands
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
#if defined(__CUDA)
  CALL stop_clock_gpu('solve_qp')
#else
  CALL stop_clock( 'solve_qp' )
#endif
  !
END SUBROUTINE
!
SUBROUTINE output_eqp_report(iteration,en1,en2,sc1)
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : qp_bands,n_bands,trev_secant,logfile
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
  REAL(DP),INTENT(IN) :: en1(n_bands,k_grid%nps)
  REAL(DP),INTENT(IN) :: en2(n_bands,k_grid%nps)
  COMPLEX(DP),INTENT(IN) :: sc1(n_bands,k_grid%nps)
  !
  ! Workspace
  !
  INTEGER :: ib, iks, ik, is, ib_index
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
  WRITE(stdout,*)
  CALL io_push_bar()
  IF( iteration >= 0 ) WRITE(stdout,'(5X,"Iter: ",I6.6)') iteration
  DO ik = 1, k_grid%np
     WRITE( stdout, '(5X,"k(",I6.6,") = (",3F12.7,") cryst. coord.")') ik, k_grid%p_cryst(1:3,ik)
     IF( lnospin ) THEN
        WRITE(stdout,'(5X,A,1X,A,1X,A)') 'band  ', '   QP en. [eV]', 'conv'
        iks = k_grid%ipis2ips(ik,1)
        DO ib = 1, n_bands
           symb(1)='  no'
           IF( (iteration /= 0) .AND. (ABS(en2(ib,iks)-en1(ib,iks)) < trev_secant) ) symb(1)=' yes'
           WRITE(stdout,'(5X,I6.6,1X,1F14.6,1X,A)') qp_bands(ib,1), en2(ib,iks)*rytoev, symb(1)
        ENDDO
     ELSE
        WRITE(stdout,'(5X,A,1X,A,1X,A,1X,A,1X,A)') 'band  ', '   QP en. [eV]', 'conv', '   QP en. [eV]', 'conv'
        DO ib = 1, n_bands
           symb(1:2)='  no'
           IF( (iteration /= 0) .AND. (ABS(en2(ib,k_grid%ipis2ips(ik,1))-en1(ib,k_grid%ipis2ips(ik,1))) &
             & < trev_secant) ) symb(1)=' yes'
           IF( (iteration /= 0) .AND. (ABS(en2(ib,k_grid%ipis2ips(ik,2))-en1(ib,k_grid%ipis2ips(ik,2))) &
             & < trev_secant) ) symb(2)=' yes'
           WRITE(stdout,'(5X,I6.6,1X,1F14.6,1X,A,1F14.6,1X,A)') qp_bands(ib,1), en2(ib,k_grid%ipis2ips(ik,1))*rytoev, &
           symb(1), en2(ib,k_grid%ipis2ips(ik,2))*rytoev, symb(2)
        ENDDO
     ENDIF
     IF(k_grid%np > 1 .AND. ik < k_grid%np) WRITE(stdout,'(5X, 33(A))') '---------------------------------'
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
     WRITE(csecitr,'(I10)') secitr
     !
     CALL json%update('exec.Q.secitr', secitr, found )
     !
     counter = 0
     DO iks = 1, k_grid%nps
        ik = k_grid%ip(iks)
        is = k_grid%is(iks)
        DO ib_index = 1, n_bands
           ib = qp_bands(ib_index,is)
           counter = counter + 1
           WRITE(ccounter, '(I10)') counter
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').ksb',(/ik,is,ib/))
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').ein('//TRIM(ADJUSTL(csecitr))//')', &
           & en1(ib_index,iks)*rytoev)
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').eout('//TRIM(ADJUSTL(csecitr))//')', &
           & en2(ib_index,iks)*rytoev)
           CALL json%add('exec.Q.en('//TRIM(ADJUSTL(ccounter))//').sc_ein('//TRIM(ADJUSTL(csecitr))//')', &
           & (/REAL(sc1(ib_index,iks)*rytoev,KIND=DP),AIMAG(sc1(ib_index,iks)*rytoev)/))
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
