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
SUBROUTINE wfreq_setup
  !-----------------------------------------------------------------------
  !
  USE mp_global,              ONLY : inter_image_comm,my_image_id,inter_pool_comm,intra_bgrp_comm,nbgrp
  USE mp,                     ONLY : mp_bcast,mp_sum
  USE westcom,                ONLY : lrwfc,iuwfc,wfreq_save_dir,wfreq_calculation,nbnd_occ,occupation,&
                                   & qp_bands,n_bands,alphapv_dfpt,n_imfreq,n_refreq,n_pdep_eigen_to_use,&
                                   & l_macropol,macropol_calculation,sigma_exx,sigma_vxcl,sigma_vxcnl,&
                                   & sigma_hf,sigma_z,sigma_eqplin,sigma_eqpsec,sigma_sc_eks,&
                                   & sigma_sc_eqplin,sigma_sc_eqpsec,sigma_diff,sigma_spectralf,&
                                   & sigma_freq,n_spectralf,l_enable_off_diagonal,ijpmap,pijmap,n_pairs,&
                                   & sigma_exx_full,sigma_vxcl_full,sigma_vxcnl_full,sigma_hf_full,&
                                   & sigma_sc_eks_full,sigma_sc_eqplin_full,sigma_corr_full,proj_c,&
                                   & qdet_dc,l_dc2025
  USE wavefunctions,          ONLY : evc
  USE buffers,                ONLY : get_buffer
  USE pwcom,                  ONLY : nbnd,nkstot,nks,npw,npwx,nspin,ngk
  USE noncollin_module,       ONLY : npol
  USE kinds,                  ONLY : DP
  USE xc_lib,                 ONLY : xclib_dft_is
  USE distribution_center,    ONLY : pert,kpt_pool,band_group,macropert,ifr,rfr,aband,occband
  USE class_idistribute,      ONLY : idistribute,IDIST_BLK
  USE types_bz_grid,          ONLY : k_grid
  USE ldaU,                   ONLY : lda_plus_u
  USE bp,                     ONLY : lelfield
  USE realus,                 ONLY : real_space
  USE control_flags,          ONLY : gamma_only
  USE wfreq_db,               ONLY : qdet_db_write_overlap
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: i,ib,jb,ipair,iks,iks_g,is,ib_index
  LOGICAL :: l_generate_plot
  LOGICAL :: l_QDET
  REAL(DP) :: this_occ,next_occ
  REAL(DP),ALLOCATABLE :: overlap_ab(:,:)
  !
  CALL do_setup()
  !
  n_bands = SIZE(qp_bands,1)
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  CALL set_npwq()
  !
  CALL set_nbndocc()
  !
  CALL my_mkdir(wfreq_save_dir)
  !
  pert = idistribute()
  CALL pert%init(n_pdep_eigen_to_use,'i','npdep',.TRUE.)
  macropert = idistribute()
  CALL macropert%init(n_pdep_eigen_to_use+3,'i','npdep+macro',.TRUE.)
  ifr = idistribute()
  CALL ifr%init(n_imfreq,'z','n_imfreq',.TRUE.)
  rfr = idistribute()
  CALL rfr%init(n_refreq,'z','n_refreq',.TRUE.)
  aband = idistribute()
  CALL aband%init(nbnd,'i','nbnd',.TRUE.)
  occband = idistribute()
  band_group = idistribute()
  !
  kpt_pool = idistribute()
  CALL kpt_pool%init(nkstot,'p','nkstot',.FALSE.,IDIST_BLK)
  !
  IF(kpt_pool%nloc /= nks) CALL errore('wfreq_setup','unexpected kpt_pool init error',1)
  IF(nbgrp > n_bands) CALL errore('wfreq_setup','nbgrp>nbnd_qp',1)
  IF(nbgrp > MINVAL(nbnd_occ)) CALL errore('wfreq_setup','nbgrp>nbnd_occ',1)
  !
  CALL set_freqlists()
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
  CASE('n','N')
     l_macropol = .FALSE.
  END SELECT
  !
  SELECT CASE(qdet_dc)
  CASE('DC2025','dc2025')
     l_dc2025 = .TRUE.
  CASE DEFAULT
     l_dc2025 = .FALSE.
  END SELECT
  !
  IF(xclib_dft_is('hybrid')) THEN
     IF(l_macropol) THEN
        IF(macropert%nlocx > nbnd) CALL errore('wfreq_setup','EXX and macropert%nlocx>nbnd',1)
     ELSE
        IF(pert%nlocx > nbnd) CALL errore('wfreq_setup','EXX and pert%nlocx>nbnd',1)
     ENDIF
  ENDIF
  !
  ! Allocate for output
  !
  ALLOCATE(sigma_exx    (n_bands,k_grid%nps))
  ALLOCATE(sigma_vxcl   (n_bands,k_grid%nps))
  ALLOCATE(sigma_vxcnl  (n_bands,k_grid%nps))
  ALLOCATE(sigma_hf     (n_bands,k_grid%nps))
  ALLOCATE(sigma_z      (n_bands,k_grid%nps))
  ALLOCATE(sigma_eqplin (n_bands,k_grid%nps))
  ALLOCATE(sigma_eqpsec (n_bands,k_grid%nps))
  ALLOCATE(sigma_diff   (n_bands,k_grid%nps))
  ALLOCATE(sigma_sc_eks (n_bands,k_grid%nps))
  sigma_exx = 0._DP
  sigma_vxcl = 0._DP
  sigma_vxcnl = 0._DP
  sigma_hf = 0._DP
  sigma_z = 0._DP
  sigma_eqplin = 0._DP
  sigma_eqpsec = 0._DP
  sigma_diff = 0._DP
  sigma_sc_eks = 0._DP
  !
  l_QDET = .FALSE.
  DO i = 1,9
     IF(wfreq_calculation(i:i) == 'H') l_QDET = .TRUE.
  ENDDO
  IF(l_QDET .OR. l_enable_off_diagonal) THEN
     n_pairs = n_bands*(n_bands+1)/2
     ALLOCATE(ijpmap(n_bands,n_bands))
     ALLOCATE(pijmap(2,n_pairs))
     ipair = 1
     DO ib = 1, n_bands
        DO jb = ib, n_bands
           ijpmap(ib,jb) = ipair
           ijpmap(jb,ib) = ipair
           pijmap(1,ipair) = ib
           pijmap(2,ipair) = jb
           ipair = ipair + 1
        ENDDO
     ENDDO
  ENDIF
  !
  IF(.NOT. l_enable_off_diagonal) THEN
     ALLOCATE(sigma_sc_eqplin (n_bands,k_grid%nps))
     ALLOCATE(sigma_sc_eqpsec (n_bands,k_grid%nps))
     sigma_sc_eqplin = 0._DP
     sigma_sc_eqpsec = 0._DP
  ELSE
     ALLOCATE(sigma_exx_full       (n_pairs,k_grid%nps))
     ALLOCATE(sigma_vxcl_full      (n_pairs,k_grid%nps))
     ALLOCATE(sigma_vxcnl_full     (n_pairs,k_grid%nps))
     ALLOCATE(sigma_hf_full        (n_pairs,k_grid%nps))
     ALLOCATE(sigma_sc_eks_full    (n_pairs,k_grid%nps))
     ALLOCATE(sigma_sc_eqplin_full (n_pairs,k_grid%nps))
     ALLOCATE(sigma_corr_full      (n_pairs,k_grid%nps))
     sigma_exx_full = 0._DP
     sigma_vxcl_full = 0._DP
     sigma_vxcnl_full = 0._DP
     sigma_hf_full = 0._DP
     sigma_sc_eks_full = 0._DP
     sigma_sc_eqplin_full = 0._DP
     sigma_corr_full = 0._DP
  ENDIF
  !
  IF(l_QDET) THEN
     !
     IF(real_space) CALL errore('wfreq_setup','QDET with real_space not supported',1)
     IF(lda_plus_u) CALL errore('wfreq_setup','QDET with lda_plus_u not supported',1)
     IF(lelfield) CALL errore('wfreq_setup','QDET with lelfield not supported',1)
     IF(.NOT. gamma_only) CALL errore('wfreq_setup','QDET requires gamma_only',1)
     !
     ! qp_bands can be sorted or unsorted, but all occupied bands must appear before empty ones
     !
     DO iks = 1, kpt_pool%nloc
        !
        iks_g = kpt_pool%l2g(iks)
        is = k_grid%is(iks_g)
        !
        DO ib_index = 1, n_bands-1
           !
           ib = qp_bands(ib_index,is)
           this_occ = occupation(ib,iks)
           !
           ib = qp_bands(ib_index+1,is)
           next_occ = occupation(ib,iks)
           !
           IF(next_occ > this_occ) &
           & CALL errore('wfreq_setup','occupied bands must be given before empty ones in qp_bands',1)
           !
        ENDDO
        !
     ENDDO
     !
     ALLOCATE(proj_c(npwx,n_bands,k_grid%nps))
     !
     proj_c(:,:,:) = 0._DP
     DO iks = 1, kpt_pool%nloc
        !
        iks_g = kpt_pool%l2g(iks)
        is = k_grid%is(iks_g)
        npw = ngk(iks)
        !
        IF(kpt_pool%nloc > 1) THEN
           IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc,0,inter_image_comm)
        ENDIF
        !
        DO ib_index = 1, n_bands
           ib = qp_bands(ib_index,is)
           proj_c(:,ib_index,iks_g) = evc(:,ib)
        ENDDO
        !
     ENDDO
     !
     CALL mp_sum(proj_c,inter_pool_comm)
     !
     !$acc enter data copyin(proj_c)
     !
     IF(nspin == 2) THEN
        !
        ALLOCATE(overlap_ab(n_bands,n_bands))
        !$acc enter data create(overlap_ab)
        !
        CALL glbrak_gamma(proj_c(:,:,1),proj_c(:,:,2),overlap_ab,npw,npwx,n_bands,n_bands,n_bands,npol)
        !
        !$acc update host(overlap_ab)
        !
        CALL mp_sum(overlap_ab,intra_bgrp_comm)
        !
        CALL qdet_db_write_overlap(overlap_ab)
        !
        !$acc exit data delete(overlap_ab)
        DEALLOCATE(overlap_ab)
        !
     ENDIF
     !
  ENDIF
  !
  l_generate_plot = .FALSE.
  DO i = 1,9
     IF(wfreq_calculation(i:i) == 'P') l_generate_plot = .TRUE.
  ENDDO
  IF(l_generate_plot) THEN
     ALLOCATE(sigma_spectralf(n_spectralf,n_bands,k_grid%nps))
     ALLOCATE(sigma_freq(n_spectralf))
     sigma_spectralf = 0._DP
     sigma_freq = 0._DP
  ENDIF
  !
  CALL wfreq_memory_report()
  !
END SUBROUTINE
