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
SUBROUTINE wfreq_setup
  !-----------------------------------------------------------------------
  !
  USE mp_global,              ONLY : nbgrp,inter_image_comm,my_image_id,mp_bcast
  USE westcom,                ONLY : nbnd_occ,alphapv_dfpt,wfreq_save_dir,n_pdep_eigen_to_use,&
                                   & n_imfreq,l_macropol,macropol_calculation,n_refreq,qp_bandrange,&
                                   & qp_bands,n_bands,wfreq_calculation,sigma_exx,sigma_vxcl,&
                                   & sigma_vxcnl,sigma_hf,sigma_z,sigma_eqplin,sigma_eqpsec,&
                                   & sigma_sc_eks,sigma_sc_eqplin,sigma_sc_eqpsec,sigma_diff,&
                                   & sigma_spectralf,sigma_freq,n_spectralf,&
                                   & l_enable_off_diagonal,ijpmap,pijmap,equalpairmap,n_pairs,&
                                   & sigma_exx_full,sigma_vxcl_full,sigma_vxcnl_full,sigma_corr_full,&
                                   & iuwfc,lrwfc,proj_c,proj_r,npwq,npwqx,fftdriver
  USE wavefunctions,          ONLY : evc,psic
  USE fft_base,               ONLY : dffts
  USE buffers,                ONLY : get_buffer
  USE fft_at_gamma,           ONLY : single_invfft_gamma
  USE pwcom,                  ONLY : nbnd,nkstot,nks,npw,npwx
  USE kinds,                  ONLY : DP
  USE xc_lib,                 ONLY : xclib_dft_is
  USE distribution_center,    ONLY : pert,macropert,ifr,rfr,aband,occband,band_group,kpt_pool
  USE class_idistribute,      ONLY : idistribute,IDIST_BLK
  USE types_bz_grid,          ONLY : k_grid
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: i,ib,jb,index,iks,ib_index
  LOGICAL :: l_generate_plot
  !
  CALL do_setup()
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  CALL set_npwq()
  !
  IF(qp_bands(1) == 0) THEN
     IF(ALLOCATED(qp_bands)) DEALLOCATE(qp_bands)
     ALLOCATE(qp_bands(qp_bandrange(2)-qp_bandrange(1)+1))
     DO i = 1, SIZE(qp_bands)
        qp_bands(i) = qp_bandrange(1)+i-1
     ENDDO
  ENDIF
  !
  n_bands = SIZE(qp_bands)
  !
  IF(qp_bands(1) > nbnd) CALL errore('wfreq_setup','qp_bands(1)>nbnd',1)
  IF(qp_bands(n_bands) > nbnd) CALL errore('wfreq_setup','qp_bands(n_bands)>nbnd',1)
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
  IF(kpt_pool%nloc /= nks) CALL errore('wfreq_setup','unexpected kpt_pool initialization error',1)
  IF(nbgrp > n_bands) CALL errore('wfreq_setup','nbgrp>nbnd_qp',1)
  IF(nbgrp > MINVAL(nbnd_occ)) CALL errore('wfreq_setup','nbgrp>nbnd_occ',1)
  !
  CALL set_freqlists()
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
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
  ALLOCATE( sigma_exx       (n_bands,k_grid%nps) )
  ALLOCATE( sigma_vxcl      (n_bands,k_grid%nps) )
  ALLOCATE( sigma_vxcnl     (n_bands,k_grid%nps) )
  ALLOCATE( sigma_hf        (n_bands,k_grid%nps) )
  ALLOCATE( sigma_z         (n_bands,k_grid%nps) )
  ALLOCATE( sigma_eqplin    (n_bands,k_grid%nps) )
  ALLOCATE( sigma_eqpsec    (n_bands,k_grid%nps) )
  ALLOCATE( sigma_sc_eks    (n_bands,k_grid%nps) )
  ALLOCATE( sigma_sc_eqplin (n_bands,k_grid%nps) )
  ALLOCATE( sigma_sc_eqpsec (n_bands,k_grid%nps) )
  ALLOCATE( sigma_diff      (n_bands,k_grid%nps) )
  IF (l_enable_off_diagonal) THEN
     n_pairs = n_bands*(n_bands+1)/2
     ALLOCATE(ijpmap(n_bands,n_bands))
     ALLOCATE(pijmap(2,n_pairs))
     ALLOCATE(equalpairmap(n_bands))
     index = 1
     DO ib = 1, n_bands
        DO jb = ib, n_bands
           ijpmap(ib,jb) = index
           ijpmap(jb,ib) = index
           pijmap(1,index) = ib
           pijmap(2,index) = jb
           IF (ib == jb) equalpairmap(ib) = index
           index = index + 1
        ENDDO 
     ENDDO
     ALLOCATE( sigma_exx_full (n_pairs,k_grid%nps) )
     ALLOCATE( sigma_vxcl_full (n_pairs,k_grid%nps) )
     ALLOCATE( sigma_vxcnl_full (n_pairs,k_grid%nps) )
     ALLOCATE( sigma_corr_full (n_pairs,k_grid%nps) )
  ENDIF
  sigma_exx = 0._DP      
  sigma_vxcl = 0._DP
  sigma_vxcnl = 0._DP
  sigma_hf = 0._DP
  sigma_z = 0._DP
  sigma_eqplin = 0._DP
  sigma_eqpsec = 0._DP
  sigma_sc_eks = 0._DP
  sigma_sc_eqplin = 0._DP
  sigma_sc_eqpsec = 0._DP
  sigma_diff = 0._DP
  IF (l_enable_off_diagonal) THEN
     sigma_exx_full = 0._DP
     sigma_vxcl_full = 0._DP
     sigma_vxcnl_full = 0._DP
     sigma_corr_full = 0._DP
  ENDIF
  !
  DO i = 1,9
     IF(wfreq_calculation(i:i) == 'H') THEN
        ALLOCATE( proj_r(dffts%nnr,n_bands,k_grid%nps) )
        ALLOCATE( proj_c(npwx,n_bands,k_grid%nps) )
        DO iks = 1, k_grid%nps 
           !
           IF(kpt_pool%nloc > 1) THEN
              IF ( my_image_id == 0 ) CALL get_buffer( evc, lrwfc, iuwfc, iks )
              CALL mp_bcast( evc, 0, inter_image_comm )
           ENDIF
           !
           DO ib_index = 1, n_bands
              !
              ib = qp_bands(ib_index)
              proj_c(:,ib_index,iks) = evc(:,ib)
              CALL single_invfft_gamma(dffts,npwq,npwqx,proj_c(1,ib_index,iks),psic,TRIM(fftdriver))
              proj_r(:,ib_index,iks) = psic(:)
              !
           END DO
           !
        END DO
        EXIT
        !
     ENDIF
  ENDDO
  !
  l_generate_plot = .FALSE.
  DO i = 1,9
     IF(wfreq_calculation(i:i) == 'P') l_generate_plot = .TRUE.
  ENDDO
  IF(l_generate_plot) THEN
     ALLOCATE(sigma_spectralf(n_spectralf,n_bands,k_grid%nps))
     ALLOCATE(sigma_freq     (n_spectralf))
     sigma_spectralf = 0._DP
     sigma_freq = 0._DP
  ENDIF
  !
END SUBROUTINE
