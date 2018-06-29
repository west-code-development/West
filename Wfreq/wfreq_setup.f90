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
SUBROUTINE wfreq_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq,npwqx,west_prefix,wfreq_save_dir,&
                                   & n_pdep_eigen_to_use,n_imfreq,nbnd_occ,l_macropol,macropol_calculation,&
                                   & n_refreq,qp_bandrange,wfreq_calculation, fftdriver
  USE westcom,                ONLY : sigma_exx,sigma_vxcl,sigma_vxcnl,sigma_hf,sigma_z,sigma_eqplin,sigma_eqpsec,sigma_sc_eks,&
                                     & sigma_sc_eqplin,sigma_sc_eqpsec,sigma_diff,sigma_spectralf,sigma_freq,n_spectralf
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE pwcom,                  ONLY : nbnd,nks
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g
  USE io_files,               ONLY : tmp_dir
  USE distribution_center,    ONLY : pert,macropert,ifr,rfr,aband
  USE class_idistribute,      ONLY : idistribute
  USE wavefunctions_module,   ONLY : evc
  USE mod_mpiio,              ONLY : set_io_comm
  USE types_bz_grid,          ONLY : k_grid
  !
  IMPLICIT NONE
  !
  REAL(DP) :: q(3)
  REAL(DP) :: qq
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: ig,i
  LOGICAL :: l_generate_plot 
  !
  CALL do_setup ( ) 
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  CALL set_npwq()
  !
  IF(qp_bandrange(1)>nbnd) CALL errore('wfreq_setup','Err: qp_bandrange(1)>nbnd', 1) 
  IF(qp_bandrange(2)>nbnd) CALL errore('wfreq_setup','Err: qp_bandrange(2)>nbnd', 1) 
  !
  CALL set_nbndocc()
  !
  CALL my_mkdir( wfreq_save_dir )
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
  !
  CALL set_freqlists( )
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
  END SELECT
  !
  CALL set_io_comm( ) ! this defines communicator between heads of each image (me_bgrp==0) 
  !
  ! Allocate for output
  !
  ALLOCATE( sigma_exx       (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_vxcl      (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_vxcnl     (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_hf        (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_z         (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_eqplin    (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_eqpsec    (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_sc_eks    (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_sc_eqplin (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_sc_eqpsec (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
  ALLOCATE( sigma_diff      (qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
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
  l_generate_plot = .FALSE.
  DO i = 1, 8
     IF( wfreq_calculation(i:i) == 'P' ) l_generate_plot = .TRUE.
  ENDDO
  IF( l_generate_plot ) THEN 
     ALLOCATE( sigma_spectralf      (n_spectralf,qp_bandrange(1):qp_bandrange(2),k_grid%nps) )
     ALLOCATE( sigma_freq           (n_spectralf) )
     sigma_spectralf = 0._DP
     sigma_freq      = 0._DP
  ENDIF 
  !
END SUBROUTINE 
