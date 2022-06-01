!
! Copyright (C) 2015-2021 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `LICENSE'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Marco Govoni
!
!-----------------------------------------------------------------------
MODULE scratch_area
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE fft_types, ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  ! COULOMB
  !
  INTEGER :: npwq
  INTEGER :: npwqx
  INTEGER :: npwq_g
  CHARACTER(LEN=6) :: fftdriver
  !
  ! DBS
  !
  REAL(DP),            ALLOCATABLE :: ev(:)
  REAL(DP),            ALLOCATABLE :: ev_distr(:)
  COMPLEX(DP),         ALLOCATABLE :: dng(:,:)
  COMPLEX(DP),         ALLOCATABLE :: dvg(:,:)
  LOGICAL,             ALLOCATABLE :: conv(:)
  !
  ! BANDS
  !
  INTEGER, ALLOCATABLE :: nbnd_occ(:)
  !
  ! Q-POINTS
  !
  INTEGER, ALLOCATABLE :: ngq(:)     ! equivalent of ngk(:) --> ex. ngq(iq) = LOCAL number of PW for (q+G) (global in iq)
  INTEGER, ALLOCATABLE :: igq_q(:,:) ! equivalent of igk_k(:,:) --> ex. igq_q(ig,iq) = map for FFT (global in iq )
  INTEGER, ALLOCATABLE :: ngq_g(:)   ! equivalent of ngk_g(:) --> ex. ngk_g(iq) = TOTAL number of PW for (q+G) (global in iq)
  !
  ! EPSILON
  !
  REAL(DP),    ALLOCATABLE :: d_epsm1_ifr(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_ifr(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_rfr(:,:,:)
  !
  ! EPSILON with q-points
  !
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_ifr_q(:,:,:,:) ! EPSILON + iq (global in iq)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_rfr_q(:,:,:,:) ! EPSILON + iq (global in iq)
  !
  ! CORRELATION
  !
  REAL(DP),    ALLOCATABLE :: d_head_ifr(:)
  COMPLEX(DP), ALLOCATABLE :: z_head_ifr(:)
  REAL(DP),    ALLOCATABLE :: d_body1_ifr(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_body1_ifr(:,:,:,:)
  REAL(DP),    ALLOCATABLE :: d_body2_ifr(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_body2_ifr(:,:,:,:,:)
  REAL(DP),    ALLOCATABLE :: d_diago(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_head_rfr(:)
  COMPLEX(DP), ALLOCATABLE :: z_body_rfr(:,:,:,:)
  !
  ! CORRELATION with q-points
  !
  COMPLEX(DP), ALLOCATABLE :: z_body1_ifr_q(:,:,:,:,:)   ! CORRELATION + iq (global in iq)
  COMPLEX(DP), ALLOCATABLE :: z_body2_ifr_q(:,:,:,:,:,:) ! CORRELATION + iq (global in iq)
  REAL(DP),    ALLOCATABLE :: d_diago_q(:,:,:,:,:)       ! CORRELATION + iq (global in iq)
  COMPLEX(DP), ALLOCATABLE :: z_body_rfr_q (:,:,:,:,:)   ! CORRELATION + iq (global in iq)
  !
  ! I/O
  !
  TYPE(fft_type_descriptor) :: dfft_io
  !
END MODULE
!
!
MODULE westin
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=512) :: outdir          ! main directory
  CHARACTER(LEN=512) :: west_prefix
  CHARACTER(LEN=512) :: qe_prefix
  CHARACTER(LEN=512) :: savedir         ! outdir/west_prefix.code.save
  CHARACTER(LEN=512) :: main_input_file ! input file
  CHARACTER(LEN=512) :: logfile         ! savedir/logfile.json
  !
END MODULE
!
!
MODULE wstat_center
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! INPUT FOR wstat_control
  !
  CHARACTER(LEN=2) :: wstat_calculation
  INTEGER :: n_pdep_basis
  INTEGER :: n_pdep_times
  INTEGER :: n_pdep_eigen
  INTEGER :: n_pdep_maxiter
  INTEGER :: n_dfpt_maxiter
  INTEGER :: n_steps_write_restart
  INTEGER :: n_pdep_restart_from_itr
  INTEGER :: n_pdep_read_from_file
  REAL(DP) :: trev_pdep
  REAL(DP) :: trev_pdep_rel
  REAL(DP) :: tr2_dfpt
  LOGICAL :: l_deflate
  LOGICAL :: l_kinetic_only
  LOGICAL :: l_minimize_exx_if_active
  LOGICAL :: l_use_ecutrho
  INTEGER, ALLOCATABLE :: qlist(:)
  !
  ! Common workspace
  !
  COMPLEX(DP) :: alphapv_dfpt
  CHARACTER(LEN=512) :: wstat_save_dir
  CHARACTER(LEN=512) :: wstat_restart_dir
  LOGICAL :: l_is_wstat_converged
  !
END MODULE
!
!
MODULE server_center
  !
  IMPLICIT NONE
  !
  ! INPUT for server_control
  !
  CHARACTER(LEN=:), ALLOCATABLE :: document
  !
END MODULE
!
!
MODULE wfreq_center
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! INPUT FOR wfreq_control
  !
  CHARACTER(LEN=8) :: wfreq_calculation
  INTEGER :: n_lanczos
  INTEGER :: n_pdep_eigen_to_use
  INTEGER :: n_imfreq
  INTEGER :: n_refreq
  INTEGER :: qp_bandrange(2)
  REAL(DP) :: ecut_imfreq
  REAL(DP) :: ecut_refreq
  REAL(DP) :: wfreq_eta
  INTEGER :: n_secant_maxiter
  REAL(DP) :: trev_secant
  LOGICAL :: l_enable_lanczos
  CHARACTER(LEN=1) :: macropol_calculation
  REAL(DP) :: exx_etot
  REAL(DP) :: o_restart_time
  REAL(DP) :: ecut_spectralf(2)
  INTEGER :: n_spectralf
  !
  ! Common workspace
  !
  CHARACTER(LEN=512) :: wfreq_save_dir
  CHARACTER(LEN=512) :: wfreq_restart_dir
  LOGICAL, PARAMETER :: l_skip_nl_part_of_hcomr = .FALSE.
  LOGICAL :: l_macropol
  !
  ! re freq
  !
  REAL(DP), ALLOCATABLE :: refreq_list(:)
  !
  ! im freq
  !
  REAL(DP), ALLOCATABLE :: imfreq_list(:)
  REAL(DP), ALLOCATABLE :: imfreq_list_integrate(:,:)
  REAL(DP), PARAMETER :: frequency_list_power = 2._DP
  !
  ! gw_etot
  !
  REAL(DP) :: dft_etot
  REAL(DP) :: dft_exc
  REAL(DP) :: gw_ecorr
  REAL(DP) :: gw_exx
  REAL(DP) :: gw_exc
  !
  ! output
  !
  REAL(DP),    ALLOCATABLE :: sigma_exx(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_vxcl(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_vxcnl(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_hf(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_z(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_eqplin(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_eqpsec(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eks(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eqplin(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eqpsec(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_diff(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_spectralf(:,:,:)
  REAL(DP),    ALLOCATABLE :: sigma_freq(:)
  !
END MODULE
!
!
MODULE westpp_center
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! INPUT FOR wfreq_control
  !
  CHARACTER(LEN=8) :: westpp_calculation
  CHARACTER(LEN=7) :: westpp_format
  INTEGER :: westpp_n_pdep_eigen_to_use
  INTEGER :: westpp_range(2)
  LOGICAL :: westpp_sign
  REAL(DP) :: westpp_r0(3)
  INTEGER :: westpp_nr
  REAL(DP) :: westpp_rmax
  REAL(DP) :: westpp_epsinfty
  !
  ! Common workspace
  !
  CHARACTER(LEN=512) :: westpp_save_dir
  !
END MODULE
!
!
MODULE wbse_init_center
  !
  USE kinds, ONLY :  DP
  !
  IMPLICIT NONE
  !
  ! INPUT FOR wbse_init
  !
  CHARACTER(LEN=1)   :: wbse_init_calculation
  CHARACTER(LEN=256) :: localization
  CHARACTER(LEN=256) :: chi_kernel
  CHARACTER(LEN=256) :: wfc_from_qbox
  CHARACTER(LEN=256) :: bisection_info      ! bisection info file name, the extension is spin channel. e.g. bisection_info='info.bis', default file = 'info.bis.1'
  !CHARACTER(LEN=256) :: qbox_bisec_wfc_filename
  INTEGER  :: which_spin_channel   ! when nspin>1,which_spin_channel determine the spin channel for func wbse_init_qboxcoupling_single_q

  REAL(DP) :: overlap_thr !overlap threshold for index_matrix calculation from wbse_init_qboxcoupling.f90(ovl_thr)
  !INTEGER :: n_pdep_eigen
  LOGICAL  :: l_use_localise_repr = .FALSE.   !flag  depend on localization
  LOGICAL  :: l_use_bisection_thr = .FALSE.   !flag  depend on localization
  !
  !LOGICAL  :: use_qbox            = .FALSE.   !control flow in wbse_init
  !LOGICAL  :: l_test_ovl          = .FALSE.   !control flow in wbse_init
  !LOGICAL  :: use_wstat_pdep      = .FALSE.   !control flow in wbse_init
  !
  LOGICAL :: l_xcchi       = .FALSE.   ! control XC CHI flow in wbse_init_qboxcoupling
  !
  ! FOR qbox_control
  !
!  INTEGER :: nrowmax
!  CHARACTER(LEN=256)  :: xml_file             ! xml file from qbox ground state calculation
!  CHARACTER(LEN=256) ::  xc                   ! xc functional
!  REAL(DP) ::            alpha_pbe0           ! alpha for PBE0 calculation
!  REAL(DP) ::            amplitude            ! amplitude for vext
!  CHARACTER(LEN=256) ::  wf_dyn               ! wavefunction update algorithm
!  REAL(DP) ::            btHF                 ! bisection threshold for HF exchange computation
!  CHARACTER(LEN=256) ::  blHF                 ! bisection levels for HF exchange computation
!  INTEGER :: nitscf
!  INTEGER :: nite
  !
  CHARACTER(LEN=512) :: wbse_init_save_dir
  CHARACTER(LEN=512) :: wbse_init_restart_dir
  !
END MODULE
!
!
MODULE wbse_center
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! INPUT FOR wbse_control
  !
  CHARACTER(LEN=1)   :: wbse_calculation !  wbse input
  CHARACTER(LEN=256) :: solver           !  wbse input
  CHARACTER(LEN=256) :: qp_correction    !  wbse input filename qp_correction
  REAL(DP) :: scissor_ope          !  wbse input

  INTEGER  :: n_liouville_times          !  wbse input
  INTEGER  :: n_liouville_eigen          !  wbse input
  INTEGER  :: n_liouville_maxiter        !  wbse input
  INTEGER  :: n_liouville_read_from_file !  wbse input

  REAL(DP) :: trev_liouville            !  wbse input
  REAL(DP) :: trev_liouville_rel        !  wbse input

  CHARACTER(LEN=3) :: ipol_input       !wbse input
  CHARACTER(LEN=1)   :: wbse_macropol_calculation !  wbse input
  !
  LOGICAL  :: l_qp_correction             !depend on qp_correction
  LOGICAL  :: l_bse_calculation           !depend on solver if BSE True if TDDFT False
  LOGICAL  :: l_diag_term_only  = .FALSE. !flag nolonger read from wbse.in, always FALSE. affect td_liouville_oper.f90
  LOGICAL  :: l_preconditioning           !flag read from wbse.in
  !
  REAL(DP) :: epsinfty              !  wbse input
  !REAL(DP) :: eps_macro            !  wbse input
  !
  !CHARACTER(LEN=256) :: wbse_diag_method !  wbse input
  CHARACTER(LEN=256) :: spin_excitation  !  wbse input
  !
  ! FOR global variables
  !
  INTEGER :: nbndval0x                 !wbse  wbse_init
  LOGICAL :: l_lanzcos     = .FALSE.   !wbse flag depend on wbse_calculation
  LOGICAL :: l_davidson    = .FALSE.   !wbse flag depend on wbse_calculation
  LOGICAL :: l_bse_triplet = .FALSE.   !wbse flag depend on spin_excitation
  REAL(DP):: sigma_c_head  = 0.0_DP
  REAL(DP):: sigma_x_head  = 0.0_DP

  LOGICAL          :: macropol_dfpt  !TODO: QUESTION  depend on wbse_macropol_calculation? macropol_calculation=N macropol_dfpt=False? affect wbse_solve_e_psi.f90
  !
  ! FOR INPUT Lanzcos diago
  !
  !CHARACTER(LEN=1) :: wlz_calculation  !wbse lanzcos calculation (no bcast)

  !INTEGER          :: n_lzstep = 0     !wbse input

  !
  ! FOR global Lanzcos diago vars
  !
  COMPLEX(DP),ALLOCATABLE :: d0psi(:,:,:,:)
  REAL(DP),ALLOCATABLE    :: alpha_store(:,:,:)
  REAL(DP),ALLOCATABLE    :: beta_store(:,:,:)
  REAL(DP),ALLOCATABLE    :: gamma_store(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zeta_store(:,:,:,:)
  !
  ! FOR Davidson method
  !
  COMPLEX(DP),ALLOCATABLE :: dng_exc(:,:,:,:)
  COMPLEX(DP),ALLOCATABLE :: dvg_exc(:,:,:,:)
  !
  CHARACTER(LEN=512) :: wbse_save_dir
  CHARACTER(LEN=512) :: wbse_restart_dir
  !
  !
END MODULE
!
!
MODULE wbsepp_center
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! INPUT FOR wbsepp_control
  !
  LOGICAL :: l_meg                 = .FALSE.       ! local flag of wbsepp determained by wbsepp_type
  LOGICAL :: l_eig_decomp          = .FALSE.       ! local flag of wbsepp determained by wbsepp_type
  LOGICAL :: l_lz_spec             = .FALSE.       ! local flag of wbsepp determained by wbsepp_type
  LOGICAL :: l_exc_plot            = .FALSE.       ! local flag of wbsepp determained by wbsepp_type
  LOGICAL :: l_exc_rho_res_plot    = .FALSE.       ! local flag of wbsepp determained by wbsepp_type
  !
  !
  INTEGER :: wbsepp_type
  !
  ! lzc part
  !
  INTEGER :: itermax, itermax0, ipol, sym_op, units, verbosity
  INTEGER :: spin_channel
  !
  CHARACTER(len=60) :: extrapolation
  REAL(dp) :: start,end,increment
  REAL(dp) :: epsil
  !
  ! exc plot part
  !
  REAL(dp) :: r0_input(3)
  INTEGER  :: iexc_plot
  !
  CHARACTER(LEN=512) :: wbsepp_save_dir
  !
END MODULE
!
!
MODULE wan_center
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: wanc(:,:)
  REAL(DP), ALLOCATABLE :: wanu(:,:)
  INTEGER :: wantot
  !
END MODULE
!
!
MODULE io_unit_numbers
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: iuwfc = 20
  INTEGER :: lrwfc
  !
END MODULE
!
!
MODULE westcom
  !
  USE scratch_area
  USE westin
  USE wstat_center
  USE server_center
  USE wfreq_center
  USE westpp_center
  USE wan_center
  USE io_unit_numbers
  USE wbse_init_center
  USE wbse_center
  USE wbsepp_center
  !
END MODULE
