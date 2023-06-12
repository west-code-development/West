!
! Copyright (C) 2015-2023 M. Govoni
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
  REAL(DP),    ALLOCATABLE :: ev(:)
  REAL(DP),    ALLOCATABLE :: ev_distr(:)
  COMPLEX(DP), ALLOCATABLE :: dng(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvg(:,:)
  LOGICAL,     ALLOCATABLE :: conv(:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dng
  ATTRIBUTES(PINNED) :: dvg
#endif
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
  REAL(DP),    ALLOCATABLE :: d_epsm1_ifr_a(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_ifr_a(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_rfr_a(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: d_epsm1_ifr
#endif
  !
  ! EPSILON with q-points
  !
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_ifr_q(:,:,:,:) ! EPSILON + iq (global in iq)
  COMPLEX(DP), ALLOCATABLE :: z_epsm1_rfr_q(:,:,:,:) ! EPSILON + iq (global in iq)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: z_epsm1_ifr_q
#endif
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
  REAL(DP),    ALLOCATABLE :: d_body1_ifr_full(:,:,:,:)
  REAL(DP),    ALLOCATABLE :: d_body2_ifr_full(:,:,:,:,:)
  REAL(DP),    ALLOCATABLE :: d_diago_full(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z_body_rfr_full(:,:,:,:)
  REAL(DP),    ALLOCATABLE :: d_head_ifr_a(:)
  COMPLEX(DP), ALLOCATABLE :: z_head_rfr_a(:)
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
  INTEGER, PARAMETER :: iuwfc = 20
  INTEGER :: lrwfc
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
  INTEGER :: n_exx_lowrank
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
  ! INPUT FOR server_control
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
  CHARACTER(LEN=9) :: wfreq_calculation
  INTEGER :: n_lanczos
  INTEGER :: n_pdep_eigen_to_use
  INTEGER :: n_imfreq
  INTEGER :: n_refreq
  INTEGER :: qp_bandrange(2)
  INTEGER, ALLOCATABLE :: qp_bands(:)
  REAL(DP) :: ecut_imfreq
  REAL(DP) :: ecut_refreq
  REAL(DP) :: wfreq_eta
  INTEGER :: n_secant_maxiter
  REAL(DP) :: trev_secant
  LOGICAL :: l_enable_lanczos
  LOGICAL :: l_enable_off_diagonal
  INTEGER :: n_pdep_eigen_off_diagonal
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
  ! qp_bands
  !
  INTEGER :: n_bands
  !
  ! off-diagonal entries mapping
  !
  INTEGER :: n_pairs
  INTEGER, ALLOCATABLE :: ijpmap(:,:)
  INTEGER, ALLOCATABLE :: pijmap(:,:)
  !
  ! downfolded Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: proj_c(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: proj_c
#endif
  COMPLEX(DP), ALLOCATABLE :: eri_w(:,:,:,:)
  LOGICAL :: l_qdet_verbose
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
  REAL(DP),    ALLOCATABLE :: sigma_z(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_eqplin(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_eqpsec(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_diff(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_exx(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_vxcl(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_vxcnl(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_hf(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eks(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eqplin(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eqpsec(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_spectralf(:,:,:)
  REAL(DP),    ALLOCATABLE :: sigma_freq(:)
  REAL(DP),    ALLOCATABLE :: sigma_exx_full(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_vxcl_full(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_vxcnl_full(:,:)
  REAL(DP),    ALLOCATABLE :: sigma_hf_full(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eks_full(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_sc_eqplin_full(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_corr_full(:,:)
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
  ! INPUT FOR westpp_control
  !
  CHARACTER(LEN=9) :: westpp_calculation
  CHARACTER(LEN=7) :: westpp_format
  INTEGER :: westpp_n_pdep_eigen_to_use
  INTEGER :: westpp_range(2)
  LOGICAL :: westpp_sign
  REAL(DP) :: westpp_r0(3)
  REAL(DP) :: westpp_box(6)
  INTEGER :: westpp_nr
  REAL(DP) :: westpp_rmax
  REAL(DP) :: westpp_epsinfty
  INTEGER :: westpp_n_liouville_to_use
  LOGICAL :: westpp_l_spin_flip
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
  CHARACTER(LEN=1) :: wbse_init_calculation
  CHARACTER(LEN=7) :: bse_method
  CHARACTER(LEN=1) :: localization
  CHARACTER(LEN=20) :: chi_kernel
  CHARACTER(LEN=512) :: wfc_from_qbox  ! wavefunction file name, extension is spin channel, e.g. 'qb_wfc.1'
  CHARACTER(LEN=512) :: bisection_info ! bisection file name, extension is spin channel, e.g. 'bis_info.1'
  REAL(DP) :: overlap_thr              ! overlap threshold for idx_matrix in wbse_init_qboxcoupling
  INTEGER :: spin_channel
  INTEGER :: n_trunc_bands
  LOGICAL :: l_local_repr
  LOGICAL :: l_pdep
  !
  ! Common workspace
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
  CHARACTER(LEN=1) :: wbse_calculation
  CHARACTER(LEN=20) :: solver
  CHARACTER(LEN=512) :: qp_correction
  REAL(DP) :: scissor_ope
  INTEGER :: n_liouville_times
  INTEGER :: n_liouville_eigen
  INTEGER :: n_liouville_maxiter
  INTEGER :: n_liouville_read_from_file
  REAL(DP) :: trev_liouville
  REAL(DP) :: trev_liouville_rel
  CHARACTER(LEN=3) :: wbse_ipol
  LOGICAL :: l_qp_correction
  LOGICAL :: l_preconditioning
  LOGICAL :: l_dipole_realspace
  LOGICAL :: l_pre_shift
  LOGICAL :: l_spin_flip
  LOGICAL :: l_spin_flip_kernel
  LOGICAL :: l_spin_flip_alda0
  LOGICAL :: l_print_spin_flip_kernel
  REAL(DP) :: spin_flip_cut1
  REAL(DP) :: wbse_epsinfty
  CHARACTER(LEN=1) :: spin_excitation
  LOGICAL :: l_forces
  INTEGER :: forces_state
  REAL(DP) :: forces_zeq_cg_tr
  REAL(DP) :: ddvxc_fd_coeff
  LOGICAL :: l_slow_tddft_k1d
  !
  ! FOR global variables
  !
  INTEGER :: nbndval0x
  LOGICAL :: l_bse ! BSE True, TDDFT False
  LOGICAL :: l_hybrid_tddft
  LOGICAL :: l_lanczos
  LOGICAL :: l_davidson
  LOGICAL :: l_bse_triplet
  LOGICAL :: l_reduce_io
  INTEGER :: n_tau
  REAL(DP) :: sigma_c_head
  REAL(DP) :: sigma_x_head
  !
  ! FOR global Lanzcos diago vars
  !
  COMPLEX(DP), ALLOCATABLE :: d0psi(:,:,:,:)
  REAL(DP),    ALLOCATABLE :: beta_store(:,:,:)
  REAL(DP),    ALLOCATABLE :: zeta_store(:,:,:,:)
  !
  ! FOR Davidson method
  !
  COMPLEX(DP), ALLOCATABLE :: dng_exc(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dvg_exc(:,:,:,:)
  !
  REAL(DP),    ALLOCATABLE :: et_qp(:,:)
  COMPLEX(DP), ALLOCATABLE :: u_matrix(:,:,:)
  REAL(DP),    ALLOCATABLE :: ovl_matrix(:,:,:)
  INTEGER,     ALLOCATABLE :: n_bse_idx(:)
  INTEGER,     ALLOCATABLE :: idx_matrix(:,:,:)
  INTEGER,     ALLOCATABLE :: tau_is_read(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: tau_all(:,:)
  REAL(DP),    ALLOCATABLE :: sf_kernel(:)
  !
  ! Common workspace
  !
  CHARACTER(LEN=512) :: wbse_save_dir
  CHARACTER(LEN=512) :: wbse_restart_dir
  !
END MODULE
!
!
MODULE occ_center
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: occupation(:,:)     ! Occupation of each band and k-point:
                                               ! 0 <= occupation(ib,iks) <= 1
                                               ! 1 <= ib <= nbnd, 1 <= iks <= nks == kpt_pool%nloc (iks NOT global)
  INTEGER,  ALLOCATABLE :: nbnd_occ(:)         ! max index of occupied bands per k-point:
                                               ! occupation(ib,iks) > 0, forall ib <= nbnd_occ(iks)
                                               ! occupation(ib,iks) = 0, forall ib  > nbnd_occ(iks)
                                               ! 1 <= ib <= nbnd, 1 <= iks <= nks == kpt_pool%nloc (iks NOT global)
  INTEGER,  ALLOCATABLE :: nbnd_occ_full(:)    ! max index of contiguous and fully occupied bands per k-point:
                                               ! occupation(ib,iks) = 1, forall ib <= nbnd_occ_full(iks)
                                               ! occupation(ib,iks) < 1, forall ib  > nbnd_occ_full(iks)
                                               ! 1 <= ib <= nbnd, 1 <= iks <= nks == kpt_pool%nloc (iks NOT global)
  LOGICAL               :: l_frac_occ          ! If .true. then occupations may be fractional
                                               ! nbnd_occ_full == nbnd_occ when l_frac_occ is .false.
  REAL(DP), PARAMETER   :: docc_thr = 0.001_DP ! Threshold for comparing occupation numbers
  REAL(DP), PARAMETER   :: de_thr = 0.001_DP   ! Threshold for comparing single-particle energies
                                               ! When two orbitals' energies differ by less than this threshold,
                                               ! they are considered degenerate. An error will be raised if two
                                               ! orbitals with different occupation numbers are degenerate in
                                               ! the fractional occupation case. See subroutine dfpt and
                                               ! compute_pt1_dpsi for details.
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
  USE wbse_init_center
  USE wbse_center
  USE occ_center
  !
END MODULE
