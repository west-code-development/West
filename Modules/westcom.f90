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
  INTEGER, ALLOCATABLE :: qp_bands(:)
  REAL(DP) :: ecut_imfreq
  REAL(DP) :: ecut_refreq
  REAL(DP) :: wfreq_eta
  INTEGER :: n_secant_maxiter
  REAL(DP) :: trev_secant
  LOGICAL :: l_enable_lanczos
  LOGICAL :: l_enable_off_diagonal
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
  INTEGER,ALLOCATABLE :: ijpmap(:,:)
  INTEGER,ALLOCATABLE :: pijmap(:,:)
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
  REAL(DP),ALLOCATABLE    :: sigma_z         (:,:)
  REAL(DP),ALLOCATABLE    :: sigma_eqplin    (:,:)
  REAL(DP),ALLOCATABLE    :: sigma_eqpsec    (:,:)
  REAL(DP),ALLOCATABLE    :: sigma_diff      (:,:)
  REAL(DP),ALLOCATABLE    :: sigma_exx       (:,:)
  REAL(DP),ALLOCATABLE    :: sigma_vxcl      (:,:)
  REAL(DP),ALLOCATABLE    :: sigma_vxcnl     (:,:)
  REAL(DP),ALLOCATABLE    :: sigma_hf        (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_sc_eks    (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_sc_eqplin (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_sc_eqpsec (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_spectralf (:,:,:)
  REAL(DP),ALLOCATABLE    :: sigma_freq      (:)
  REAL(DP),ALLOCATABLE    :: sigma_exx_full  (:,:) 
  REAL(DP),ALLOCATABLE    :: sigma_vxcl_full (:,:) 
  REAL(DP),ALLOCATABLE    :: sigma_vxcnl_full(:,:) 
  REAL(DP),ALLOCATABLE    :: sigma_hf_full   (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_sc_eks_full (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_sc_eqplin_full (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_corr_full (:,:)
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
  CHARACTER(LEN=8) :: westpp_calculation
  CHARACTER(LEN=7) :: westpp_format
  INTEGER :: westpp_n_pdep_eigen_to_use
  INTEGER :: westpp_range(2)
  LOGICAL :: westpp_sign
  REAL(DP) :: westpp_r0(3)
  REAL(DP) :: westpp_box(6)
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
  REAL(DP)              :: docc_thr = 0.001_DP ! Threshold for comparing occupation numbers
  REAL(DP)              :: de_thr = 0.001_DP   ! Threshold for comparing single-particle energies
                                               ! When two orbitals' energies differ by less than this threshold,
                                               ! they are considered degenerate. An error will be raised if two
                                               ! orbitals with different occupation numbers are degenerate in
                                               ! the fractional occupation case. See subroutine dfpt and
                                               ! compute_pt1_dpsi for details.
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
  USE occ_center
  USE io_unit_numbers
  !
END MODULE
