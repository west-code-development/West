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
MODULE wbsecom
  !
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  ! INPUT FOR wbse_init
  !
  CHARACTER(LEN=1)   :: wbse_init_calculation
  CHARACTER(LEN=256) :: chi_kernel
  INTEGER  :: which_spin_channel
  REAL(DP) :: overlap_thr
  LOGICAL  :: l_use_localise_repr = .FALSE. 
  LOGICAL  :: l_test_ovl          = .FALSE. 
  LOGICAL  :: use_wstat_pdep      = .FALSE.
  LOGICAL  :: l_use_bisection_thr = .FALSE.
  !
  ! INPUT FOR wbse_control
  !
  INTEGER  :: n_plep_basis
  INTEGER  :: n_plep_times
  INTEGER  :: n_plep_eigen
  INTEGER  :: n_plep_maxiter
  INTEGER  :: n_plep_read_from_file
  ! 
  LOGICAL  :: l_qp_correction
  LOGICAL  :: l_bse_calculation
  LOGICAL  :: l_diag_term_only 
  LOGICAL  :: l_preconditioning 
  !
  REAL(DP) :: trev_plep
  REAL(DP) :: trev_plep_rel
  REAL(DP) :: scissor_ope 
  REAL(DP) :: eps_macro
  !  
  CHARACTER(LEN=1)   :: wbse_calculation
  CHARACTER(LEN=256) :: which_bse_method
  CHARACTER(LEN=256) :: wbse_diag_method
  CHARACTER(LEN=256) :: spin_excitation
  !
  ! FOR global variables
  !
  INTEGER :: nbndval0x
  LOGICAL :: l_lanzcos     = .FALSE.
  LOGICAL :: l_davidson    = .FALSE. 
  LOGICAL :: l_bse_triplet = .FALSE. 
  LOGICAL :: l_xcchi       = .FALSE.
  REAL(DP):: mac_isz
  !
  ! FOR INPUT Lanzcos diago
  !  
  CHARACTER(LEN=1) :: wlz_calculation
  CHARACTER(LEN=3) :: ipol_input 
  INTEGER          :: n_lzstep = 0
  LOGICAL          :: macropol_dfpt = .FALSE.
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
END MODULE
