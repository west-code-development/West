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
MODULE bse_module
  !
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  LOGICAL                  :: bse_calc = .FALSE.
  LOGICAL                  :: l_wannier_repr   = .FALSE.
  INTEGER                  :: ngm_g_max
  REAL(DP)                 :: ovl_thr
  REAL(DP),    ALLOCATABLE :: et_qp(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_ks(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: u_matrix(:,:,:)
  REAL(DP),    ALLOCATABLE :: ovl_matrix(:,:,:)
  !
  ! pdep method
  !
  COMPLEX(DP), ALLOCATABLE :: kernel_kd1(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: kernel_kd2(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: fock_term_g(:,:,:)
  REAL(DP),    ALLOCATABLE :: ev_loc(:)
  LOGICAL                  :: epsilon_rpa = .TRUE.
  LOGICAL                  :: wstat_epsilon_rpa = .TRUE.
  LOGICAL                  :: wstat_use_ecutrho = .FALSE.
  INTEGER                  :: wstat_n_pdep_eigen = 10
  !
  ! bse parallel
  !
  INTEGER,     ALLOCATABLE :: size_index_matrix_lz(:)
  REAL(DP),    ALLOCATABLE :: index_matrix_lz(:,:,:)
  !
END MODULE
