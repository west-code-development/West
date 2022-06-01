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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
MODULE bse_module
  !
  USE kinds,                ONLY : DP
  !USE class_idistribute,    ONLY : idistribute
  !
  IMPLICIT NONE
  !
  ! general vars
  !
  LOGICAL,     PUBLIC :: bse_calc = .false.
  LOGICAL,     PUBLIC :: l_wannier_repr   = .false.
  INTEGER,     PUBLIC :: ngm_g_max
  REAL(DP),    PUBLIC :: ovl_thr
  REAL(DP),    PUBLIC, ALLOCATABLE :: et_qp(:,:)                   !qp correction
  COMPLEX(DP), PUBLIC, ALLOCATABLE :: evc_ks(:,:,:)
  COMPLEX(DP), PUBLIC, ALLOCATABLE :: u_matrix(:,:,:)
  REAL(DP),    PUBLIC, ALLOCATABLE :: ovl_matrix(:,:,:)
  !
  ! pdep method
  !
  COMPLEX(DP), PUBLIC, ALLOCATABLE :: kernel_kd1(:,:,:)
  COMPLEX(DP), PUBLIC, ALLOCATABLE :: kernel_kd2(:,:,:)
  COMPLEX(DP), PUBLIC, ALLOCATABLE :: fock_term_g(:,:,:)
  REAL(DP),    PUBLIC, ALLOCATABLE :: ev_loc(:)
  LOGICAL,     PUBLIC :: epsilon_rpa = .true.
  LOGICAL,     PUBLIC :: wstat_epsilon_rpa = .true.
  LOGICAL,     PUBLIC :: wstat_use_ecutrho = .false.
  INTEGER,     PUBLIC :: wstat_n_pdep_eigen = 10
  !
  ! bse parallel
  !
  INTEGER, ALLOCATABLE    :: size_index_matrix_lz(:)
  REAL(DP),ALLOCATABLE    :: index_matrix_lz(:,:,:)
  !TYPE(idistribute)       :: bseparal
  !
ENDMODULE bse_module
