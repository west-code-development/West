!
! Copyright (C) 2015-2024 M. Govoni
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
!----------------------------------------------------------------------------
SUBROUTINE solve_hf ( )
  !----------------------------------------------------------------------------
  !
  USE control_flags,   ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  IF( gamma_only ) THEN
    CALL solve_hf_gamma( )
  ELSE
    CALL solve_hf_k( )
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_hf_gamma( )
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies.
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE westcom,              ONLY : n_bands,sigma_exx,sigma_vxcl,sigma_vxcnl,sigma_hf,&
                                 & l_enable_off_diagonal,sigma_exx_full,sigma_vxcl_full,&
                                 & sigma_vxcnl_full,sigma_hf_full
  USE io_push,              ONLY : io_push_title
  USE wfreq_io,             ONLY : write_hf
  !
  IMPLICIT NONE
  !
  CALL start_clock('solve_hf')
  !
  CALL io_push_title('Hartree-Fock Exact E(X)change')
  !
  ! Get SIGMA Vxc
  !
  CALL calc_vxc(sigma_vxcl,sigma_vxcnl)
  !
  ! Get SIGMA EXX
  !
  CALL calc_exx2(sigma_exx, .FALSE.)
  !
  ! Get SIGMA X
  !
  sigma_hf(:,:) = sigma_exx(:,:) - sigma_vxcl(:,:) - sigma_vxcnl(:,:)
  IF (l_enable_off_diagonal) &
  & sigma_hf_full(:,:) = sigma_exx_full(:,:) - sigma_vxcl_full(:,:) - sigma_vxcnl_full(:,:)
  !
  CALL write_hf(sigma_hf,n_bands)
  !
  CALL stop_clock( 'solve_hf' )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_hf_k( )
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies.
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE westcom,              ONLY : n_bands,sigma_exx,sigma_vxcl,sigma_vxcnl,sigma_hf
  USE io_push,              ONLY : io_push_title
  USE wfreq_io,             ONLY : write_hf
  !
  IMPLICIT NONE
  !
  CALL start_clock('solve_hf')
  !
  CALL io_push_title('Hartree-Fock Exact E(X)change')
  !
  ! Get SIGMA Vxc
  !
  CALL calc_vxc(sigma_vxcl,sigma_vxcnl)
  !
  ! Get SIGMA EXX
  !
  CALL calc_exx2(sigma_exx,.FALSE.)
  !
  ! Get SIGMA X
  !
  sigma_hf(:,:) = sigma_exx(:,:) - sigma_vxcl(:,:) - sigma_vxcnl(:,:)
  !
  CALL write_hf(sigma_hf,n_bands)
  !
  CALL stop_clock( 'solve_hf' )
  !
END SUBROUTINE
