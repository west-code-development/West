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
!----------------------------------------------------------------------------
SUBROUTINE solve_hf ()
  !----------------------------------------------------------------------------
  !
  USE control_flags,   ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  IF( gamma_only ) THEN
    CALL solve_hf_gamma()
  ELSE
    CALL solve_hf_k( )
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE solve_hf_gamma()
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies.
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : qp_bands,n_bands,sigma_exx,sigma_vxcl,sigma_vxcnl,&
                                 & sigma_hf,l_enable_off_diagonal,sigma_exx_full,&
                                 & sigma_vxcl_full,sigma_vxcnl_full,sigma_hf_full
  USE mp_world,             ONLY : mpime,root
  USE pwcom,                ONLY : et
  USE io_push,              ONLY : io_push_title
  USE constants,            ONLY : rytoev
  USE west_io,              ONLY : serial_table_output
  USE wfreq_io,             ONLY : writeout_solvehf
  USE types_bz_grid,        ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: myglobk
  INTEGER :: ib, iks
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: sigma_exx_all_occupied(:,:)
  !
  CALL start_clock('solve_hf')
  !
  CALL io_push_title("Hartree-Fock Exact E(X)change")
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
  IF (l_enable_off_diagonal) sigma_hf_full(:,:) = sigma_exx_full(:,:) &
  & - sigma_vxcl_full(:,:) - sigma_vxcnl_full(:,:)
  !
  CALL writeout_solvehf( sigma_hf(1,1), n_bands, k_grid%nps )
  !
  ! Output it per k-point
  !
  ALLOCATE(out_tab(n_bands,6))
  DO iks=1,k_grid%nps
     DO ib = 1, n_bands
        out_tab( ib, 1) = REAL( qp_bands(ib), KIND=DP)
        out_tab( ib, 2) = et(qp_bands(ib),iks) * rytoev
        out_tab( ib, 3) = sigma_exx(ib,iks) * rytoev
        out_tab( ib, 4) = sigma_vxcl(ib,iks) * rytoev
        out_tab( ib, 5) = sigma_vxcnl(ib,iks) * rytoev
        out_tab( ib, 6) = ( et(qp_bands(ib),iks) + sigma_hf(ib,iks) ) * rytoev
     ENDDO
     WRITE(myglobk,'(i5.5)') iks
     !
     CALL serial_table_output(mpime==root,'ehf_K'//myglobk,out_tab,&
     & n_bands,6,&
     & (/'      band','    E0[eV]','    Sx[eV]','  Vxcl[eV]',' Vxcnl[eV]','   EHF[eV]'/))
  ENDDO
  DEALLOCATE(out_tab)
  !
  !DEALLOCATE( sigma_exx  )
  !DEALLOCATE( sigma_vxcl )
  !DEALLOCATE( sigma_vxcnl)
  !DEALLOCATE( sigma_hf   )
  !
  CALL stop_clock( "solve_hf" )
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
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : qp_bands,n_bands,sigma_exx,&
                                   & sigma_vxcl,sigma_vxcnl,sigma_hf
  USE mp_world,             ONLY : mpime,root
  USE pwcom,                ONLY : et
  USE io_push,              ONLY : io_push_title
  USE constants,            ONLY : rytoev
  USE west_io,              ONLY : serial_table_output
  USE wfreq_io,             ONLY : writeout_solvehf
  USE types_bz_grid,        ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: myglobk
  INTEGER :: ib, iks
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: sigma_exx_all_occupied(:,:)
  !
  CALL start_clock('solve_hf')
  !
  CALL io_push_title("Hartree-Fock Exact E(X)change")
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
  CALL writeout_solvehf( sigma_hf(1,1), n_bands, k_grid%nps  )
  !
  ! Output it per k-point
  !
  ALLOCATE(out_tab(n_bands,6))
  DO iks=1,k_grid%nps
     DO ib = 1, n_bands
        out_tab( ib, 1) = REAL( qp_bands(ib), KIND=DP)
        out_tab( ib, 2) = et(ib,iks) * rytoev
        out_tab( ib, 3) = sigma_exx(ib,iks) * rytoev
        out_tab( ib, 4) = sigma_vxcl(ib,iks) * rytoev
        out_tab( ib, 5) = sigma_vxcnl(ib,iks) * rytoev
        out_tab( ib, 6) = ( et(ib,iks) + sigma_hf(ib,iks) ) * rytoev
     ENDDO
     WRITE(myglobk,'(i5.5)') iks
     !
     CALL serial_table_output(mpime==root,'ehf_K'//myglobk,out_tab,&
     & n_bands,6,&
     & (/'      band','    E0[eV]','    Sx[eV]','  Vxcl[eV]',' Vxcnl[eV]','   EHF[eV]'/))
  ENDDO
  DEALLOCATE(out_tab)
  !
  !DEALLOCATE( sigma_exx  )
  !DEALLOCATE( sigma_vxcl )
  !DEALLOCATE( sigma_vxcnl)
  !DEALLOCATE( sigma_hf   )
  !
  CALL stop_clock( "solve_hf" )
  !
END SUBROUTINE
