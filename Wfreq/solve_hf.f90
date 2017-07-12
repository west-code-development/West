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
!-----------------------------------------------------------------------
SUBROUTINE solve_hf()
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies. 
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE kinds,                ONLY : DP 
  USE westcom,              ONLY : sqvc,west_prefix,n_pdep_eigen_to_use,n_lanczos,qp_bandrange,iks_l2g,&
                                 & n_secant_maxiter,n_imfreq,nbnd_occ,trev_secant,l_enable_gwetot,exx_etot, &
                                 & sigma_exx, sigma_vxcl, sigma_vxcnl, sigma_hf
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm, &
                                 & root_bgrp,me_bgrp,me_image,root_image
  USE mp_world,             ONLY : nproc,mpime,root,world_comm
  USE mp,                   ONLY : mp_bcast,mp_barrier
  USE io_global,            ONLY : stdout, ionode
  USE pwcom,                ONLY : et,nks,current_spin,isk,xk,nbnd,lsda,g2kin,current_k
  USE wavefunctions_module, ONLY : evc,psic,psic_nc
  USE io_files,             ONLY : tmp_dir
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE noncollin_module,     ONLY : noncolin,npol 
  USE constants,            ONLY : rytoev
  !USE west_io,              ONLY : serial_table_output
  USE distribution_center,  ONLY : pert
  USE funct,                ONLY : get_exx_fraction,dft_is_hybrid
  USE klist,                ONLY : wk
  USE wfreq_io,             ONLY : writeout_solvehf 
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
  CALL calc_exx2(sigma_exx, qp_bandrange(1), qp_bandrange(2))
  !
  ! Get SIGMA X
  !
  sigma_hf(:,:) = sigma_exx(:,:) - sigma_vxcl(:,:) - sigma_vxcnl(:,:)
  !
  CALL writeout_solvehf( sigma_hf(qp_bandrange(1),1), qp_bandrange(2)-qp_bandrange(1)+1, nks  )
! !
! ! Output it per k-point
! !
! ALLOCATE(out_tab(qp_bandrange(2)-qp_bandrange(1)+1,6))
! DO iks=1,nks
!    DO ib = qp_bandrange(1), qp_bandrange(2)
!       out_tab( ib - qp_bandrange(1) + 1, 1) = REAL( ib, KIND=DP) 
!       out_tab( ib - qp_bandrange(1) + 1, 2) = et(ib,iks) * rytoev
!       out_tab( ib - qp_bandrange(1) + 1, 3) = sigma_exx(ib,iks) * rytoev
!       out_tab( ib - qp_bandrange(1) + 1, 4) = sigma_vxcl(ib,iks) * rytoev
!       out_tab( ib - qp_bandrange(1) + 1, 5) = sigma_vxcnl(ib,iks) * rytoev
!       out_tab( ib - qp_bandrange(1) + 1, 6) = ( et(ib,iks) + sigma_hf(ib,iks) ) * rytoev
!    ENDDO
!    WRITE(myglobk,'(i5.5)') iks_l2g(iks)
!    !
!    CALL serial_table_output(mpime==root,4000,'ehf_K'//myglobk,out_tab,&
!    & qp_bandrange(2)-qp_bandrange(1)+1,6,&
!    & (/'      band','    E0[eV]','    Sx[eV]','  Vxcl[eV]',' Vxcnl[eV]','   EHF[eV]'/))
! ENDDO
! DEALLOCATE(out_tab)
  !
  ! Compute the exact exchange energy (used to calculate total GW energy) 
  !
  IF( l_enable_gwetot) THEN
     !
     nbndval = MIN( MAXVAL( nbnd_occ(:) ), nbnd ) 
     ALLOCATE(sigma_exx_all_occupied(nbndval,nks))
     !
     CALL calc_exx2( sigma_exx_all_occupied, 1, nbndval ) 
     !
     exx_etot = 0._DP
     DO iks = 1, nks 
        DO ib = 1, nbnd_occ(iks) 
           exx_etot = exx_etot + sigma_exx_all_occupied( ib, iks) * wk(iks) / 2._DP
        ENDDO
     ENDDO
     !
     DEALLOCATE( sigma_exx_all_occupied )  
     !
  ENDIF
  !
  !DEALLOCATE( sigma_exx  )
  !DEALLOCATE( sigma_vxcl )
  !DEALLOCATE( sigma_vxcnl)
  !DEALLOCATE( sigma_hf   )
  !
  CALL stop_clock( "solve_hf" )
  !
END SUBROUTINE
