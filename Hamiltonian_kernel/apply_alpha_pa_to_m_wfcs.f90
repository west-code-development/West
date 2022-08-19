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
!-----------------------------------------------------------------------
SUBROUTINE apply_alpha_pa_to_m_wfcs(m,f,alpha)
!------------------------------------------------------------------------
  !
  ! | f_i > = alpha * P_a | f_i >            forall i = 1:n_bands
  ! P_a = sum_{i from 1 to n_bands} |proj_c_i><proj_c_i| 
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : npw,npwx
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE wavefunctions,        ONLY : evc
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : npol
  USE westcom,              ONLY : n_bands,proj_c
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m
  COMPLEX(DP), INTENT(IN) :: alpha
  COMPLEX(DP), INTENT(INOUT) :: f(npwx*npol,m)
  !
  ! Workspace
  !
  REAL(DP) :: alpha_r
  REAL(DP), ALLOCATABLE :: ps_r(:,:)
  COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
  !
  CALL start_clock ('alphapa')
  !
  ! ps_{ij} = < westcom/proj_c_i | f_j >
  !
  IF( gamma_only ) THEN
     !
     ALLOCATE( ps_r(n_bands,m) )
     ps_r = 0.0_DP
     alpha_r = REAL(alpha,KIND=DP)
     !
     CALL glbrak_gamma( proj_c, f, ps_r, npw, npwx, n_bands, m, n_bands, npol)
     CALL mp_sum(ps_r,intra_bgrp_comm)
     !
     CALL DGEMM('N','N',2*npwx*npol,m,n_bands,alpha_r,proj_c,2*npwx*npol,ps_r,n_bands,0.0_DP,f,2*npwx*npol)
     !
     DEALLOCATE( ps_r )
     !
  ELSE
     !
     ALLOCATE( ps_c(n_bands,m) )
     ps_c = (0.0_DP,0.0_DP)
     !
     CALL glbrak_k( proj_c, f, ps_c, npw, npwx, n_bands, m, n_bands, npol)
     CALL mp_sum(ps_c,intra_bgrp_comm)
     !
     CALL ZGEMM('N','N',npwx*npol,m,n_bands,alpha,proj_c,npwx*npol,ps_c,n_bands,(0.0_DP,0.0_DP),f,npwx*npol)
     !
     DEALLOCATE( ps_c )
     !
  ENDIF
  !
  CALL stop_clock ('alphapa')
  !
END SUBROUTINE
