!
! Copyright (C) 2015-2023 M. Govoni
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
SUBROUTINE apply_alpha_pa_to_m_wfcs(iks,m,f,alpha)
  !-----------------------------------------------------------------------
  !
  ! | f_i > = alpha * P_a | f_i >            forall i = 1:n_bands
  ! P_a = sum_{i from 1 to n_bands} |proj_c_i><proj_c_i|
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : npw,npwx
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : npol
  USE westcom,              ONLY : n_bands,proj_c
#if defined(__CUDA)
  USE west_gpu,             ONLY : ps_r,ps_c
#if defined(__NCCL)
  USE west_gpu,             ONLY : gpu_sum,gpu_intra_bgrp_comm
#endif
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: iks,m
  COMPLEX(DP), INTENT(IN) :: alpha
  COMPLEX(DP), INTENT(INOUT) :: f(npwx*npol,m)
#if defined(__CUDA)
  ATTRIBUTES(DEVICE) :: f
#endif
  !
  ! Workspace
  !
  REAL(DP) :: alpha_r
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: ps_r(:,:)
  COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
#endif
  !
#if defined(__CUDA)
  CALL start_clock_gpu('alphapa')
#else
  CALL start_clock('alphapa')
#endif
  !
  ! ps_{ij} = < westcom/proj_c_i | f_j >
  !
  IF( gamma_only ) THEN
     !
     alpha_r = REAL(alpha,KIND=DP)
     !
#if !defined(__CUDA)
     ALLOCATE( ps_r(n_bands,m) )
     ps_r = 0.0_DP
#endif
     !
     !$acc host_data use_device(proj_c,ps_r)
     CALL glbrak_gamma( proj_c(:,:,iks), f, ps_r, npw, npwx, n_bands, m, n_bands, npol)
     !$acc end host_data
     !
#if defined(__NCCL)
     CALL gpu_sum(ps_r,n_bands*m,gpu_intra_bgrp_comm)
#else
     !$acc host_data use_device(ps_r)
     CALL mp_sum(ps_r,intra_bgrp_comm)
     !$acc end host_data
#endif
     !
     !$acc host_data use_device(proj_c,ps_r)
     CALL DGEMM('N','N',2*npwx*npol,m,n_bands,alpha_r,proj_c(1,1,iks),2*npwx*npol,ps_r,n_bands,0.0_DP,f,2*npwx*npol)
     !$acc end host_data
     !
#if !defined(__CUDA)
     DEALLOCATE( ps_r )
#endif
     !
  ELSE
     !
#if !defined(__CUDA)
     ALLOCATE( ps_c(n_bands,m) )
     ps_c = (0.0_DP,0.0_DP)
#endif
     !
     !$acc host_data use_device(proj_c,ps_c)
     CALL glbrak_k( proj_c(:,:,iks), f, ps_c, npw, npwx, n_bands, m, n_bands, npol)
     !$acc end host_data
     !
#if defined(__NCCL)
     CALL gpu_sum(ps_c,n_bands*m,gpu_intra_bgrp_comm)
#else
     !$acc host_data use_device(ps_c)
     CALL mp_sum(ps_c,intra_bgrp_comm)
     !$acc end host_data
#endif
     !
     !$acc host_data use_device(proj_c,ps_c)
     CALL ZGEMM('N','N',npwx*npol,m,n_bands,alpha,proj_c(1,1,iks),npwx*npol,ps_c,n_bands,(0.0_DP,0.0_DP),f,npwx*npol)
     !$acc end host_data
     !
#if !defined(__CUDA)
     DEALLOCATE( ps_c )
#endif
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('alphapa')
#else
  CALL stop_clock('alphapa')
#endif
  !
END SUBROUTINE
