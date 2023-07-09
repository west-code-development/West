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
SUBROUTINE apply_alpha_pc_to_m_wfcs(nbndval,m,f,alpha)
  !-----------------------------------------------------------------------
  !
  !       | f_i > =   alpha * P_c | f_i >                              forall i = 1:m
  ! i.e.  | f_i > = - alpha * P_v | f_i >  + alpha | f_i >             forall i = 1:m
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : npw,npwx
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : npol
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : evc=>evc_d
  USE west_gpu,             ONLY : ps_r,ps_c
#if defined(__NCCL)
  USE west_gpu,             ONLY : gpu_sum,gpu_intra_bgrp_comm
#endif
  USE cublas
#else
  USE wavefunctions,        ONLY : evc
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nbndval,m
  COMPLEX(DP), INTENT(IN) :: alpha
  COMPLEX(DP), INTENT(INOUT) :: f(npwx*npol,m)
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
  CALL start_clock_gpu('alphapc')
#else
  CALL start_clock('alphapc')
#endif
  !
  ! ps = < evc | f >
  !
  IF( gamma_only ) THEN
     !
     alpha_r = REAL(alpha,KIND=DP)
     !
#if !defined(__CUDA)
     ALLOCATE( ps_r(nbndval,m) )
     ps_r = 0.0_DP
#endif
     !
     !$acc host_data use_device(f,ps_r)
     CALL glbrak_gamma( evc, f, ps_r, npw, npwx, nbndval, m, nbndval, npol)
     !$acc end host_data
     !
#if defined(__NCCL)
     CALL gpu_sum(ps_r,nbndval*m,gpu_intra_bgrp_comm)
#else
     !$acc host_data use_device(ps_r)
     CALL mp_sum(ps_r,intra_bgrp_comm)
     !$acc end host_data
#endif
     !
     !$acc host_data use_device(ps_r,f)
     CALL DGEMM('N','N',2*npwx*npol,m,nbndval,-alpha_r,evc,2*npwx*npol,ps_r,nbndval,alpha_r,f,2*npwx*npol)
     !$acc end host_data
     !
#if !defined(__CUDA)
     DEALLOCATE( ps_r )
#endif
     !
  ELSE
     !
#if !defined(__CUDA)
     ALLOCATE( ps_c(nbndval,m) )
     ps_c = (0.0_DP,0.0_DP)
#endif
     !
     !$acc host_data use_device(f,ps_c)
     CALL glbrak_k( evc, f, ps_c, npw, npwx, nbndval, m, nbndval, npol)
     !$acc end host_data
     !
#if defined(__NCCL)
     CALL gpu_sum(ps_c,nbndval*m,gpu_intra_bgrp_comm)
#else
     !$acc host_data use_device(ps_c)
     CALL mp_sum(ps_c,intra_bgrp_comm)
     !$acc end host_data
#endif
     !
     !$acc host_data use_device(ps_c,f)
     CALL ZGEMM('N','N',npwx*npol,m,nbndval,-alpha,evc,npwx*npol,ps_c,nbndval,alpha,f,npwx*npol)
     !$acc end host_data
     !
#if !defined(__CUDA)
     DEALLOCATE( ps_c )
#endif
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('alphapc')
#else
  CALL stop_clock('alphapc')
#endif
  !
END SUBROUTINE
