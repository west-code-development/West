!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE apply_hqp_to_m_wfcs(iks,m,f,g)
  !-----------------------------------------------------------------------
  !
  ! g = S|evc><evc|f>*(et_qp - et - Delta) + Delta|f>
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : npw,npwx,nbnd,et
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : npol
  USE westcom,              ONLY : et_qp
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : evc=>evc_d
  USE west_gpu,             ONLY : ps_r,ps_c
  USE cublas
#else
  USE wavefunctions,        ONLY : evc
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: iks, m
  COMPLEX(DP), INTENT(IN) :: f(npwx*npol,m)
  COMPLEX(DP), INTENT(INOUT) :: g(npwx*npol,m)
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,ig
  REAL(DP) :: delta
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: ps_r(:,:)
  COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
#endif
  !
#if defined(_CUDA)
  CALL start_clock_gpu('hqp')
#else
  CALL start_clock('hqp')
#endif
  !
  !$acc enter data copyin(et,et_qp)
  !
  delta = et_qp(nbnd,iks)-et(nbnd,iks)
  !
  ! ps = < evc | f >
  !
  IF(gamma_only) THEN
     !
     !$acc host_data use_device(f,ps_r)
#if defined(__CUDA)
     CALL glbrak_gamma_gpu(evc,f,ps_r,npw,npwx,nbnd,m,nbnd,npol)
#else
     ALLOCATE(ps_r(nbnd,m))
     ps_r = 0.0_DP
     !
     CALL glbrak_gamma(evc,f,ps_r,npw,npwx,nbnd,m,nbnd,npol)
#endif
     !
     CALL mp_sum(ps_r,intra_bgrp_comm)
     !$acc end host_data
     !
     !$acc parallel loop collapse(2) present(ps_r,et_qp,et)
     DO ibnd = 1,m
        DO jbnd = 1,nbnd
           ps_r(jbnd,ibnd) = ps_r(jbnd,ibnd)*(et_qp(jbnd,iks)-et(jbnd,iks)-delta)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(ps_r,g)
     CALL DGEMM('N','N',2*npwx*npol,m,nbnd,1.0_DP,evc,2*npwx*npol,ps_r,nbnd,1.0_DP,g,2*npwx*npol)
     !$acc end host_data
     !
  ELSE
     !
     !$acc host_data use_device(f,ps_c)
#if defined(__CUDA)
     CALL glbrak_k_gpu(evc,f,ps_c,npw,npwx,nbnd,m,nbnd,npol)
#else
     ALLOCATE(ps_c(nbnd,m))
     ps_c = (0.0_DP,0.0_DP)
     !
     CALL glbrak_k(evc,f,ps_c,npw,npwx,nbnd,m,nbnd,npol)
#endif
     !
     CALL mp_sum(ps_c,intra_bgrp_comm)
     !$acc end host_data
     !
     !$acc parallel loop collapse(2) present(ps_c,et_qp,et)
     DO ibnd = 1,m
        DO jbnd = 1,nbnd
           ps_c(jbnd,ibnd) = ps_c(jbnd,ibnd)*(et_qp(jbnd,iks)-et(jbnd,iks)-delta)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(ps_c,g)
     CALL ZGEMM('N','N',npwx*npol,m,nbnd,(1.0_DP,0.0_DP),evc,npwx*npol,ps_c,nbnd,(1.0_DP,0.0_DP),g,npwx*npol)
     !$acc end host_data
     !
#if !defined(__CUDA)
     DEALLOCATE(ps_c)
#endif
     !
  ENDIF
  !
  !$acc parallel loop collapse(2) present(g,f)
  DO ibnd = 1,m
     DO ig = 1,npwx*npol
        g(ig,ibnd) = g(ig,ibnd)+delta*f(ig,ibnd)
     ENDDO
  ENDDO
  !$acc end parallel
  !
  !$acc exit data delete(et,et_qp)
  !
#if defined(_CUDA)
  CALL stop_clock_gpu('hqp')
#else
  CALL stop_clock('hqp')
#endif
  !
END SUBROUTINE
