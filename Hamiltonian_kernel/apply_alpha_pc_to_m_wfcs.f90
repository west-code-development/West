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
SUBROUTINE apply_alpha_pc_to_m_wfcs(nbndval,m,f,alpha)
!------------------------------------------------------------------------
  !
  !       | f_i > =   alpha * P_c | f_i >                              forall i = 1:m
  ! i.e.  | f_i > = - alpha * P_v | f_i >  + alpha | f_i >             forall i = 1:m
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : npw,npwx
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE wavefunctions,        ONLY : evc
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : npol
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
  REAL(DP), ALLOCATABLE :: ps_r(:,:)
  COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
  !
  CALL start_clock ('alphapc')
  !
  ! ps = < evc | f >
  !
  IF( gamma_only ) THEN
     !
     ALLOCATE( ps_r(nbndval,m) )
     ps_r = 0.0_DP
     alpha_r = REAL(alpha,KIND=DP)
     !
     CALL glbrak_gamma( evc, f, ps_r, npw, npwx, nbndval, m, nbndval, npol)
     CALL mp_sum(ps_r,intra_bgrp_comm)
     !
     CALL DGEMM('N','N',2*npwx*npol,m,nbndval,-alpha_r,evc,2*npwx*npol,ps_r,nbndval,alpha_r,f,2*npwx*npol)
     !
     DEALLOCATE( ps_r )
     !
  ELSE
     !
     ALLOCATE( ps_c(nbndval,m) )
     ps_c = (0.0_DP,0.0_DP)
     !
     CALL glbrak_k( evc, f, ps_c, npw, npwx, nbndval, m, nbndval, npol)
     CALL mp_sum(ps_c,intra_bgrp_comm)
     !
     CALL ZGEMM('N','N',npwx*npol,m,nbndval,-alpha,evc,npwx*npol,ps_c,nbndval,alpha,f,npwx*npol)
     !
     DEALLOCATE( ps_c )
     !
  ENDIF
  !
  CALL stop_clock ('alphapc')
  !
END SUBROUTINE
