!
! Copyright (C) 2015-2025 M. Govoni
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
SUBROUTINE glbrak_gamma(a,b,c,ng,ngx,na,nb,ldc,np)
  !-----------------------------------------------------------------------
  !
  ! g -> gspace
  ! l -> local
  !
  ! c_ij = < a_i | b_j >    ( i = 1, na ; j = 1, nb )
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ng,ngx,np
  INTEGER, INTENT(IN) :: na,nb,ldc
  REAL(DP), INTENT(IN) :: a(2*ngx,na) ! complex passed as real for cublas compatibility
  REAL(DP), INTENT(IN) :: b(2*ngx,nb) ! same as a
  REAL(DP), INTENT(OUT) :: c(ldc,nb)
  !
  !$acc host_data use_device(a,b,c)
  CALL DGEMM('C','N',na,nb,2*ng,2._DP,a,2*ngx,b,2*ngx,0._DP,c,ldc)
  IF(gstart == 2) CALL DGER(na,nb,-1._DP,a,2*ngx,b,2*ngx,c,ldc)
  !$acc end host_data
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE glbrak_k(a,b,c,ng,ngx,na,nb,ldc,np)
  !-----------------------------------------------------------------------
  !
  ! g -> gspace
  ! l -> local
  !
  ! c_ij = < a_i | b_j >    ( i = 1, na ; j = 1, nb )
  !
  USE kinds,                ONLY : DP
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ng,ngx,np
  INTEGER, INTENT(IN) :: na,nb,ldc
  COMPLEX(DP), INTENT(IN) :: a(ngx*np,na)
  COMPLEX(DP), INTENT(IN) :: b(ngx*np,nb)
  COMPLEX(DP), INTENT(OUT) :: c(ldc,nb)
  !
  ! Workspace
  !
  COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
  COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
  !
  !$acc host_data use_device(a,b,c)
  CALL ZGEMM('C','N',na,nb,ng,one,a,ngx*np,b,ngx*np,zero,c,ldc)
  IF(np == 2) CALL ZGEMM('C','N',na,nb,ng,one,a(ngx+1,1),ngx*np,b(1+ngx,1),ngx*np,one,c,ldc)
  !$acc end host_data
  !
END SUBROUTINE
