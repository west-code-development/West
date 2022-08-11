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
SUBROUTINE glbrak_gamma( a, b, c, ng, ngx, na, nb, ldc, np)
!------------------------------------------------------------------------
  !
  ! g -> gspace
  ! l -> local
  !
  ! c_ij = < a_i | b_j >    ( i = 1, na ; j = 1, nb )
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ng, ngx, np
  INTEGER, INTENT(IN) :: na,nb,ldc
  COMPLEX(DP), INTENT(IN) :: a(ngx,na)
  COMPLEX(DP), INTENT(IN) :: b(ngx,nb)
  REAL(DP), INTENT(OUT) :: c(ldc,nb)
  !
  CALL DGEMM( 'C', 'N', na, nb, 2*ng, 2._DP, a, 2*ngx, b, 2*ngx, 0._DP, c, ldc )
  IF (gstart == 2 ) CALL DGER( na, nb, -1._DP, a, 2*ngx, b, 2*ngx, c, ldc )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE glbrak_k( a, b, c, ng, ngx, na, nb, ldc, np )
!------------------------------------------------------------------------
  !
  ! g -> gspace
  ! l -> local
  !
  ! c_ij = < a_i | b_j >    ( i = 1, na ; j = 1, nb )
  !
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ng, ngx, np
  INTEGER, INTENT(IN) :: na,nb,ldc
  COMPLEX(DP), INTENT(IN) :: a(ngx*np,na)
  COMPLEX(DP), INTENT(IN) :: b(ngx*np,nb)
  COMPLEX(DP), INTENT(OUT) :: c(ldc,nb)
  !
  ! Workspace
  !
  COMPLEX(DP),PARAMETER :: zero=(0._DP,0._DP), one=(1._DP,0._DP)
  !
  CALL ZGEMM('C','N',na,nb,ng,one,a(1,1),ngx*np,b(1,1),ngx*np,zero,c,ldc)
  IF(np==2) THEN
     CALL ZGEMM('C','N',na,nb,ng,one,a(ngx+1,1),ngx*np,b(1+ngx,1),ngx*np,one,c,ldc)
  ENDIF
  !
END SUBROUTINE
!
#if defined(__CUDA)
!-----------------------------------------------------------------------
SUBROUTINE glbrak_gamma_gpu(a_d,b_d,c_d,ng,ngx,na,nb,ldc,np)
!------------------------------------------------------------------------
  !
  ! g -> gspace
  ! l -> local
  !
  ! c_ij = < a_i | b_j >    ( i = 1, na ; j = 1, nb )
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE cublas
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ng,ngx,np
  INTEGER, INTENT(IN) :: na,nb,ldc
  REAL(DP), DEVICE, INTENT(IN) :: a_d(2*ngx,na) ! COMPLEX passed as REAL
  REAL(DP), DEVICE, INTENT(IN) :: b_d(2*ngx,nb) ! COMPLEX passed as REAL
  REAL(DP), DEVICE, INTENT(OUT) :: c_d(ldc,nb)
  !
  CALL DGEMM('C','N',na,nb,2*ng,2._DP,a_d,2*ngx,b_d,2*ngx,0._DP,c_d,ldc)
  !
  IF(gstart == 2) THEN
     CALL DGER(na,nb,-1._DP,a_d,2*ngx,b_d,2*ngx,c_d,ldc)
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE glbrak_k_gpu(a_d,b_d,c_d,ng,ngx,na,nb,ldc,np)
!------------------------------------------------------------------------
  !
  ! g -> gspace
  ! l -> local
  !
  ! c_ij = < a_i | b_j >    ( i = 1, na ; j = 1, nb )
  !
  USE kinds,                ONLY : DP
  USE cublas
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ng,ngx,np
  INTEGER, INTENT(IN) :: na,nb,ldc
  COMPLEX(DP), DEVICE, INTENT(IN) :: a_d(ngx*np,na)
  COMPLEX(DP), DEVICE, INTENT(IN) :: b_d(ngx*np,nb)
  COMPLEX(DP), DEVICE, INTENT(OUT) :: c_d(ldc,nb)
  !
  ! Workspace
  !
  COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
  COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
  !
  CALL ZGEMM('C','N',na,nb,ng,one,a_d,ngx*np,b_d,ngx*np,zero,c_d,ldc)
  !
  IF(np == 2) THEN
     CALL ZGEMM('C','N',na,nb,ng,one,a_d(ngx+1,1),ngx*np,b_d(ngx+1,1),ngx*np,one,c_d,ldc)
  ENDIF
  !
END SUBROUTINE
#endif
