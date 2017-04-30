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
  CALL DGEMM( 'C', 'N', na, nb, 2*ng, 2.0_DP, a, 2*ngx, b, 2*ngx, 0.0_DP, c, ldc )
  IF (gstart == 2 ) CALL DGER( na, nb, -1.0_DP, a, 2*ngx, b, 2*ngx, c, ldc )
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
  CALL ZGEMM('C','N',na,nb,ng,one,a(1,1)    ,ngx*np,b(1,1)    ,ngx*np,zero,c,ldc)
  IF(np==2) THEN 
     CALL ZGEMM('C','N',na,nb,ng,one,a(ngx+1,1),ngx*np,b(1+ngx,1),ngx*np,one ,c,ldc)
  ENDIF
  !
END SUBROUTINE 
