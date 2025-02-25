!
! Copyright (C) 2015-2025 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `LICENSE'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Marco Govoni
!
!-----------------------------------------------------------------------
MODULE conversions
   !-----------------------------------------------------------------------
   USE kinds, ONLY : DP
   IMPLICIT NONE
   CONTAINS
    !
    FUNCTION ltoa(l) RESULT(res)
       IMPLICIT NONE
       CHARACTER(:),ALLOCATABLE :: res
       LOGICAL,INTENT(IN) :: l
       CHARACTER(4) :: t="true"
       CHARACTER(5) :: f="false"
       IF( l ) THEN
          res = t
       ELSE
          res =f
       ENDIF
    END FUNCTION
    !
    FUNCTION itoa(i) RESULT(res)
       IMPLICIT NONE
       CHARACTER(:),ALLOCATABLE :: res
       INTEGER,INTENT(IN) :: i
       CHARACTER(RANGE(i)+2) :: tmp
       WRITE(tmp,'(I0)') i
       res = TRIM(tmp)
    END FUNCTION
    !
    FUNCTION dtoa(d) RESULT(res)
       IMPLICIT NONE
       CHARACTER(:),ALLOCATABLE :: res
       REAL(DP),INTENT(IN) :: d
       CHARACTER(14) :: tmp
       WRITE(tmp,'(ES14.6)') d
       res = TRIM(tmp)
    END FUNCTION
END MODULE
