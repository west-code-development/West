!
! Copyright (C) 2015-2017 M. Govoni 
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
PURE FUNCTION to_upper_case (str) RESULT (string)
!-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   CHARACTER(*), INTENT(IN) :: str
   CHARACTER(LEN(str))      :: string
   !
   INTEGER :: ic, i
   !
   CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'
   !
   ! Capitalize each letter if it is low
   string = str
   DO i = 1, LEN_TRIM(str)
      ic = INDEX(low, str(i:i))
      IF (ic > 0) string(i:i) = cap(ic:ic)
   ENDDO
   !
END FUNCTION
!
!-----------------------------------------------------------------------
PURE FUNCTION to_lower_case (str) RESULT (string)
!-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   CHARACTER(*), INTENT(IN) :: str
   CHARACTER(LEN=LEN(str)), INTENT(OUT) :: string
   !
   INTEGER :: ic, i
   !
   CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'
   !
   ! Lower case each letter if it is capital 
   string = str
   DO i = 1, LEN_TRIM(str)
      ic = INDEX(cap, str(i:i))
      IF (ic > 0) string(i:i) = low(ic:ic)
   ENDDO
   !
END FUNCTION
