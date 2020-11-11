!
! Copyright (C) 2015 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file: 
! He Ma
!----------------------------------------------------------------------------
MODULE cpp_wrappers
  !----------------------------------------------------------------------------
  !
  USE iso_c_binding
  !
  IMPLICIT NONE
  !
  !----------------------------------------------------------------------------
  INTERFACE 
    !----------------------------------------------------------------------------
    !
    SUBROUTINE c_sleep(n) bind(C) 
      IMPORT
      INTEGER(KIND=C_INT), INTENT(IN)   ::  n
    END SUBROUTINE
    !
    SUBROUTINE c_wait_for_file(max_seconds, file_exists, lockfilename) bind(C)
      IMPORT
      INTEGER(KIND=C_INT), INTENT(IN)    :: max_seconds
      LOGICAL(KIND=C_BOOL), INTENT(OUT)  :: file_exists
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lockfilename(*)
    END SUBROUTINE
    !
  END INTERFACE
  !
  !----------------------------------------------------------------------------
END MODULE