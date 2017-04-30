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
SUBROUTINE report_dynamical_memory() 
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE iso_c_binding,        ONLY : c_int
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER(KIND=c_int) :: kilobytes
  !
  CALL memstat(kilobytes)
  IF(kilobytes>0) THEN
     WRITE(stdout,"(5X, 'per-process dynamical memory ',f7.1,' Mb')") kilobytes/1000.0
     WRITE(stdout,*) ' '
  ENDIF
  !
END SUBROUTINE
