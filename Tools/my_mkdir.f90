!
! Copyright (C) 2015-2024 M. Govoni
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
SUBROUTINE my_mkdir(dirname)
  !------------------------------------------------------------------------
  !
  USE mp,           ONLY : mp_barrier
  USE mp_world,     ONLY : mpime,root,world_comm
  USE forpy_mod,    ONLY : call_py_noret,import_py,module_py,tuple,tuple_create,dict,dict_create
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  CHARACTER(LEN=*), INTENT(IN) :: dirname
  !
  ! Workspace
  !
  INTEGER :: ierr
  TYPE(tuple) :: args
  TYPE(dict) :: kwargs
  TYPE(module_py) :: pymod
  !
  IF(mpime == root) THEN
     !
     ierr = import_py(pymod, "west_utils")
     !
     ierr = tuple_create(args, 1)
     ierr = args%setitem(0, TRIM(ADJUSTL(dirname)))
     ierr = dict_create(kwargs)
     ierr = call_py_noret(pymod, "my_mkdir", args, kwargs)
     !
     CALL kwargs%destroy
     CALL args%destroy
     CALL pymod%destroy
     !
  ENDIF
  !
  ! BARRIER
  !
  CALL mp_barrier(world_comm)
  !
END SUBROUTINE
