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
SUBROUTINE exx_ungo()
  !-----------------------------------------------------------------------
  !
  USE exx,                    ONLY : deallocate_exx
  USE xc_lib,                 ONLY : xclib_dft_is,stop_exx
  USE command_line_options,   ONLY : command_line
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: is_westpp
  LOGICAL, EXTERNAL :: matches
  !
  is_westpp = matches('westpp.x',command_line)
  !
  IF(xclib_dft_is('hybrid') .AND. .NOT. is_westpp) THEN
     CALL stop_exx
     CALL deallocate_exx
  ENDIF
  !
END SUBROUTINE
