!
! Copyright (C) 2015-2022 M. Govoni
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
  !
  IMPLICIT NONE
  !
  IF( xclib_dft_is('hybrid') ) THEN
     CALL stop_exx
     CALL deallocate_exx
  ENDIF
  !
END SUBROUTINE
