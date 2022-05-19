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
SUBROUTINE destroy_pw_arrays( )
  !-----------------------------------------------------------------------
  !
  USE mp_global,              ONLY : my_image_id
  USE buffers,                ONLY : close_buffer
  USE control_flags,          ONLY : gamma_only
  USE xc_lib,                 ONLY : xclib_dft_is
  USE klist,                  ONLY : nks
  USE westcom,                ONLY : iuwfc
  !
  IMPLICIT NONE
  !
  CALL start_clock('des_pw_ar')
  !
  IF(my_image_id == 0) THEN
     IF(.NOT. (gamma_only .AND. nks == 1 .AND. .NOT. xclib_dft_is('hybrid'))) THEN
        CALL close_buffer(iuwfc,'DELETE')
     ENDIF
  ENDIF
  !
  CALL stop_clock('des_pw_ar')
  !
END SUBROUTINE
