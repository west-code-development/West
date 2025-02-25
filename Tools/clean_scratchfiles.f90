!
! Copyright (C) 2015-2025 M. Govoni
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
SUBROUTINE clean_scratchfiles()
  !-----------------------------------------------------------------------
  !
  USE buffers,       ONLY : close_buffer
  USE mp_global,     ONLY : my_image_id
  USE westcom,       ONLY : iuwfc
  !
  IMPLICIT NONE
  !
  IF (my_image_id==0) THEN
     CALL close_buffer ( iuwfc, 'DELETE' )
  ENDIF
  !
END SUBROUTINE
