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
SUBROUTINE clean_scratchfiles()
  !-----------------------------------------------------------------------
  !
  USE buffers,       ONLY : close_buffer
!  USE io_files,      ONLY : iunigk
  USE mp_world,      ONLY : world_comm
  USE mp_global,     ONLY : my_image_id
  USE mp,            ONLY : mp_barrier 
  USE westcom,       ONLY : iuwfc
  !
  IMPLICIT NONE
  !
  LOGICAL :: opnd
  !
  IF (my_image_id==0) THEN
     CALL close_buffer ( iuwfc, 'DELETE' )
!     INQUIRE( UNIT = iunigk, OPENED = opnd )
!     IF ( opnd ) CLOSE( UNIT = iunigk, STATUS = 'DELETE' )
  ENDIF
  !
  CALL mp_barrier(world_comm)
  !
END SUBROUTINE 
