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
SUBROUTINE my_mkdir( dirname )
  !------------------------------------------------------------------------
  !
  USE wrappers,  ONLY : f_mkdir_safe
  USE mp,        ONLY : mp_barrier,mp_bcast
  USE mp_world,  ONLY : mpime, root, world_comm
  USE io_files,  ONLY : check_writable
  !
  ! I/O
  !
  CHARACTER(LEN=*), INTENT(IN) :: dirname
  !
  ! Workspace
  !
  INTEGER                    :: ierr
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  ! BARRIER
  !
  !
  IF ( mpime == root ) ierr = f_mkdir_safe( TRIM( dirname ) )
  CALL mp_bcast ( ierr, root, world_comm )
  !
  CALL errore( 'create_directory', &
               'unable to create directory ' // TRIM( dirname ), ierr )
  !
  ! ... check whether the scratch directory is writable
  !
  IF ( mpime == root ) ierr = check_writable ( dirname, mpime )
  CALL mp_bcast( ierr, root, world_comm )
  !
  CALL errore( 'create_directory:', &
               TRIM( dirname ) // ' non existent or non writable', ierr )
  CALL mp_barrier( world_comm )
  !
  RETURN
  !
END SUBROUTINE
!
!
!
SUBROUTINE my_rmdir( dirname )
  !
  USE mp_world,             ONLY : root,mpime,world_comm
  USE mp,                   ONLY : mp_barrier,mp_bcast
  USE wrappers,             ONLY : f_rmdir
  !
  IMPLICIT NONE
  ! 
  ! I/O
  ! 
  CHARACTER(LEN=*), INTENT(IN) :: dirname
  !
  ! Workspace
  !
  INTEGER :: ierr,ip
  CHARACTER(6) :: my_label
  CHARACTER(30) :: fname
  !
  ! BARRIER
  !
  CALL mp_barrier( world_comm )
  !
  ! ... clear the directory
  !
  IF(mpime==root) THEN
     ierr =  f_rmdir( TRIM( dirname ) )
  ENDIF
  !
  CALL mp_bcast( ierr, root, world_comm)
  !
  CALL errore( 'rm_directory', 'cannot rm dir', ierr )
  !
  RETURN
  !
END SUBROUTINE
