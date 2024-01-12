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
SUBROUTINE parse_command_arguments( )
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : main_input_file !, main_output_file
  !
  IMPLICIT NONE
  !
  INTEGER :: nargs, iiarg
  LOGICAL :: ifound !, ofound
  CHARACTER(LEN=512) :: string
  !
  ifound = .FALSE.
  !ofound = .FALSE.
  !
  nargs = command_argument_count()
  string = ' '
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL get_command_argument( iiarg, string )
     !
     IF ( .NOT. ifound .AND. TRIM( string ) == '-i' ) THEN
        CALL get_command_argument( ( iiarg + 1 ) , main_input_file )
        ifound =.TRUE.
     ENDIF
     !
     !IF ( .NOT. ofound .AND. TRIM( string ) == '-o' ) THEN
     !   CALL get_command_argument( ( iiarg + 1 ) , main_output_file )
     !   ofound =.TRUE.
     !ENDIF
     !
     !IF( ifound .AND. ofound ) EXIT
     IF( ifound ) EXIT
     !
  ENDDO
  !
  IF( .NOT. ifound ) CALL errore('parse_cmd','Cannot find input file (-i)',1)
  !IF( .NOT. ofound ) CALL errore('parse_cmd','Cannot find output file (-o)',1)
  !
END SUBROUTINE
