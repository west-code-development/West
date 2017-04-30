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
MODULE bar
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP 
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  SAVE
  !
  TYPE :: bar_type
     INTEGER :: tot_load
     INTEGER :: counter
     INTEGER :: printing_number
     REAL(DP) :: time_start
     REAL(DP) :: time_now
  END TYPE bar_type
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE start_bar_type(b,label,load)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bar_type) :: b
      INTEGER,INTENT(IN) :: load
      CHARACTER(*),INTENT(IN) :: label
      REAL(DP) :: time
      REAL(DP),EXTERNAL :: get_clock
      !
      CALL start_clock( label )
      time = get_clock( label )
      ! 
      b%tot_load = load
      b%counter = 0
      b%printing_number = 1
      b%time_start = time
      b%time_now = time
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE stop_bar_type(b,label)
    !-----------------------------------------------------------------------
      !
      USE io_push,   ONLY : io_push_bar
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bar_type) :: b
      CHARACTER(*),INTENT(IN) :: label
      !
      CALL io_push_bar()
      CALL stop_clock( label )
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE update_bar_type(b,label,add)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      TYPE(bar_type) :: b
      INTEGER,INTENT(IN) :: add
      CHARACTER(*),INTENT(IN) :: label
      REAL(DP) :: time
      REAL(DP),EXTERNAL :: get_clock
      !
      CALL start_clock( label )
      time = get_clock( label )
      ! 
      b%counter = b%counter + add
      b%time_now = time
      !
      CALL eventually_print_out_the_bar(b)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eventually_print_out_the_bar(b)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      TYPE(bar_type) :: b
      REAL(DP) :: work_percent
      INTEGER :: ihas
      CHARACTER(LEN=20) :: hashes_perc
      CHARACTER(LEN=20),EXTERNAL :: human_readable_time 
      REAL(DP) :: elapsed_time, expected_time
      !
      work_percent=REAL(b%counter*100,DP)/REAL(b%tot_load,DP)
      !
      IF( INT(work_percent/5.0_DP) .GE. b%printing_number ) THEN
         !
         hashes_perc=""
         DO ihas=1,INT(work_percent/5.0_DP)
            hashes_perc=TRIM(hashes_perc)//"#"
         ENDDO
         !
         elapsed_time=b%time_now-b%time_start
         expected_time=elapsed_time * 100.0_DP / (5.0_DP* INT(work_percent/5.0_DP) )
         WRITE(stdout, "(5X, 'in progress... |', a20, '| ', i3.3, '% :',a20,' (E)', a20, ' (X)')") &
                hashes_perc, INT(work_percent/5.0_DP)*5, &
                TRIM( human_readable_time(elapsed_time)), TRIM( human_readable_time(expected_time))
         !
         b%printing_number = INT(work_percent/5.0_DP)+1
         !
      ENDIF 
      !
    END SUBROUTINE
    !
END MODULE
