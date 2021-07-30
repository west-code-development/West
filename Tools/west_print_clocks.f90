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
SUBROUTINE west_print_clocks( )
  !----------------------------------------------------------------------------
  !
  USE json_module,   ONLY : json_file
  USE kinds,         ONLY : DP
  USE mytime,        ONLY : nclock, clock_label, cputime, walltime, notrunning, &
                            called, t0cpu, t0wall, f_tcpu, f_wall
  USE westcom,       ONLY : logfile
  USE mp_world,      ONLY : mpime, root
  !
  IMPLICIT NONE
  !
  INTEGER :: n
  REAL(DP) :: elapsed_cpu_time, elapsed_wall_time
  TYPE(json_file) :: json
  INTEGER :: iunit, nmax
  CHARACTER(20), EXTERNAL :: human_readable_time
  !
  IF( mpime == root ) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
  !
  DO n = 1, nclock
     !
     CALL print_this_clock( n )
     !
     IF ( t0cpu(n) == notrunning ) THEN
        !
        ! ... clock stopped, print the stored value for the cpu time
        !
        elapsed_cpu_time = cputime(n)
        elapsed_wall_time= walltime(n)
        !
     ELSE
        !
        ! ... clock not stopped, print the current value of the cpu time
        !
        elapsed_cpu_time   = cputime(n) + f_tcpu() - t0cpu(n)
        elapsed_wall_time  = walltime(n) + f_wall() - t0wall(n)
        called(n)  = called(n) + 1
        !
     ENDIF
     !
     nmax = called(n)
     !
     IF( mpime == root ) THEN
        CALL json%add('timing.'//TRIM(clock_label(n))//'.cpu:sec',elapsed_cpu_time)
        CALL json%add('timing.'//TRIM(clock_label(n))//'.cpu:hum',TRIM(human_readable_time(elapsed_cpu_time)))
        CALL json%add('timing.'//TRIM(clock_label(n))//'.wall:sec',elapsed_wall_time)
        CALL json%add('timing.'//TRIM(clock_label(n))//'.wall:hum',TRIM(human_readable_time(elapsed_wall_time)))
        CALL json%add('timing.'//TRIM(clock_label(n))//'.nocalls',nmax)
     ENDIF
     !
  ENDDO
  !
  IF( mpime == root ) THEN
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE get_clock_called( label, ncalls )
  !----------------------------------------------------------------------------
  !
  USE mytime,        ONLY : nclock, clock_label, called
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: label
  INTEGER, INTENT(OUT) :: ncalls
  !
  INTEGER :: n
  !
  ncalls = 0
  !
  DO n = 1, nclock
     !
     IF( label /= clock_label(n) ) CYCLE
     !
     ncalls = called(n)
     EXIT
     !
  ENDDO
  !
END SUBROUTINE
