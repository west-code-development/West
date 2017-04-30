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
FUNCTION human_readable_time(time)
  !-----------------------------------------------------------------------
  !
  ! ... Given a time in seconds, the result is :
  ! ...  9999d-23h-59m-59.9s
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE 
  ! 
  ! I/O
  !
  REAL(DP), INTENT(IN) :: time
  CHARACTER(20) :: human_readable_time
  !
  ! Workspace
  !
  CHARACTER(20) :: ds_temp,s_temp,m_temp,h_temp,d_temp
  REAL(DP) :: seconds
  INTEGER :: minutes,hours,days 
  !
  ! If < 0.1 s
  !
  IF(time<0.1_DP) THEN
     human_readable_time="< 00.1s"
     RETURN
  ENDIF
  !
  ! If seconds
  !
  IF(time<60.0_DP) THEN
     seconds=time
     WRITE(s_temp,'(i2.2)') INT(seconds)
     WRITE(ds_temp,'(i1.1)') INT(seconds*10.0 - 10.0*INT(seconds))
     human_readable_time=TRIM(ADJUSTL(s_temp))//"."//TRIM(ADJUSTL(ds_temp))//"s"
     RETURN
  ENDIF
  !
  ! If minutes
  !
  IF(time<3600.0_DP) THEN
     minutes=INT(time/60.0_DP)
     WRITE(m_temp,'(i2.2)') minutes
     seconds=time-minutes*60.0_DP
     WRITE(s_temp,'(i2.2)') INT(seconds)
     WRITE(ds_temp,'(i1.1)') INT(seconds*10.0_DP - 10.0_DP*INT(seconds))
     human_readable_time=TRIM(ADJUSTL(m_temp))//"m-"//TRIM(ADJUSTL(s_temp))//"."//TRIM(ADJUSTL(ds_temp))//"s"
     RETURN
  ENDIF
  !
  ! If hours
  !
  IF(time<86400.0_DP) THEN
     hours=INT(time/3600.0_DP)
     WRITE(h_temp,'(i2.2)') hours
     minutes=INT((time-hours*3600.0_DP)/60.0_DP)
     WRITE(m_temp,'(i2.2)') minutes
     seconds=time-hours*3600.0_DP-minutes*60.0_DP
     WRITE(s_temp,'(i2.2)') INT(seconds)
     WRITE(ds_temp,'(i1.1)') INT(seconds*10.0_DP - 10.0_DP*INT(seconds))
     human_readable_time=TRIM(ADJUSTL(h_temp))//"h-"//TRIM(ADJUSTL(m_temp))//"m-"//TRIM(ADJUSTL(s_temp))//&
     "."//TRIM(ADJUSTL(ds_temp))//"s"
     RETURN
  ENDIF
  !
  ! If days
  !
  days=INT(time/86400.0_DP)
  WRITE(d_temp,'(i5)') days
  hours=INT((time-days*86400.0_DP)/3600.0_DP)
  WRITE(h_temp,'(i2.2)') hours
  minutes=INT((time-days*86400.0_DP-hours*3600.0_DP)/60.0_DP)
  WRITE(m_temp,'(i2.2)') minutes
  seconds=time-days*86400.0_DP-hours*3600.0_DP-minutes*60.0_DP
  WRITE(s_temp,'(i2.2)') INT(seconds)
  WRITE(ds_temp,'(i1.1)') INT(seconds*10.0_DP - 10.0_DP*INT(seconds))
  human_readable_time=&
  TRIM(ADJUSTL(d_temp))//"d-"//TRIM(ADJUSTL(h_temp))//"h-"//TRIM(ADJUSTL(m_temp))//"m-"//TRIM(ADJUSTL(s_temp))//&
  "."//TRIM(ADJUSTL(ds_temp))//"s"
  !
END FUNCTION
