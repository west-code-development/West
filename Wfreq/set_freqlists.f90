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
SUBROUTINE set_freqlists()
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine sets both the real and im freq lists 
  !
  USE kinds,                ONLY : DP 
  USE westcom,              ONLY : ecut_imfreq,imfreq_list,n_imfreq,frequency_list_power,ecut_refreq,refreq_list,n_refreq
  USE distribution_center,  ONLY : ifr,rfr,aband
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: ifreq,glob_ifreq
  !
  REAL(DP) :: p,y
  !
  ! Im Freq
  !
  IF(ALLOCATED(imfreq_list)) DEALLOCATE(imfreq_list)
  ALLOCATE(imfreq_list(ifr%nloc))
  !
  ! freq_n = p * ( y_n ** (-1/alpha) -1 )
  !    y_n = 1 - (i-1)/N
  !      p = Ec / ( N**(1/alpha)-1 )
  !
  p = ecut_imfreq / ( REAL( n_imfreq, KIND=DP)**(1._DP/frequency_list_power) - 1._DP ) 
  !
  DO ifreq = 1, ifr%nloc
     glob_ifreq = ifr%l2g(ifreq)
     !
     y = 1._DP - REAL( glob_ifreq-1, KIND=DP) / REAL(n_imfreq, KIND=DP)
     imfreq_list(ifreq) = p * ( y**(-1._DP/frequency_list_power) - 1._DP )
     !
  ENDDO
  !
  ! Re Freq
  !
  IF(ALLOCATED(refreq_list)) DEALLOCATE(refreq_list)
  ALLOCATE(refreq_list(rfr%nloc))
  !
  DO ifreq = 1, rfr%nloc
     glob_ifreq = rfr%l2g(ifreq)
     !
     refreq_list(ifreq) = ecut_refreq / REAL(n_refreq-1,KIND=DP) * REAL(glob_ifreq-1,KIND=DP) 
     !
  ENDDO
  !
END SUBROUTINE
