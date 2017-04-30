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
SUBROUTINE set_nbndocc( )
  !-----------------------------------------------------------------------
  !
  ! current_spin is needed
  !
  USE kinds,                  ONLY : DP
  USE pwcom,                  ONLY : nbnd,nelup,neldw,isk,nelec,nspin,lsda,nks
  USE constants,              ONLY : degspin
  USE noncollin_module,       ONLY : noncolin,npol
  USE io_global,              ONLY : stdout
  USE westcom,                ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  INTEGER :: spin,iks
  !
  ! Calculate NBNDVAL
  !
  IF(ALLOCATED(nbnd_occ)) DEALLOCATE(nbnd_occ)
  ALLOCATE( nbnd_occ(nks) )
  !
  IF(lsda) THEN  
     !
     DO iks = 1, nks
        spin = isk(iks)
        !
        SELECT CASE(spin)
        CASE(1)
           nbnd_occ(iks) = NINT (nelup) 
        CASE(2)
           nbnd_occ(iks) = NINT (neldw)
        END SELECT
        !
     ENDDO
     !
  ELSEIF(noncolin) THEN
     !
     nbnd_occ(:) = NINT( nelec )
     !
  ELSE
     !
     nbnd_occ(:) = NINT( nelec ) / degspin
     !
  ENDIF
  !
END SUBROUTINE 
