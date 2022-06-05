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
SUBROUTINE set_nbndocc()
  !-----------------------------------------------------------------------
  !
  ! current_spin is needed
  !
  USE kinds,                  ONLY : DP
  USE pwcom,                  ONLY : nelup,neldw,isk,nelec,lsda,nks,&
                                   & nbnd,lgauss,ltetra,wg,wk,tfixed_occ
  USE constants,              ONLY : degspin
  USE noncollin_module,       ONLY : noncolin
  USE westcom,                ONLY : nbnd_occ,l_frac_occ,occupation,&
                                   & nbnd_occ_full
  USE types_bz_grid,          ONLY : k_grid
  USE control_flags,          ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER :: spin,iks,ibnd
  REAL(DP) :: occ
  !
  ! Calculate NBNDVAL
  !
  IF(ALLOCATED(nbnd_occ)) DEALLOCATE(nbnd_occ)
  ALLOCATE(nbnd_occ(nks))
  !
  IF(lsda) THEN
     !
     DO iks = 1,nks
        !
        spin = isk(iks)
        !
        SELECT CASE(spin)
        CASE(1)
           nbnd_occ(iks) = NINT(nelup)
        CASE(2)
           nbnd_occ(iks) = NINT(neldw)
        END SELECT
        !
     ENDDO
     !
  ELSEIF(noncolin) THEN
     !
     nbnd_occ(:) = NINT(nelec)
     !
  ELSE
     !
     nbnd_occ(:) = NINT(nelec)/degspin
     !
  ENDIF
  !
  IF (ltetra) CALL errore("set_nbndocc", "tetrahedral occupation not implemented", 1)
  !
  IF (tfixed_occ .or. lgauss) THEN
     !
     l_frac_occ = .true.
     !
     IF(ALLOCATED(occupation)) DEALLOCATE(occupation)
     IF(ALLOCATED(nbnd_occ_full)) DEALLOCATE(nbnd_occ_full)
     IF(ALLOCATED(nbnd_occ)) DEALLOCATE(nbnd_occ)
     ALLOCATE( occupation(nbnd, k_grid%nps) )
     ALLOCATE( nbnd_occ_full(k_grid%nps) )
     ALLOCATE( nbnd_occ(k_grid%nps) )
     !
     occupation = 0._DP
     nbnd_occ_full = 0
     nbnd_occ = 0
     !
     DO iks = 1, k_grid%nps
        DO ibnd = 1, nbnd
           !
           occ = wg(ibnd,iks) / wk(iks)
           occupation(ibnd,iks) = occ
           !
           IF (occ > 0.999) nbnd_occ_full(iks) = ibnd
           IF (occ > 0.001) nbnd_occ(iks) = ibnd
           !
        ENDDO
     ENDDO
     !
  ELSE
     !
     l_frac_occ = .false.
     !
  ENDIF
  !
  IF (l_frac_occ .and. .not. gamma_only) THEN
     CALL errore("set_nbndocc", "fraction occupation only implemented for gamma-only case", 1)
  ENDIF
  !
END SUBROUTINE 
