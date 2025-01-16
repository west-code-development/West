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
SUBROUTINE set_nbndocc()
  !-----------------------------------------------------------------------
  !
  ! current_spin is needed
  !
  USE kinds,                  ONLY : DP
  USE pwcom,                  ONLY : nbnd,nelup,neldw,isk,nelec,lsda,nks,lgauss,ltetra,wg,wk,&
                                   & tfixed_occ
  USE constants,              ONLY : degspin
  USE noncollin_module,       ONLY : noncolin
  USE westcom,                ONLY : nbndval0x,nbnd_occ,l_frac_occ,occupation,nbnd_occ_full,docc_thr
  USE control_flags,          ONLY : gamma_only
  USE mp_global,              ONLY : inter_pool_comm
  USE mp,                     ONLY : mp_max
  !
  IMPLICIT NONE
  !
  INTEGER :: spin,iks,ibnd
  !
  !
  IF(ltetra) CALL errore("set_nbndocc", "tetrahedral occupation not implemented", 1)
  !
  ! Determine if occupations are fractional
  !
  l_frac_occ = tfixed_occ .OR. lgauss
  !
  IF(l_frac_occ .AND. .NOT. gamma_only) THEN
     CALL errore("set_nbndocc", "fraction occupation only implemented for gamma-only case", 1)
  ENDIF
  !
  IF(ALLOCATED(occupation)) DEALLOCATE(occupation)
  IF(ALLOCATED(nbnd_occ)) DEALLOCATE(nbnd_occ)
  IF(ALLOCATED(nbnd_occ_full)) DEALLOCATE(nbnd_occ_full)
  !
  ! 2nd dimension of occupation, nbnd_occ, and nbnd_occ_full distributed across pools
  !
  ALLOCATE(occupation(nbnd,nks))
  ALLOCATE(nbnd_occ(nks))
  ALLOCATE(nbnd_occ_full(nks))
  !
  occupation(:,:) = 0._DP
  nbnd_occ_full(:) = 0
  nbnd_occ(:) = 0
  !
  IF(l_frac_occ) THEN
     !
     ! FRACTIONAL OCCUPATIONS
     !
     DO iks = 1,nks
        DO ibnd = 1,nbnd
           !
           occupation(ibnd,iks) = wg(ibnd,iks)/wk(iks)
           !
           IF(occupation(ibnd,iks) > docc_thr) nbnd_occ(iks) = MAX(nbnd_occ(iks),ibnd)
           IF(occupation(ibnd,iks) > (1-docc_thr)) nbnd_occ_full(iks) = MAX(nbnd_occ_full(iks),ibnd)
           !
        ENDDO
     ENDDO
     !
  ELSE
     !
     ! WHOLE OCCUPATIONS
     !
     IF(lsda) THEN
        !
        ! Collinear spin
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
        ! Noncollinear spin
        !
        nbnd_occ(:) = NINT(nelec)
        !
     ELSE
        !
        ! No spin
        !
        nbnd_occ(:) = NINT(nelec)/degspin
        !
     ENDIF
     !
     DO iks = 1,nks
        occupation(1:nbnd_occ(iks),iks) = 1._DP
     ENDDO
     nbnd_occ_full(:) = nbnd_occ
     !
  ENDIF
  !
  IF(MAXVAL(nbnd_occ) == 0) CALL errore("set_nbndocc", "nbnd_occ was not set", 1)
  !
  nbndval0x = MAXVAL(nbnd_occ)
  !
  CALL mp_max(nbndval0x,inter_pool_comm)
  !
END SUBROUTINE
