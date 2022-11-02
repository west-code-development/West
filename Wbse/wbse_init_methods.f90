!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE wbse_init_methods()
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : isk,nks,npw,ngk
  USE wavefunctions,        ONLY : evc
  USE westcom,              ONLY : lrwfc,iuwfc,ev,dvg,n_pdep_eigen,npwqx,nbnd_occ,&
                                 & wbse_init_calculation,which_spin_channel
  USE lsda_mod,             ONLY : nspin
  USE pdep_db,              ONLY : pdep_db_read
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE buffers,              ONLY : get_buffer
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : pert
  !
  IMPLICIT NONE
  !
  !LOGICAL, INTENT(IN) :: ff_activate
  !
  INTEGER :: iks, current_spin
  INTEGER :: iq, nkq, ikq
  REAL(DP):: xq(3)
  LOGICAL :: l_restart_calc, spin_resolve
  !
  SELECT CASE(wbse_init_calculation)
  CASE('r','R')
     !
     ! RESTART
     !
     l_restart_calc = .True.
     !
  CASE('s','S')
     !
     ! FROM SCRATCH
     !
     l_restart_calc = .False.
     !
  CASE DEFAULT
     CALL errore('Wbse_init', 'Wrong wbse_init_calculation', 1)
  END SELECT
  !
  !IF(.NOT. ff_activate) THEN
  !   !
  !   ! Activate band group parallel PDEP
  !   !
  !   pert = idistribute()
  !   CALL pert%init(n_pdep_eigen, 'B','nvecx',.TRUE.)
  !   !
  !ENDIF
  !
  ALLOCATE(dvg(npwqx,pert%nlocx))
  ALLOCATE(ev(n_pdep_eigen))
  !
  spin_resolve = (which_spin_channel > 0) .AND. (nspin > 1)
  !
  nkq = 1
  !
  DO iq = 1, nkq
     !
     xq(:) = 0._DP
     !
     !IF(.NOT. ff_activate) THEN
     !   CALL pdep_db_read(n_pdep_eigen)
     !ENDIF
     !
     DO iks = 1, nks
        !
        current_spin = isk(iks)
        ikq = 1 !grid_ikq (iq,iks)
        !
        ! ... Number of G vectors for PW expansion of wfs at k
        !
        npw = ngk(iks)
        !
        ! ... read in GS wavefunctions iks
        !
        IF(nkq > 1) THEN
           !IF (my_image_id==0) CALL get_buffer (evc, lrwfc, iuwfc, iks)
           !IF (my_image_id==0) CALL get_buffer (evq, lrwfc, iuwfc, ikq)
           !CALL mp_bcast(evc,0,inter_image_comm)
           !CALL mp_bcast(evq,0,inter_image_comm)
        ELSE
           IF(nks > 1) THEN
              IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
              CALL mp_bcast(evc,0,inter_image_comm)
           ENDIF
           !
           IF(spin_resolve) THEN
              IF(current_spin == which_spin_channel) THEN
                 !IF(ff_activate) THEN
                    CALL wbse_init_qboxcoupling_single_q(iks,iq,xq,current_spin,nbnd_occ(iks),l_restart_calc)
                 !ELSE
                 !   CALL wbse_init_pdep_single_q(iks,iq,xq,current_spin,nbnd_occ(iks),n_pdep_eigen,dvg,ev,l_restart_calc)
                 !ENDIF
              ENDIF
           ELSE
              !IF(ff_activate) THEN
                 CALL wbse_init_qboxcoupling_single_q(iks,iq,xq,current_spin,nbnd_occ(iks),l_restart_calc)
              !ELSE
              !   CALL wbse_init_pdep_single_q(iks,iq,xq,current_spin,nbnd_occ(iks),n_pdep_eigen,dvg,ev,l_restart_calc)
              !ENDIF
           ENDIF
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE
