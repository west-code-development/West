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
SUBROUTINE do_resp()
  !
  USE kinds,                 ONLY : DP
  USE io_push,               ONLY : io_push_title
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE wavefunctions,         ONLY : psic,evc
  USE pwcom,                 ONLY : npw,npwx,igk_k,current_k,ngk,wg
  USE control_flags,         ONLY : gamma_only
  USE mp,                    ONLY : mp_bcast
  USE mp_global,             ONLY : my_image_id,inter_image_comm
  USE buffers,               ONLY : get_buffer
  USE westcom,               ONLY : iuwfc,lrwfc,nbnd_occ,dvg_exc,westpp_n_liouville_to_use,&
                                  & westpp_range,westpp_save_dir
  USE fft_at_gamma,          ONLY : double_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE plep_db,               ONLY : plep_db_read
  USE distribution_center,   ONLY : pert
  USE class_idistribute,     ONLY : idistribute
  USE types_bz_grid,         ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ibnd,nbndval,iks,iexc,lexc
  INTEGER :: barra_load
  REAL(DP) :: w1
  REAL(DP), ALLOCATABLE :: rho(:)
  COMPLEX(DP), ALLOCATABLE :: psic_aux(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  !
  ! ... DISTRIBUTE
  !
  pert = idistribute()
  CALL pert%init(westpp_n_liouville_to_use,'i','nvec',.TRUE.)
  !
  ! READ EIGENVALUES AND VECTORS FROM OUTPUT
  !
  CALL plep_db_read(westpp_n_liouville_to_use)
  !
  ALLOCATE(rho(dffts%nnr))
  IF(.NOT. gamma_only) ALLOCATE(psic_aux(dffts%nnr))
  !
  CALL io_push_title('Charge Density Res(P)onse')
  !
  barra_load = 0
  DO lexc = 1,pert%nloc
     iexc = pert%l2g(lexc)
     IF(iexc < westpp_range(1) .OR. iexc > westpp_range(2)) CYCLE
     barra_load = barra_load+1
  ENDDO
  !
  CALL start_bar_type(barra,'westpp',barra_load*k_grid%nps*SUM(nbnd_occ))
  !
  DO lexc = 1,pert%nloc
     !
     ! local -> global
     !
     iexc = pert%l2g(lexc)
     IF(iexc < westpp_range(1) .OR. iexc > westpp_range(2)) CYCLE
     !
     DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
        !
        ! ... Set k-point, spin, kinetic energy, needed by Hpsi
        !
        current_k = iks
        npw = ngk(iks)
        nbndval = nbnd_occ(iks)
        !
        ! ... read in wavefunctions from the previous iteration
        !
        IF(k_grid%nps > 1) THEN
           IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc,0,inter_image_comm)
        ENDIF
        !
        rho(:) = 0._DP
        !
        DO ibnd = 1, nbndval
           !
           w1 = wg(ibnd,iks)/omega
           !
           IF(gamma_only) THEN
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),dvg_exc(:,ibnd,iks,lexc),psic,'Wave')
              !
              rho(:) = rho(:) + w1 * REAL(psic(:),KIND=DP)*AIMAG(psic(:))
              !
           ELSE
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,current_k))
              CALL single_invfft_k(dffts,npw,npwx,dvg_exc(:,ibnd,iks,lexc),psic_aux,'Wave',igk_k(:,current_k))
              !
              rho(:) = rho(:) + w1 * CONJG(psic(:)) * psic_aux(:)
              !
           ENDIF
           !
           CALL update_bar_type(barra,'westpp',1)
           !
        ENDDO
        !
        WRITE(fname,'(a,i6.6,a,i6.6)') TRIM(westpp_save_dir)//'/respK',1,'E',iexc
        CALL dump_r(rho,TRIM(fname))
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'westpp')
  !
  DEALLOCATE(rho)
  IF(.NOT. gamma_only) DEALLOCATE(psic_aux)
  !
END SUBROUTINE
