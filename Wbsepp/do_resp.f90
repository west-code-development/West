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
  USE pwcom,                 ONLY : npw,npwx,igk_k,current_k,ngk,nks,wg
  USE control_flags,         ONLY : gamma_only
  USE mp,                    ONLY : mp_bcast
  USE mp_global,             ONLY : my_image_id,inter_image_comm
  USE buffers,               ONLY : get_buffer
  USE westcom,               ONLY : iuwfc,lrwfc,nbnd_occ,iexc_plot,westpp_format,dvg_exc,&
                                  & n_liouville_read_from_file,wbsepp_save_dir
  USE fft_at_gamma,          ONLY : double_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE plep_db,               ONLY : plep_db_read
  USE distribution_center,   ONLY : pert
  USE class_idistribute,     ONLY : idistribute
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ibnd,nbndval,iks
  REAL(DP) :: w1
  REAL(DP), ALLOCATABLE :: rho(:)
  COMPLEX(DP), ALLOCATABLE :: psic_aux(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  !
  westpp_format = 'C'
  !
  ! ... DISTRIBUTE
  !
  pert = idistribute()
  CALL pert%init(n_liouville_read_from_file,'i','nvec',.TRUE.)
  !
  ! READ EIGENVALUES AND VECTORS FROM OUTPUT
  !
  CALL plep_db_read(n_liouville_read_from_file)
  !
  ALLOCATE(rho(dffts%nnr))
  IF(.NOT. gamma_only) ALLOCATE(psic_aux(dffts%nnr))
  !
  rho(:) = 0._DP
  !
  CALL io_push_title('Charge Density (R)esponse')
  !
  CALL start_bar_type(barra,'wbsepp',SUM(nbnd_occ))
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     DO ibnd = 1, nbndval
        !
        w1 = wg(ibnd,iks)/omega
        !
        IF(gamma_only) THEN
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),dvg_exc(:,ibnd,iks,iexc_plot),psic,'Wave')
           !
           rho(:) = rho(:) + w1 * REAL(psic(:),KIND=DP)*AIMAG(psic(:))
           !
        ELSE
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,current_k))
           CALL single_invfft_k(dffts,npw,npwx,dvg_exc(:,ibnd,iks,iexc_plot),psic_aux,'Wave',igk_k(:,current_k))
           !
           rho(:) = rho(:) + w1 * CONJG(psic(:)) * psic_aux(:)
           !
        ENDIF
        !
        CALL update_bar_type(barra,'wbsepp',1)
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'wbsepp')
  !
  WRITE(fname,'(a,i6.6,a,i6.6)') TRIM(wbsepp_save_dir)//'/respK',1,'E',iexc_plot
  IF(my_image_id == 0) CALL dump_r(rho,TRIM(fname))
  !
  DEALLOCATE(rho)
  IF(.NOT. gamma_only) DEALLOCATE(psic_aux)
  !
END SUBROUTINE
