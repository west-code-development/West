!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE wbse_calc_dens(devc, drho)
  !---------------------------------------------------------------------
  !
  ! This subroutine calculates the response charge density
  ! from linear response orbitals and ground state orbitals.
  !
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE fft_base,               ONLY : dffts
  USE lsda_mod,               ONLY : nspin,lsda
  USE wavefunctions,          ONLY : psic,evc
  USE noncollin_module,       ONLY : npol
  USE pwcom,                  ONLY : npw,npwx,igk_k,current_k,nks,current_spin,isk,wg,ngk
  USE control_flags,          ONLY : gamma_only
  USE mp,                     ONLY : mp_sum,mp_bcast
  USE mp_global,              ONLY : my_image_id,inter_image_comm
  USE buffers,                ONLY : get_buffer
  USE westcom,                ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,l_lanczos
  USE fft_at_gamma,           ONLY : double_invfft_gamma
  USE fft_at_k,               ONLY : single_invfft_k
  USE distribution_center,    ONLY : aband
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: devc(npwx*npol,nbndval0x,nks)
  COMPLEX(DP), INTENT(OUT) :: drho(dffts%nnr,nspin)
  !
  ! Workspace
  !
  INTEGER :: ir, ibnd, iks, nbndval, lbnd
  REAL(DP) :: w1
  REAL(DP), ALLOCATABLE :: tmp_r(:)
  COMPLEX(DP), ALLOCATABLE :: tmp_c(:)
  !
  CALL start_clock('calc_dens')
  !
  IF(gamma_only) THEN
     ALLOCATE(tmp_r(dffts%nnr))
     !
     tmp_r(:) = 0._DP
  ELSE
     ALLOCATE(tmp_c(dffts%nnr))
     !
     drho(:,:) = (0._DP,0._DP)
  ENDIF
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     nbndval = nbnd_occ(iks)
     !
     ! ... Set k-point and spin
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read GS wavefunctions
     !
     IF(nks > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     IF(gamma_only) THEN
        !
        ! double bands @ gamma
        !
        DO lbnd = 1, aband%nloc
           !
           ibnd = aband%l2g(lbnd)
           IF(ibnd < 1 .OR. ibnd > nbndval) CYCLE
           !
           w1 = wg(ibnd,iks)/omega
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),devc(:,ibnd,iks),psic,'Wave')
           !
           DO ir = 1, dffts%nnr
              tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
           ENDDO
           !
        ENDDO
        !
        drho(:,current_spin) = CMPLX(tmp_r,KIND=DP)
        !
     ELSE
        !
        ! only single bands
        !
        DO lbnd = 1, aband%nloc
           !
           ibnd = aband%l2g(lbnd)
           IF(ibnd < 1 .OR. ibnd > nbndval) CYCLE
           !
           w1 = wg(ibnd,iks)/omega
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,current_k))
           CALL single_invfft_k(dffts,npw,npwx,devc(:,ibnd,iks),tmp_c,'Wave',igk_k(:,current_k))
           !
           DO ir = 1, dffts%nnr
              drho(ir,current_spin) = drho(ir,current_spin) + w1*CONJG(psic(ir))*tmp_c(ir)
           ENDDO
           !
        ENDDO
        !
        IF(npol == 2) THEN
           !
           DO lbnd = 1, aband%nloc
              !
              ibnd = aband%l2g(lbnd)
              IF(ibnd < 1 .OR. ibnd > nbndval) CYCLE
              !
              w1 = wg(ibnd,iks)/omega
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1:npwx*2,ibnd),psic,'Wave',igk_k(:,current_k))
              CALL single_invfft_k(dffts,npw,npwx,devc(npwx+1:npwx*2,ibnd,iks),tmp_c,'Wave',igk_k(:,current_k))
              !
              DO ir = 1, dffts%nnr
                 drho(ir,current_spin) = drho(ir,current_spin) + w1 * CONJG(psic(ir))*tmp_c(ir)
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  IF(l_lanczos) THEN
     CALL mp_sum(drho,inter_image_comm)
  ENDIF
  !
  IF(gamma_only) THEN
     DEALLOCATE(tmp_r)
  ELSE
     DEALLOCATE(tmp_c)
  ENDIF
  !
  CALL stop_clock('calc_dens')
  !
END SUBROUTINE
