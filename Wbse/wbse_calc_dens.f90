!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE wbse_calc_dens( devc, drho )
  !---------------------------------------------------------------------
  !
  ! This subroutine calculates the response charge density 
  ! from linear response orbitals and ground state orbitals.
  !
  !
  USE kinds,                  ONLY : dp
  USE io_global,              ONLY : stdout
  USE cell_base,              ONLY : omega
  USE fft_base,               ONLY : dffts,dfftp
  USE lsda_mod,               ONLY : nspin,lsda
  USE wavefunctions_module,   ONLY : psic, evc
  USE noncollin_module,       ONLY : npol
  USE pwcom,                  ONLY : npw,npwx,igk_k,current_k,ngk,nks,current_spin,isk,wg 
  USE control_flags,          ONLY : gamma_only
  USE mp,                     ONLY : mp_sum, mp_bcast
  USE mp_global,              ONLY : inter_pool_comm, intra_bgrp_comm,&
                                     inter_bgrp_comm, my_image_id, inter_image_comm 
  USE buffers,                ONLY : get_buffer
  USE wbsecom,                ONLY : nbndval0x, l_lanzcos
  USE westcom,                ONLY : iuwfc,lrwfc,nbnd_occ
  USE fft_at_gamma,           ONLY : single_invfft_gamma, double_invfft_gamma
  USE fft_at_k,               ONLY : single_fwfft_k,single_invfft_k
  USE distribution_center,    ONLY : aband
  !
  IMPLICIT NONE
  !
  COMPLEX(dp), INTENT(IN) :: devc(npwx*npol, nbndval0x, nks)
  COMPLEX(dp), INTENT(OUT):: drho(dffts%nnr, nspin)
  !
  ! Local variables
  !
  INTEGER       :: ir, ibnd, iks, nbndval
  INTEGER       :: nbvalloc, il1
  REAL(kind=dp) :: w1, prod
  COMPLEX(dp), ALLOCATABLE :: psic_aux(:)
  !
  CALL start_clock('wbse_calc_dens')
  !
  drho(:,:) =  (0.0_DP,0.0_DP)
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     nbndval = nbnd_occ(iks)
     !
     nbvalloc = 0
     DO il1 = 1, aband%nloc
        ibnd = aband%l2g(il1)
        IF( ibnd < 1 .OR. ibnd > nbndval ) CYCLE
        nbvalloc = nbvalloc + 1
     ENDDO
     !
     ! ... Set k-point and spin
     !
     current_k = iks
     ! 
     IF ( lsda ) current_spin = isk(iks)
     !
     ! ... read in GS wavefunctions from the dir
     !
     IF (nks>1) THEN
        !
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        CALL mp_bcast(evc,0,inter_image_comm)
        !
     ENDIF
     !
     IF (gamma_only) THEN
        !
        ! double bands @ gamma
        !
        DO il1 = 1, nbvalloc
           !
           ibnd = aband%l2g(il1)
           !
           w1 = wg(ibnd, iks)/omega
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),devc(1,ibnd,iks),psic,'Wave')
           !
           DO ir=1, dffts%nnr
              !
              prod =  REAL( psic(ir),KIND=DP) * DIMAG( psic(ir)) 
              !
              drho(ir,current_spin) = drho(ir,current_spin) + w1 * CMPLX( prod, 0.0_DP, KIND=DP)
              !
           ENDDO
           !
        ENDDO
        !
     ELSE
        !
        ALLOCATE (psic_aux(dffts%nnr))
        !
        ! only single bands
        !
        DO il1 = 1, nbvalloc
           !
           ibnd = aband%l2g(il1)
           ! 
           w1 = wg(ibnd, iks)/omega
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,current_k))
           !
           CALL single_invfft_k(dffts,npw,npwx,devc(1,ibnd,iks),psic_aux,'Wave',igk_k(1,current_k))
           ! 
           DO ir=1, dffts%nnr
              !
              drho(ir,current_spin) = drho(ir,current_spin) + w1 * DCONJG(psic(ir))* psic_aux(ir) 
              !
           ENDDO
           !
        ENDDO
        !
        IF (npol==2) THEN
           !
           DO il1 = 1, nbvalloc 
              !
              ibnd = aband%l2g(il1) 
              !
              w1 = wg(ibnd, iks)/omega
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ibnd),psic,'Wave',igk_k(1,current_k))
              !
              CALL single_invfft_k(dffts,npw,npwx,devc(npwx+1,ibnd,iks),psic_aux,'Wave',igk_k(1,current_k))
              !
              DO ir=1, dffts%nnr
                 !
                 drho(ir,current_spin) = drho(ir,current_spin) + w1 * DCONJG(psic(ir))* psic_aux(ir)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
        DEALLOCATE ( psic_aux )
        !
     ENDIF
     ! 
  ENDDO
  !
  IF (l_lanzcos) THEN
     !
     CALL mp_sum (drho(:,:),inter_image_comm)
     !
  ENDIF
  !
  CALL stop_clock('wbse_calc_dens')
  !
  RETURN
  !
ENDSUBROUTINE wbse_calc_dens
