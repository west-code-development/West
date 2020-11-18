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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
SUBROUTINE wbsepp_plot_charged_density_res_exc()
  !
  ! ... This pp reads eig-values and -vectors from davidson diago
  ! ... and project them on KS empty eig-vectors
  !
  USE kinds,                  ONLY : dp
  USE io_global,              ONLY : stdout
  USE cell_base,              ONLY : omega, at, alat
  USE fft_base,               ONLY : dffts,dfftp
  USE lsda_mod,               ONLY : nspin,lsda
  USE wavefunctions_module,   ONLY : psic, evc
  USE noncollin_module,       ONLY : npol
  USE pwcom,                  ONLY : npw,npwx,igk_k,current_k,ngk,nks,current_spin,isk
  USE wvfct,                  ONLY : wg
  USE control_flags,          ONLY : gamma_only
  USE mp,                     ONLY : mp_sum, mp_bcast, mp_min
  USE mp_global,              ONLY : me_bgrp, inter_pool_comm, intra_bgrp_comm,&
                                     inter_bgrp_comm, my_image_id, inter_image_comm
  USE buffers,                ONLY : get_buffer
  USE westcom,                ONLY : iuwfc,lrwfc,nbnd_occ,ev
  USE fft_at_gamma,           ONLY : single_invfft_gamma, double_invfft_gamma
  USE fft_at_k,               ONLY : single_fwfft_k,single_invfft_k
  !wbsecom combined into westcom
  !USE wbsecom,                ONLY : dvg_exc, n_plep_read_from_file
  USE westcom,              ONLY : iexc_plot, r0_input
  USE plep_db,                ONLY : plep_db_read
  USE distribution_center,    ONLY : pert
  USE class_idistribute,      ONLY : idistribute
  USE westcom,                ONLY : westpp_format, dvg_exc, n_plep_read_from_file
  USE ions_base,              ONLY : nat,tau
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER     :: ibnd, nbndval, nvec,count
  INTEGER     :: ir, i, j, k, ip, iks, index0, ir_end
  REAL(DP)    :: drmin_g
  REAL(DP)    :: rcoeff, w1
  REAL(DP)    :: inv_nr1,inv_nr2,inv_nr3
  COMPLEX(DP) :: zcoeff
  REAL(DP), ALLOCATABLE :: r(:,:),dr(:), rho_out(:,:), segno(:)
  COMPLEX(DP), ALLOCATABLE :: psic_aux(:)
  INTEGER, PARAMETER :: n_ipol = 3
  CHARACTER(LEN=6)   :: my_label
  CHARACTER(LEN=256) :: fname
  !
  INTEGER  :: iexc
  !
  INTEGER  :: n_point, ni
  REAL(DP) :: dstep, ri1, ri0, dx,dy,dz, drr, summ, summ0, vi1, vi0
  !
  iexc = iexc_plot
  !
  ! ... DISTRIBUTE nvec
  !
  nvec = n_plep_read_from_file
  pert = idistribute()
  CALL pert%init(nvec,'i','nvec',.TRUE.)
  CALL wbse_memory_report()
  !
  ! READ EIGENVALUES AND VECTORS FROM OUTPUT
  !
  CALL plep_db_read( n_plep_read_from_file )
  !
  ALLOCATE(rho_out(dfftp%nnr,nspin))
  !
  rho_out(:,:) = 0.0_DP
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF (nks>1) THEN
        !
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        !
        CALL mp_bcast(evc,0,inter_image_comm)
        !
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     !
     DO ibnd = 1, nbndval
        !
        w1 = wg(ibnd, iks)/omega
        !
        IF (gamma_only) THEN
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),dvg_exc(1,ibnd,iks,iexc),psic,'Wave')
           !
           rho_out(:, current_spin) = rho_out(:, current_spin) + w1 * DBLE(psic(:))*AIMAG(psic(:))
           !
        ELSE
           !
           ALLOCATE (psic_aux(dffts%nnr))
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,dvg_exc(1,ibnd,iks,iexc),psic_aux,'Wave',igk_k(1,current_k))
           !
           rho_out(:, current_spin) = rho_out(:, current_spin) + w1 * CONJG(psic(:)) * psic_aux(:)
           !
           DEALLOCATE (psic_aux)
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  summ0 = 0.0_DP
  DO ir = 1, dffts%nnr
     summ0 = summ0 + rho_out(ir,1)
  ENDDO
  !
  summ0 = summ0*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  !
  CALL mp_sum(summ0, intra_bgrp_comm)
  !
  !
  WRITE (stdout,*) "Plot of charge density response of the exciton states (Ry unit): ", iexc, ev(iexc)
  WRITE (stdout,*) "summ0", summ0
  !
  westpp_format = 'C'
  !
  WRITE (my_label,'(i6.6)') iexc
  fname = './file_plot_exc_'//TRIM(my_label)
  WRITE (stdout,*) "Write to file : ", fname
  IF(my_image_id==0) CALL dump_r( rho_out(:,1), fname)
  !
  DEALLOCATE( ev )
  DEALLOCATE( dvg_exc )
  DEALLOCATE( rho_out )
  !
  RETURN
  !
ENDSUBROUTINE
