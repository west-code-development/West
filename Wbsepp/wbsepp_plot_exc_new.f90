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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
SUBROUTINE wbsepp_plot_exc()
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
                                   & inter_bgrp_comm, my_image_id, inter_image_comm
  USE buffers,                ONLY : get_buffer
  USE westcom,                ONLY : iuwfc,lrwfc,nbnd_occ,ev,dvg_exc,n_plep_read_from_file,&
                                   & iexc_plot,r0_input
  USE fft_at_gamma,           ONLY : single_invfft_gamma, double_invfft_gamma
  USE fft_at_k,               ONLY : single_fwfft_k,single_invfft_k
  USE plep_db,                ONLY : plep_db_read
  USE distribution_center,    ONLY : pert
  USE class_idistribute,      ONLY : idistribute
  USE westcom,                ONLY : westpp_format
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
  COMPLEX(DP), ALLOCATABLE :: rho_aux(:,:)
  INTEGER, PARAMETER :: n_ipol = 3
  CHARACTER(LEN=6)   :: my_label
  CHARACTER(LEN=256) :: fname
  !
  INTEGER  :: iexc
  REAL(DP) :: r0_input_aux(3)
  !
  INTEGER  :: n_point, ni
  REAL(DP) :: dstep, ri1, ri0, dx,dy,dz, drr, summ, summ0, vi1, vi0
  !
  iexc = iexc_plot
  r0_input_aux(1) =  r0_input(1)
  r0_input_aux(2) =  r0_input(2)
  r0_input_aux(3) =  r0_input(3)
  !
write(stdout,*) at(1,1)
  DO ip = 1, n_ipol
     !
     IF (r0_input_aux(ip) .LT. 0) r0_input_aux(ip) = r0_input_aux(ip) + at(ip,ip)
     !
     IF (r0_input_aux(ip) .GT. at(ip,ip)) r0_input_aux(ip) = r0_input_aux(ip) - at(ip,ip)
     !
  ENDDO
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
  ALLOCATE(r(dfftp%nnr,n_ipol))
  ALLOCATE(dr(dfftp%nnr))
  ALLOCATE(segno(dfftp%nnr))
  ALLOCATE(rho_aux(dfftp%nnr,nspin))
  ALLOCATE(rho_out(dfftp%nnr,nspin))
  !
  r(:,:) = 0.0d0
  !
  ! Calculate r
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
#if defined (__MPI)
  index0 = dfftp%nr1x*dfftp%nr2x*SUM(dfftp%npp(1:me_bgrp))
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
  index0 = 0
  ir_end = dfftp%nnr
#endif
  !
write (stdout,*) ir_end
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     !
     DO ip = 1, n_ipol
        !
        r(ir,ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                   DBLE( j )*inv_nr2*at(ip,2) + &
                   DBLE( k )*inv_nr3*at(ip,3)
        !
     ENDDO
     !
  ENDDO
  !
  dr(:) = 10000000
  !
  DO ir = 1, dffts%nnr
     !
     dr(ir) = DSQRT((r(ir,1) - r0_input_aux(1))**2 + &
                    (r(ir,2) - r0_input_aux(2))**2 + &
                    (r(ir,3) - r0_input_aux(3))**2)
     !
  ENDDO
  !
  drmin_g = MINVAL(dr(:))
  !
  CALL mp_min(drmin_g, intra_bgrp_comm)
  !
  rho_aux(:,:) = 0.0_DP
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
        w1 = wg(ibnd, iks)
        !
        IF (gamma_only) THEN
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),dvg_exc(1,ibnd,iks,iexc),psic,'Wave')
           !
           rcoeff = 0.0_DP
           !
           DO ir = 1, dffts%nnr
              !
              IF (dr(ir) == drmin_g) THEN
                 !
                 rcoeff = DBLE( psic(ir) )
                 !
              ENDIF
              !
           ENDDO
           !
           CALL mp_sum(rcoeff, intra_bgrp_comm)
           !
           segno(:) = SIGN( 1.d0, AIMAG(psic(:)) )
           !
           rho_aux(:, current_spin) = rho_aux(:, current_spin) + w1 * CMPLX(rcoeff*AIMAG(psic(:)), 0.0_DP)
           !
        ELSE
           !
           ALLOCATE (psic_aux(dffts%nnr))
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,dvg_exc(1,ibnd,iks,iexc),psic_aux,'Wave',igk_k(1,current_k))
           !
           zcoeff = 0.0_DP
           !
           DO ir = 1, dffts%nnr
              !
              IF (dr(ir) == drmin_g) THEN
                 !
                 zcoeff = psic(ir)
                 !
              ENDIF
              !
           ENDDO
           !
           CALL mp_sum(zcoeff, intra_bgrp_comm)
           !
           rho_aux(:, current_spin) = rho_aux(:, current_spin) + w1 * zcoeff * psic_aux(:)
           !
           DEALLOCATE (psic_aux)
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  rho_out(:,:) = 0.0_DP
  rho_out(:, current_spin) = (DBLE(rho_aux(:, current_spin))**2 + AIMAG(rho_aux(:, current_spin))**2)/omega
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
  rho_out(:,1) = rho_out(:,1)/summ0
  !
  WRITE (stdout,*) "Plot of exciton states (Ry unit): ", iexc, ev(iexc)
  WRITE (stdout,*) "Hole state is fixed at (alat unit): ", r0_input_aux(1), r0_input_aux(2), r0_input_aux(3)
  !
  westpp_format = 'C'
  !
  WRITE (my_label,'(i6.6)') iexc
  fname = './file_plot_exc_'//TRIM(my_label)
  WRITE (stdout,*) "Write to file : ", fname
  IF(my_image_id==0) CALL dump_r( rho_out(:,1), fname)
  !
  summ = 0.0_DP
  DO ir = 1, dffts%nnr
     summ = summ + rho_out(ir,1)
  ENDDO
  !
  summ = summ*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  !
  CALL mp_sum(summ, intra_bgrp_comm)
  !
  WRITE (stdout,*) "Check normalization: ", summ
  !
  ! plot distribution density
  !
  WRITE (stdout,*) "Plot of distribution elecrtron "
  !
  dstep = 0.01
  n_point = int(at(1,1)*0.5/dstep)
  !
  DO ni = 1, n_point
     !
     ri1 = (ni-1)*dstep
     ri0 =  ni*dstep
     !
     vi1 = (4.0/3.0)*3.1415*ri1*ri1*ri1*alat*alat*alat
     vi0 = (4.0/3.0)*3.1415*ri0*ri0*ri0*alat*alat*alat
     !
     summ = 0.0_DP
     summ0 = 0.0
     DO ir = 1, dffts%nnr
        !
        dx = r(ir,1) - r0_input_aux(1)
        IF (dx .gt.  at(1,1) * 0.5) dx = dx - at(1,1)
        IF (dx .le. -at(1,1) * 0.5) dx = dx + at(1,1)
        !
        dy = r(ir,2) - r0_input_aux(2)
        IF (dy .gt.  at(2,2) * 0.5) dy = dy - at(2,2)
        IF (dy .le. -at(2,2) * 0.5) dy = dy + at(2,2)
        !
        dz = r(ir,3) - r0_input_aux(3)
        IF (dz .gt.  at(3,3) * 0.5) dz = dz - at(3,3)
        IF (dz .le. -at(3,3) * 0.5) dz = dz + at(3,3)
        !
        drr = SQRT(dx*dx + dy*dy + dz*dz)
        !
        IF ((ri1 .lt. drr).and.(drr .le. ri0)) THEN
           !
           summ0 = summ0 + 1
           summ = summ + rho_out(ir,1)
           !
        ENDIF
        !
     ENDDO
     !
     CALL mp_sum(summ, intra_bgrp_comm)
     CALL mp_sum(summ0, intra_bgrp_comm)
     !
     WRITE(stdout, "(i5, 5x, f12.6, 5x, f12.6, 5x, 2f12.6)") ni, ri1*alat, (vi0-vi1), summ0, summ
     !
  ENDDO
  !
  DEALLOCATE( ev )
  DEALLOCATE( dvg_exc )
  DEALLOCATE( r,dr,rho_aux,segno,rho_out)
  !
  RETURN
  !
ENDSUBROUTINE
