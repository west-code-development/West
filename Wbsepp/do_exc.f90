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
SUBROUTINE do_exc()
  !
  USE kinds,                 ONLY : DP
  USE io_push,               ONLY : io_push_title
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE cell_base,             ONLY : omega,at,alat
  USE fft_base,              ONLY : dffts
  USE wavefunctions,         ONLY : psic,evc
  USE pwcom,                 ONLY : npw,npwx,igk_k,current_k,ngk,nks,wg
  USE control_flags,         ONLY : gamma_only
  USE mp,                    ONLY : mp_sum,mp_bcast,mp_min
  USE mp_global,             ONLY : me_bgrp,intra_bgrp_comm,my_image_id,inter_image_comm
  USE buffers,               ONLY : get_buffer
  USE westcom,               ONLY : iuwfc,lrwfc,nbnd_occ,dvg_exc,n_liouville_read_from_file,&
                                  & iexc_plot,wbsepp_r0,westpp_format,wbsepp_save_dir
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
  INTEGER :: ibnd,nbndval,i,j,k,ir,ip,iks,index0,ir_end
  REAL(DP) :: drmin_g,rcoeff,w1,summ0,inv_nr1,inv_nr2,inv_nr3
  REAL(DP) :: r0(3)
  COMPLEX(DP) :: zcoeff
  REAL(DP), ALLOCATABLE :: r(:,:),dr(:),rho(:)
  COMPLEX(DP), ALLOCATABLE :: psic_aux(:),rho_aux(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: n_ipol = 3
  !
  westpp_format = 'C'
  !
  r0(1) = wbsepp_r0(1)/alat
  r0(2) = wbsepp_r0(2)/alat
  r0(3) = wbsepp_r0(3)/alat
  !
  DO ip = 1, n_ipol
     !
     DO WHILE(r0(ip) < 0)
        r0(ip) = r0(ip) + at(ip,ip)
     ENDDO
     !
     DO WHILE(r0(ip) >= at(ip,ip))
        r0(ip) = r0(ip) - at(ip,ip)
     ENDDO
     !
  ENDDO
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
  ALLOCATE(r(dffts%nnr,n_ipol))
  ALLOCATE(dr(dffts%nnr))
  ALLOCATE(rho(dffts%nnr))
  ALLOCATE(rho_aux(dffts%nnr))
  !
  r(:,:) = 0._DP
  !
  ! Calculate r
  !
  inv_nr1 = 1._DP / REAL(dffts%nr1,KIND=DP)
  inv_nr2 = 1._DP / REAL(dffts%nr2,KIND=DP)
  inv_nr3 = 1._DP / REAL(dffts%nr3,KIND=DP)
  !
  index0 = dffts%nr1x*dffts%nr2x*SUM(dffts%nr3p(1:me_bgrp))
  ir_end = MIN(dffts%nnr,dffts%nr1x*dffts%nr2x*dffts%nr3p(me_bgrp+1))
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dffts%nr1x*dffts%nr2x)
     i = i - (dffts%nr1x*dffts%nr2x)*k
     j = i / dffts%nr1x
     i = i - dffts%nr1x*j
     !
     DO ip = 1, n_ipol
        r(ir,ip) = REAL(i,KIND=DP)*inv_nr1*at(ip,1) + &
                   REAL(j,KIND=DP)*inv_nr2*at(ip,2) + &
                   REAL(k,KIND=DP)*inv_nr3*at(ip,3)
     ENDDO
     !
  ENDDO
  !
  dr(:) = 10000000
  !
  DO ir = 1, dffts%nnr
     dr(ir) = SQRT((r(ir,1) - r0(1))**2 + &
                   (r(ir,2) - r0(2))**2 + &
                   (r(ir,3) - r0(3))**2)
  ENDDO
  !
  drmin_g = MINVAL(dr)
  !
  CALL mp_min(drmin_g,intra_bgrp_comm)
  !
  rho_aux(:) = 0._DP
  !
  CALL io_push_title('(E)xciton State')
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
        w1 = wg(ibnd,iks)
        !
        IF(gamma_only) THEN
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),dvg_exc(:,ibnd,iks,iexc_plot),psic,'Wave')
           !
           rcoeff = 0._DP
           !
           DO ir = 1, dffts%nnr
              IF(dr(ir) == drmin_g) THEN
                 rcoeff = REAL(psic(ir),KIND=DP)
              ENDIF
           ENDDO
           !
           CALL mp_sum(rcoeff,intra_bgrp_comm)
           !
           rho_aux(:) = rho_aux(:) + w1 * CMPLX(rcoeff*AIMAG(psic(:)),KIND=DP)
           !
        ELSE
           !
           ALLOCATE(psic_aux(dffts%nnr))
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,current_k))
           CALL single_invfft_k(dffts,npw,npwx,dvg_exc(:,ibnd,iks,iexc_plot),psic_aux,'Wave',igk_k(:,current_k))
           !
           zcoeff = 0._DP
           !
           DO ir = 1, dffts%nnr
              IF(dr(ir) == drmin_g) THEN
                 zcoeff = psic(ir)
              ENDIF
           ENDDO
           !
           CALL mp_sum(zcoeff,intra_bgrp_comm)
           !
           rho_aux(:) = rho_aux(:) + w1 * zcoeff * psic_aux(:)
           !
           DEALLOCATE(psic_aux)
           !
        ENDIF
        !
        CALL update_bar_type(barra,'wbsepp',1)
        !
     ENDDO
     !
     rho(:) = (REAL(rho_aux(:),KIND=DP)**2 + AIMAG(rho_aux(:))**2)/omega
     summ0 = SUM(rho)
     summ0 = summ0*omega/(dffts%nr1*dffts%nr2*dffts%nr3)
     !
     CALL mp_sum(summ0,intra_bgrp_comm)
     !
     rho(:) = rho/summ0
     !
     WRITE(fname,'(a,i6.6,a,i6.6)') TRIM(wbsepp_save_dir)//'/excK',iks,'E',iexc_plot
     IF(my_image_id == 0) CALL dump_r(rho,TRIM(fname))
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'wbsepp')
  !
  DEALLOCATE(r)
  DEALLOCATE(dr)
  DEALLOCATE(rho)
  DEALLOCATE(rho_aux)
  !
END SUBROUTINE
