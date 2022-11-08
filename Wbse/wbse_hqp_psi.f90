!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE bse_hqp_psi(iks, nbvalloc, psi, dpsi)
  !
  ! evc is qp wfcs, et_qp is qp eigenvalues
  ! evc_ks is ks wfc, et is ks wfcs
  !
  ! This routine computes:
  ! dpsi = S|evc><evc|psi>*(et_qp - et - Delta) + Delta|psi>
  !
  USE kinds,               ONLY : DP
  USE noncollin_module,    ONLY : npol
  USE wavefunctions,       ONLY : evc
  USE mp_global,           ONLY : intra_bgrp_comm
  USE mp,                  ONLY : mp_sum
  USE control_flags,       ONLY : gamma_only
  USE pwcom,               ONLY : npw,npwx,nbnd,et
  USE westcom,             ONLY : et_qp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iks, nbvalloc
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,nbvalloc)
  COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbvalloc)
  !
  ! local variables
  !
  INTEGER :: ibnd
  REAL(DP) :: delta
  REAL(DP), ALLOCATABLE :: ps_r(:,:)
  COMPLEX(DP), ALLOCATABLE :: ps_c(:,:)
  !
  delta = et_qp(nbnd,iks) - et(nbnd,iks)
  !
  IF(gamma_only) THEN
     ALLOCATE(ps_r(nbnd,nbvalloc))
     ps_r = 0._DP
  ENDIF
  ALLOCATE(ps_c(nbnd,nbvalloc))
  ps_c = (0._DP,0._DP)
  !
  ! ps = < evc | f >
  !
  IF(gamma_only) THEN
     CALL glbrak_gamma(evc,psi,ps_r,npw,npwx,nbnd,nbvalloc,nbnd,npol)
     CALL mp_sum(ps_r,intra_bgrp_comm)
     ps_c(:,:) = CMPLX(ps_r,KIND=DP)
  ELSE
     CALL glbrak_k(evc,psi,ps_c,npw,npwx,nbnd,nbvalloc,nbnd,npol)
     CALL mp_sum(ps_c,intra_bgrp_comm)
  ENDIF
  !
  DO ibnd = 1, nbnd
     ps_c(ibnd,:) =  ps_c(ibnd,:) * (et_qp(ibnd,iks) - et(ibnd,iks) - delta)
  ENDDO
  !
  CALL ZGEMM('N','N',npwx*npol,nbvalloc,nbnd,(1._DP,0._DP),evc,npwx*npol,&
       & ps_c,nbnd,(1._DP,0._DP),dpsi,npwx*npol)
  !
  dpsi(:,:) = dpsi + delta * psi
  !
  IF(gamma_only) DEALLOCATE(ps_r)
  DEALLOCATE(ps_c)
  !
END SUBROUTINE
!
SUBROUTINE read_qp_eigs()
  !
  USE io_global,     ONLY : ionode,ionode_id
  USE mp,            ONLY : mp_bcast,mp_barrier
  USE mp_world,      ONLY : world_comm
  USE westcom,       ONLY : et_qp,qp_correction
  USE pwcom,         ONLY : nks,nbnd
  USE lsda_mod,      ONLY : lsda
  !
  IMPLICIT NONE
  !
  INTEGER :: iun
  INTEGER :: ibnd, nqp_eigs
  INTEGER :: num_k_points,ik
  CHARACTER(LEN=3) :: my_ik
  CHARACTER(LEN=256) :: file_qp
  !
  ! read eigenvalues from file
  !
  CALL mp_barrier(world_comm)
  !
  ALLOCATE(et_qp(nbnd,nks))
  !
  IF(lsda) THEN
     num_k_points = nks/2
  ELSE
     num_k_points = nks
  ENDIF
  !
  DO ik = 1,num_k_points
     !
     WRITE(my_ik,'(i1)') ik
     !
     file_qp = TRIM(qp_correction)//'.'//TRIM(my_ik)
     !
     IF(ionode) THEN
        OPEN(NEWUNIT=iun,FILE=TRIM(file_qp),FORM='FORMATTED',STATUS='OLD')
        READ(iun,*) nqp_eigs
     ENDIF
     !
     CALL mp_bcast(nqp_eigs,ionode_id,world_comm)
     !
     IF(nqp_eigs /= nbnd) CALL errore ('read_qp_eigs', 'nqp_eigs /= nbnd', 1)
     !
     IF(ionode) THEN
        !
        DO ibnd = 1,nqp_eigs
           IF(lsda) THEN
              READ(iun,*) et_qp(ibnd,ik),et_qp(ibnd,ik+num_k_points)
           ELSE
              READ(iun,*) et_qp(ibnd,ik)
           ENDIF
        ENDDO
        !
        CLOSE(iun)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_bcast(et_qp,ionode_id,world_comm)
  !
END SUBROUTINE
!
SUBROUTINE read_ks_wfc()
  !
  USE mp,            ONLY : mp_bcast,mp_barrier
  USE mp_world,      ONLY : world_comm
  USE westcom,       ONLY : evc_ks
  USE pwcom,         ONLY : npwx,nks,nbnd
  USE plep_io,       ONLY : plep_read_G_and_distribute_wfc
  !
  IMPLICIT NONE
  !
  INTEGER :: iks
  CHARACTER(LEN=3) :: my_ik
  CHARACTER(LEN=256) :: fname
  !
  ! read eigenvalues from file
  !
  CALL mp_barrier(world_comm)
  !
  ALLOCATE(evc_ks(npwx,nbnd,nks))
  !
  DO iks = 1, nks
     !
     WRITE(my_ik,'(i1)') iks
     !
     fname = './ks_wfc_tmp_'//TRIM(my_ik)//'.dat'
     CALL plep_read_G_and_distribute_wfc(fname,evc_ks(:,:,iks),nbnd)
     !
  ENDDO
  !
END SUBROUTINE
