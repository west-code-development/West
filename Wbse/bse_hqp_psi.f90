!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE bse_hqp_psi(iks, current_spin, nbvalloc, psi, dpsi)
  !
  ! evc is qp wfcs, et_qp is qp eigenvalues
  ! evc_ks is ks wfc, et is ks wfcs
  !
  ! This routine compute:
  ! dpsi = S|evc><evc|psi>*(et_qp - Delta)
  !      - S|evc_ks><evc_ks|psi>*et + Delta|psi>
  !
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE noncollin_module,    ONLY : npol
  USE gvect,               ONLY : gstart
  USE pwcom,               ONLY : npwx
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: iks, current_spin, nbvalloc
  COMPLEX(DP), INTENT(IN)    :: psi(npwx*npol,nbvalloc)
  COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbvalloc)
  !
  IF (.false.) THEN
     !
     CALL bse_hqp_psi_scf(iks, current_spin, nbvalloc, psi, dpsi)
     !
  ELSE
     !
     CALL bse_hqp_psi_nscf(iks, current_spin, nbvalloc, psi, dpsi)
     !
  ENDIF
  !
  RETURN
  !
ENDSUBROUTINE


SUBROUTINE bse_hqp_psi_scf(iks, current_spin, nbvalloc, psi, dpsi)
  !
  ! evc is qp wfcs, et_qp is qp eigenvalues
  ! evc_ks is ks wfc, et is ks wfcs
  !
  ! This routine compute:
  ! dpsi = S|evc><evc|psi>*(et_qp - Delta)
  !      - S|evc_ks><evc_ks|psi>*et + Delta|psi>
  !
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE noncollin_module,    ONLY : noncolin, npol
  USE wavefunctions       ,ONLY : evc
  USE mp_bands,            ONLY : intra_bgrp_comm
  USE mp,                  ONLY : mp_sum
  USE control_flags,       ONLY : gamma_only
  USE gvect,               ONLY : gstart
  USE pwcom,               ONLY : npw, npwx, nbnd, et
  USE bse_module,          ONLY : et_qp, evc_ks
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: iks, current_spin, nbvalloc
  COMPLEX(DP), INTENT(IN)    :: psi(npwx*npol,nbvalloc)
  COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbvalloc)
  !
  ! local variables
  !
  INTEGER  :: ibnd
  REAL(DP) :: delta
  REAL(DP),    ALLOCATABLE   :: ps_r(:,:), ps_r_ks(:,:)
  COMPLEX(DP), ALLOCATABLE   :: ps_c(:,:), ps_c_ks(:,:)
  !
  delta = et_qp(nbnd,iks) - et(nbnd,iks)
  !
  IF( gamma_only ) THEN
     ALLOCATE( ps_r(nbnd, nbvalloc) )
     ALLOCATE( ps_r_ks(nbnd, nbvalloc) )
     ps_r = 0.0_DP
     ps_r_ks = 0.0_DP
  ENDIF
  !
  ALLOCATE( ps_c(nbnd, nbvalloc) )
  ALLOCATE( ps_c_ks(nbnd, nbvalloc) )
  ps_c = 0.0_DP
  ps_c_ks = 0.0_DP
  !
  ! ps = < evc | f >
  ! ps_ks = < evc_ks | f >
  !
  IF( gamma_only ) THEN
     !
     CALL glbrak_gamma ( evc, psi, ps_r, npw, npwx, nbnd, nbvalloc, nbnd, npol)
     CALL glbrak_gamma ( evc_ks(:,:,iks), psi, ps_r_ks, npw, npwx, nbnd, nbvalloc, nbnd, npol)
     CALL mp_sum(ps_r,intra_bgrp_comm)
     CALL mp_sum(ps_r_ks,intra_bgrp_comm)
     ps_c(:,:) = CMPLX (ps_r(:,:),0.0_DP, KIND=DP)
     ps_c_ks(:,:) = CMPLX (ps_r_ks(:,:),0.0_DP, KIND=DP)
     !
  ELSE
     !
     CALL glbrak_k ( evc, psi, ps_c, npw, npwx, nbnd, nbvalloc, nbnd, npol)
     CALL glbrak_k ( evc_ks(:,:,iks), psi, ps_c_ks, npw, npwx, nbnd, nbvalloc, nbnd, npol)
     CALL mp_sum(ps_c,intra_bgrp_comm)
     CALL mp_sum(ps_c_ks,intra_bgrp_comm)
     !
  ENDIF
  !
  DO ibnd = 1, nbnd
     !
     ps_c(ibnd,:) =  ps_c(ibnd,:) * (et_qp(ibnd,iks) - delta)
     ps_c_ks(ibnd,:) =  -ps_c_ks(ibnd,:) * et(ibnd,iks)
     !
  ENDDO
  !
  CALL ZGEMM('N','N',npwx*npol,nbvalloc,nbnd,(1.0_DP,0.0_DP), &
             evc,npwx*npol,ps_c,nbnd,(1.0_DP,0.0_DP),dpsi,npwx*npol)
  CALL ZGEMM('N','N',npwx*npol,nbvalloc,nbnd,(1.0_DP,0.0_DP), &
             evc_ks(:,:,iks),npwx*npol,ps_c_ks,nbnd,(1.0_DP,0.0_DP),dpsi,npwx*npol)
  !
  dpsi(:,:) = dpsi(:,:) + delta * psi(:,:)
  !
  IF( gamma_only ) THEN
     DEALLOCATE( ps_r )
     DEALLOCATE( ps_r_ks )
  ENDIF
  !
  DEALLOCATE( ps_c )
  DEALLOCATE( ps_c_ks )
  !
  RETURN
  !
ENDSUBROUTINE bse_hqp_psi_scf
!
!
!
SUBROUTINE bse_hqp_psi_nscf(iks, current_spin, nbvalloc, psi, dpsi)
  !
  ! This routine compute:
  ! dpsi = S|evc><evc|psi>*(et_qp - et - Delta) + Delta|psi>
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE noncollin_module,    ONLY : noncolin, npol
  USE wavefunctions,       ONLY : evc
  USE mp_bands,            ONLY : intra_bgrp_comm
  USE mp,                  ONLY : mp_sum
  USE control_flags,       ONLY : gamma_only
  USE gvect,               ONLY : gstart
  USE pwcom,               ONLY : npw, npwx, nbnd, et
  USE bse_module,          ONLY : et_qp
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: iks, current_spin, nbvalloc
  COMPLEX(DP), INTENT(IN)    :: psi(npwx*npol,nbvalloc)
  COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbvalloc)
  !
  ! local variables
  !
  INTEGER  :: ibnd
  REAL(DP) :: delta
  REAL(DP),    ALLOCATABLE   :: ps_r(:,:)
  COMPLEX(DP), ALLOCATABLE   :: ps_c(:,:)
  !
  delta = et_qp(nbnd,iks) - et(nbnd,iks)
  !
  IF( gamma_only ) THEN
     ALLOCATE( ps_r(nbnd, nbvalloc) )
     ps_r = 0.0_DP
  ENDIF
  !
  ALLOCATE( ps_c(nbnd, nbvalloc) )
  ps_c = (0.0_DP, 0.0_DP)
  !
  ! ps = < evc | f >
  !
  IF( gamma_only ) THEN
     !
     CALL glbrak_gamma ( evc, psi, ps_r, npw, npwx, nbnd, nbvalloc, nbnd, npol)
     CALL mp_sum(ps_r,intra_bgrp_comm)
     ps_c(:,:) = CMPLX (ps_r(:,:), 0.0_DP, KIND=DP)
     !
  ELSE
     !
     CALL glbrak_k ( evc, psi, ps_c, npw, npwx, nbnd, nbvalloc, nbnd, npol)
     CALL mp_sum(ps_c,intra_bgrp_comm)
     !
  ENDIF
  !
  DO ibnd = 1, nbnd
     !
     ps_c(ibnd,:) =  ps_c(ibnd,:) * (et_qp(ibnd,iks) - et(ibnd,iks) - delta)
     !
  ENDDO
  !
  CALL ZGEMM('N','N',npwx*npol,nbvalloc,nbnd,(1.0_DP,0.0_DP), &
              evc,npwx*npol,ps_c,nbnd,(1.0_DP,0.0_DP),dpsi,npwx*npol)
  !
  dpsi(:,:) = dpsi(:,:) + delta * psi(:,:)
  !
  IF( gamma_only ) THEN
     DEALLOCATE( ps_r )
  ENDIF
  !
  DEALLOCATE( ps_c )
  !
  RETURN
  !
ENDSUBROUTINE bse_hqp_psi_nscf
!
!
!
SUBROUTINE read_qp_eigs()
  !
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE mp_world,      ONLY : world_comm
  USE bse_module,    ONLY : et_qp
  USE pwcom,         ONLY : nks, nbnd, et
  USE lsda_mod,      ONLY : lsda
  USE westcom,       ONLY : qp_correction
  !
  IMPLICIT NONE
  !
  INTEGER :: ibnd, nqp_eigs
  INTEGER :: num_k_points,ik
  CHARACTER(LEN=3) :: my_ik
  CHARACTER(LEN=256) :: file_qp
  !
  ! read eigenvalues from file
  !
  call mp_barrier(world_comm)
  !
  allocate(et_qp(nbnd,nks))
  !
  IF ( lsda ) THEN
     !
     num_k_points = nks / 2
     !
  ELSE
     !
     num_k_points = nks
     !
  ENDIF
  !
  DO ik = 1, num_k_points
     !
     WRITE(my_ik,'(i1)') ik
     !
     file_qp = TRIM(qp_correction)// "." //TRIM(my_ik)
     !file_qp = 'qp_eigs_'//TRIM(my_ik)//'.dat'
     !
     if (ionode) then
        !
        open(unit = 99,file =TRIM(file_qp),form = 'formatted',status = 'old')
        !
        read(99,*) nqp_eigs
        !
     endif
     !
     call mp_bcast (nqp_eigs,ionode_id,world_comm )
     !
     if (nqp_eigs .ne. nbnd) then
        !
        call errore ('read_qp_eigs', 'number qp eigenvalues not equal nbnd in pw.in', nqp_eigs)
        !
     endif
     !
     if (ionode) then
        !
        do ibnd = 1, nqp_eigs
           !
           IF ( lsda ) THEN
              read (99, * ) et_qp(ibnd,ik), et_qp(ibnd,ik+num_k_points)
           ELSE
              read (99, * ) et_qp(ibnd,ik)
           ENDIF
           !
        enddo
        !
        close (99)
        !
     endif
     !
  ENDDO
  !
  call mp_bcast(et_qp,ionode_id,world_comm)
  !
  return
  !
ENDSUBROUTINE read_qp_eigs
!
!
!
SUBROUTINE read_ks_wfc()
  !
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE mp_world,      ONLY : world_comm
  USE bse_module,    ONLY : evc_ks
  USE pwcom,         ONLY : npwx,nks,nbnd
  USE lsda_mod,      ONLY : lsda
  USE plep_io,       ONLY : plep_read_G_and_distribute_wfc
  !
  IMPLICIT NONE
  !`
  INTEGER            :: iks
  CHARACTER(LEN=3)   :: my_ik
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
     fname = "./ks_wfc_tmp_"//TRIM( my_ik )//'.dat'
     CALL plep_read_G_and_distribute_wfc(fname,evc_ks(:,:,iks),nbnd)
     !
  ENDDO
  !
  RETURN
  !
ENDSUBROUTINE read_ks_wfc
