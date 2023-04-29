!
! Copyright (C) 2015-2023 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_setup()
  !-----------------------------------------------------------------------
  !
  USE westcom,              ONLY : localization,l_local_repr,macropol_calculation,l_macropol,&
                                 & solver,l_bse,wbse_calculation,l_davidson,l_lanczos,&
                                 & qp_correction,l_qp_correction,spin_excitation,l_bse_triplet,&
                                 & wstat_calculation,n_pdep_times,n_pdep_eigen,n_pdep_basis,&
                                 & n_pdep_maxiter,n_pdep_read_from_file,trev_pdep_rel,trev_pdep,&
                                 & n_liouville_times,n_liouville_eigen,n_liouville_maxiter,&
                                 & n_liouville_read_from_file,trev_liouville_rel,trev_liouville,&
                                 & alphapv_dfpt,l_use_ecutrho,wbse_save_dir,l_hybrid_tddft
  USE kinds,                ONLY : DP
  USE types_coulomb,        ONLY : pot3D
  USE wbse_dv,              ONLY : wbse_dv_setup
  USE xc_lib,               ONLY : xclib_dft_is
  USE exx_base,             ONLY : exxdiv_treatment,erfc_scrlen
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  SELECT CASE(TRIM(localization))
  CASE('N','n')
     l_local_repr = .FALSE.
  CASE('B','b','W','w')
     l_local_repr = .TRUE.
  END SELECT
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
  END SELECT
  !
  SELECT CASE(TRIM(solver))
  CASE('BSE','bse')
     l_bse = .TRUE.
  CASE('TDDFT','tddft')
     l_bse = .FALSE.
  END SELECT
  !
  ! ground state hybrid DFT + TDDFT -> TD-hybrid-DFT
  !
  IF((.NOT. l_bse) .AND. xclib_dft_is('hybrid')) THEN
     l_hybrid_tddft = .TRUE.
  ELSE
     l_hybrid_tddft = .FALSE.
  ENDIF
  !
  SELECT CASE(wbse_calculation)
  CASE('D','d')
     l_davidson = .TRUE.
     l_lanczos = .FALSE.
  CASE('L','l')
     l_lanczos = .TRUE.
     l_davidson = .FALSE.
  END SELECT
  !
  IF(TRIM(qp_correction) == '') THEN
     l_qp_correction = .FALSE.
  ELSE
     l_qp_correction = .TRUE.
  ENDIF
  !
  SELECT CASE(spin_excitation)
  CASE('S','s','singlet')
     l_bse_triplet = .FALSE.
  CASE('T','t','triplet')
     l_bse_triplet = .TRUE.
  END SELECT
  !
  SELECT CASE(wbse_calculation)
  CASE('l','d','r','R')
     wstat_calculation = 'R'
  CASE('L','D','s','S')
     wstat_calculation = 'S'
  CASE DEFAULT
     wstat_calculation = wbse_calculation
  END SELECT
  !
  n_pdep_times = n_liouville_times
  n_pdep_eigen = n_liouville_eigen
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  n_pdep_maxiter = n_liouville_maxiter
  n_pdep_read_from_file = n_liouville_read_from_file
  trev_pdep_rel = trev_liouville_rel
  trev_pdep = trev_liouville
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  l_use_ecutrho = .FALSE.
  !
  CALL set_npwq()
  !
  IF(l_hybrid_tddft) THEN
     !
     IF(erfc_scrlen > 0._DP) THEN
        !
        ! HSE functional, mya = 1._DP, myb = -1._DP, mymu = erfc_scrlen
        !
        CALL pot3D%init('Rho',.FALSE.,exxdiv_treatment,mya=1._DP,myb=-1._DP,mymu=erfc_scrlen)
        !
     ELSE
        !
        ! PBE0 functional, mya = 1._DP, myb = 0._DP, mymu = 1._DP to avoid divergence
        !
        CALL pot3D%init('Rho',.FALSE.,exxdiv_treatment,mya=1._DP,myb=0._DP,mymu=1._DP)
        !
     ENDIF
     !
  ELSE
     !
     CALL pot3D%init('Rho',.FALSE.,'gb')
     !
  ENDIF
  !
  CALL pot3D%print_divergence()
  !
  CALL set_nbndocc()
  !
  CALL wbse_dv_setup(l_bse)
  !
  CALL my_mkdir(wbse_save_dir)
  !
  IF(l_qp_correction) THEN
     CALL read_qp_eigs()
  ENDIF
  !
  ! read ovl_matrix and u_matrix, and compute macroscopic term, if any
  !
  IF(l_bse .OR. l_hybrid_tddft) CALL bse_start()
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE bse_start()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE pwcom,                ONLY : isk,nks,npwx
  USE westcom,              ONLY : l_reduce_io,tau_is_read,tau_all,n_tau,nbnd_occ,nbndval0x,&
                                 & n_trunc_bands,sigma_c_head,sigma_x_head,wbse_epsinfty,&
                                 & l_local_repr,overlap_thr,u_matrix,ovl_matrix,n_bse_idx,idx_matrix
  USE lsda_mod,             ONLY : nspin
  USE constants,            ONLY : e2,pi
  USE cell_base,            ONLY : omega
  USE types_coulomb,        ONLY : pot3D
  USE wbse_io,              ONLY : read_umatrix_and_omatrix
  USE distribution_center,  ONLY : aband
  USE class_idistribute,    ONLY : idistribute
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iks,is
  INTEGER :: lbnd,ibnd,jbnd,my_ibnd,nbnd_do,nbndval
  INTEGER :: ipair,do_idx
  REAL(DP) :: ovl_value
  COMPLEX(DP), ALLOCATABLE :: u_tmp(:,:)
  !
  ! the divergence term in Fock potential
  !
  sigma_x_head = pot3D%div
  !
  ! compute macroscopic term, it needs macroscopic dielectric constant from input
  !
  sigma_c_head = ((1._DP/wbse_epsinfty) - 1._DP) * (2._DP*e2/pi) * ((6._DP*pi*pi/omega)**(1._DP/3._DP))
  !
  WRITE(stdout,'(/,5X,"Macroscopic dielectric constant correction:",f9.5)') sigma_c_head
  !
  nbnd_do = nbndval0x-n_trunc_bands
  !
  aband = idistribute()
  !
  CALL aband%init(nbnd_do,'b','nbndval',.FALSE.)
  !
  ! allocate and read unitary matrix and overlap matrix, if any
  !
  IF(l_local_repr) THEN
     !
     ALLOCATE(u_tmp(nbnd_do,nbnd_do))
     ALLOCATE(u_matrix(aband%nloc,nbnd_do,nspin))
     ALLOCATE(ovl_matrix(nbnd_do,nbnd_do,nspin))
     !
     DO is = 1,nspin
        !
        CALL read_umatrix_and_omatrix(nbnd_do,is,u_tmp,ovl_matrix(:,:,is))
        !
        DO lbnd = 1,aband%nloc
           ibnd = aband%l2g(lbnd)
           u_matrix(lbnd,:,is) = u_tmp(ibnd,:)
        ENDDO
        !
     ENDDO
     !
     DEALLOCATE(u_tmp)
     !
  ENDIF
  !
  ALLOCATE(idx_matrix(nbnd_do*nbnd_do,2,nspin))
  !
  idx_matrix(:,:,:) = 0
  !
  ALLOCATE(n_bse_idx(nspin))
  !
  n_bse_idx(:) = 0
  !
  DO iks = 1,nks
     !
     nbndval = nbnd_occ(iks)
     is = isk(iks)
     do_idx = 0
     !
     DO ibnd = 1,nbndval-n_trunc_bands
        DO jbnd = 1,nbndval-n_trunc_bands
           IF(l_local_repr) THEN
              ovl_value = ovl_matrix(ibnd,jbnd,is)
           ELSE
              ovl_value = 0._DP
           ENDIF
           !
           IF(l_local_repr) THEN
              IF(ovl_value >= overlap_thr) THEN
                 do_idx = do_idx + 1
                 idx_matrix(do_idx,1,is) = ibnd+n_trunc_bands
                 idx_matrix(do_idx,2,is) = jbnd+n_trunc_bands
              ENDIF
           ELSE
              do_idx = do_idx + 1
              idx_matrix(do_idx,1,is) = ibnd+n_trunc_bands
              idx_matrix(do_idx,2,is) = jbnd+n_trunc_bands
           ENDIF
        ENDDO
     ENDDO
     !
     n_bse_idx(is) = do_idx
     !
  ENDDO
  !
  IF(l_reduce_io) THEN
     !
     ! I/O is reduced by reading tau only once
     !
     ALLOCATE(tau_is_read(nbnd_do,nbnd_do,nspin))
     !
     tau_is_read(:,:,:) = 0
     n_tau = 0
     !
     DO iks = 1,nks
        !
        is = isk(iks)
        do_idx = n_bse_idx(is)
        !
        ! count number of tau needed by my band group
        !
        DO lbnd = 1,aband%nloc
           !
           my_ibnd = aband%l2g(lbnd)
           !
           DO ipair = 1, do_idx
              !
              ibnd = idx_matrix(ipair,1,is)-n_trunc_bands
              jbnd = idx_matrix(ipair,2,is)-n_trunc_bands
              !
              IF(ibnd == my_ibnd) tau_is_read(ibnd,jbnd,is) = 1
              !
           ENDDO
           !
        ENDDO
        !
        DO jbnd = 1,nbnd_do
           DO ibnd = jbnd,nbnd_do
              IF(tau_is_read(ibnd,jbnd,is) == 1 .OR. tau_is_read(jbnd,ibnd,is) == 1) n_tau = n_tau+1
           ENDDO
        ENDDO
        !
     ENDDO
     !
     ALLOCATE(tau_all(npwx,n_tau))
     !
     ! reset counter for actual reading in wbse_io
     !
     tau_is_read(:,:,:) = 0
     n_tau = 0
     !
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE read_qp_eigs()
  !-----------------------------------------------------------------------
  !
  ! read qp eigenvalues from JSON file
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : RYTOEV
  USE io_global,            ONLY : ionode
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : intra_image_comm
  USE westcom,              ONLY : et_qp,qp_correction
  USE pwcom,                ONLY : nspin,nbnd
  USE json_module,          ONLY : json_file
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: found
  INTEGER :: is,nspin_tmp
  CHARACTER(LEN=6) :: labels
  REAL(DP),ALLOCATABLE :: rvals(:)
  TYPE(json_file) :: json
  !
  ALLOCATE(et_qp(nbnd,nspin))
  !
  IF(ionode) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(qp_correction))
     !
     CALL json%get('system.electron.nspin',nspin_tmp,found)
     IF(.NOT. found) CALL errore('read_qp_eigs','nspin not found',1)
     IF(nspin_tmp /= nspin) CALL errore('read_qp_eigs','nspin mismatch',1)
     !
     DO is = 1,nspin
        !
        WRITE(labels,'(i6.6)') is
        !
        CALL json%get('output.Q.K'//labels//'.eqpSec',rvals,found)
        IF(.NOT. found) CALL errore('read_qp_eigs','eqpSec not found',1)
        IF(SIZE(rvals) /= nbnd) CALL errore('read_qp_eigs','nbnd mismatch',1)
        et_qp(:,is) = rvals/RYTOEV
        !
     ENDDO
     !
     CALL json%destroy()
     !
  ENDIF
  !
  CALL mp_bcast(et_qp,0,intra_image_comm)
  !
END SUBROUTINE
