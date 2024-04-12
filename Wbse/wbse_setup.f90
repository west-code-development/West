!
! Copyright (C) 2015-2024 M. Govoni
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
  USE westcom,              ONLY : localization,l_local_repr,solver,l_bse,wbse_calculation,&
                                 & l_davidson,l_lanczos,qp_correction,l_qp_correction,&
                                 & spin_excitation,l_bse_triplet,wstat_calculation,n_pdep_times,&
                                 & n_pdep_eigen,n_pdep_basis,n_pdep_maxiter,n_pdep_read_from_file,&
                                 & trev_pdep_rel,trev_pdep,n_liouville_times,n_liouville_eigen,&
                                 & n_liouville_maxiter,n_liouville_read_from_file,&
                                 & trev_liouville_rel,trev_liouville,alphapv_dfpt,l_use_ecutrho,&
                                 & wbse_save_dir,l_hybrid_tddft,l_spin_flip,l_spin_flip_kernel,&
                                 & do_inexact_krylov
  USE kinds,                ONLY : DP
  USE types_coulomb,        ONLY : pot3D
  USE wbse_dv,              ONLY : wbse_dv_setup,wbse_sf_kernel_setup
  USE xc_lib,               ONLY : xclib_dft_is
  USE exx_base,             ONLY : erfc_scrlen
  USE pwcom,                ONLY : nkstot,nks,nspin
  USE distribution_center,  ONLY : kpt_pool
  USE class_idistribute,    ONLY : idistribute,IDIST_BLK
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  SELECT CASE(localization)
  CASE('N','n')
     l_local_repr = .FALSE.
  CASE('B','b','W','w')
     l_local_repr = .TRUE.
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
  l_bse_triplet = .FALSE.
  IF(nspin == 1 .AND. l_bse .AND. (spin_excitation == 'T' .OR. spin_excitation == 't')) &
  & l_bse_triplet = .TRUE.
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
        CALL pot3D%init('Rho',.FALSE.,'gb',mya=1._DP,myb=-1._DP,mymu=erfc_scrlen)
        !
     ELSE
        !
        ! PBE0 functional, mya = 1._DP, myb = 0._DP, mymu = 1._DP to avoid divergence
        !
        CALL pot3D%init('Rho',.FALSE.,'gb',mya=1._DP,myb=0._DP,mymu=1._DP)
        !
     ENDIF
     !
  ELSE
     !
     CALL pot3D%init('Rho',.FALSE.,'gb')
     !
  ENDIF
  !
  !$acc enter data copyin(pot3D)
  !$acc enter data copyin(pot3D%sqvc)
  !
  CALL pot3D%print_divergence()
  !
  CALL set_nbndocc()
  !
  CALL wbse_dv_setup(l_bse)
  !
  IF(l_spin_flip .AND. l_spin_flip_kernel) CALL wbse_sf_kernel_setup()
  !
  CALL my_mkdir(wbse_save_dir)
  !
  kpt_pool = idistribute()
  CALL kpt_pool%init(nkstot,'p','nkstot',.FALSE.,IDIST_BLK)
  !
  IF(kpt_pool%nloc /= nks) CALL errore('wbse_init_setup','unexpected kpt_pool init error',1)
  !
  IF(l_qp_correction) CALL read_qp_eigs()
  !
  ! read ovl_matrix and u_matrix, and compute macroscopic term, if any
  !
  IF(l_bse .OR. l_hybrid_tddft) CALL bse_start()
  !
  do_inexact_krylov = .FALSE.
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE bse_start()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE pwcom,                ONLY : npwx
  USE westcom,              ONLY : l_reduce_io,tau_is_read,tau_all,n_tau,nbnd_occ,nbndval0x,&
                                 & n_trunc_bands,sigma_c_head,sigma_x_head,wbse_epsinfty,&
                                 & l_local_repr,overlap_thr,u_matrix,ovl_matrix,n_bse_idx,idx_matrix
  USE constants,            ONLY : e2,pi
  USE cell_base,            ONLY : omega
  USE types_coulomb,        ONLY : pot3D
  USE wbse_io,              ONLY : read_umatrix_and_omatrix
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE class_idistribute,    ONLY : idistribute,IDIST_BLK
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iks,is,is_g
  INTEGER :: lbnd,ibnd,jbnd,my_ibnd,nbnd_do,nbndval
  INTEGER :: ipair,do_idx
  REAL(DP) :: ovl_value
  !
  ! the divergence term in Fock potential
  !
  sigma_x_head = pot3D%compute_divergence('gb')
  !
  ! compute macroscopic term, it needs macroscopic dielectric constant from input
  !
  sigma_c_head = pot3D%compute_divergence('default')
  sigma_c_head = sigma_c_head * ((1._DP/wbse_epsinfty) - 1._DP)
  !
  WRITE(stdout,'(/,5X,"Macroscopic dielectric constant correction:",f9.5)') sigma_c_head
  !
  nbnd_do = nbndval0x-n_trunc_bands
  !
  band_group = idistribute()
  CALL band_group%init(nbnd_do,'b','nbndval',.FALSE.,IDIST_BLK)
  !
  ! allocate and read unitary matrix and overlap matrix, if any
  !
  IF(l_local_repr) THEN
     !
     ALLOCATE(u_matrix(nbnd_do,nbnd_do,kpt_pool%nloc))
     ALLOCATE(ovl_matrix(nbnd_do,nbnd_do,kpt_pool%nloc))
     !
     DO is = 1,kpt_pool%nloc
        !
        is_g = kpt_pool%l2g(is)
        !
        CALL read_umatrix_and_omatrix(nbnd_do,is_g,u_matrix(:,:,is),ovl_matrix(:,:,is))
        !
     ENDDO
     !
  ENDIF
  !
  ALLOCATE(idx_matrix(nbnd_do*nbnd_do,2,kpt_pool%nloc))
  !
  idx_matrix(:,:,:) = 0
  !
  ALLOCATE(n_bse_idx(kpt_pool%nloc))
  !
  n_bse_idx(:) = 0
  !
  DO iks = 1,kpt_pool%nloc
     !
     nbndval = nbnd_occ(iks)
     do_idx = 0
     !
     DO ibnd = 1,nbndval-n_trunc_bands
        DO jbnd = 1,nbndval-n_trunc_bands
           IF(l_local_repr) THEN
              ovl_value = ovl_matrix(ibnd,jbnd,iks)
           ELSE
              ovl_value = 0._DP
           ENDIF
           !
           IF(l_local_repr) THEN
              IF(ovl_value >= overlap_thr) THEN
                 do_idx = do_idx+1
                 idx_matrix(do_idx,1,iks) = ibnd+n_trunc_bands
                 idx_matrix(do_idx,2,iks) = jbnd+n_trunc_bands
              ENDIF
           ELSE
              do_idx = do_idx+1
              idx_matrix(do_idx,1,iks) = ibnd+n_trunc_bands
              idx_matrix(do_idx,2,iks) = jbnd+n_trunc_bands
           ENDIF
        ENDDO
     ENDDO
     !
     n_bse_idx(iks) = do_idx
     !
  ENDDO
  !
  IF(l_reduce_io) THEN
     !
     ! I/O is reduced by reading tau only once
     !
     ALLOCATE(tau_is_read(nbnd_do,nbnd_do,kpt_pool%nloc))
     !
     tau_is_read(:,:,:) = 0
     n_tau = 0
     !
     DO iks = 1,kpt_pool%nloc
        !
        do_idx = n_bse_idx(iks)
        !
        ! count number of tau needed by my band group
        !
        DO lbnd = 1,band_group%nloc
           !
           my_ibnd = band_group%l2g(lbnd)
           !
           DO ipair = 1, do_idx
              !
              ibnd = idx_matrix(ipair,1,iks)-n_trunc_bands
              jbnd = idx_matrix(ipair,2,iks)-n_trunc_bands
              !
              IF(ibnd == my_ibnd) tau_is_read(ibnd,jbnd,iks) = 1
              !
           ENDDO
           !
        ENDDO
        !
        DO jbnd = 1,nbnd_do
           DO ibnd = jbnd,nbnd_do
              IF(tau_is_read(ibnd,jbnd,iks) == 1 .OR. tau_is_read(jbnd,ibnd,iks) == 1) n_tau = n_tau+1
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
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : intra_pool_comm,my_bgrp_id,me_bgrp
  USE wvfct,                ONLY : et
  USE westcom,              ONLY : nbnd_occ,et_qp,delta_qp,qp_correction
  USE pwcom,                ONLY : nspin,nbnd
  USE json_module,          ONLY : json_file
  USE distribution_center,  ONLY : kpt_pool
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: found
  INTEGER :: is,is_g,ib,nspin_tmp,qp_bands_start,qp_bands_end,qp_bands_mid,nbndval
  CHARACTER(LEN=6) :: labels
  REAL(DP) :: delta
  INTEGER,ALLOCATABLE :: ivals(:)
  REAL(DP),ALLOCATABLE :: rvals(:)
  TYPE(json_file) :: json
  !
  ALLOCATE(et_qp(nbnd,kpt_pool%nloc))
  ALLOCATE(delta_qp(kpt_pool%nloc))
  !
  IF(my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(qp_correction))
     !
     CALL json%get('system.electron.nspin',nspin_tmp,found)
     IF(.NOT. found) CALL errore('read_qp_eigs','nspin not found',1)
     IF(nspin_tmp /= nspin) CALL errore('read_qp_eigs','nspin mismatch',1)
     !
     CALL json%get('input.wfreq_control.qp_bandrange',ivals,found)
     IF(.NOT. found) CALL errore('read_qp_eigs','qp_bandrange not found',1)
     qp_bands_start = ivals(1)
     qp_bands_end = ivals(2)
     !
     DO is = 1,kpt_pool%nloc
        !
        is_g = kpt_pool%l2g(is)
        nbndval = nbnd_occ(is)
        !
        WRITE(labels,'(i6.6)') is_g
        !
        CALL json%get('output.Q.K'//labels//'.eqpSec',rvals,found)
        IF(.NOT. found) CALL errore('read_qp_eigs','eqpSec not found',1)
        et_qp(qp_bands_start:qp_bands_end,is) = rvals/RYTOEV
        !
        qp_bands_mid = (qp_bands_start+nbndval)/2
        delta = 0._DP
        !
        ! delta_qp for valence bands is averaged between qp_bands_start and qp_bands_mid
        !
        ! band index: 1 ... qp_bands_start ... qp_bands_mid ... nbndval ... qp_bands_end ... nbnd
        !                   |<--------- average --------->|
        !
        DO ib = qp_bands_start,qp_bands_mid
           delta = delta+et_qp(ib,is)-et(ib,is)
        ENDDO
        !
        delta = delta/(qp_bands_mid-qp_bands_start+1)
        !
        DO ib = 1,qp_bands_start-1
           et_qp(ib,is) = et(ib,is)+delta
        ENDDO
        !
        qp_bands_mid = (qp_bands_end+nbndval+1)/2
        delta = 0._DP
        !
        ! delta_qp for conduction bands is averaged between qp_bands_mid and qp_bands_end
        !
        ! band index: 1 ... qp_bands_start ... nbndval ... qp_bands_mid ... qp_bands_end ... nbnd
        !                                                  |<-------- average -------->|
        !
        DO ib = qp_bands_mid,qp_bands_end
           delta = delta+et_qp(ib,is)-et(ib,is)
        ENDDO
        !
        delta = delta/(qp_bands_end-qp_bands_mid+1)
        delta_qp(is) = delta
        !
        DO ib = qp_bands_end+1,nbnd
           et_qp(ib,is) = et(ib,is)+delta
        ENDDO
        !
     ENDDO
     !
     CALL json%destroy()
     !
  ENDIF
  !
  CALL mp_bcast(et_qp,0,intra_pool_comm)
  CALL mp_bcast(delta_qp,0,intra_pool_comm)
  !
END SUBROUTINE
