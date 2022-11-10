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
!-----------------------------------------------------------------------
SUBROUTINE wbse_setup()
  !-----------------------------------------------------------------------
  !
  USE westcom,          ONLY : localization,l_use_localise_repr,l_use_bisection_thr,&
                             & macropol_calculation,l_macropol,solver,l_bse_calculation,&
                             & wbse_calculation,l_davidson,l_lanczos,qp_correction,l_qp_correction,&
                             & spin_excitation,l_bse_triplet,wstat_calculation,n_pdep_times,&
                             & n_pdep_eigen,n_pdep_basis,n_pdep_maxiter,n_pdep_read_from_file,&
                             & trev_pdep_rel,trev_pdep,n_liouville_times,n_liouville_eigen,&
                             & n_liouville_maxiter,n_liouville_read_from_file,trev_liouville_rel,&
                             & trev_liouville,alphapv_dfpt,l_use_ecutrho,nbndval0x,nbnd_occ,&
                             & wbse_save_dir
  USE kinds,            ONLY : DP
  USE types_coulomb,    ONLY : pot3D
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  CALL do_setup()
  !
  SELECT CASE(TRIM(localization))
  CASE('N','n')
     l_use_localise_repr = .FALSE.
     l_use_bisection_thr = .FALSE.
  CASE('B','b')
     l_use_localise_repr = .TRUE.
     l_use_bisection_thr = .TRUE.
  END SELECT
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
  END SELECT
  !
  SELECT CASE(TRIM(solver))
  CASE('BSE','bse')
     l_bse_calculation = .TRUE.
  CASE('TDDFT','tddft')
     l_bse_calculation = .FALSE.
  END SELECT
  !
  SELECT CASE(wbse_calculation)
  CASE('D','d')
     l_davidson = .TRUE.
  CASE('L','l')
     l_lanczos  = .TRUE.
  END SELECT
  !
  IF (TRIM(qp_correction) == '') THEN
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
  CALL pot3D%init('Rho',.FALSE.,'gb')
  !
  CALL set_nbndocc()
  !
  nbndval0x = MAXVAL(nbnd_occ(:))
  !
  CALL wbse_dv_setup(l_bse_calculation)
  !
  CALL my_mkdir(wbse_save_dir)
  !
  IF (l_qp_correction) THEN
     CALL read_qp_eigs()
  ENDIF
  !
  ! read ovl_matrix and u_matrix, and compute macroscopic term, if any
  !
  IF (l_bse_calculation) THEN
     CALL bse_start()
  ENDIF
  !
END SUBROUTINE
!
SUBROUTINE bse_start()
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE pwcom,            ONLY : isk,nks,npwx
  USE westcom,          ONLY : l_reduce_io,tau_is_read,tau_all,n_tau,nbnd_occ,nbndval0x,&
                             & sigma_c_head,sigma_x_head,epsinfty,l_use_localise_repr,overlap_thr,&
                             & u_matrix,ovl_matrix,size_index_matrix_lz,index_matrix_lz
  USE lsda_mod,         ONLY : nspin
  USE constants,        ONLY : e2,pi
  USE cell_base,        ONLY : omega
  USE types_coulomb,    ONLY : pot3D
  USE wbse_io,          ONLY : read_umatrix_and_omatrix
  !
  IMPLICIT NONE
  !
  INTEGER :: do_index, nbndval, tmp_size, is
  INTEGER :: ibnd, jbnd, iks, current_spin
  REAL(DP) :: ovl_value
  !
  ! compute the divergence term in Fock potential, using F-G method
  !
  CALL pot3D%init('Rho',.FALSE.,'gb')
  !
  sigma_x_head = pot3D%div
  !
  ! compute macroscopic term, it needs macroscopic dielectric constant
  ! from input.
  !
  sigma_c_head = ((1._DP/epsinfty) - 1._DP) * (2._DP*e2/pi) * ((6._DP*pi*pi/omega)**(1._DP/3._DP))
  !
  WRITE(stdout,'(/,5X,"Macroscopic dielectric constant correction:",f9.5)') sigma_c_head
  !
  ! allocate and read unitary matrix and overlap matrix, if any
  !
  IF(l_use_localise_repr) THEN
     !
     ALLOCATE(u_matrix(nbndval0x,nbndval0x,nspin))
     ALLOCATE(ovl_matrix(nbndval0x,nbndval0x,nspin))
     !
     DO is = 1,nspin
        CALL read_umatrix_and_omatrix(nbndval0x,is,u_matrix(:,:,is),ovl_matrix(:,:,is))
     ENDDO
     !
  ENDIF
  !
  !IF(.NOT. use_wstat_pdep) THEN
     !
     ! Using coupling approach, single k and q
     ! define an index_matrix_lz, for bse_kernel paralel
     !
     tmp_size = nbndval0x*nbndval0x
     ALLOCATE(index_matrix_lz(tmp_size,2,nspin))
     !
     index_matrix_lz(:,:,:) = 0
     !
     ALLOCATE(size_index_matrix_lz(nspin))
     !
     size_index_matrix_lz(:) = 0
     !
     DO iks = 1,nks
        !
        nbndval = nbnd_occ(iks)
        current_spin = isk(iks)
        do_index = 0
        !
        DO ibnd = 1,nbndval
           DO jbnd = 1,nbndval
              IF(l_use_localise_repr) THEN
                 ovl_value = ovl_matrix(ibnd,jbnd,current_spin)
              ELSE
                 ovl_value = 0._DP
              ENDIF
              !
              IF(l_use_localise_repr) THEN
                 IF(ovl_value >= overlap_thr) THEN
                    do_index = do_index + 1
                    index_matrix_lz(do_index,1,current_spin) = ibnd
                    index_matrix_lz(do_index,2,current_spin) = jbnd
                 ENDIF
              ELSE
                 do_index = do_index + 1
                 index_matrix_lz(do_index,1,current_spin) = ibnd
                 index_matrix_lz(do_index,2,current_spin) = jbnd
              ENDIF
           ENDDO
        ENDDO
        !
        size_index_matrix_lz(current_spin) = do_index
        !
     ENDDO
     !
!  ELSE
!     !
!     ! Using coupling approach, single k and q
!     ! define an index_matrix_lz, for bse_kernel paralel
!     !
!     tmp_size = nbndval0x*nbndval0x*n_pdep_eigen
!     ALLOCATE(index_matrix_lz(tmp_size,3,nspin))
!     !
!     index_matrix_lz(:,:,:) = 0
!     !
!     ALLOCATE(size_index_matrix_lz(nspin))
!     !
!     size_index_matrix_lz(:) = 0
!     !
!     DO iks = 1,nks
!        !
!        nbndval = nbnd_occ(iks)
!        current_spin = isk(iks)
!        do_index = 0
!        !
!        DO alnd = 1,n_pdep_eigen
!           DO ibnd = 1,nbndval
!              DO jbnd = 1,nbndval
!                 IF(l_use_localise_repr) THEN
!                    ovl_value = ovl_matrix(ibnd,jbnd,current_spin)
!                 ELSE
!                    ovl_value = 0._DP
!                 ENDIF
!                 !
!                 IF(ovl_value >= overlap_thr) THEN
!                    IF (gamma_only) THEN
!                       do_index = do_index + 1
!                       index_matrix_lz(do_index,1,current_spin) = ibnd
!                       index_matrix_lz(do_index,2,current_spin) = jbnd
!                       index_matrix_lz(do_index,3,current_spin) = alnd
!                    ENDIF
!                 ELSE
!                    do_index = do_index + 1
!                    index_matrix_lz(do_index,1,current_spin) = ibnd
!                    index_matrix_lz(do_index,2,current_spin) = jbnd
!                    index_matrix_lz(do_index,3,current_spin) = alnd
!                 ENDIF
!              ENDDO
!           ENDDO
!        ENDDO
!        !
!        size_index_matrix_lz(current_spin) = do_index
!        !
!     ENDDO
!     !
!  ENDIF
  !
  IF(l_reduce_io) THEN
     !
     do_index = SUM(size_index_matrix_lz)
     !
     ALLOCATE(tau_is_read(nbndval0x,nbndval0x,nspin))
     ALLOCATE(tau_all(npwx,do_index))
     !
     tau_is_read(:,:,:) = 0
     n_tau = 0
     !
  ENDIF
  !
END SUBROUTINE
