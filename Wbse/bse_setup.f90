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
SUBROUTINE bse_init()
     !
     USE kinds,              ONLY : DP
     USE io_global,          ONLY : stdout
     USE bse_module,         ONLY : u_matrix,ovl_matrix,l_wannier_repr,ngm_g_max,ovl_thr,&
                                    size_index_matrix_lz,index_matrix_lz
     USE pwcom,              ONLY : isk,nks
     USE westcom,            ONLY : nbnd_occ,nbndval0x,sigma_c_head,sigma_x_head,epsinfty
     USE lsda_mod,           ONLY : nspin
     USE mp,                 ONLY : mp_max
     USE gvect,              ONLY : ngm,ig_l2g
     USE mp_global,          ONLY : intra_bgrp_comm
     USE constants,          ONLY : e2,pi
     USE cell_base,          ONLY : omega
     USE types_coulomb,      ONLY : pot3D
     !
     IMPLICIT NONE
     !
     INTEGER     :: do_index, nbndval, tmp_size, is
     INTEGER     :: ibnd, jbnd, iks, current_spin
     REAL(DP)    :: ovl_value
     !
     ! compute the divergence term in Fock potential, using F-G method
     !
     CALL pot3D%init('Rho',.FALSE.,'gb')
     CALL pot3D%print_divergence()
     sigma_x_head = pot3D%div
     !
     ! compute macroscopic term, it needs macroscopic dielectric constant
     ! from input.
     !
     sigma_c_head = ((1.0_DP/epsinfty) - 1.0_DP)*(2.0_DP*e2/pi) * &
                    ((6.0_DP*pi*pi/omega)**(1.0_DP/3.0_DP))
     !
     WRITE(stdout,'(/,5X,"Macroscopic dielectric constant correction:", f9.5)') sigma_c_head
     !
     ! set max ngm_g for reading pot file
     !
     ngm_g_max = MAXVAL(ig_l2g(1:ngm))
     CALL mp_max(ngm_g_max,intra_bgrp_comm)
     !
     ! allocate and read unitary matrix and overlap matrix, if any
     !
     IF (l_wannier_repr) THEN
        !
        ALLOCATE (u_matrix (nbndval0x, nbndval0x, nspin))
        ALLOCATE (ovl_matrix (nbndval0x, nbndval0x, nspin))
        !
        DO is = 1, nspin
           !
           CALL read_umatrix_and_omatrix(nbndval0x,is,u_matrix(:,:,is),ovl_matrix(:,:,is))
           !
        ENDDO
        !
     ENDIF
     !
     !IF (.NOT.use_wstat_pdep) THEN
        !
        ! Using Coupling approach, single k and q
        ! define an index_matrix_lz, for bse_kernel paralel
        !
        tmp_size = nbndval0x*nbndval0x
        ALLOCATE (index_matrix_lz(tmp_size, 2, nspin))
        !
        index_matrix_lz(:,:,:) = 0.0
        !
        ALLOCATE(size_index_matrix_lz(nspin))
        !
        size_index_matrix_lz(:) = 0
        !
        DO iks = 1, nks
           !
           nbndval = nbnd_occ(iks)
           !
           current_spin = isk(iks)
           !
           do_index = 0
           !
           DO ibnd = 1, nbndval, 1
              !
              DO jbnd = 1, nbndval, 1
                 !
                 IF (l_wannier_repr) THEN
                    !
                    ovl_value = ovl_matrix(ibnd,jbnd,current_spin)
                    !
                 ELSE
                    !
                    ovl_value = 0.0_DP
                    !
                 ENDIF
                 !
                 IF (l_wannier_repr) THEN
                    !
                    IF (ovl_value >= ovl_thr) THEN
                       !
                       do_index = do_index + 1
                       !
                       index_matrix_lz(do_index, 1, current_spin) = ibnd
                       index_matrix_lz(do_index, 2, current_spin) = jbnd
                       !
                    ENDIF
                    !
                 ELSE
                    !
                    do_index = do_index + 1
                    !
                    index_matrix_lz(do_index, 1, current_spin) = ibnd
                    index_matrix_lz(do_index, 2, current_spin) = jbnd
                    !
                 ENDIF
                 !
              ENDDO
              !
           ENDDO
           !
           size_index_matrix_lz(current_spin) = do_index
           !
        ENDDO
        !
     !ELSE
        !
        ! Using Coupling approach, single k and q
        ! define an index_matrix_lz, for bse_kernel paralel
        !
!        tmp_size = nbndval0x*nbndval0x*n_pdep_eigen
!        ALLOCATE (index_matrix_lz(tmp_size, 3, nspin))
!        !
!        index_matrix_lz(:,:,:) = 0.0
!        !
!        ALLOCATE(size_index_matrix_lz(nspin))
!        !
!        size_index_matrix_lz(:) = 0
!        !
!        DO iks = 1, nks
!           !
!           nbndval = nbnd_occ(iks)
!           !
!           current_spin = isk(iks)
!           !
!           do_index = 0
!           !
!           DO alnd = 1, n_pdep_eigen
!              !
!              DO ibnd = 1, nbndval
!                 !
!                 DO jbnd = 1, nbndval
!                    !
!                    IF (l_wannier_repr) THEN
!                       !
!                       ovl_value = ovl_matrix(ibnd,jbnd,current_spin)
!                       !
!                    ELSE
!                       !
!                       ovl_value = 0.0_DP
!                       !
!                    ENDIF
!                    !
!                    IF (ovl_value >= ovl_thr ) THEN
!                       !
!                       IF (gamma_only) THEN
!                          !
!                          do_index = do_index + 1
!                          !
!                          index_matrix_lz(do_index, 1, current_spin) = ibnd
!                          index_matrix_lz(do_index, 2, current_spin) = jbnd
!                          index_matrix_lz(do_index, 3, current_spin) = alnd
!                          !
!                       ENDIF
!                       !
!                    ELSE
!                       !
!                       do_index = do_index + 1
!                       !
!                       index_matrix_lz(do_index, 1, current_spin) = ibnd
!                       index_matrix_lz(do_index, 2, current_spin) = jbnd
!                       index_matrix_lz(do_index, 3, current_spin) = alnd
!                       !
!                    ENDIF
!                    !
!                 ENDDO
!                 !
!              ENDDO
!              !
!           ENDDO
!           !
!           size_index_matrix_lz(current_spin) = do_index
!           !
!        ENDDO
!        !
!     ENDIF
     !
     ! allocate and read pdep eigen-pots and -values from WSTAT output, if any
     !
     !IF (bse_pdel) THEN
        !
        !
     !ENDIF
     !
END SUBROUTINE
!
SUBROUTINE wbse_clear()
     !
     USE eqv,                ONLY: dmuxc
     USE bse_module,         ONLY: u_matrix,ovl_matrix,et_qp,evc_ks
     USE bse_module,         ONLY: index_matrix_lz, size_index_matrix_lz
     !
     IMPLICIT NONE
     !
     IF (ALLOCATED(u_matrix)) DEALLOCATE (u_matrix)
     IF (ALLOCATED(ovl_matrix)) DEALLOCATE (ovl_matrix)
     IF (ALLOCATED(dmuxc)) DEALLOCATE(dmuxc)
     IF (ALLOCATED(et_qp)) DEALLOCATE(et_qp)
     IF (ALLOCATED(evc_ks)) DEALLOCATE(evc_ks)
     IF (ALLOCATED(index_matrix_lz)) DEALLOCATE(index_matrix_lz)
     IF (ALLOCATED(size_index_matrix_lz)) DEALLOCATE(size_index_matrix_lz)
     !
END SUBROUTINE
