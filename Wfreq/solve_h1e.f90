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
!-----------------------------------------------------------------------
SUBROUTINE solve_h1e()
  !-----------------------------------------------------------------------
  !
  USE westcom,              ONLY : n_bands,n_pairs,proj_c,h1e,iuwfc,lrwfc,&
                                 & sigma_exx_full,sigma_corr_full,wfreq_save_dir,&
                                 & sigma_exx
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : inter_image_comm,my_image_id,mp_bcast
  USE io_push,              ONLY : io_push_title, io_push_bar
  USE ldaU,                 ONLY : lda_plus_u
  USE bp,                   ONLY : lelfield
  USE realus,               ONLY : real_space
  USE control_flags,        ONLY : gamma_only
  USE wfreq_db,             ONLY : qdet_db_write_h1e
  USE mp_world,             ONLY : mpime,root
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), h1e_tmp(:,:)
  INTEGER  :: s,i,iunit
  REAL(DP)  :: ry_to_ha = 0.5_DP
  !
  ! Compute 1-electron integrals h1e(i,j) = <i|H|j> where H is the 1-electron Hamiltonian
  !
  CALL io_push_title("h1e")
  !
  ! You need to modify this file if you plan to use things like
  ! Hubbard U, electric field, meta-GGA, etc. that are not considered here.
  
  ! Initial checks
  IF ( real_space ) CALL errore("solve_h1e", "real_space not supported", 1)
  IF ( lda_plus_u ) CALL errore("solve_h1e", "lda_plus_u not supported", 1)
  IF ( lelfield ) CALL errore("solve_h1e", "lelfield not supported", 1)
  IF ( .NOT. gamma_only ) CALL errore("solve_h1e", "only support gamma-point case", 1)
  !
  ALLOCATE( psi(npwx,n_bands,nspin), hpsi(npwx,n_bands,nspin) )
  ALLOCATE( h1e(n_pairs,nspin), h1e_tmp(n_pairs,nspin) )
  !
  psi = 0._DP
  h1e = 0._DP
  h1e_tmp = 0._DP
  !
  DO s = 1, nspin
    !
    psi(:,:,s) = proj_c(:,:,s)
    !
  ENDDO
  !
  CALL compute_hks(psi, hpsi, h1e_tmp)
  ! H1e = H^{KS}
  h1e = h1e_tmp
  !
  CALL compute_vxc(psi, hpsi, h1e_tmp)
  ! H1e = H^{KS} - V_{xc}
  h1e = h1e - h1e_tmp
  !
  CALL compute_vxx(psi, hpsi, h1e_tmp)
  ! H1e = H^{KS} - V_{xc} - V_{xx}
  h1e = h1e - h1e_tmp
  ! H1e = H^{KS} - V_{xc} - V_{xx} + \Sigma^{x}
  h1e = h1e + sigma_exx_full
  !
  CALL calc_exx2(sigma_exx, .TRUE.)
  ! H1e = H^{KS} - V_{xc} - V_{xx} + \Sigma^{x} - \Sigma^{x}_{dc}
  h1e = h1e - sigma_exx_full
  ! H1e = H^{KS} - V_{xc} - V_{xx} + \Sigma^{x} - \Sigma^{x}_{dc} + \Sigma^{c}
  h1e = h1e + REAL(sigma_corr_full)
  ! Call to solve_qp with
  ! l_secant = .FALSE.
  ! l_generate_plot = .FALSE.
  ! l_QDET = .TRUE. 
  CALL solve_qp( .FALSE., .FALSE., .TRUE. )
  ! H1e = H^{KS} - V_{xc} - V_{xx} + \Sigma^{x} - \Sigma^{x}_{dc} + \Sigma^{c} - \Sigma^{c}_{dc}
  h1e = h1e - REAL(sigma_corr_full)
  ! write H1e to JSON file
  CALL qdet_db_write_h1e(h1e)
  !
  CALL io_push_bar()
  !
  DEALLOCATE( psi, hpsi )
  DEALLOCATE( h1e, h1e_tmp )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------  
SUBROUTINE compute_hks(psi, hpsi, h1e_tmp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx,nbnd
  USE westcom,              ONLY : n_bands,n_pairs
  USE uspp,                 ONLY : nkb
  USE becmod,               ONLY : becp,allocate_bec_type,is_allocated_bec_type,deallocate_bec_type
  USE wvfct,                ONLY : current_k
  USE lsda_mod,             ONLY : current_spin
  USE wavefunctions, ONLY : evc
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: psi(npwx,n_bands,nspin)
  COMPLEX(DP), INTENT(OUT) :: hpsi(npwx,n_bands,nspin)
  REAL(DP), INTENT(OUT)    :: h1e_tmp(n_pairs,nspin)
  INTEGER  :: i, s
  !
  ! Computes matrix elements of the Kohn-Sham Hamiltonian in the pair-basis of Kohn-Sham eigenstates. 
  !
  hpsi = 0._DP
  !
  IF ( is_allocated_bec_type( becp ) )  CALL deallocate_bec_type( becp )
  CALL allocate_bec_type(nkb, nbnd, becp)
  !
  DO s = 1, nspin
     !
     current_k = s
     current_spin = s
     !
     CALL h_psi(npwx, npw, n_bands, psi(:,:,s), hpsi(:,:,s))
     !
  ENDDO
  !
  ! compute integrals from psi and hpsi and transform into pair-basis.
  CALL compute_integrals(psi, hpsi, h1e_tmp)
  !
  CALL deallocate_bec_type(becp)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------  
SUBROUTINE compute_vxc(psi, hpsi, h1e_tmp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx
  USE westcom,              ONLY : n_bands,n_pairs
  USE scf,                  ONLY : scf_type,rho,rho_core,rhog_core, &
                                 & create_scf_type,destroy_scf_type
  USE xc_lib,               ONLY : exx_is_active
  USE exx,                  ONLY : vexx,vexxace_gamma,use_ace
  USE wvfct,                ONLY : current_k
  USE lsda_mod,             ONLY : current_spin
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: psi(npwx,n_bands,nspin)
  COMPLEX(DP), INTENT(OUT) :: hpsi(npwx,n_bands,nspin)
  REAL(DP), INTENT(OUT)    :: h1e_tmp(n_pairs,nspin)
  INTEGER  :: i, s
  TYPE(scf_type) :: vhxc
  REAL(DP) :: ehart, charge, etxc, vtxc, eth, etotefield, ee
  !
  ! Computes Matrix elements of V_{xc} in the pair-basis.
  !
  CALL create_scf_type(vhxc, do_not_allocate_becsum = .true.)
  !
  ! XC (semilocal)
  !
  hpsi = 0._DP
  vhxc%of_r = 0._DP
  !
  CALL v_xc(rho, rho_core, rhog_core, etxc, vtxc, vhxc%of_r)
  !
  DO s = 1, nspin
    !
    CALL vloc_psi_gamma(npwx, npw, n_bands, psi(:,:,s), vhxc%of_r(:,s), hpsi(:,:,s))
    !
  ENDDO
  !
  CALL compute_integrals(psi, hpsi, h1e_tmp)
  !
  CALL destroy_scf_type(vhxc)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------  
SUBROUTINE compute_vxx(psi, hpsi, h1e_tmp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx
  USE westcom,              ONLY : n_bands,n_pairs
  USE scf,                  ONLY : scf_type,rho,rho_core,rhog_core, &
                                 & create_scf_type,destroy_scf_type
  USE xc_lib,               ONLY : exx_is_active
  USE exx,                  ONLY : vexx,vexxace_gamma,use_ace
  USE wvfct,                ONLY : current_k
  USE lsda_mod,             ONLY : current_spin
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: psi(npwx,n_bands,nspin)
  COMPLEX(DP), INTENT(OUT) :: hpsi(npwx,n_bands,nspin)
  REAL(DP), INTENT(OUT)    :: h1e_tmp(n_pairs,nspin)
  INTEGER  :: i, s
  TYPE(scf_type) :: vhxc
  REAL(DP) :: ehart, charge, etxc, vtxc, eth, etotefield, ee
  !
  CALL create_scf_type(vhxc, do_not_allocate_becsum = .true.)
  !
  ! EXX
  !
  hpsi = 0._DP
  !
  IF ( exx_is_active() ) THEN
     !
     DO s = 1, nspin
        !
        current_k = s
        current_spin = s
        !
        IF ( use_ace) THEN
           CALL vexxace_gamma(npwx, n_bands, psi(:,:,s), ee, hpsi(:,:,s))
        ELSE
           CALL vexx(npwx, npw, n_bands, psi(:,:,s), hpsi(:,:,s))
        END IF
     END DO
     !
  END IF
  !
  CALL compute_integrals(psi, hpsi, h1e_tmp)
  !
  CALL destroy_scf_type(vhxc)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------  
SUBROUTINE compute_integrals(psi, hpsi, h1e)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx
  USE westcom,              ONLY : n_bands,n_pairs,ijpmap,wfreq_save_dir
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE mp_world,             ONLY : mpime,root  
  USE gvect,                ONLY : gstart
  USE noncollin_module,     ONLY : npol
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: psi(npwx,n_bands,nspin)
  COMPLEX(DP), INTENT(OUT) :: hpsi(npwx,n_bands,nspin)
  REAL(DP), INTENT(OUT)    :: h1e(n_pairs,nspin)
  !
  REAL(DP),ALLOCATABLE     :: h1e_tmp(:,:,:)
  INTEGER   :: i, j, s, ib_index, jb_Index
  INTEGER   :: iunit
  REAL(DP)  :: ry_to_ha = 0.5_DP
  REAL(DP), EXTERNAL    :: DDOT
  !
  ALLOCATE( h1e_tmp(n_bands,n_bands,nspin) )
  h1e_tmp = 0._DP
  !
  ! Following h_psi subroutine, set imaginary part of G=0 component of hpsi to zero
  !
  IF ( gstart == 2 ) THEN
     hpsi(1,:,:) = CMPLX( DBLE( hpsi(1,:,:) ), 0.D0 ,kind=DP)
  ENDIF
  !
  DO s = 1, nspin
     !
     ! Following regterg subroutine, compute overlap integrals between psi and hpsi
     ! Hamiltonian matrix (or its components) is assumed to be real and symmetric
     !
     CALL glbrak_gamma( psi(:,:,s), hpsi(:,:,s), h1e_tmp(:,:,s), npw, npwx, n_bands, n_bands, n_bands, npol)
     !
  ENDDO
  !
  CALL mp_sum(h1e_tmp, intra_bgrp_comm)
  !
  DO s = 1, nspin    
    DO ib_index = 1, n_bands
       DO jb_index = 1, n_bands
          IF (ib_index > jb_index) CYCLE
          h1e(ijpmap(ib_index,jb_index),s) = h1e_tmp(ib_index,jb_index,s)
       ENDDO
    ENDDO
  ENDDO
  !
  DEALLOCATE(h1e_tmp)
  !
END SUBROUTINE
