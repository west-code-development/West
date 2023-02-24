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
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE solve_h1e()
  !-----------------------------------------------------------------------
  !
  USE westcom,              ONLY : l_enable_off_diagonal,n_bands,n_pairs,qp_bands,ijpmap,&
                                 & sigma_hf_full,sigma_corr_full,sigma_exx_full,sigma_hf,&
                                 & sigma_sc_eqpsec,sigma_exx
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,et
  USE io_push,              ONLY : io_push_title
  USE wfreq_db,             ONLY : qdet_db_write_h1e
  !
  IMPLICIT NONE
  !
  INTEGER :: is,ib
  REAL(DP), ALLOCATABLE :: h1e(:,:)
  REAL(DP), ALLOCATABLE :: h1e_diag(:,:)
  !
  ! Compute 1-electron integrals h1e(i,j) = <i|H|j> where H is the 1-electron Hamiltonian
  !
  CALL io_push_title("h1e")
  !
  ! You need to modify this file if you plan to use things like
  ! Hubbard U, electric field, meta-GGA, etc. that are not considered here.
  ! If any of these are requested, the code stops in wfreq_setup.
  !
  ALLOCATE(h1e(n_pairs,nspin))
  h1e(:,:) = 0._DP
  !
  IF(l_enable_off_diagonal) THEN
     DO is = 1, nspin
        DO ib = 1, n_bands
           h1e(ijpmap(ib,ib),is) = et(qp_bands(ib),is)
        ENDDO
     ENDDO
  ELSE
     ALLOCATE(h1e_diag(n_bands,nspin))
     h1e_diag(:,:) = 0._DP
     !
     DO is = 1, nspin
        DO ib = 1, n_bands
           h1e_diag(ib,is) = et(qp_bands(ib),is)
        ENDDO
     ENDDO
  ENDIF
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c}
  !
  IF(l_enable_off_diagonal) THEN
     h1e(:,:) = h1e + sigma_hf_full + REAL(sigma_corr_full,KIND=DP)
  ELSE
     h1e_diag(:,:) = h1e_diag + sigma_hf + REAL(sigma_sc_eqpsec,KIND=DP)
  ENDIF
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc}
  !
  IF(l_enable_off_diagonal) THEN
     CALL compute_hartree_double_counting(h1e)
  ELSE
     CALL compute_hartree_double_counting_diag(h1e_diag)
  ENDIF
  !
  CALL calc_exx2(sigma_exx, .TRUE.)
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc} - \Sigma^{x}_{dc}
  !
  IF(l_enable_off_diagonal) THEN
     h1e(:,:) = h1e - sigma_exx_full
  ELSE
     h1e_diag(:,:) = h1e_diag - sigma_exx
  ENDIF
  !
  ! Call solve_qp with
  ! l_secant = .FALSE.
  ! l_generate_plot = .FALSE.
  ! l_QDET = .TRUE.
  !
  CALL solve_qp(.FALSE., .FALSE., .TRUE.)
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc} - \Sigma^{x}_{dc} - \Sigma^{c}_{dc}
  !
  IF(l_enable_off_diagonal) THEN
     h1e(:,:) = h1e - REAL(sigma_corr_full,KIND=DP)
  ELSE
     h1e_diag(:,:) = h1e_diag - REAL(sigma_sc_eqpsec,KIND=DP)
  ENDIF
  !
  IF(.NOT. l_enable_off_diagonal) THEN
     DO is = 1, nspin
        DO ib = 1, n_bands
           h1e(ijpmap(ib,ib),is) = h1e_diag(ib,is)
        ENDDO
     ENDDO
  ENDIF
  !
  CALL qdet_db_write_h1e(h1e)
  !
  DEALLOCATE(h1e)
  IF(.NOT. l_enable_off_diagonal) DEALLOCATE(h1e_diag)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_hartree_double_counting(h1e_tmp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin
  USE westcom,              ONLY : n_bands,qp_bands,n_pairs,ijpmap,eri_w,occupation
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: h1e_tmp(n_pairs,nspin)
  !
  INTEGER :: is1, is2, ipair, jpair, jb
  REAL(DP) :: prefactor
  !
  ! if nspin=1, a factor of 2 has to be applied to the occupation
  !
  IF(nspin == 1) THEN
     prefactor = 2._DP
  ELSE
     prefactor = 1._DP
  ENDIF
  !
  DO is1 = 1, nspin
     DO ipair = 1, n_pairs
        DO is2 = 1, nspin
           DO jb = 1, n_bands
              jpair = ijpmap(jb,jb)
              h1e_tmp(ipair,is1) = h1e_tmp(ipair,is1) &
              & - prefactor*eri_w(ipair,jpair,is1,is2)*occupation(qp_bands(jb),is2)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_hartree_double_counting_diag(h1e_diag_tmp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin
  USE westcom,              ONLY : n_bands,qp_bands,ijpmap,eri_w,occupation
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: h1e_diag_tmp(n_bands,nspin)
  !
  INTEGER :: is1, is2, ipair, jpair, ib, jb
  REAL(DP) :: prefactor
  !
  ! if nspin=1, a factor of 2 has to be applied to the occupation
  !
  IF(nspin == 1) THEN
     prefactor = 2._DP
  ELSE
     prefactor = 1._DP
  ENDIF
  !
  DO is1 = 1, nspin
     DO ib = 1, n_bands
        ipair = ijpmap(ib,ib)
        DO is2 = 1, nspin
           DO jb = 1, n_bands
              jpair = ijpmap(jb,jb)
              h1e_diag_tmp(ib,is1) = h1e_diag_tmp(ib,is1) &
              & - prefactor*eri_w(ipair,jpair,is1,is2)*occupation(qp_bands(jb),is2)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
