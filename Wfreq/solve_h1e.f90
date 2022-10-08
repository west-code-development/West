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
SUBROUTINE solve_h1e()
  !-----------------------------------------------------------------------
  !
  USE westcom,              ONLY : n_bands,n_pairs,h1e,sigma_exx_full,sigma_corr_full,sigma_exx,&
                                 & sigma_hf_full,qp_bands,ijpmap
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,et
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE wfreq_db,             ONLY : qdet_db_write_h1e
  !
  IMPLICIT NONE
  !
  INTEGER :: is,ib
  !
  ! Compute 1-electron integrals h1e(i,j) = <i|H|j> where H is the 1-electron Hamiltonian
  !
  CALL io_push_title("h1e")
  !
  ! You need to modify this file if you plan to use things like
  ! Hubbard U, electric field, meta-GGA, etc. that are not considered here.
  ! If any of these are requested, the code stops in wfreq_setup.
  !
  ALLOCATE( h1e(n_pairs,nspin) )
  !
  h1e(:,:) = 0._DP
  !
  DO is = 1, nspin
     DO ib = 1, n_bands
        h1e(ijpmap(ib,ib),is) = et(qp_bands(ib),is)
     ENDDO
  ENDDO
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c}
  !
  h1e(:,:) = h1e + sigma_hf_full + REAL(sigma_corr_full,KIND=DP)
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc}
  !
  CALL compute_hartree_double_counting(h1e)
  !
  CALL calc_exx2(sigma_exx, .TRUE.)
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc} - \Sigma^{x}_{dc}
  !
  h1e(:,:) = h1e - sigma_exx_full
  !
  ! Call solve_qp with
  ! l_secant = .FALSE.
  ! l_generate_plot = .FALSE.
  ! l_QDET = .TRUE.
  !
  CALL solve_qp( .FALSE., .FALSE., .TRUE. )
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc} - \Sigma^{x}_{dc} - \Sigma^{c}_{dc}
  !
  h1e(:,:) = h1e - REAL(sigma_corr_full,KIND=DP)
  !
  CALL qdet_db_write_h1e(h1e)
  !
  CALL io_push_bar()
  !
  DEALLOCATE( h1e )
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
  REAL(DP), INTENT(INOUT) :: h1e_tmp(n_pairs,nspin)
  !
  INTEGER :: is1, is2, ipair, jpair, ib, ib_index
  REAL(DP) :: prefactor
  !
  ! if nspin=1, a factor of 2 has to be applied to the occupation
  !
  IF (nspin == 1) THEN
     prefactor = 2._DP
  ELSE
     prefactor = 1._DP
  ENDIF
  !
  DO is1 = 1, nspin
     DO ipair = 1, n_pairs
        DO is2 = 1, nspin
           DO ib = 1, n_bands
              jpair = ijpmap(ib,ib)
              ib_index = qp_bands(ib)
              h1e_tmp(ipair,is1) = h1e_tmp(ipair,is1) &
              & - prefactor*eri_w(ipair,jpair,is1,is2)*occupation(ib_index,is2)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
