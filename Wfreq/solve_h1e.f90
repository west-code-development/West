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
  USE westcom,              ONLY : n_bands,n_pairs,h1e,iuwfc,&
                                 & sigma_exx_full,sigma_corr_full,&
                                 & sigma_exx, sigma_hf_full, qp_bands, ijpmap
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,et
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : mp_bcast
  USE io_push,              ONLY : io_push_title, io_push_bar
  USE ldaU,                 ONLY : lda_plus_u
  USE bp,                   ONLY : lelfield
  USE realus,               ONLY : real_space
  USE control_flags,        ONLY : gamma_only
  USE wfreq_db,             ONLY : qdet_db_write_h1e
  !
  IMPLICIT NONE
  !
  INTEGER  :: s, iband
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
  ALLOCATE( h1e(n_pairs,nspin))
  !
  h1e = 0._DP
  !
  !
  DO s = 1, nspin
    DO iband = 1, n_bands
      h1e(ijpmap(iband, iband),s) = et(qp_bands(iband), s)  
    ENDDO
  ENDDO 
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c}
  h1e = h1e + sigma_hf_full + REAL(sigma_corr_full)
  !
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc}
  CALL compute_hartree_double_counting(h1e)
  !
  CALL calc_exx2(sigma_exx, .TRUE.)
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc} - \Sigma^{x}_{dc}
  h1e = h1e - sigma_exx_full
  ! Call to solve_qp with
  ! l_secant = .FALSE.
  ! l_generate_plot = .FALSE.
  ! l_QDET = .TRUE. 
  CALL solve_qp( .FALSE., .FALSE., .TRUE. )
  ! H1e = H^{KS} - V_{xc} + \Sigma^{x} + \Sigma^{c} - V^{H}_{dc} - \Sigma^{x}_{dc} - \Sigma^{c}_{dc}
  h1e = h1e - REAL(sigma_corr_full)
  ! write H1e to JSON file
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
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin
  USE westcom,              ONLY : n_bands, n_pairs, ijpmap, eri_w, occupation, qp_bands
  
  REAL(DP), INTENT(INOUT)    :: h1e_tmp(n_pairs,nspin)
  !
  INTEGER                  :: s1, s2, ipair, jpair, iband, band_index
  REAL(DP) :: prefactor
  
  ! if nspin=1, a factor of 2 has to be applied to the occupation
  IF (nspin == 1) THEN
    prefactor = 2.0_DP
  ELSE
    prefactor = 1.0_DP
  ENDIF

  DO s1 = 1, nspin
    DO ipair = 1, n_pairs
      DO s2= 1, nspin
        DO iband = 1, n_bands
          jpair = ijpmap(iband, iband)
          band_index = qp_bands(iband)
          h1e_tmp(ipair, s1) = h1e_tmp(ipair, s1) - prefactor*eri_w(ipair,jpair,s1,s2)*occupation(band_index,s2)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE
