!
! Copyright (C) 2015-2016 M. Govoni 
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
SUBROUTINE exx_go()
  !-----------------------------------------------------------------------
  !
  USE io_global,              ONLY : stdout
  USE kinds,                  ONLY : DP
  USE funct,                  ONLY : dft_is_hybrid,init_dft_exxrpa,stop_exx, &
                                     &get_screening_parameter,get_exx_fraction,start_exx,get_gau_parameter
  USE fft_base,               ONLY : dfftp,dffts 
  USE exx,                    ONLY : x_gamma_extrapolation,exxdiv_treatment,exx_grid_init,exx_div_check,exx_divergence,&
                                     &deallocate_exx,exxinit,vexx,exx_restart,erfc_scrlen,exxdiv,exxalfa,gau_scrlen,&
                                     &exxdiv_treatment, ecutfock
  USE gvecw,                  ONLY : ecutwfc
  USE wvfct,                  ONLY : nbnd, npwx
  USE noncollin_module,       ONLY : npol
  USE io_files,               ONLY : nwordwfc, iunwfc, tmp_dir, wfc_dir
  USE buffers,                ONLY : open_buffer,get_buffer,close_buffer
  USE control_flags,          ONLY : io_level
  USE westcom,                ONLY : l_minimize_exx_if_active
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: exst
  !
  !
  IF( dft_is_hybrid() ) THEN
     exxdiv_treatment='gb'
     !exx_nwordwfc=2*dffts%nnr
     erfc_scrlen = get_screening_parameter()
     gau_scrlen = get_gau_parameter()
     exxalfa = get_exx_fraction()
     IF( l_minimize_exx_if_active ) ecutfock = ecutwfc
!     CALL exx_restart(.true.)
!
     WRITE(stdout, '(7X,"** WARNING : EXX-alpha            = ", f14.6)') exxalfa
     WRITE(stdout, '(7X,"** WARNING : EXX-erfc_scrlen      = ", f14.6)') erfc_scrlen
     WRITE(stdout, '(7X,"** WARNING : EXX-gau_scrlen       = ", f14.6)') gau_scrlen
     WRITE(stdout, '(7X,"** WARNING : EXX-ecutfock         = ", f14.6)') ecutfock
     WRITE(stdout, '(7X,"** WARNING : EXX-exxdiv_treatment = ", a14  )') exxdiv_treatment
     !
     wfc_dir = tmp_dir
     nwordwfc = nbnd * npwx * npol
     io_level = 1
     CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
     !
     CALL start_exx()
     CALL weights()
     CALL exx_grid_init
     CALL exx_div_check
     CALL exxinit
     exxdiv = exx_divergence() 
     WRITE(stdout, '(7X,"** WARNING : EXX-exxdiv           = ", f14.6)') exxdiv
     !
     CALL close_buffer( iunwfc, 'KEEP' )
     !
  ENDIF
  !
END SUBROUTINE 
