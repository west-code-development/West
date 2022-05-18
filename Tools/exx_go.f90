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
SUBROUTINE exx_go()
  !-----------------------------------------------------------------------
  !
  USE io_global,              ONLY : stdout
  USE xc_lib,                 ONLY : xclib_dft_is,get_screening_parameter,xclib_get_exx_fraction,&
                                     get_gau_parameter,start_exx
  USE exx,                    ONLY : exxinit,ecutfock,exxalfa,use_ace
  USE exx_base,               ONLY : exxdiv_treatment,exx_grid_init,exx_div_check,exxdiv,&
                                     exx_divergence,exx_mp_init,gau_scrlen,erfc_scrlen
  USE gvecw,                  ONLY : ecutwfc
  USE wvfct,                  ONLY : nbnd,npwx
  USE noncollin_module,       ONLY : npol
  USE io_files,               ONLY : nwordwfc,iunwfc,tmp_dir,wfc_dir
  USE buffers,                ONLY : open_buffer,close_buffer
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
  IF( xclib_dft_is('hybrid') ) THEN
     exxdiv_treatment='gb'
     erfc_scrlen = get_screening_parameter()
     gau_scrlen = get_gau_parameter()
     exxalfa = xclib_get_exx_fraction()
     use_ace = .FALSE.
     IF( l_minimize_exx_if_active ) THEN
        ecutfock = ecutwfc
     ELSE
        ecutfock = ecutwfc*4
     ENDIF
     !
     WRITE(stdout, '(7X,"** WARNING : EXX-use_ace          = ", l1)') use_ace
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
     CALL exx_grid_init()
     ! exx_mp_init necessary when k points are used
     CALL exx_mp_init()
     CALL exx_div_check()
     CALL exxinit(DoLoc=.FALSE.)
     exxdiv = exx_divergence()
     WRITE(stdout, '(7X,"** WARNING : EXX-exxdiv           = ", f14.6)') exxdiv
     CALL close_buffer ( iunwfc, 'KEEP' )
     !
  ENDIF
  !
END SUBROUTINE
