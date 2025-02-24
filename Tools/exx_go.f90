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
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE exx_go()
  !-----------------------------------------------------------------------
  !
  ! See also PW/src/setup.f90: setup_exx
  !
  USE io_global,              ONLY : stdout
  USE xc_lib,                 ONLY : xclib_dft_is,get_screening_parameter,xclib_get_exx_fraction,&
                                   & get_gau_parameter,start_exx
  USE exx,                    ONLY : exxinit,ecutfock,exxalfa,use_ace,nbndproj,xi
  USE exx_base,               ONLY : exxdiv_treatment,exx_grid_init,exx_div_check,exxdiv,&
                                   & exx_divergence,exx_mp_init,gau_scrlen,erfc_scrlen
  USE gvecw,                  ONLY : ecutwfc
  USE pwcom,                  ONLY : nbnd,npwx,nks
  USE noncollin_module,       ONLY : npol
  USE io_files,               ONLY : nwordwfc,iunwfc,tmp_dir,wfc_dir
  USE buffers,                ONLY : open_buffer,close_buffer
  USE control_flags,          ONLY : io_level
  USE westcom,                ONLY : l_minimize_exx_if_active,n_exx_lowrank,westpp_l_compute_tdm,&
                                   & westpp_l_spin_flip,westpp_l_dipole_realspace
  USE mp_global,              ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE mp_exx,                 ONLY : mp_start_exx
  USE mp,                     ONLY : mp_bcast
  USE command_line_options,   ONLY : ntg_,command_line
#if defined(__CUDA)
  USE exx,                    ONLY : xi_d
#endif
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: is_westpp
  LOGICAL :: is_wbse_init
  LOGICAL :: do_exxinit
  LOGICAL :: exst
  LOGICAL, EXTERNAL :: matches
  !
  ! Disable band parallelization (egrp) in vexx, as WEST handles band parallelization separately
  !
  CALL mp_start_exx(1,ntg_,intra_bgrp_comm)
  !
  is_westpp = matches('westpp.x',command_line)
  is_wbse_init = matches('wbse_init.x',command_line)
  !
  ! Initialize EXX only if calling h_psi
  !
  IF(is_westpp .OR. is_wbse_init) THEN
     do_exxinit = .FALSE.
     IF(is_westpp .AND. westpp_l_compute_tdm .AND. (.NOT. westpp_l_spin_flip) &
     & .AND. (.NOT. westpp_l_dipole_realspace)) do_exxinit = .TRUE.
  ELSE
     do_exxinit = .TRUE.
  ENDIF
  !
  IF(xclib_dft_is('hybrid')) THEN
     !
     exxdiv_treatment = 'gb'
     erfc_scrlen = get_screening_parameter()
     gau_scrlen = get_gau_parameter()
     exxalfa = xclib_get_exx_fraction()
     IF(n_exx_lowrank > 0) THEN
        use_ace = .TRUE.
     ELSE
        use_ace = .FALSE.
     ENDIF
     IF(l_minimize_exx_if_active) THEN
        ecutfock = ecutwfc
     ELSE
        ecutfock = ecutwfc*4
     ENDIF
     !
     WRITE(stdout,'(7X,"** WARNING : EXX use_ace          = ",L1)') use_ace
     WRITE(stdout,'(7X,"** WARNING : EXX alpha            = ",F14.6)') exxalfa
     WRITE(stdout,'(7X,"** WARNING : EXX erfc_scrlen      = ",F14.6)') erfc_scrlen
     WRITE(stdout,'(7X,"** WARNING : EXX gau_scrlen       = ",F14.6)') gau_scrlen
     WRITE(stdout,'(7X,"** WARNING : EXX ecutfock         = ",F14.6)') ecutfock
     WRITE(stdout,'(7X,"** WARNING : EXX exxdiv_treatment = ",A)') TRIM(exxdiv_treatment)
     !
     wfc_dir = tmp_dir
     nwordwfc = nbnd*npwx*npol
     io_level = 1
     !
     CALL start_exx()
     CALL exx_grid_init()
     ! exx_mp_init necessary when k points are used
     CALL exx_mp_init()
     IF(do_exxinit) THEN
        IF(use_ace) THEN
           nbndproj = n_exx_lowrank
           ALLOCATE(xi(npwx*npol,nbndproj,nks))
#if defined(__CUDA)
           ALLOCATE(xi_d(npwx*npol,nbndproj))
#endif
           IF(my_image_id == 0) CALL aceinit0()
           CALL mp_bcast(xi,0,inter_image_comm)
           nbndproj = n_exx_lowrank
#if defined (__CUDA)
           IF(nks == 1) xi_d(:,:) = xi(:,:,1)
#endif
        ELSE
           CALL open_buffer(iunwfc,'wfc',nwordwfc,io_level,exst)
           CALL exxinit(DoLoc=.FALSE.)
           CALL close_buffer(iunwfc,'KEEP')
        ENDIF
     ENDIF
     CALL exx_div_check()
     exxdiv = exx_divergence()
     WRITE(stdout,'(7X,"** WARNING : EXX-exxdiv           = ",F14.6)') exxdiv
     !
  ENDIF
  !
END SUBROUTINE
