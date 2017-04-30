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
SUBROUTINE init_pw_arrays(ncalbec)
  !-----------------------------------------------------------------------
  !
  USE pwcom
  USE kinds,                  ONLY : DP
  USE scf,                    ONLY : rho, rho_core, v, vltot, vrs, kedtau
  USE uspp,                   ONLY : vkb
  USE uspp_param,             ONLY : upf
  USE io_files,               ONLY : prefix, iunres, diropn, tmp_dir, iunwfc, nwordwfc
  USE io_global,              ONLY : ionode_id,ionode
  USE funct,                  ONLY : dft_is_gradient, dmxc
  USE dfunct,                 ONLY : newd
  USE fft_base,               ONLY : dfftp
  USE mp,                     ONLY : mp_barrier,mp_sum,mp_bcast,mp_min,mp_max
  USE mp_global,              ONLY : intra_bgrp_comm,inter_image_comm,nproc_bgrp,my_image_id
  USE wavefunctions_module,   ONLY : evc
  USE klist,                  ONLY : nelec, nks, xk, ngk, igk_k
  USE constants,              ONLY : degspin
  USE io_global,              ONLY : stdout
  USE io_files,               ONLY : prefix, iunpun
  USE becmod,                 ONLY : becp,allocate_bec_type
  USE uspp,                   ONLY : vkb, nkb
  USE control_flags,          ONLY : io_level
  USE noncollin_module,       ONLY : noncolin,npol
  USE buffers,                ONLY : open_buffer,get_buffer,close_buffer
  USE funct,                  ONLY : dft_is_hybrid,init_dft_exxrpa,stop_exx
  USE exx,                    ONLY : x_gamma_extrapolation,exxdiv_treatment,exx_grid_init,exx_div_check,&
                                     &deallocate_exx,exxinit,vexx
  USE westcom,                ONLY : iuwfc,lrwfc
  USE gvecw,                  ONLY : gcutw
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ncalbec
  !
  INTEGER :: i, l, nt, kpoint, ik
  CHARACTER(256) :: filint
!  REAL(DP) :: rhotot
  LOGICAL :: exst,exst_mem
  !
  CALL start_clock('init_pw_ar')
  !
  !  sum self-consistent part (vr) and local part (vltot) of potential
  !
  CALL hinit0()
  !
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r, dfftp%nnr,nspin,doublegrid)
  !
  !  allocate memory for gradient corrections (if needed)
  !
!  IF ( dft_is_gradient() ) THEN
!     ALLOCATE  ( dvxc_rr(dfftp%nnr,nspin,nspin))
!     ALLOCATE  ( dvxc_sr(dfftp%nnr,nspin,nspin))
!     ALLOCATE  ( dvxc_ss(dfftp%nnr,nspin,nspin))
!     ALLOCATE  ( dvxc_s (dfftp%nnr,nspin,nspin))
!     ALLOCATE  ( grho   (3, dfftp%nnr, nspin))
!  ENDIF
  !
  CALL allocate_bec_type ( nkb, ncalbec, becp )
  !
  !  local potential
  !
!M  CALL init_vloc
  !
!M  CALL init_us_1
  !
  CALL newd
  !
  !  derivative of the xc potential
  !
!  dmuxc(:) = 0.d0
!  DO i = 1,dfftp%nnr
!     rhotot = rho%of_r(i,current_spin)+rho_core(i)
!     IF ( rhotot> 1.d-30 ) dmuxc(i)= dmxc( rhotot)
!     IF ( rhotot<-1.d-30 ) dmuxc(i)=-dmxc(-rhotot)
!  ENDDO
  !
  !  initialize data needed for gradient corrections
  !
!  CALL cg_setupdgc
  !
  !  read wave functions and calculate indices
  !
!M  CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )
!M  REWIND( iunigk )
!M  DO ik = 1, nks
!M     CALL gk_sort( xk(1,ik), ngm, g, gcutw, npw, igk, g2kin )
!M     ngk(ik) = npw
!M     igk_k(1:npw,ik)= igk(1:npw)
!M     IF ( nks > 1 ) WRITE( iunigk ) igk
!M  ENDDO
  !
!  kpoint=1
!  CALL gk_sort (xk(1,kpoint),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)
!  g2kin=g2kin*tpiba2
  !
  !  Kleinman-Bylander PPs
  !
!  CALL init_us_2 (npw, igk, xk(1,kpoint), vkb)
  !
  CALL exx_go()
  !
  iunres=88
  !
  !iuwfc = 20
  lrwfc = nbnd * npwx * npol
  IF(my_image_id==0) THEN
     ! CALL open_buffer ( iuwfc, 'wfc', lrwfc, io_level, exst )
     CALL open_buffer ( iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir )
     CALL get_buffer( evc, lrwfc, iuwfc, 1 )
  ENDIF
  CALL mp_bcast(evc,0,inter_image_comm)
  !
  !   open the wavefunction file (already existing)
  !
  CALL stop_clock('init_pw_ar')
  !
END SUBROUTINE 
