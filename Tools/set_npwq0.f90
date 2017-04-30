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
SUBROUTINE set_npwq0()
  !-----------------------------------------------------------------------
  !
  USE westcom,   ONLY : npwq0,npwq0x,npwq0_g,l_use_ecutrho,fftdriver
  USE mp,        ONLY : mp_max
  USE mp_global, ONLY : intra_bgrp_comm
  USE gvect,     ONLY : ig_l2g,ngm,ngmx
  USE pwcom,     ONLY : npw,npwx 
  !
  IMPLICIT NONE
  !
  IF( l_use_ecutrho ) THEN 
     npwq0     = ngm
     npwq0x    = ngmx
     fftdriver = 'Dense'
  ELSE
     npwq0     = npw
     npwq0x    = npwx
     fftdriver = 'Wave'
  ENDIF
  ! 
  !ALLOCATE(q0ig_l2g(npwq0))
  !q0ig_l2g(1:npwq0) = ig_l2g(1:npwq0)
  !npwq0_g=MAXVAL(q0ig_l2g(1:npwq0))
  npwq0_g=MAXVAL(ig_l2g(1:npwq0))
  CALL mp_max(npwq0_g,intra_bgrp_comm)
  !
END SUBROUTINE 
