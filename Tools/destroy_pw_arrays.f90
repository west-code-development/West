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
SUBROUTINE destroy_pw_arrays( )
  !-----------------------------------------------------------------------
  !
  USE pwcom
  USE kinds,                  ONLY : DP
  USE scf,                    ONLY : rho, rho_core, v, vltot, vrs, kedtau
  USE uspp,                   ONLY : vkb
  USE uspp_param,             ONLY : upf
  USE io_files,               ONLY : prefix, iunres, diropn, tmp_dir
  USE io_global,              ONLY : ionode_id,ionode
  USE funct,                  ONLY : dft_is_gradient, dmxc
  USE dfunct,                 ONLY : newd
  USE fft_base,               ONLY : dfftp
  USE mp,                     ONLY : mp_barrier,mp_sum,mp_bcast,mp_min,mp_max
  USE mp_global,              ONLY : my_image_id
  USE wavefunctions_module,   ONLY : evc
  USE klist,                  ONLY : nelec
  USE constants,              ONLY : degspin
  USE io_global,              ONLY : stdout
  USE io_files,               ONLY : prefix, iunpun
  USE becmod,                 ONLY : becp,allocate_bec_type
  USE uspp,                   ONLY : vkb, nkb
  USE control_flags,          ONLY : io_level
  USE noncollin_module,       ONLY : noncolin,npol
  USE buffers,                ONLY : close_buffer
  USE westcom,                ONLY : iuwfc
  !
  IMPLICIT NONE
  !
  CALL start_clock('des_pw_ar')
  !
  IF(my_image_id==0) CALL close_buffer(iuwfc,'delete')
  !
  CALL stop_clock('des_pw_ar')
  !
END SUBROUTINE 
