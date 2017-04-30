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
SUBROUTINE set_iks_l2g
  !-----------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE pwcom,                  ONLY : nks,nkstot
  USE westcom,                ONLY : iks_l2g
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : inter_pool_comm,npool,my_pool_id
  !
  IMPLICIT NONE
  !
  INTEGER :: iks, ip, my_offset
  INTEGER :: my_nks(0:npool-1)
  !
  ALLOCATE( iks_l2g(nkstot) )
  !
  my_nks=0
  my_nks(my_pool_id) = nks
  CALL mp_sum( my_nks, inter_pool_comm )
  !
  my_offset = 0 
  DO ip = 0, my_pool_id-1
     my_offset = my_offset + my_nks(ip)
  ENDDO
  !
  DO iks = 1, nks
     iks_l2g(iks) = my_offset + iks 
  ENDDO
  !
END SUBROUTINE 
