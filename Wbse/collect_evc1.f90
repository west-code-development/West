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
! Victor Yu
!
SUBROUTINE collect_evc1(evc1)
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_bgrp_comm
  USE pwcom,                ONLY : npwx
  USE westcom,              ONLY : evc1_all,nbndval0x,n_trunc_bands
  USE distribution_center,  ONLY : kpt_pool,band_group
#if defined(__NCCL)
  USE west_gpu,             ONLY : gpu_sum,gpu_inter_bgrp_comm
#else
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_bgrp_comm
#endif
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  !
  INTEGER :: iks,lbnd,ibnd,ig,nbnd_do
  INTEGER :: kpt_pool_nloc,band_group_nloc,band_group_myoffset
  !
  nbnd_do = nbndval0x-n_trunc_bands
  kpt_pool_nloc = kpt_pool%nloc
  band_group_nloc = band_group%nloc
  band_group_myoffset = band_group%myoffset
  !
  IF(.NOT. ALLOCATED(evc1_all)) THEN
     ALLOCATE(evc1_all(npwx,nbnd_do,kpt_pool_nloc))
     !$acc enter data create(evc1_all)
  ENDIF
  !
  !$acc kernels present(evc1_all)
  evc1_all(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  !$acc parallel loop collapse(3) present(evc1_all,evc1)
  DO iks = 1,kpt_pool_nloc
     DO lbnd = 1,band_group_nloc
        DO ig = 1,npwx
           ibnd = band_group_myoffset+lbnd
           evc1_all(ig,ibnd,iks) = evc1(ig,lbnd,iks)
        ENDDO
     ENDDO
  ENDDO
  !$acc end parallel
  !
#if defined(__NCCL)
  CALL gpu_sum(evc1_all,npwx*nbnd_do*kpt_pool_nloc,gpu_inter_bgrp_comm)
#else
  !$acc update host(evc1_all)
  CALL mp_sum(evc1_all,inter_bgrp_comm)
  !$acc update device(evc1_all)
#endif
  !
END SUBROUTINE
