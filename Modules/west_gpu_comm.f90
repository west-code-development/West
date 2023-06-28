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
! Marco Govoni
!
MODULE west_gpu_comm
  !
  USE kinds,                  ONLY : DP
#if defined(__NCCL)
  USE nccl
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: gpu_comm_start
  PUBLIC :: gpu_comm_end
  PUBLIC :: gpu_sum
  PUBLIC :: gpu_inter_image_comm
  PUBLIC :: gpu_inter_pool_comm
  PUBLIC :: gpu_inter_bgrp_comm
  PUBLIC :: gpu_intra_bgrp_comm
  !
  TYPE(ncclComm) :: gpu_inter_image_comm
  TYPE(ncclComm) :: gpu_inter_pool_comm
  TYPE(ncclComm) :: gpu_inter_bgrp_comm
  TYPE(ncclComm) :: gpu_intra_bgrp_comm
  !
  INTERFACE gpu_sum
    MODULE PROCEDURE sum_r8,sum_c16
  END INTERFACE
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE gpu_comm_start()
    !-----------------------------------------------------------------------
    !
    USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_pool_comm,npool,&
                                   & my_pool_id,inter_bgrp_comm,nbgrp,my_bgrp_id,intra_bgrp_comm,&
                                   & nproc_bgrp,me_bgrp
    USE mp,                   ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    ! Workspace
    !
    TYPE(ncclUniqueId) :: gpu_inter_image_id
    TYPE(ncclUniqueId) :: gpu_inter_pool_id
    TYPE(ncclUniqueId) :: gpu_inter_bgrp_id
    TYPE(ncclUniqueId) :: gpu_intra_bgrp_id
    TYPE(ncclResult) :: nccl_result
    !
    IF(my_image_id == 0) nccl_result = ncclGetUniqueId(gpu_inter_image_id)
    IF(my_pool_id == 0) nccl_result = ncclGetUniqueId(gpu_inter_pool_id)
    IF(my_bgrp_id == 0) nccl_result = ncclGetUniqueId(gpu_inter_bgrp_id)
    IF(me_bgrp == 0) nccl_result = ncclGetUniqueId(gpu_intra_bgrp_id)
    !
    CALL mp_bcast(gpu_inter_image_id%internal,0,inter_image_comm)
    CALL mp_bcast(gpu_inter_pool_id%internal,0,inter_pool_comm)
    CALL mp_bcast(gpu_inter_bgrp_id%internal,0,inter_bgrp_comm)
    CALL mp_bcast(gpu_intra_bgrp_id%internal,0,intra_bgrp_comm)
    !
    nccl_result = ncclCommInitRank(gpu_inter_image_comm,nimage,gpu_inter_image_id,my_image_id)
    nccl_result = ncclCommInitRank(gpu_inter_pool_comm,npool,gpu_inter_pool_id,my_pool_id)
    nccl_result = ncclCommInitRank(gpu_inter_bgrp_comm,nbgrp,gpu_inter_bgrp_id,my_bgrp_id)
    nccl_result = ncclCommInitRank(gpu_intra_bgrp_comm,nproc_bgrp,gpu_intra_bgrp_id,me_bgrp)
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE gpu_comm_end()
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! Workspace
    !
    TYPE(ncclResult) :: nccl_result
    !
    nccl_result = ncclCommDestroy(gpu_inter_image_comm)
    nccl_result = ncclCommDestroy(gpu_inter_pool_comm)
    nccl_result = ncclCommDestroy(gpu_inter_bgrp_comm)
    nccl_result = ncclCommDestroy(gpu_intra_bgrp_comm)
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE sum_r8(buf,length,comm)
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    REAL(DP), INTENT(INOUT) :: buf(*)
    INTEGER, INTENT(IN) :: length
    TYPE(ncclComm), INTENT(IN) :: comm
    !
    ! Workspace
    !
    TYPE(ncclResult) :: nccl_result
    !
    !$acc host_data use_device(buf)
    nccl_result = ncclAllReduce(buf,buf,length,ncclDouble,ncclSum,comm,0)
    !$acc end host_data
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE sum_c16(buf,length,comm)
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    COMPLEX(DP), INTENT(INOUT) :: buf(*)
    INTEGER, INTENT(IN) :: length
    TYPE(ncclComm), INTENT(IN) :: comm
    !
    ! Workspace
    !
    TYPE(ncclResult) :: nccl_result
    !
    ! No native support for complex data type
    !
    !$acc host_data use_device(buf)
    nccl_result = ncclAllReduce(C_DEVLOC(buf),C_DEVLOC(buf),length*2,ncclDouble,ncclSum,comm,0)
    !$acc end host_data
    !
  END SUBROUTINE
  !
#endif
END MODULE
