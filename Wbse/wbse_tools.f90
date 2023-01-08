!
! Copyright (C) 2015-2022 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
MODULE wbse_tools
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTERFACE wbse_build_hr
     MODULE PROCEDURE build_hr_real
  END INTERFACE
  !
  INTERFACE wbse_update_with_vr_distr
     MODULE PROCEDURE update_with_vr_distr_real
  END INTERFACE
  !
  INTERFACE wbse_refresh_with_vr_distr
     MODULE PROCEDURE refresh_with_vr_distr_real
  END INTERFACE
  !
  INTERFACE apply_preconditioning_dvg
     MODULE PROCEDURE preconditioner_complex
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE build_hr_real(ag,bg,l2_s,l2_e,c_distr,g_e)
      !------------------------------------------------------------------------
      !
      !  c_distr = < ag | bg >
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm,&
                                     & inter_bgrp_comm,my_bgrp_id,inter_pool_comm,my_pool_id
      USE mp,                   ONLY : mp_sum,mp_bcast
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,nbndval0x
      USE gvect,                ONLY : gstart
      USE west_mp,              ONLY : mp_circular_shift_left_c16_4d
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: bg(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      REAL(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_e
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      REAL(DP):: reduce
      !
#if defined(__CUDA)
      CALL start_clock_gpu('build_hr')
#else
      CALL start_clock('build_hr')
#endif
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         IF(l2_e >= l2_s) THEN
            !
            !$acc enter data create(ag,c_distr(1:pert%nglob,l2_s:l2_e)) copyin(bg)
            !
            !$acc kernels present(c_distr(1:pert%nglob,l2_s:l2_e))
            c_distr(1:pert%nglob,l2_s:l2_e) = 0._DP
            !$acc end kernels
            !
         ENDIF
         !
         DO icycl = 0,nimage-1
            !
            IF(l2_e >= l2_s) THEN
               !
               idx = MOD(my_image_id+icycl,nimage)
               nloc = pert%nglob/nimage
               IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
               !
               !$acc update device(ag) wait
               !
               DO il1 = 1,nloc
                  !
                  ig1 = pert%l2g(il1,idx)
                  IF(ig1 < 1 .OR. ig1 > g_e) CYCLE
                  !
                  !$acc parallel async present(ag,bg,c_distr(1:pert%nglob,l2_s:l2_e),nbnd_occ,ngk)
                  !$acc loop
                  DO il2 = l2_s,l2_e
                     !
                     reduce = 0._DP
                     !
                     !$acc loop seq
                     DO iks = 1,nks
                        !
                        nbndval = nbnd_occ(iks)
                        npw = ngk(iks)
                        !
                        !$acc loop collapse(2) reduction(+:reduce)
                        DO ibnd = 1,nbndval
                           DO il3 = 1,npw
                              reduce = reduce+REAL(ag(il3,ibnd,iks,il1),KIND=DP)*REAL(bg(il3,ibnd,iks,il2),KIND=DP) &
                              & +AIMAG(ag(il3,ibnd,iks,il1))*AIMAG(bg(il3,ibnd,iks,il2))
                           ENDDO
                        ENDDO
                        !
                        reduce = reduce*2._DP
                        !
                        IF(gstart == 2) THEN
                           !$acc loop reduction(+:reduce)
                           DO ibnd = 1,nbndval
                              reduce = reduce-REAL(ag(1,ibnd,iks,il1),KIND=DP)*REAL(bg(1,ibnd,iks,il2),KIND=DP)
                           ENDDO
                        ENDIF
                        !
                     ENDDO
                     !
                     c_distr(ig1,il2) = reduce
                     !
                  ENDDO
                  !$acc end parallel
                  !
               ENDDO
               !
            ENDIF
            !
            ! Cycle ag
            !
            CALL mp_circular_shift_left_c16_4d(ag,icycl,inter_image_comm)
            !
         ENDDO
         !
         IF(l2_e >= l2_s) THEN
            !
            !$acc update host(c_distr(1:pert%nglob,l2_s:l2_e)) wait
            !$acc exit data delete(ag,bg,c_distr(1:pert%nglob,l2_s:l2_e))
            !
            CALL mp_sum(c_distr(:,l2_s:l2_e),intra_bgrp_comm)
            !
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(c_distr,0,inter_bgrp_comm)
      CALL mp_bcast(c_distr,0,inter_pool_comm)
      !
#if defined(__CUDA)
      CALL stop_clock_gpu('build_hr')
#else
      CALL stop_clock('build_hr')
#endif
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_real(ag,bg,nselect,n,lda,vr_distr,ew)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,my_bgrp_id,my_pool_id
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,nbndval0x
      USE west_mp,              ONLY : mp_circular_shift_left_c16_4d
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: bg(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      REAL(DP) :: dconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
#if defined(__CUDA)
      CALL start_clock_gpu('update_vr')
#else
      CALL start_clock('update_vr')
#endif
      !
      ! ag, bg only needed by pool 0 and band group 0 in the next step
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         ALLOCATE(hg(npwx,nbndval0x,nks,pert%nlocx))
         !
         !$acc enter data create(ag,bg,hg)
         !
         !$acc kernels present(hg)
         hg(:,:,:,:) = 0._DP
         !$acc end kernels
         !
         DO icycl = 0,nimage-1
            !
            idx = MOD(my_image_id+icycl,nimage)
            nloc = pert%nglob/nimage
            IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
            !
            !$acc update device(ag) wait
            !
            DO il1 = 1,nloc
               !
               ig1 = pert%l2g(il1,idx)
               IF(ig1 < 1 .OR. ig1 > n) CYCLE
               !
               DO il2 = 1,pert%nloc
                  !
                  ig2 = pert%l2g(il2)
                  IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  DO iks = 1,nks
                     !
                     nbndval = nbnd_occ(iks)
                     npw = ngk(iks)
                     !
                     !$acc parallel loop collapse(2) async present(hg,ag)
                     DO ibnd = 1,nbndval
                        DO il3 = 1,npw
                           hg(il3,ibnd,iks,il2) = dconst*ag(il3,ibnd,iks,il1)+hg(il3,ibnd,iks,il2)
                        ENDDO
                     ENDDO
                     !$acc end parallel
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            ! Cycle ag
            !
            CALL mp_circular_shift_left_c16_4d(ag,icycl,inter_image_comm)
            !
         ENDDO
         !
         !$acc update device(ag) wait
         !
         DO il2 = 1,pert%nloc
            !
            ig2 = pert%l2g(il2)
            IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
            !
            dconst = -ew(ig2)
            !
            DO iks = 1,nks
               !
               nbndval = nbnd_occ(iks)
               npw = ngk(iks)
               !
               !$acc parallel loop collapse(2) async present(ag,hg)
               DO ibnd = 1,nbndval
                  DO il3 = 1,npw
                     ag(:,ibnd,iks,il2) = dconst*hg(:,ibnd,iks,il2)
                  ENDDO
               ENDDO
               !$acc end parallel
               !
            ENDDO
            !
         ENDDO
         !
         !$acc wait
         !$acc kernels present(hg)
         hg(:,:,:,:) = 0._DP
         !$acc end kernels
         !
         DO icycl = 0,nimage-1
            !
            idx = MOD(my_image_id+icycl,nimage)
            nloc = pert%nglob/nimage
            IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
            !
            !$acc update device(bg) wait
            !
            DO il1 = 1,nloc
               !
               ig1 = pert%l2g(il1,idx)
               IF(ig1 < 1 .OR. ig1 > n) CYCLE
               !
               DO il2 = 1,pert%nloc
                  !
                  ig2 = pert%l2g(il2)
                  IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  DO iks = 1,nks
                     !
                     nbndval = nbnd_occ(iks)
                     npw = ngk(iks)
                     !
                     !$acc parallel loop collapse(2) async present(hg,bg)
                     DO ibnd = 1,nbndval
                        DO il3 = 1,npw
                           hg(il3,ibnd,iks,il2) = dconst*bg(il3,ibnd,iks,il1)+hg(il3,ibnd,iks,il2)
                        ENDDO
                     ENDDO
                     !$acc end parallel
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            ! Cycle bg
            !
            CALL mp_circular_shift_left_c16_4d(bg,icycl,inter_image_comm)
            !
         ENDDO
         !
         !$acc wait
         !
         DO il2 = 1,pert%nloc
            !
            ig2 = pert%l2g(il2)
            IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
            !
            DO iks = 1,nks
               !
               nbndval = nbnd_occ(iks)
               npw = ngk(iks)
               !
               !$acc parallel loop collapse(2) async present(ag,hg)
               DO ibnd = 1,nbndval
                  DO il3 = 1,npw
                     ag(il3,ibnd,iks,il2) = ag(il3,ibnd,iks,il2)+hg(il3,ibnd,iks,il2)
                  ENDDO
               ENDDO
               !$acc end parallel
               !
            ENDDO
            !
         ENDDO
         !
         !$acc update host(ag) wait
         !$acc exit data delete(ag,bg,hg)
         DEALLOCATE(hg)
         !
      ENDIF
      !
#if defined(__CUDA)
      CALL stop_clock_gpu('update_vr')
#else
      CALL stop_clock('update_vr')
#endif
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE refresh_with_vr_distr_real(ag,nselect,n,lda,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_bgrp_comm,&
                                     & my_bgrp_id,inter_pool_comm,my_pool_id
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,nbndval0x
      USE west_mp,              ONLY : mp_circular_shift_left_c16_4d
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      REAL(DP) :: dconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
#if defined(__CUDA)
      CALL start_clock_gpu('refresh_vr')
#else
      CALL start_clock('refresh_vr')
#endif
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         ALLOCATE(hg(npwx,nbndval0x,nks,pert%nlocx))
         !
         !$acc enter data create(ag,hg)
         !
         !$acc kernels present(hg)
         hg(:,:,:,:) = 0._DP
         !$acc end kernels
         !
         DO icycl = 0,nimage-1
            !
            idx = MOD(my_image_id+icycl,nimage)
            nloc = pert%nglob/nimage
            IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
            !
            !$acc update device(ag) wait
            !
            DO il1 = 1,nloc
               !
               ig1 = pert%l2g(il1,idx)
               IF(ig1 < 1 .OR. ig1 > n) CYCLE
               !
               DO il2 = 1,pert%nloc
                  !
                  ig2 = pert%l2g(il2)
                  IF(ig2 > nselect) CYCLE
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  DO iks = 1,nks
                     !
                     nbndval = nbnd_occ(iks)
                     npw = ngk(iks)
                     !
                     !$acc parallel loop collapse(2) async present(hg,ag)
                     DO ibnd = 1,nbndval
                        DO il3 = 1,npw
                           hg(il3,ibnd,iks,il2) = dconst*ag(il3,ibnd,iks,il1)+hg(il3,ibnd,iks,il2)
                        ENDDO
                     ENDDO
                     !$acc end parallel
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            ! Cycle ag
            !
            CALL mp_circular_shift_left_c16_4d(ag,icycl,inter_image_comm)
            !
         ENDDO
         !
         !$acc wait
         !
         DO il2 = 1,pert%nloc
            !
            ig2 = pert%l2g(il2)
            !
            IF(ig2 > nselect) THEN
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  !$acc parallel loop collapse(2) async present(ag)
                  DO ibnd = 1,nbndval
                     DO il3 = 1,npw
                        ag(il3,ibnd,iks,il2) = 0._DP
                     ENDDO
                  ENDDO
                  !$acc end parallel
                  !
               ENDDO
            ELSE
               DO iks  = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  !$acc parallel loop collapse(2) async present(ag,hg)
                  DO ibnd = 1,nbndval
                     DO il3 = 1,npw
                        ag(il3,ibnd,iks,il2) = hg(il3,ibnd,iks,il2)
                     ENDDO
                  ENDDO
                  !$acc end parallel
                  !
               ENDDO
            ENDIF
            !
         ENDDO
         !
         !$acc update host(ag) wait
         !$acc exit data delete(ag,hg)
         DEALLOCATE(hg)
         !
      ENDIF
      !
      CALL mp_bcast(ag,0,inter_bgrp_comm)
      CALL mp_bcast(ag,0,inter_pool_comm)
      !
#if defined(__CUDA)
      CALL stop_clock_gpu('refresh_vr')
#else
      CALL stop_clock('refresh_vr')
#endif
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE preconditioner_complex(ag,nselect,n,turn_shift)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : inter_image_comm
      USE mp,                   ONLY : mp_max
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx
      USE westcom,              ONLY : nbnd_occ,nbndval0x
      USE wvfct,                ONLY : g2kin,et
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n
      LOGICAL,INTENT(IN) :: turn_shift
      !
      ! Workspace
      !
      INTEGER :: il1,ig1,ig,ibnd,nbndval
      INTEGER :: iks,current_k
      INTEGER :: mloc,mstart,max_mloc
      REAL(DP):: temp,minimum
      !
      minimum = 0.001_DP
      !
      CALL start_clock('precd_ag')
      !
      mloc = 0
      mstart = 1
      DO il1 = 1,pert%nloc
         ig1 = pert%l2g(il1)
         IF(ig1 <= n .OR. ig1 > n+nselect) CYCLE
         IF(mloc == 0) mstart = il1
         mloc = mloc + 1
      ENDDO
      !
      ! Apply Liouville operator
      !
      max_mloc = mloc
      CALL mp_max(max_mloc,inter_image_comm)
      !
      DO il1 = mstart,mstart+max_mloc-1
         !
         ig1 = pert%l2g(il1)
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            current_k = iks
            !
            CALL g2_kin(iks)
            !
            IF(.NOT. (ig1 <= n .OR. ig1 > n+nselect)) THEN
               DO ibnd = 1,nbndval
                  DO ig = 1,npwx
                     IF(turn_shift) THEN
                        temp = (g2kin(ig) - et(ibnd,iks))
                     ELSE
                        temp = g2kin(ig)
                     ENDIF
                     !
                     IF(ABS(temp) < minimum) temp = SIGN(minimum,temp)
                     ag(ig,ibnd,iks,il1) = ag(ig,ibnd,iks,il1)/temp
                  ENDDO
               ENDDO
            ENDIF
            !
         ENDDO
         !
      ENDDO
      !
      CALL stop_clock('precd_ag')
      !
    END SUBROUTINE
    !
END MODULE
