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
  INTERFACE wbse_precondition_dvg
     MODULE PROCEDURE precondition_dvg_complex
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE build_hr_real(ag,bg,l2_s,l2_e,c_distr,g_e,sf)
      !------------------------------------------------------------------------
      !
      !  c_distr = < ag | bg >
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_pool_comm,&
                                     & inter_bgrp_comm,intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      USE pwcom,                ONLY : npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
      USE gvect,                ONLY : gstart
      USE west_mp,              ONLY : west_mp_circ_shift
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: bg(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      REAL(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_e
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,lbnd,ibnd,iks,iks_do,nbndval
      INTEGER :: l1_e
      INTEGER :: icycl,idx,nloc
      INTEGER :: pert_nglob,kpt_pool_nloc
      REAL(DP):: reduce
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      INTEGER,PARAMETER :: flks(2) = [2,1]
      !
#if defined(__CUDA)
      CALL start_clock_gpu('build_hr')
#else
      CALL start_clock('build_hr')
#endif
      !
      pert_nglob = pert%nglob
      kpt_pool_nloc = kpt_pool%nloc
      !
      ALLOCATE(nbnd_loc(kpt_pool%nloc))
      !
      DO iks = 1,kpt_pool%nloc
         !
         IF(sf) THEN
            iks_do = flks(iks)
         ELSE
            iks_do = iks
         ENDIF
         !
         nbndval = nbnd_occ(iks_do)
         !
         nbnd_loc(iks) = 0
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
         ENDDO
         !
      ENDDO
      !
      !$acc enter data copyin(nbnd_loc)
      !
      IF(l2_e >= l2_s) THEN
         !
         !$acc enter data create(c_distr(1:pert_nglob,l2_s:l2_e))
         !
         !$acc kernels present(c_distr(1:pert_nglob,l2_s:l2_e))
         c_distr(1:pert_nglob,l2_s:l2_e) = 0._DP
         !$acc end kernels
         !
      ENDIF
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
         !
         l1_e = 0
         DO il1 = nloc,1,-1
            ig1 = pert%l2g(il1,idx)
            IF(ig1 <= g_e) THEN
               l1_e = il1
               EXIT
            ENDIF
         ENDDO
         !
         IF(l1_e > 0 .AND. l2_e >= l2_s) THEN
            !
            !$acc parallel vector_length(1024) present(ag,bg,c_distr(1:pert_nglob,l2_s:l2_e),nbnd_loc,ngk)
            !$acc loop collapse(2)
            DO il1 = 1,l1_e
               DO il2 = l2_s,l2_e
                  !
                  ! ig1 = pert%l2g(il1,idx)
                  !
                  ig1 = nimage*(il1-1)+idx+1
                  !
                  reduce = 0._DP
                  !
                  !$acc loop seq
                  DO iks = 1,kpt_pool_nloc
                     !
                     nbndval = nbnd_loc(iks)
                     npw = ngk(iks)
                     !
                     !$acc loop collapse(2) reduction(+:reduce)
                     DO lbnd = 1,nbndval
                        DO il3 = 1,npw
                           reduce = reduce+2._DP*REAL(ag(il3,lbnd,iks,il1),KIND=DP)*REAL(bg(il3,lbnd,iks,il2),KIND=DP) &
                           & +2._DP*AIMAG(ag(il3,lbnd,iks,il1))*AIMAG(bg(il3,lbnd,iks,il2))
                        ENDDO
                     ENDDO
                     !
                     IF(gstart == 2) THEN
                        !$acc loop reduction(+:reduce)
                        DO lbnd = 1,nbndval
                           reduce = reduce-REAL(ag(1,lbnd,iks,il1),KIND=DP)*REAL(bg(1,lbnd,iks,il2),KIND=DP)
                        ENDDO
                     ENDIF
                     !
                  ENDDO
                  !
                  c_distr(ig1,il2) = reduce
                  !
               ENDDO
            ENDDO
            !$acc end parallel
            !
         ENDIF
         !
         ! Cycle ag
         !
         IF(nimage > 1) THEN
            !
            CALL west_mp_circ_shift(ag,icycl,inter_image_comm)
            !
            IF(l2_e >= l2_s) THEN
               !$acc update device(ag)
            ENDIF
            !
         ENDIF
         !
      ENDDO
      !
      IF(l2_e >= l2_s) THEN
         !
         !$acc update host(c_distr(1:pert_nglob,l2_s:l2_e))
         !$acc exit data delete(c_distr(1:pert_nglob,l2_s:l2_e))
         !
         CALL mp_sum(c_distr(:,l2_s:l2_e),intra_bgrp_comm)
         CALL mp_sum(c_distr(:,l2_s:l2_e),inter_bgrp_comm)
         CALL mp_sum(c_distr(:,l2_s:l2_e),inter_pool_comm)
         !
      ENDIF
      !
      !$acc exit data delete(nbnd_loc)
      DEALLOCATE(nbnd_loc)
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
    SUBROUTINE update_with_vr_distr_real(ag,bg,nselect,n,lda,vr_distr,ew,sf)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      USE pwcom,                ONLY : npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
      USE west_mp,              ONLY : west_mp_circ_shift
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: bg(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2,lbnd,ibnd,iks,iks_do,nbndval
      INTEGER :: l1_e,l2_s,l2_e
      INTEGER :: icycl,idx,nloc
      INTEGER :: kpt_pool_nloc
      REAL(DP) :: dconst
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      INTEGER,PARAMETER :: flks(2) = [2,1]
      !
#if defined(__CUDA)
      CALL start_clock_gpu('update_vr')
#else
      CALL start_clock('update_vr')
#endif
      !
      kpt_pool_nloc = kpt_pool%nloc
      !
      l2_s = 0
      DO il2 = 1,pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2 > n) THEN
            l2_s = il2
            EXIT
         ENDIF
      ENDDO
      !
      l2_e = 0
      DO il2 = pert%nloc,1,-1
         ig2 = pert%l2g(il2)
         IF(ig2 <= n+nselect) THEN
            l2_e = il2
            EXIT
         ENDIF
      ENDDO
      !
      ALLOCATE(nbnd_loc(kpt_pool%nloc))
      !
      DO iks = 1,kpt_pool%nloc
         !
         IF(sf) THEN
            iks_do = flks(iks)
         ELSE
            iks_do = iks
         ENDIF
         !
         nbndval = nbnd_occ(iks_do)
         !
         nbnd_loc(iks) = 0
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
         ENDDO
         !
      ENDDO
      !
      ALLOCATE(hg(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx))
      !
      !$acc enter data create(hg) copyin(vr_distr,ew,nbnd_loc)
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
         l1_e = 0
         DO il1 = nloc,1,-1
            ig1 = pert%l2g(il1,idx)
            IF(ig1 <= n) THEN
               l1_e = il1
               EXIT
            ENDIF
         ENDDO
         !
         IF(l1_e > 0 .AND. l2_s > 0 .AND. l2_e >= l2_s) THEN
            !
            !$acc parallel vector_length(1024) present(vr_distr,nbnd_loc,ngk,hg,ag)
            !$acc loop seq
            DO il1 = 1,l1_e
               !
               ! ig1 = pert%l2g(il1,idx)
               !
               ig1 = nimage*(il1-1)+idx+1
               !
               !$acc loop
               DO il2 = l2_s,l2_e
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc loop seq
                  DO iks = 1,kpt_pool_nloc
                     !
                     nbndval = nbnd_loc(iks)
                     npw = ngk(iks)
                     !
                     !$acc loop collapse(2)
                     DO lbnd = 1,nbndval
                        DO il3 = 1,npw
                           hg(il3,lbnd,iks,il2) = dconst*ag(il3,lbnd,iks,il1)+hg(il3,lbnd,iks,il2)
                        ENDDO
                     ENDDO
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !$acc end parallel
            !
         ENDIF
         !
         ! Cycle ag
         !
         IF(nimage > 1) THEN
            !
            CALL west_mp_circ_shift(ag,icycl,inter_image_comm)
            !
            !$acc update device(ag)
            !
         ENDIF
         !
      ENDDO
      !
      IF(l2_s > 0 .AND. l2_e >= l2_s) THEN
         !
         !$acc parallel vector_length(1024) present(ew,nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = l2_s,l2_e
            !
            ! ig2 = pert%l2g(il2)
            !
            ig2 = nimage*(il2-1)+my_image_id+1
            !
            dconst = -ew(ig2)
            !
            !$acc loop seq
            DO iks = 1,kpt_pool_nloc
               !
               nbndval = nbnd_loc(iks)
               npw = ngk(iks)
               !
               !$acc loop collapse(2)
               DO lbnd = 1,nbndval
                  DO il3 = 1,npw
                     ag(il3,lbnd,iks,il2) = dconst*hg(il3,lbnd,iks,il2)
                  ENDDO
               ENDDO
               !
            ENDDO
            !
         ENDDO
         !$acc end parallel
         !
      ENDIF
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
         l1_e = 0
         DO il1 = nloc,1,-1
            ig1 = pert%l2g(il1,idx)
            IF(ig1 <= n) THEN
               l1_e = il1
               EXIT
            ENDIF
         ENDDO
         !
         IF(l1_e > 0 .AND. l2_s > 0 .AND. l2_e >= l2_s) THEN
            !
            !$acc parallel vector_length(1024) present(vr_distr,nbnd_loc,ngk,hg,bg)
            !$acc loop seq
            DO il1 = 1,l1_e
               !
               ! ig1 = pert%l2g(il1,idx)
               !
               ig1 = nimage*(il1-1)+idx+1
               !
               !$acc loop
               DO il2 = l2_s,l2_e
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc loop seq
                  DO iks = 1,kpt_pool_nloc
                     !
                     nbndval = nbnd_loc(iks)
                     npw = ngk(iks)
                     !
                     !$acc loop collapse(2)
                     DO lbnd = 1,nbndval
                        DO il3 = 1,npw
                           hg(il3,lbnd,iks,il2) = dconst*bg(il3,lbnd,iks,il1)+hg(il3,lbnd,iks,il2)
                        ENDDO
                     ENDDO
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !$acc end parallel
            !
         ENDIF
         !
         ! Cycle bg
         !
         IF(nimage > 1) THEN
            !
            CALL west_mp_circ_shift(bg,icycl,inter_image_comm)
            !
            !$acc update device(bg)
            !
         ENDIF
         !
      ENDDO
      !
      IF(l2_s > 0 .AND. l2_e >= l2_s) THEN
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = l2_s,l2_e
            !$acc loop seq
            DO iks = 1,kpt_pool_nloc
               !
               nbndval = nbnd_loc(iks)
               npw = ngk(iks)
               !
               !$acc loop collapse(2)
               DO lbnd = 1,nbndval
                  DO il3 = 1,npw
                     ag(il3,lbnd,iks,il2) = ag(il3,lbnd,iks,il2)+hg(il3,lbnd,iks,il2)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !$acc end parallel
         !
      ENDIF
      !
      !$acc exit data delete(hg,vr_distr,ew,nbnd_loc)
      DEALLOCATE(hg)
      DEALLOCATE(nbnd_loc)
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
    SUBROUTINE refresh_with_vr_distr_real(ag,nselect,n,lda,vr_distr,sf)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      USE pwcom,                ONLY : npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
      USE west_mp,              ONLY : west_mp_circ_shift
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2,lbnd,ibnd,iks,iks_do,nbndval
      INTEGER :: l1_e,l2_s,l2_e
      INTEGER :: icycl,idx,nloc
      INTEGER :: pert_nloc,kpt_pool_nloc
      REAL(DP) :: dconst
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      INTEGER,PARAMETER :: flks(2) = [2,1]
      !
#if defined(__CUDA)
      CALL start_clock_gpu('refresh_vr')
#else
      CALL start_clock('refresh_vr')
#endif
      !
      pert_nloc = pert%nloc
      kpt_pool_nloc = kpt_pool%nloc
      !
      l2_s = 1
      l2_e = 0
      DO il2 = pert%nloc,1,-1
         ig2 = pert%l2g(il2)
         IF(ig2 <= nselect) THEN
            l2_e = il2
            EXIT
         ENDIF
      ENDDO
      !
      ALLOCATE(nbnd_loc(kpt_pool%nloc))
      !
      DO iks = 1,kpt_pool%nloc
         !
         IF(sf) THEN
            iks_do = flks(iks)
         ELSE
            iks_do = iks
         ENDIF
         !
         nbndval = nbnd_occ(iks_do)
         !
         nbnd_loc(iks) = 0
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
         ENDDO
         !
      ENDDO
      !
      ALLOCATE(hg(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx))
      !
      !$acc enter data create(hg) copyin(vr_distr,nbnd_loc)
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
         l1_e = 0
         DO il1 = nloc,1,-1
            ig1 = pert%l2g(il1,idx)
            IF(ig1 <= n) THEN
               l1_e = il1
               EXIT
            ENDIF
         ENDDO
         !
         IF(l1_e > 0 .AND. l2_e >= l2_s) THEN
            !
            !$acc parallel vector_length(1024) present(vr_distr,nbnd_loc,ngk,hg,ag)
            !$acc loop seq
            DO il1 = 1,l1_e
               !
               ! ig1 = pert%l2g(il1,idx)
               !
               ig1 = nimage*(il1-1)+idx+1
               !
               !$acc loop
               DO il2 = l2_s,l2_e
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc loop seq
                  DO iks = 1,kpt_pool_nloc
                     !
                     nbndval = nbnd_loc(iks)
                     npw = ngk(iks)
                     !
                     !$acc loop collapse(2)
                     DO lbnd = 1,nbndval
                        DO il3 = 1,npw
                           hg(il3,lbnd,iks,il2) = dconst*ag(il3,lbnd,iks,il1)+hg(il3,lbnd,iks,il2)
                        ENDDO
                     ENDDO
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !$acc end parallel
            !
         ENDIF
         !
         ! Cycle ag
         !
         IF(nimage > 1) THEN
            !
            CALL west_mp_circ_shift(ag,icycl,inter_image_comm)
            !
            !$acc update device(ag)
            !
         ENDIF
         !
      ENDDO
      !
      l2_s = 1
      l2_e = 0
      DO il2 = pert%nloc,1,-1
         ig2 = pert%l2g(il2)
         IF(ig2 <= nselect) THEN
            l2_e = il2
            EXIT
         ENDIF
      ENDDO
      !
      IF(l2_e > 0) THEN
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = l2_s,l2_e
            !$acc loop seq
            DO iks = 1,kpt_pool_nloc
               !
               nbndval = nbnd_loc(iks)
               npw = ngk(iks)
               !
               !$acc loop collapse(2)
               DO lbnd = 1,nbndval
                  DO il3 = 1,npw
                     ag(il3,lbnd,iks,il2) = hg(il3,lbnd,iks,il2)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !$acc end parallel
         !
      ENDIF
      !
      l2_s = 0
      l2_e = pert_nloc
      DO il2 = 1,pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2 > nselect) THEN
            l2_s = il2
            EXIT
         ENDIF
      ENDDO
      !
      IF(l2_s > 0) THEN
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = l2_s,l2_e
            !$acc loop seq
            DO iks = 1,kpt_pool_nloc
               !
               nbndval = nbnd_loc(iks)
               npw = ngk(iks)
               !
               !$acc loop collapse(2)
               DO lbnd = 1,nbndval
                  DO il3 = 1,npw
                     ag(il3,lbnd,iks,il2) = 0._DP
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !$acc end parallel
         !
      ENDIF
      !
      !$acc exit data delete(hg,vr_distr,nbnd_loc)
      DEALLOCATE(hg)
      DEALLOCATE(nbnd_loc)
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
    SUBROUTINE precondition_dvg_complex(ag,nselect,n,turn_shift,sf)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE distribution_center,  ONLY : pert,kpt_pool,band_group
      USE pwcom,                ONLY : npwx
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
      USE wvfct,                ONLY : g2kin,et
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n
      LOGICAL,INTENT(IN) :: turn_shift
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,ig1,ig,lbnd,ibnd,iks,iks_do,nbndval
      INTEGER :: l1_s,l1_e
      INTEGER :: kpt_pool_nloc,band_group_myoffset
      REAL(DP):: tmp,tmp_abs,tmp_sgn
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      REAL(DP),ALLOCATABLE :: g2kin_save(:,:)
      INTEGER,PARAMETER :: flks(2) = [2,1]
      REAL(DP),PARAMETER :: minimum = 1._DP
      !
#if defined(__CUDA)
      CALL start_clock_gpu('precd_ag')
#else
      CALL start_clock('precd_ag')
#endif
      !
      kpt_pool_nloc = kpt_pool%nloc
      band_group_myoffset = band_group%myoffset
      !
      l1_s = 0
      DO il1 = 1,pert%nloc
         ig1 = pert%l2g(il1)
         IF(ig1 > n) THEN
            l1_s = il1
            EXIT
         ENDIF
      ENDDO
      !
      l1_e = 0
      DO il1 = pert%nloc,1,-1
         ig1 = pert%l2g(il1)
         IF(ig1 <= n+nselect) THEN
            l1_e = il1
            EXIT
         ENDIF
      ENDDO
      !
      IF(l1_s > 0 .AND. l1_e >= l1_s) THEN
         !
         ALLOCATE(g2kin_save(npwx,kpt_pool%nloc))
         !$acc enter data create(g2kin_save)
         ALLOCATE(nbnd_loc(kpt_pool%nloc))
         !
         DO iks = 1,kpt_pool%nloc
            !
            CALL g2_kin(iks)
            !
            !$acc kernels present(g2kin_save,g2kin)
            g2kin_save(:,iks) = g2kin
            !$acc end kernels
            !
            IF(sf) THEN
               iks_do = flks(iks)
            ELSE
               iks_do = iks
            ENDIF
            !
            nbndval = nbnd_occ(iks_do)
            !
            nbnd_loc(iks) = 0
            DO lbnd = 1,band_group%nloc
               ibnd = band_group%l2g(lbnd)+n_trunc_bands
               IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
            ENDDO
            !
         ENDDO
         !
         !$acc enter data copyin(nbnd_loc)
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,g2kin_save,et,ag)
         !$acc loop
         DO il1 = l1_s,l1_e
            !$acc loop seq
            DO iks = 1,kpt_pool_nloc
               !
               nbndval = nbnd_loc(iks)
               !
               !$acc loop collapse(2)
               DO lbnd = 1,nbndval
                  DO ig = 1,npwx
                     !
                     ! ibnd = band_group%l2g(lbnd)
                     !
                     ibnd = band_group_myoffset+lbnd
                     !
                     IF(turn_shift) THEN
                        tmp = g2kin_save(ig,iks)-et(ibnd+n_trunc_bands,iks_do)
                     ELSE
                        tmp = g2kin_save(ig,iks)
                     ENDIF
                     !
                     ! Same as the following line but without thread divergence
                     ! IF(ABS(tmp) < minimum) tmp = SIGN(minimum,tmp)
                     !
                     tmp_abs = MAX(ABS(tmp),minimum)
                     tmp_sgn = SIGN(1._DP,tmp)
                     tmp = tmp_sgn*tmp_abs
                     !
                     ag(ig,lbnd,iks,il1) = ag(ig,lbnd,iks,il1)/tmp
                     !
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !$acc end parallel
         !
         !$acc exit data delete(g2kin_save,nbnd_loc)
         DEALLOCATE(g2kin_save)
         DEALLOCATE(nbnd_loc)
         !
      ENDIF
      !
#if defined(__CUDA)
      CALL stop_clock_gpu('precd_ag')
#else
      CALL stop_clock('precd_ag')
#endif
      !
    END SUBROUTINE
    !
END MODULE
