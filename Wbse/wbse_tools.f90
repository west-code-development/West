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
                                     & my_pool_id,inter_bgrp_comm,intra_bgrp_comm
      USE mp,                   ONLY : mp_sum,mp_bcast
      USE distribution_center,  ONLY : pert,aband
      USE pwcom,                ONLY : nks,npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
      USE gvect,                ONLY : gstart
      USE west_mp,              ONLY : mp_circular_shift_left_c16_4d
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,aband%nlocx,nks,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: bg(npwx,aband%nlocx,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      REAL(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_e
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,lbnd,ibnd,iks,nbndval,iks_do
      INTEGER, DIMENSION(2), PARAMETER :: flks = (/ 2, 1 /)
      INTEGER :: il1_end
      INTEGER :: icycl,idx,nloc
      INTEGER :: pert_nglob
      REAL(DP):: reduce
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      !
#if defined(__CUDA)
      CALL start_clock_gpu('build_hr')
#else
      CALL start_clock('build_hr')
#endif
      !
      pert_nglob = pert%nglob
      !
      IF(my_pool_id == 0) THEN
         !
         ALLOCATE(nbnd_loc(nks))
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            !
            nbnd_loc(iks) = 0
            DO lbnd = 1,aband%nloc
               ibnd = aband%l2g(lbnd)+n_trunc_bands
               IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
            ENDDO
            !
         ENDDO
         !
         !$acc enter data copyin(nbnd_loc)
         !
         IF(l2_e >= l2_s) THEN
            !
            !$acc enter data create(c_distr(1:pert_nglob,l2_s:l2_e)) copyin(ag,bg)
            !
            !$acc kernels present(c_distr(1:pert_nglob,l2_s:l2_e))
            c_distr(1:pert_nglob,l2_s:l2_e) = 0._DP
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
               DO il1 = 1,nloc
                  ig1 = pert%l2g(il1,idx)
                  IF(ig1 > g_e) EXIT
                  il1_end = il1
               ENDDO
               !
               !$acc parallel vector_length(1024) present(ag,bg,c_distr(1:pert_nglob,l2_s:l2_e),nbnd_loc,ngk)
               !$acc loop collapse(2)
               DO il1 = 1,il1_end
                  DO il2 = l2_s,l2_e
                     !
                     ! ig1 = pert%l2g(il1,idx)
                     !
                     ig1 = nimage*(il1-1)+idx+1
                     !
                     reduce = 0._DP
                     !
                     !$acc loop seq
                     DO iks = 1,nks
                        !
                        IF(sf) THEN
                           iks_do = flks(iks)
                        ELSE
                           iks_do = iks
                        ENDIF
                        !
                        nbndval = nbnd_loc(iks_do)
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
                        !reduce = reduce*2._DP
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
            CALL mp_circular_shift_left_c16_4d(ag,icycl,inter_image_comm)
            !
            IF(l2_e >= l2_s) THEN
               !$acc update device(ag)
            ENDIF
            !
         ENDDO
         !
         IF(l2_e >= l2_s) THEN
            !
            !$acc update host(c_distr(1:pert_nglob,l2_s:l2_e))
            !$acc exit data delete(ag,bg,c_distr(1:pert_nglob,l2_s:l2_e))
            !
            CALL mp_sum(c_distr(:,l2_s:l2_e),intra_bgrp_comm)
            CALL mp_sum(c_distr(:,l2_s:l2_e),inter_bgrp_comm)
            !
         ENDIF
         !
         !$acc exit data delete(nbnd_loc)
         DEALLOCATE(nbnd_loc)
         !
      ENDIF
      !
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
    SUBROUTINE update_with_vr_distr_real(ag,bg,nselect,n,lda,vr_distr,ew,sf)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,my_pool_id
      USE distribution_center,  ONLY : pert,aband
      USE pwcom,                ONLY : nks,npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
      USE west_mp,              ONLY : mp_circular_shift_left_c16_4d
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,aband%nlocx,nks,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: bg(npwx,aband%nlocx,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2,lbnd,ibnd,iks,nbndval,iks_do
      INTEGER, DIMENSION(2), PARAMETER :: flks = (/ 2, 1 /)
      INTEGER :: il1_end,il2_start,il2_end
      INTEGER :: icycl,idx,nloc
      REAL(DP) :: dconst
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
#if defined(__CUDA)
      CALL start_clock_gpu('update_vr')
#else
      CALL start_clock('update_vr')
#endif
      !
      ! ag, bg only needed by pool 0 in the next step
      !
      IF(my_pool_id == 0) THEN
         !
         il2_start = 1
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 <= n) CYCLE
            il2_start = il2
            EXIT
         ENDDO
         !
         il2_end = 1
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 > n+nselect) EXIT
            il2_end = il2
         ENDDO
         !
         ALLOCATE(nbnd_loc(nks))
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            !
            nbnd_loc(iks) = 0
            DO lbnd = 1,aband%nloc
               ibnd = aband%l2g(lbnd)+n_trunc_bands
               IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
            ENDDO
            !
         ENDDO
         !
         ALLOCATE(hg(npwx,aband%nlocx,nks,pert%nlocx))
         !
         !$acc enter data create(hg) copyin(ag,bg,vr_distr,ew,nbnd_loc)
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
            DO il1 = 1,nloc
               ig1 = pert%l2g(il1,idx)
               IF(ig1 > n) EXIT
               il1_end = il1
            ENDDO
            !
            !$acc parallel vector_length(1024) present(vr_distr,nbnd_loc,ngk,hg,ag)
            !$acc loop seq
            DO il1 = 1,il1_end
               !
               ! ig1 = pert%l2g(il1,idx)
               !
               ig1 = nimage*(il1-1)+idx+1
               !
               !$acc loop
               DO il2 = il2_start,il2_end
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc loop seq
                  DO iks = 1,nks
                     !
                     IF(sf) THEN
                        iks_do = flks(iks)
                     ELSE
                        iks_do = iks
                     ENDIF
                     !
                     nbndval = nbnd_loc(iks_do)
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
            ! Cycle ag
            !
            CALL mp_circular_shift_left_c16_4d(ag,icycl,inter_image_comm)
            !
            !$acc update device(ag)
            !
         ENDDO
         !
         !$acc parallel vector_length(1024) present(ew,nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = il2_start,il2_end
            !
            ! ig2 = pert%l2g(il2)
            !
            ig2 = nimage*(il2-1)+my_image_id+1
            !
            dconst = -ew(ig2)
            !
            !$acc loop seq
            DO iks = 1,nks
               !
               IF(sf) THEN
                  iks_do = flks(iks)
               ELSE
                  iks_do = iks
               ENDIF
               !
               nbndval = nbnd_loc(iks_do)
               npw = ngk(iks)
               !
               !$acc loop collapse(2)
               DO lbnd = 1,nbndval
                  DO il3 = 1,npw
                     ag(:,lbnd,iks,il2) = dconst*hg(:,lbnd,iks,il2)
                  ENDDO
               ENDDO
               !
            ENDDO
            !
         ENDDO
         !$acc end parallel
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
            DO il1 = 1,nloc
               ig1 = pert%l2g(il1,idx)
               IF(ig1 > n) EXIT
               il1_end = il1
            ENDDO
            !
            !$acc parallel vector_length(1024) present(vr_distr,nbnd_loc,ngk,hg,bg)
            !$acc loop seq
            DO il1 = 1,il1_end
               !$acc loop
               DO il2 = il2_start,il2_end
                  !
                  ! ig1 = pert%l2g(il1,idx)
                  !
                  ig1 = nimage*(il1-1)+idx+1
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc loop seq
                  DO iks = 1,nks
                     !
                     IF(sf) THEN
                        iks_do = flks(iks)
                     ELSE
                        iks_do = iks
                     ENDIF
                     !
                     nbndval = nbnd_loc(iks_do)
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
            ENDDO
            !$acc end parallel
            !
            ! Cycle bg
            !
            CALL mp_circular_shift_left_c16_4d(bg,icycl,inter_image_comm)
            !
            !$acc update device(bg)
            !
         ENDDO
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = il2_start,il2_end
            !$acc loop seq
            DO iks = 1,nks
               !
               IF(sf) THEN
                  iks_do = flks(iks)
               ELSE
                  iks_do = iks
               ENDIF
               !
               nbndval = nbnd_loc(iks_do)
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
         !$acc exit data delete(bg,hg,vr_distr,ew,nbnd_loc) copyout(ag)
         DEALLOCATE(hg)
         DEALLOCATE(nbnd_loc)
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
    SUBROUTINE refresh_with_vr_distr_real(ag,nselect,n,lda,vr_distr,sf)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_pool_comm,&
                                     & my_pool_id
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : pert,aband
      USE pwcom,                ONLY : nks,npwx,npw,ngk
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
      USE west_mp,              ONLY : mp_circular_shift_left_c16_4d
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,aband%nlocx,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2,lbnd,ibnd,iks,nbndval,iks_do
      INTEGER, DIMENSION(2), PARAMETER :: flks = (/ 2, 1 /)
      INTEGER :: il1_end,il2_start,il2_end
      INTEGER :: icycl,idx,nloc
      INTEGER :: pert_nloc
      REAL(DP) :: dconst
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
#if defined(__CUDA)
      CALL start_clock_gpu('refresh_vr')
#else
      CALL start_clock('refresh_vr')
#endif
      !
      pert_nloc = pert%nloc
      !
      IF(my_pool_id == 0) THEN
         !
         il2_start = 1
         il2_end = 1
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 > nselect) EXIT
            il2_end = il2
         ENDDO
         !
         ALLOCATE(nbnd_loc(nks))
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            !
            nbnd_loc(iks) = 0
            DO lbnd = 1,aband%nloc
               ibnd = aband%l2g(lbnd)+n_trunc_bands
               IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
            ENDDO
            !
         ENDDO
         !
         ALLOCATE(hg(npwx,aband%nlocx,nks,pert%nlocx))
         !
         !$acc enter data create(hg) copyin(ag,vr_distr,nbnd_loc)
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
            DO il1 = 1,nloc
               ig1 = pert%l2g(il1,idx)
               IF(ig1 > n) EXIT
               il1_end = il1
            ENDDO
            !
            !$acc parallel vector_length(1024) present(vr_distr,nbnd_loc,ngk,hg,ag)
            !$acc loop seq
            DO il1 = 1,il1_end
               !
               ! ig1 = pert%l2g(il1,idx)
               !
               ig1 = nimage*(il1-1)+idx+1
               !
               !$acc loop
               DO il2 = il2_start,il2_end
                  !
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc loop seq
                  DO iks = 1,nks
                     !
                     IF(sf) THEN
                        iks_do = flks(iks)
                     ELSE
                        iks_do = iks
                     ENDIF
                     !
                     nbndval = nbnd_loc(iks_do)
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
            ! Cycle ag
            !
            CALL mp_circular_shift_left_c16_4d(ag,icycl,inter_image_comm)
            !
            !$acc update device(ag)
            !
         ENDDO
         !
         il2_start = 1
         il2_end = 1
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 > nselect) EXIT
            il2_end = il2
         ENDDO
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = il2_start,il2_end
            !$acc loop seq
            DO iks = 1,nks
               !
               IF(sf) THEN
                  iks_do = flks(iks)
               ELSE
                  iks_do = iks
               ENDIF
               !
               nbndval = nbnd_loc(iks_do)
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
         il2_start = 1
         il2_end = pert_nloc
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 <= nselect) CYCLE
            il2_start = il2
            EXIT
         ENDDO
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,ngk,ag,hg)
         !$acc loop
         DO il2 = il2_start,il2_end
            !$acc loop seq
            DO iks = 1,nks
               !
               IF(sf) THEN
                  iks_do = flks(iks)
               ELSE
                  iks_do = iks
               ENDIF
               !
               nbndval = nbnd_loc(iks_do)
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
         !$acc exit data delete(hg,vr_distr,nbnd_loc) copyout(ag)
         DEALLOCATE(hg)
         DEALLOCATE(nbnd_loc)
         !
      ENDIF
      !
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
    SUBROUTINE precondition_dvg_complex(ag,nselect,n,turn_shift,sf)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : inter_pool_comm,my_pool_id,nbgrp,my_bgrp_id
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : pert,aband
      USE pwcom,                ONLY : nks,npwx
      USE westcom,              ONLY : nbnd_occ,n_trunc_bands
#if defined(__CUDA)
      USE wvfct_gpum,           ONLY : g2kin=>g2kin_d,et=>et_d
#else
      USE wvfct,                ONLY : g2kin,et
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,aband%nlocx,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n
      LOGICAL,INTENT(IN) :: turn_shift
      LOGICAL,INTENT(IN) :: sf
      !
      ! Workspace
      !
      INTEGER :: il1,ig1,ig,lbnd,ibnd,nbndval,iks_do
      INTEGER, DIMENSION(2), PARAMETER :: flks = (/ 2, 1 /)
      INTEGER :: iks
      INTEGER :: il1_start,il1_end
      REAL(DP):: tmp,tmp_abs,tmp_sgn
      INTEGER,ALLOCATABLE :: nbnd_loc(:)
      REAL(DP),ALLOCATABLE :: g2kin_save(:,:)
      !REAL(DP),PARAMETER :: minimum = 0.01_DP
      REAL(DP),PARAMETER :: minimum = 1._DP
      !
#if defined(__CUDA)
      CALL start_clock_gpu('precd_ag')
#else
      CALL start_clock('precd_ag')
#endif
      !
      IF(my_pool_id == 0) THEN
         !
         ALLOCATE(g2kin_save(npwx,nks))
         !$acc enter data create(g2kin_save)
         ALLOCATE(nbnd_loc(nks))
         !
         DO iks = 1,nks
            !
#if defined(__CUDA)
            CALL g2_kin_gpu(iks)
#else
            CALL g2_kin(iks)
#endif
            !
            !$acc kernels present(g2kin_save)
            g2kin_save(:,iks) = g2kin
            !$acc end kernels
            !
            nbndval = nbnd_occ(iks)
            !
            nbnd_loc(iks) = 0
            DO lbnd = 1,aband%nloc
               ibnd = aband%l2g(lbnd)+n_trunc_bands
               IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_loc(iks) = nbnd_loc(iks)+1
            ENDDO
            !
         ENDDO
         !
         !$acc enter data copyin(ag,nbnd_loc)
         !
         il1_start = 1
         DO il1 = 1,pert%nloc
            ig1 = pert%l2g(il1)
            IF(ig1 <= n) CYCLE
            il1_start = il1
            EXIT
         ENDDO
         !
         il1_end = 1
         DO il1 = 1,pert%nloc
            ig1 = pert%l2g(il1)
            IF(ig1 > n+nselect) EXIT
            il1_end = il1
         ENDDO
         !
         !$acc parallel vector_length(1024) present(nbnd_loc,g2kin_save,ag)
         !$acc loop
         DO il1 = il1_start,il1_end
            !$acc loop seq
            DO iks = 1,nks
               !
               IF(sf) THEN
                  iks_do = flks(iks)
               ELSE
                  iks_do = iks
               ENDIF
               !
               nbndval = nbnd_loc(iks_do)
               !
               !$acc loop collapse(2)
               DO lbnd = 1,nbndval
                  DO ig = 1,npwx
                     !
                     ! ibnd = aband%l2g(lbnd)
                     !
                     ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
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
         !$acc exit data delete(g2kin_save,nbnd_loc) copyout(ag)
         DEALLOCATE(g2kin_save)
         DEALLOCATE(nbnd_loc)
         !
      ENDIF
      !
      CALL mp_bcast(ag,0,inter_pool_comm)
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
