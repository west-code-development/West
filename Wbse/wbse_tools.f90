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
! Marco Govoni
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
     MODULE PROCEDURE build_hr_real, build_hr_complex
  END INTERFACE
  !
  INTERFACE wbse_update_with_vr_distr
     MODULE PROCEDURE update_with_vr_distr_real, update_with_vr_distr_complex
  END INTERFACE
  !
  INTERFACE wbse_refresh_with_vr_distr
     MODULE PROCEDURE refresh_with_vr_distr_real, refresh_with_vr_distr_complex
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
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
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
      INTEGER :: il1,il2,ig1,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      REAL(DP):: reduce
      REAL(DP),EXTERNAL :: DDOT
      !
      CALL start_clock('build_hr')
      !
      ! Initialize to zero
      !
      c_distr(:,l2_s:l2_e) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
         !
         DO il1 = 1,nloc
            !
            ig1 = pert%l2g(il1,idx)
            IF(ig1 < 1 .OR. ig1 > g_e) CYCLE
            !
            DO il2 = l2_s,l2_e
               !
               reduce = 0._DP
               !
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     reduce = reduce + 2._DP * DDOT(2*npw,ag(1,ibnd,iks,il1),1,bg(1,ibnd,iks,il2),1)
                     IF(gstart == 2) reduce = reduce - REAL(ag(1,ibnd,iks,il1),KIND=DP)*REAL(bg(1,ibnd,iks,il2),KIND=DP)
                  ENDDO
                  !
               ENDDO
               !
               c_distr(ig1,il2) = reduce
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
      CALL mp_sum(c_distr(:,l2_s:l2_e),intra_bgrp_comm)
      !
      CALL stop_clock('build_hr')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE build_hr_complex(ag,bg,l2_s,l2_e,c_distr,g_e)
      !------------------------------------------------------------------------
      !
      !  c_distr = < ag | bg >
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
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
      COMPLEX(DP),INTENT(IN) :: bg(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      COMPLEX(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_e
      !
      ! Workspace
      !
      INTEGER :: il1,il2,ig1,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      COMPLEX(DP),EXTERNAL :: ZDOTC
      !
      CALL start_clock('build_hr')
      !
      ! Initialize to zero
      !
      c_distr(:,l2_s:l2_e) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
         !
         DO il1 = 1,nloc
            !
            ig1 = pert%l2g(il1,idx)
            IF(ig1 < 1 .OR. ig1 > g_e) CYCLE
            !
            DO il2 = l2_s,l2_e
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     c_distr(ig1,il2) = c_distr(ig1,il2) + ZDOTC(npw,ag(1,ibnd,iks,il1),1,bg(1,ibnd,iks,il2),1)
                  ENDDO
                  !
               ENDDO
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
      CALL mp_sum(c_distr(:,l2_s:l2_e),intra_bgrp_comm)
      !
      CALL stop_clock('build_hr')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_real(ag,bg,nselect,n,lda,vr_distr,ew)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id
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
      INTEGER :: il1,il2,ig1,ig2,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      COMPLEX(DP) :: zconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
      CALL start_clock('update_vr')
      !
      ALLOCATE(hg(npwx,nbndval0x,nks,pert%nlocx))
      hg(:,:,:,:) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
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
               zconst = CMPLX(vr_distr(ig1,il2),KIND=DP)
               !
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     CALL ZAXPY(npw,zconst,ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
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
      DO il2 = 1,pert%nloc
         !
         ig2 = pert%l2g(il2)
         IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd = 1,nbndval
               ag(:,ibnd,iks,il2) = -ew(ig2) * hg(:,ibnd,iks,il2)
            ENDDO
            !
         ENDDO
         !
      ENDDO
      !
      hg(:,:,:,:) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
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
               zconst = CMPLX(vr_distr(ig1,il2),KIND=DP)
               !
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     CALL ZAXPY(npw,zconst,bg(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
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
      DO il2 = 1,pert%nloc
         !
         ig2 = pert%l2g(il2)
         IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd = 1,nbndval
               ag(:,ibnd,iks,il2) = ag(:,ibnd,iks,il2) + hg(:,ibnd,iks,il2)
            ENDDO
            !
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE(hg)
      !
      CALL stop_clock('update_vr')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_complex(ag,bg,nselect,n,lda,vr_distr,ew)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id
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
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,ig1,ig2,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
      CALL start_clock('update_vr')
      !
      ALLOCATE(hg(npwx,nbndval0x,nks,pert%nlocx))
      hg(:,:,:,:) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
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
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
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
      DO il2 = 1,pert%nloc
         !
         ig2 = pert%l2g(il2)
         IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd = 1,nbndval
               ag(:,ibnd,iks,il2) = -ew(ig2) * hg(:,ibnd,iks,il2)
            ENDDO
            !
         ENDDO
         !
      ENDDO
      !
      hg(:,:,:,:) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
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
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),bg(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
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
      DO il2 = 1,pert%nloc
         !
         ig2 = pert%l2g(il2)
         IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
         !
         DO iks = 1,nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd = 1,nbndval
               ag(:,ibnd,iks,il2) = ag(:,ibnd,iks,il2) + hg(:,ibnd,iks,il2)
            ENDDO
            !
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE(hg)
      !
      CALL stop_clock('update_vr')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE refresh_with_vr_distr_real(ag,nselect,n,lda,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id
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
      INTEGER :: il1,il2,ig1,ig2,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      COMPLEX(DP) :: zconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
      CALL start_clock('refresh_vr')
      !
      ALLOCATE(hg(npwx,nbndval0x,nks,pert%nlocx))
      hg(:,:,:,:) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
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
               zconst = CMPLX(vr_distr(ig1,il2),KIND=DP)
               !
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     CALL ZAXPY(npw,zconst,ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
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
      DO il2 = 1,pert%nloc
         !
         ig2 = pert%l2g(il2)
         !
         IF(ig2 > nselect) THEN
            DO iks = 1,nks
               !
               nbndval = nbnd_occ(iks)
               !
               DO ibnd = 1,nbndval
                  ag(:,ibnd,iks,il2) = 0._DP
               ENDDO
               !
            ENDDO
         ELSE
            DO iks  = 1,nks
               !
               nbndval = nbnd_occ(iks)
               !
               DO ibnd = 1,nbndval
                  ag(:,ibnd,iks,il2) = hg(:,ibnd,iks,il2)
               ENDDO
               !
            ENDDO
         ENDIF
         !
      ENDDO
      !
      DEALLOCATE(hg)
      !
      CALL stop_clock('refresh_vr')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE refresh_with_vr_distr_complex(ag,nselect,n,lda,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id
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
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,ig1,ig2,ibnd,iks,nbndval
      INTEGER :: icycl,idx,nloc
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      !
      CALL start_clock('refresh_vr')
      !
      ALLOCATE(hg(npwx,nbndval0x,nks,pert%nlocx))
      hg(:,:,:,:) = 0._DP
      !
      DO icycl = 0,nimage-1
         !
         idx = MOD(my_image_id+icycl,nimage)
         nloc = pert%nglob/nimage
         IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
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
               DO iks = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  npw = ngk(iks)
                  !
                  DO ibnd = 1,nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle ag
         !
         CALL mp_circular_shift_left_c16_4d(ag,icycl,inter_image_comm)
         !
      ENDDO
      !
      DO il2 = 1,pert%nloc
         !
         ig2 = pert%l2g(il2)
         !
         IF(ig2 > nselect) THEN
           DO iks = 1,nks
              !
              nbndval = nbnd_occ(iks)
              !
              DO ibnd = 1,nbndval
                 ag(:,ibnd,iks,il2) = 0._DP
              ENDDO
              !
           ENDDO
         ELSE
           DO iks = 1,nks
              !
              nbndval = nbnd_occ(iks)
              !
              DO ibnd = 1,nbndval
                 ag(:,ibnd,iks,il2) = hg(:,ibnd,iks,il2)
              ENDDO
              !
           ENDDO
         ENDIF
         !
      ENDDO
      !
      DEALLOCATE(hg)
      !
      CALL stop_clock('refresh_vr')
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
