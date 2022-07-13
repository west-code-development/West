!
! Copyright (C) 2015-2021 M. Govoni
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
MODULE wstat_tools
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTERFACE diagox
     MODULE PROCEDURE diagox_dsy, diagox_zhe
  END INTERFACE
  !
  INTERFACE serial_diagox
     MODULE PROCEDURE serial_diagox_dsy, serial_diagox_zhe
  END INTERFACE
  !
  INTERFACE build_hr
     MODULE PROCEDURE build_hr_real, build_hr_complex
  END INTERFACE
  !
  INTERFACE redistribute_vr_distr
     MODULE PROCEDURE redistribute_vr_distr_real, redistribute_vr_distr_complex
  END INTERFACE
  !
  INTERFACE update_with_vr_distr
     MODULE PROCEDURE update_with_vr_distr_real, update_with_vr_distr_complex
  END INTERFACE
  !
  INTERFACE refresh_with_vr_distr
     MODULE PROCEDURE refresh_with_vr_distr_real, refresh_with_vr_distr_complex
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE diagox_dsy(nbase,nvec,hr_distr,nvecx,ew,vr_distr)
      !------------------------------------------------------------------------
      !
      USE io_global,             ONLY : stdout
      USE io_push,               ONLY : io_push_bar
      USE mp_world,              ONLY : nproc
      USE distribution_center,   ONLY : pert
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: nbase,nvec,nvecx
      REAL(DP),INTENT(IN) :: hr_distr(nvecx,pert%nlocx)
      REAL(DP),INTENT(OUT) :: ew(nvecx)
      REAL(DP),INTENT(OUT) :: vr_distr(nvecx,pert%nlocx)
      !
      REAL(DP) :: time_spent(2)
      REAL(DP),EXTERNAL :: GET_CLOCK
      CHARACTER(20),EXTERNAL :: human_readable_time
      LOGICAL :: l_parallel
      INTEGER :: npur,npuc
      CHARACTER(LEN=8) :: aux_label_npur
      CHARACTER(LEN=8) :: aux_label_npuc
      !
#if defined(__SCALAPACK)
#if defined(__CUDA)
      IF(nproc > 3 .AND. nbase > 8000) THEN
#else
      IF(nproc > 3 .AND. nbase > 2000) THEN
#endif
         l_parallel = .TRUE.
      ELSE
         l_parallel = .FALSE.
      ENDIF
#else
      l_parallel = .FALSE.
#endif
      !
      CALL start_clock('diagox')
      !
      IF(l_parallel) THEN
         !
         time_spent(1) = get_clock('diagox')
         !
#if defined(__ELPA) && defined(__CUDA)
         npur = MAX(FLOOR(SQRT(REAL(MIN(nproc,nbase/64),KIND=DP))),2)
#else
         npur = MAX(FLOOR(SQRT(REAL(MIN(nproc,nbase/8),KIND=DP))),2)
#endif
         npuc = npur
         !
#if defined(__SCALAPACK)
         CALL parallel_distributed_diago_dsy(nvec,nbase,nvecx,hr_distr,vr_distr,ew(1:nvec),npur,npuc,pert)
#endif
         !
         time_spent(2) = get_clock('diagox')
         WRITE(aux_label_npur,'(i8)') npur
         WRITE(aux_label_npuc,'(i8)') npuc
#if defined(__ELPA)
#if defined(__CUDA)
         WRITE(stdout,"(5x,'p-DIAGOX done in ',a,' with ELPA (GPU), grid ',a)") &
#else
         WRITE(stdout,"(5x,'p-DIAGOX done in ',a,' with ELPA, grid ',a)") &
#endif
#else
         WRITE(stdout,"(5x,'p-DIAGOX done in ',a,' with ScaLAPACK, grid ',a)") &
#endif
           & TRIM(ADJUSTL(human_readable_time(time_spent(2)-time_spent(1)))), &
           & '('//TRIM(ADJUSTL(aux_label_npur))//'x'//TRIM(ADJUSTL(aux_label_npuc))//')'
         !
      ELSE
         !
         time_spent(1) = get_clock('diagox')
         !
         CALL serial_diagox_dsy(nvec,nbase,nvecx,hr_distr,ew(1:nvec),vr_distr)
         !
         time_spent(2) = get_clock('diagox')
#if defined(__CUDA)
         WRITE(stdout,"(5x,'s-DIAGOX done in ',a,' with cuSOLVER (GPU)')") &
#else
         WRITE(stdout,"(5x,'s-DIAGOX done in ',a,' with LAPACK')") &
#endif
           & TRIM(ADJUSTL(human_readable_time(time_spent(2)-time_spent(1))))
         !
      ENDIF
      !
      CALL stop_clock('diagox')
      !
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE diagox_zhe(nbase,nvec,hr_distr,nvecx,ew,vr_distr)
      !------------------------------------------------------------------------
      !
      USE io_global,             ONLY : stdout
      USE io_push,               ONLY : io_push_bar
      USE mp_world,              ONLY : nproc
      USE distribution_center,   ONLY : pert
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: nbase,nvec,nvecx
      COMPLEX(DP),INTENT(IN) :: hr_distr(nvecx,pert%nlocx)
      REAL(DP),INTENT(OUT) :: ew(nvecx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(nvecx,pert%nlocx)
      !
      REAL(DP) :: time_spent(2)
      REAL(DP),EXTERNAL :: GET_CLOCK
      CHARACTER(20),EXTERNAL :: human_readable_time
      LOGICAL :: l_parallel
      INTEGER :: npur,npuc
      CHARACTER(LEN=8) :: aux_label_npur
      CHARACTER(LEN=8) :: aux_label_npuc
      !
#if defined(__SCALAPACK)
#if defined(__CUDA)
      IF(nproc > 3 .AND. nbase > 8000) THEN
#else
      IF(nproc > 3 .AND. nbase > 2000) THEN
#endif
         l_parallel = .TRUE.
      ELSE
         l_parallel = .FALSE.
      ENDIF
#else
      l_parallel = .FALSE.
#endif
      !
      CALL start_clock('diagox')
      !
      IF(l_parallel) THEN
         !
         time_spent(1) = get_clock('diagox')
         !
#if defined(__ELPA) && defined(__CUDA)
         npur = MAX(FLOOR(SQRT(REAL(MIN(nproc,nbase/64),KIND=DP))),2)
#else
         npur = MAX(FLOOR(SQRT(REAL(MIN(nproc,nbase/8),KIND=DP))),2)
#endif
         npuc = npur
         !
#if defined(__SCALAPACK)
         CALL parallel_distributed_diago_zhe(nvec,nbase,nvecx,hr_distr,vr_distr,ew(1:nvec),npur,npuc,pert)
#endif
         !
         time_spent(2) = get_clock('diagox')
         WRITE(aux_label_npur,'(i8)') npur
         WRITE(aux_label_npuc,'(i8)') npuc
#if defined(__ELPA)
#if defined(__CUDA)
         WRITE(stdout,"(5x,'p-DIAGOX done in ',a,' with ELPA (GPU), grid ',a)") &
#else
         WRITE(stdout,"(5x,'p-DIAGOX done in ',a,' with ELPA, grid ',a)") &
#endif
#else
         WRITE(stdout,"(5x,'p-DIAGOX done in ',a,' with ScaLAPACK, grid ',a)") &
#endif
           & TRIM(ADJUSTL(human_readable_time(time_spent(2)-time_spent(1)))), &
           & '('//TRIM(ADJUSTL(aux_label_npur))//'x'//TRIM(ADJUSTL(aux_label_npuc))//')'
         !
      ELSE
         !
         time_spent(1) = get_clock('diagox')
         !
         CALL serial_diagox_zhe(nvec,nbase,nvecx,hr_distr,ew(1:nvec),vr_distr)
         !
         time_spent(2) = get_clock('diagox')
#if defined(__CUDA)
         WRITE(stdout,"(5x,'s-DIAGOX done in ',a,' with cuSOLVER (GPU)')") &
#else
         WRITE(stdout,"(5x,'s-DIAGOX done in ',a,' with LAPACK')") &
#endif
           & TRIM(ADJUSTL(human_readable_time(time_spent(2)-time_spent(1))))
         !
      ENDIF
      !
      CALL stop_clock('diagox')
      !
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE serial_diagox_dsy(nselect,n,lda,hr_distr,e,vr_distr)
      !------------------------------------------------------------------------
      !
      ! Diagox -- serial
      !   nselect : number of wanted ev
      !   n       : actual dimension of a
      !   lda     : leading dimension of a
      !   a       : matrix to be diagox
      !   e       : eigenval(1:nselsect), even though it is defined 1:lda
      !   z       : unitary trans.
      !
      USE mp,                    ONLY : mp_bcast,mp_sum
      USE mp_global,             ONLY : me_bgrp,root_bgrp,intra_bgrp_comm,inter_image_comm
      USE linear_algebra_kernel, ONLY : matdiago_dsy
      USE distribution_center,   ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nselect,lda,n
      REAL(DP),INTENT(IN) :: hr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(OUT) :: e(nselect)
      REAL(DP),INTENT(OUT) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      REAL(DP),ALLOCATABLE :: zz(:,:),ee(:)
      INTEGER :: i,il2,ig1,ig2
      !
      ALLOCATE(zz(n,n))
      ALLOCATE(ee(n))
      !
      ! Create zz from hr_distr
      !
      zz = 0._DP
      DO il2 = 1,pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2 > n) CYCLE
         DO ig1 = 1,pert%nglob
            IF(ig1 > n) CYCLE
            zz(ig1,ig2) = hr_distr(ig1,il2)
         ENDDO
      ENDDO
      CALL mp_sum(zz,inter_image_comm)
      ee = 0._DP
      !
      IF(me_bgrp == root_bgrp) THEN
         CALL matdiago_dsy(n,zz,ee,.FALSE.)
      ENDIF
      !
      CALL mp_bcast(ee,root_bgrp,intra_bgrp_comm)
      CALL mp_bcast(zz,root_bgrp,intra_bgrp_comm)
      !
      DO il2 = 1,pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2 > nselect) CYCLE
         DO ig1 = 1,pert%nglob
            IF(ig1 > n) CYCLE
            vr_distr(ig1,il2) = zz(ig1,ig2)
         ENDDO
      ENDDO
      DO i = 1,nselect
         e(i) = ee(i)
      ENDDO
      !
      DEALLOCATE(zz)
      DEALLOCATE(ee)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE serial_diagox_zhe(nselect,n,lda,hr_distr,e,vr_distr)
      !------------------------------------------------------------------------
      !
      ! Diagox -- serial
      !   nselect : number of wanted ev
      !   n       : actual dimension of a
      !   lda     : leading dimension of a
      !   a       : matrix to be diagox
      !   e       : eigenval(1:nselsect), even though it is defined 1:lda
      !   z       : unitary trans.
      !
      USE mp,                    ONLY : mp_bcast,mp_sum
      USE mp_global,             ONLY : me_bgrp,root_bgrp,intra_bgrp_comm,inter_image_comm
      USE linear_algebra_kernel, ONLY : matdiago_zhe
      USE distribution_center,   ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nselect,lda,n
      COMPLEX(DP),INTENT(IN) :: hr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(OUT) :: e(nselect)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: zz(:,:)
      REAL(DP),ALLOCATABLE :: ee(:)
      INTEGER :: i,il2,ig1,ig2
      !
      ALLOCATE(zz(n,n))
      ALLOCATE(ee(n))
      !
      ! Create zz from hr_distr
      !
      zz = 0._DP
      DO il2 = 1,pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2 > n) CYCLE
         DO ig1 = 1,pert%nglob
            IF(ig1 > n) CYCLE
            zz(ig1,ig2) = hr_distr(ig1,il2)
         ENDDO
      ENDDO
      CALL mp_sum(zz,inter_image_comm)
      ee = 0._DP
      !
      IF(me_bgrp == root_bgrp) THEN
         CALL matdiago_zhe(n,zz,ee,.FALSE.)
      ENDIF
      !
      CALL mp_bcast(ee,root_bgrp,intra_bgrp_comm)
      CALL mp_bcast(zz,root_bgrp,intra_bgrp_comm)
      !
      DO il2 = 1,pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2 > nselect) CYCLE
         DO ig1 = 1,pert%nglob
            IF(ig1 > n) CYCLE
            vr_distr(ig1,il2) = zz(ig1,ig2)
         ENDDO
      ENDDO
      DO i = 1,nselect
         e(i) = ee(i)
      ENDDO
      !
      DEALLOCATE(zz)
      DEALLOCATE(ee)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE build_hr_real(ag,bg,l2_s,l2_e,c_distr,g_e)
      !------------------------------------------------------------------------
      !
      !  c_distr = < ag | bg >
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm,&
                                     & inter_bgrp_comm,my_bgrp_id,inter_pool_comm,my_pool_id
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_bcast
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq,npwqx
      USE gvect,                ONLY : gstart
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwqx,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: bg(npwqx,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      REAL(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_e
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,l3_s,ig1
      INTEGER :: icycl,idx,nloc,pert_nglob,pert_nlocx
      REAL(DP) :: reduce
#if !defined(__CUDA)
      REAL(DP),EXTERNAL :: DDOT
#endif
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
            pert_nglob = pert%nglob
            pert_nlocx = pert%nlocx
            !
            !$acc enter data create(ag,c_distr(1:pert_nglob,l2_s:l2_e)) copyin(bg)
            !
            !$acc parallel loop collapse(2)
            DO il2 = l2_s,l2_e
               DO il1 = 1,pert_nglob
                  c_distr(il1,il2) = 0._DP
               ENDDO
            ENDDO
            !$acc end parallel
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
#if defined(__CUDA)
                  IF(gstart == 1) THEN
                     l3_s = 1
                  ELSE
                     l3_s = 2
                  ENDIF
                  !
                  !$acc parallel async
                  !$acc loop
                  DO il2 = l2_s,l2_e
                     reduce = 0._DP
                     !$acc loop reduction(+:reduce)
                     DO il3 = l3_s,npwq
                        reduce = reduce+REAL(ag(il3,il1),KIND=DP)*REAL(bg(il3,il2),KIND=DP) &
                        & +AIMAG(ag(il3,il1))*AIMAG(bg(il3,il2))
                     ENDDO
                     c_distr(ig1,il2) = 2._DP*reduce
                  ENDDO
                  !$acc end parallel
#else
                  DO il2 = l2_s,l2_e
                     !
                     IF(gstart == 1) THEN
                        c_distr(ig1,il2) = 2._DP*DDOT(2*npwq,ag(1,il1),1,bg(1,il2),1)
                     ELSE
                        c_distr(ig1,il2) = 2._DP*DDOT(2*npwq-2,ag(2,il1),1,bg(2,il2),1)
                     ENDIF
                     !
                  ENDDO
#endif
                  !
               ENDDO
               !
            ENDIF
            !
            ! Cycle ag
            !
            CALL mp_circular_shift_left(ag,icycl,inter_image_comm)
            !
         ENDDO
         !
         IF(l2_e >= l2_s) THEN
            !
            !$acc update host(c_distr(1:pert_nglob,l2_s:l2_e)) wait
            !$acc exit data delete(ag,bg,c_distr(1:pert_nglob,l2_s:l2_e))
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
    SUBROUTINE build_hr_complex(ag,bg,l2_s,l2_e,c_distr,g_e)
      !------------------------------------------------------------------------
      !
      !  c_distr = < ag | bg >
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,intra_bgrp_comm,&
                                     & inter_bgrp_comm,my_bgrp_id,inter_pool_comm,my_pool_id
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_bcast
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq,npwqx
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwqx,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: bg(npwqx,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      COMPLEX(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_e
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1
      INTEGER :: icycl,idx,nloc,pert_nglob,pert_nlocx
      COMPLEX(DP) :: reduce
#if !defined(__CUDA)
      COMPLEX(DP),EXTERNAL :: ZDOTC
#endif
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
            pert_nglob = pert%nglob
            pert_nlocx = pert%nlocx
            !
            !$acc enter data create(ag,c_distr(1:pert_nglob,l2_s:l2_e)) copyin(bg)
            !
            !$acc parallel loop collapse(2)
            DO il2 = l2_s,l2_e
               DO il1 = 1,pert_nglob
                  c_distr(il1,il2) = 0._DP
               ENDDO
            ENDDO
            !$acc end parallel
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
#if defined(__CUDA)
                  !$acc parallel async
                  !$acc loop
                  DO il2 = l2_s,l2_e
                     reduce = 0._DP
                     !$acc loop reduction(+:reduce)
                     DO il3 = 1,npwq
                        reduce = reduce+CONJG(ag(il3,il1))*bg(il3,il2)
                     ENDDO
                     c_distr(ig1,il2) = reduce
                  ENDDO
                  !$acc end parallel
#else
                  DO il2 = l2_s,l2_e
                     !
                     c_distr(ig1,il2) = ZDOTC(npwq,ag(1,il1),1,bg(1,il2),1)
                     !
                  ENDDO
#endif
                  !
               ENDDO
               !
            ENDIF
            !
            ! Cycle ag
            !
            CALL mp_circular_shift_left(ag,icycl,inter_image_comm)
            !
         ENDDO
         !
         IF(l2_e >= l2_s) THEN
            !
            !$acc update host(c_distr(1:pert_nglob,l2_s:l2_e)) wait
            !$acc exit data delete(ag,bg,c_distr(1:pert_nglob,l2_s:l2_e))
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
    SUBROUTINE redistribute_vr_distr_real(nselect,n,lda,vr_distr,ishift)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : me_bgrp,my_bgrp_id,nbgrp,nproc_bgrp,intra_bgrp_comm,&
                                     & my_image_id,my_pool_id,npool
      USE mp_world,             ONLY : world_comm,nproc
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : pert
      USE sort_tools,           ONLY : heapsort
      USE west_mp,              ONLY : mp_alltoallv
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(INOUT) :: vr_distr(lda,pert%nlocx)
      INTEGER,INTENT(IN) :: ishift(lda)
      !
      ! Workspace
      !
      INTEGER :: i_col,n_col,j_loc,j_loc2,j_glob,i_proc
      INTEGER :: this_dest,this_sour
      INTEGER,ALLOCATABLE :: idx_send(:)
      INTEGER,ALLOCATABLE :: idx_recv(:)
      INTEGER,ALLOCATABLE :: send_count(:)
      INTEGER,ALLOCATABLE :: recv_count(:)
      INTEGER,ALLOCATABLE :: send_displ(:)
      INTEGER,ALLOCATABLE :: recv_displ(:)
      INTEGER,ALLOCATABLE :: dest(:)
      INTEGER,ALLOCATABLE :: swap(:)
      INTEGER,ALLOCATABLE :: tmp_i(:)
      REAL(DP),ALLOCATABLE :: tmp_r(:,:)
      REAL(DP),ALLOCATABLE :: val_send(:,:)
      REAL(DP),ALLOCATABLE :: val_recv(:,:)
      !
      CALL start_clock('redistr_vr')
      !
      ! vr_distr only needed by pool 0 and band group 0 in the next step
      !
      ALLOCATE(send_count(nproc))
      ALLOCATE(recv_count(nproc))
      !
      send_count = 0
      recv_count = 0
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
         n_col = 0
         !
         DO j_glob = n+1,n+nselect
            !
            IF(ishift(j_glob) /= 0) THEN
               CALL pert%g2l(ishift(j_glob),j_loc,this_sour)
               !
               IF(this_sour == my_image_id) THEN
                  n_col = n_col+1
               ENDIF
            ENDIF
            !
         ENDDO
         !
         ALLOCATE(val_send(lda,n_col))
         ALLOCATE(idx_send(n_col))
         ALLOCATE(dest(n_col))
         ALLOCATE(swap(n_col))
         !
         n_col = 0
         !
         DO j_glob = n+1,n+nselect
            !
            IF(ishift(j_glob) /= 0) THEN
               CALL pert%g2l(ishift(j_glob),j_loc,this_sour)
               !
               IF(this_sour == my_image_id) THEN
                  n_col = n_col+1
                  !
                  CALL pert%g2l(j_glob,j_loc2,this_dest)
                  !
                  val_send(:,n_col) = vr_distr(:,j_loc)
                  idx_send(n_col) = j_loc2
                  dest(n_col) = this_dest*npool*nbgrp*nproc_bgrp
                  send_count(dest(n_col)+1) = send_count(dest(n_col)+1)+1
               ENDIF
            ENDIF
            !
         ENDDO
         !
         CALL heapsort(n_col,dest,swap)
         !
         DEALLOCATE(dest)
         ALLOCATE(tmp_r(lda,n_col))
         ALLOCATE(tmp_i(n_col))
         !
         DO i_col = 1,n_col
            tmp_r(:,i_col) = val_send(:,swap(i_col))
            tmp_i(i_col) = idx_send(swap(i_col))
         ENDDO
         !
         val_send = tmp_r
         idx_send = tmp_i
         !
         DEALLOCATE(swap)
         DEALLOCATE(tmp_r)
         DEALLOCATE(tmp_i)
         !
         n_col = 0
         !
         DO j_loc = 1,pert%nloc
            !
            j_glob = pert%l2g(j_loc)
            !
            IF(j_glob > n .AND. j_glob <= n+nselect .AND. ishift(j_glob) /= 0) THEN
               n_col = n_col+1
               !
               CALL pert%g2l(ishift(j_glob),j_loc2,this_sour)
               !
               this_sour = this_sour*npool*nbgrp*nproc_bgrp
               recv_count(this_sour+1) = recv_count(this_sour+1)+1
            ENDIF
            !
         ENDDO
         !
         ALLOCATE(val_recv(lda,n_col))
         ALLOCATE(idx_recv(n_col))
         !
      ELSE
         ALLOCATE(val_send(1,1))
         ALLOCATE(idx_send(1))
         ALLOCATE(val_recv(1,1))
         ALLOCATE(idx_recv(1))
      ENDIF
      !
      ALLOCATE(send_displ(nproc))
      ALLOCATE(recv_displ(nproc))
      !
      send_displ = 0
      recv_displ = 0
      !
      DO i_proc = 2,nproc
         send_displ(i_proc) = SUM(send_count(1:i_proc-1))
         recv_displ(i_proc) = SUM(recv_count(1:i_proc-1))
      ENDDO
      !
      CALL mp_alltoallv(idx_send,send_count,send_displ,idx_recv,recv_count,recv_displ,world_comm)
      !
      send_count = send_count*lda
      recv_count = recv_count*lda
      send_displ = send_displ*lda
      recv_displ = recv_displ*lda
      !
      CALL mp_alltoallv(val_send,send_count,send_displ,val_recv,recv_count,recv_displ,world_comm)
      !
      DEALLOCATE(val_send)
      DEALLOCATE(idx_send)
      DEALLOCATE(send_count)
      DEALLOCATE(recv_count)
      DEALLOCATE(send_displ)
      DEALLOCATE(recv_displ)
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         IF(me_bgrp == 0) THEN
            DO i_col = 1,n_col
               vr_distr(:,idx_recv(i_col)) = val_recv(:,i_col)
            ENDDO
         ENDIF
         !
         CALL mp_bcast(vr_distr,0,intra_bgrp_comm)
      ENDIF
      !
      DEALLOCATE(val_recv)
      DEALLOCATE(idx_recv)
      !
      CALL stop_clock('redistr_vr')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE redistribute_vr_distr_complex(nselect,n,lda,vr_distr,ishift)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : me_bgrp,my_bgrp_id,nbgrp,nproc_bgrp,intra_bgrp_comm,&
                                     & my_image_id,my_pool_id,npool
      USE mp_world,             ONLY : world_comm,nproc
      USE mp,                   ONLY : mp_bcast
      USE distribution_center,  ONLY : pert
      USE sort_tools,           ONLY : heapsort
      USE west_mp,              ONLY : mp_alltoallv
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nselect,n,lda
      COMPLEX(DP),INTENT(INOUT) :: vr_distr(lda,pert%nlocx)
      INTEGER,INTENT(IN) :: ishift(lda)
      !
      ! Workspace
      !
      INTEGER :: i_col,n_col,j_loc,j_loc2,j_glob,i_proc
      INTEGER :: this_dest,this_sour
      INTEGER,ALLOCATABLE :: idx_send(:)
      INTEGER,ALLOCATABLE :: idx_recv(:)
      INTEGER,ALLOCATABLE :: send_count(:)
      INTEGER,ALLOCATABLE :: recv_count(:)
      INTEGER,ALLOCATABLE :: send_displ(:)
      INTEGER,ALLOCATABLE :: recv_displ(:)
      INTEGER,ALLOCATABLE :: dest(:)
      INTEGER,ALLOCATABLE :: swap(:)
      INTEGER,ALLOCATABLE :: tmp_i(:)
      COMPLEX(DP),ALLOCATABLE :: tmp_c(:,:)
      COMPLEX(DP),ALLOCATABLE :: val_send(:,:)
      COMPLEX(DP),ALLOCATABLE :: val_recv(:,:)
      !
      CALL start_clock('redistr_vr')
      !
      ! vr_distr only needed by pool 0 and band group 0 in the next step
      !
      ALLOCATE(send_count(nproc))
      ALLOCATE(recv_count(nproc))
      !
      send_count = 0
      recv_count = 0
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
         n_col = 0
         !
         DO j_glob = n+1,n+nselect
            !
            IF(ishift(j_glob) /= 0) THEN
               CALL pert%g2l(ishift(j_glob),j_loc,this_sour)
               !
               IF(this_sour == my_image_id) THEN
                  n_col = n_col+1
               ENDIF
            ENDIF
            !
         ENDDO
         !
         ALLOCATE(val_send(lda,n_col))
         ALLOCATE(idx_send(n_col))
         ALLOCATE(dest(n_col))
         ALLOCATE(swap(n_col))
         !
         n_col = 0
         !
         DO j_glob = n+1,n+nselect
            !
            IF(ishift(j_glob) /= 0) THEN
               CALL pert%g2l(ishift(j_glob),j_loc,this_sour)
               !
               IF(this_sour == my_image_id) THEN
                  n_col = n_col+1
                  !
                  CALL pert%g2l(j_glob,j_loc2,this_dest)
                  !
                  val_send(:,n_col) = vr_distr(:,j_loc)
                  idx_send(n_col) = j_loc2
                  dest(n_col) = this_dest*npool*nbgrp*nproc_bgrp
                  send_count(dest(n_col)+1) = send_count(dest(n_col)+1)+1
               ENDIF
            ENDIF
            !
         ENDDO
         !
         CALL heapsort(n_col,dest,swap)
         !
         DEALLOCATE(dest)
         ALLOCATE(tmp_c(lda,n_col))
         ALLOCATE(tmp_i(n_col))
         !
         DO i_col = 1,n_col
            tmp_c(:,i_col) = val_send(:,swap(i_col))
            tmp_i(i_col) = idx_send(swap(i_col))
         ENDDO
         !
         val_send = tmp_c
         idx_send = tmp_i
         !
         DEALLOCATE(swap)
         DEALLOCATE(tmp_c)
         DEALLOCATE(tmp_i)
         !
         n_col = 0
         !
         DO j_loc = 1,pert%nloc
            !
            j_glob = pert%l2g(j_loc)
            !
            IF(j_glob > n .AND. j_glob <= n+nselect .AND. ishift(j_glob) /= 0) THEN
               n_col = n_col+1
               !
               CALL pert%g2l(ishift(j_glob),j_loc2,this_sour)
               !
               this_sour = this_sour*npool*nbgrp*nproc_bgrp
               recv_count(this_sour+1) = recv_count(this_sour+1)+1
            ENDIF
            !
         ENDDO
         !
         ALLOCATE(val_recv(lda,n_col))
         ALLOCATE(idx_recv(n_col))
         !
      ELSE
         ALLOCATE(val_send(1,1))
         ALLOCATE(idx_send(1))
         ALLOCATE(val_recv(1,1))
         ALLOCATE(idx_recv(1))
      ENDIF
      !
      ALLOCATE(send_displ(nproc))
      ALLOCATE(recv_displ(nproc))
      !
      send_displ = 0
      recv_displ = 0
      !
      DO i_proc = 2,nproc
         send_displ(i_proc) = SUM(send_count(1:i_proc-1))
         recv_displ(i_proc) = SUM(recv_count(1:i_proc-1))
      ENDDO
      !
      CALL mp_alltoallv(idx_send,send_count,send_displ,idx_recv,recv_count,recv_displ,world_comm)
      !
      send_count = send_count*lda
      recv_count = recv_count*lda
      send_displ = send_displ*lda
      recv_displ = recv_displ*lda
      !
      CALL mp_alltoallv(val_send,send_count,send_displ,val_recv,recv_count,recv_displ,world_comm)
      !
      DEALLOCATE(val_send)
      DEALLOCATE(idx_send)
      DEALLOCATE(send_count)
      DEALLOCATE(recv_count)
      DEALLOCATE(send_displ)
      DEALLOCATE(recv_displ)
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         IF(me_bgrp == 0) THEN
            DO i_col = 1,n_col
               vr_distr(:,idx_recv(i_col)) = val_recv(:,i_col)
            ENDDO
         ENDIF
         !
         CALL mp_bcast(vr_distr,0,intra_bgrp_comm)
      ENDIF
      !
      DEALLOCATE(val_recv)
      DEALLOCATE(idx_recv)
      !
      CALL stop_clock('redistr_vr')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_real(ag,bg,nselect,n,lda,vr_distr,ew)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_bgrp_comm,&
                                     & my_bgrp_id,inter_pool_comm,my_pool_id
      USE mp,                   ONLY : mp_circular_shift_left,mp_bcast
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq,npwqx
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwqx,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: bg(npwqx,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2
      INTEGER :: icycl,idx,nloc,pert_nlocx
      REAL(DP) :: dconst
      COMPLEX(DP) :: zconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
      COMPLEX(DP),ALLOCATABLE :: hg2(:,:)
#if defined(__CUDA)
      ATTRIBUTES(PINNED) :: hg,hg2
#endif
      !
#if defined(__CUDA)
      CALL start_clock_gpu('update_vr')
#else
      CALL start_clock('update_vr')
#endif
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         pert_nlocx = pert%nlocx
         !
         ALLOCATE(hg(npwqx,pert%nlocx))
         ALLOCATE(hg2(npwqx,pert%nlocx))
         !
         !$acc enter data create(ag,bg,hg,hg2)
         !
         !$acc parallel loop collapse(2)
         DO il2 = 1,pert_nlocx
            DO il1 = 1,npwqx
               hg(il1,il2) = 0._DP
               hg2(il1,il2) = 0._DP
            ENDDO
         ENDDO
         !$acc end parallel
         !
         DO icycl = 0,nimage-1
            !
            idx = MOD(my_image_id+icycl,nimage)
            nloc = pert%nglob/nimage
            IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
            !
            !$acc update device(ag,bg) wait
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
#if defined(__CUDA)
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc parallel loop async
                  DO il3 = 1,npwq
                     hg(il3,il2) = dconst*ag(il3,il1)+hg(il3,il2)
                     hg2(il3,il2) = dconst*bg(il3,il1)+hg2(il3,il2)
                  ENDDO
                  !$acc end parallel
#else
                  zconst = vr_distr(ig1,il2)
                  CALL ZAXPY(npwq,zconst,ag(1,il1),1,hg(1,il2),1)
                  CALL ZAXPY(npwq,zconst,bg(1,il1),1,hg2(1,il2),1)
#endif
                  !
               ENDDO
            ENDDO
            !
            ! Cycle ag, bg
            !
            CALL mp_circular_shift_left(ag,icycl,inter_image_comm)
            CALL mp_circular_shift_left(bg,icycl+nimage,inter_image_comm)
            !
         ENDDO
         !
         !$acc update host(hg,hg2) wait
         !$acc exit data delete(ag,bg,hg,hg2)
         !
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
            ag(:,il2) = -ew(ig2)*hg(:,il2)+hg2(:,il2)
         ENDDO
         !
         DEALLOCATE(hg)
         DEALLOCATE(hg2)
         !
      ENDIF
      !
      CALL mp_bcast(ag,0,inter_bgrp_comm)
      CALL mp_bcast(ag,0,inter_pool_comm)
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
    SUBROUTINE update_with_vr_distr_complex(ag,bg,nselect,n,lda,vr_distr,ew)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_bgrp_comm,&
                                     & my_bgrp_id,inter_pool_comm,my_pool_id
      USE mp,                   ONLY : mp_circular_shift_left,mp_bcast
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq,npwqx
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwqx,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: bg(npwqx,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2
      INTEGER :: icycl,idx,nloc,pert_nlocx
      COMPLEX(DP) :: zconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
      COMPLEX(DP),ALLOCATABLE :: hg2(:,:)
#if defined(__CUDA)
      ATTRIBUTES(PINNED) :: hg,hg2
#endif
      !
#if defined(__CUDA)
      CALL start_clock_gpu('update_vr')
#else
      CALL start_clock('update_vr')
#endif
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         pert_nlocx = pert%nlocx
         !
         ALLOCATE(hg(npwqx,pert%nlocx))
         ALLOCATE(hg2(npwqx,pert%nlocx))
         !
         !$acc enter data create(ag,bg,hg,hg2)
         !
         !$acc parallel loop collapse(2)
         DO il2 = 1,pert_nlocx
            DO il1 = 1,npwqx
               hg(il1,il2) = 0._DP
               hg2(il1,il2) = 0._DP
            ENDDO
         ENDDO
         !$acc end parallel
         !
         DO icycl = 0,nimage-1
            !
            idx = MOD(my_image_id+icycl,nimage)
            nloc = pert%nglob/nimage
            IF(idx < MOD(pert%nglob,nimage)) nloc = nloc+1
            !
            !$acc update device(ag,bg) wait
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
#if defined(__CUDA)
                  zconst = vr_distr(ig1,il2)
                  !
                  !$acc parallel loop async
                  DO il3 = 1,npwq
                     hg(il3,il2) = zconst*ag(il3,il1)+hg(il3,il2)
                     hg2(il3,il2) = zconst*bg(il3,il1)+hg2(il3,il2)
                  ENDDO
                  !$acc end parallel
#else
                  CALL ZAXPY(npwq,vr_distr(ig1,il2),ag(1,il1),1,hg(1,il2),1)
                  CALL ZAXPY(npwq,vr_distr(ig1,il2),bg(1,il1),1,hg2(1,il2),1)
#endif
                  !
               ENDDO
            ENDDO
            !
            ! Cycle ag, bg
            !
            CALL mp_circular_shift_left(ag,icycl,inter_image_comm)
            CALL mp_circular_shift_left(bg,icycl+nimage,inter_image_comm)
            !
         ENDDO
         !
         !$acc update host(hg,hg2) wait
         !$acc exit data delete(ag,bg,hg,hg2)
         !
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 <= n .OR. ig2 > n+nselect) CYCLE
            ag(:,il2) = -ew(ig2)*hg(:,il2)+hg2(:,il2)
         ENDDO
         !
         DEALLOCATE(hg)
         DEALLOCATE(hg2)
         !
      ENDIF
      !
      CALL mp_bcast(ag,0,inter_bgrp_comm)
      CALL mp_bcast(ag,0,inter_pool_comm)
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
      USE mp,                   ONLY : mp_circular_shift_left,mp_bcast
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq,npwqx
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwqx,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2
      INTEGER :: icycl,idx,nloc,pert_nlocx
      REAL(DP) :: dconst
      COMPLEX(DP) :: zconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
#if defined(__CUDA)
      ATTRIBUTES(PINNED) :: hg
#endif
      !
#if defined(__CUDA)
      CALL start_clock_gpu('refresh_vr')
#else
      CALL start_clock('refresh_vr')
#endif
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         pert_nlocx = pert%nlocx
         !
         ALLOCATE(hg(npwqx,pert%nlocx))
         !
         !$acc enter data create(ag,hg)
         !
         !$acc parallel loop collapse(2)
         DO il2 = 1,pert_nlocx
            DO il1 = 1,npwqx
               hg(il1,il2) = 0._DP
            ENDDO
         ENDDO
         !$acc end parallel
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
#if defined(__CUDA)
                  dconst = vr_distr(ig1,il2)
                  !
                  !$acc parallel loop async
                  DO il3 = 1,npwq
                     hg(il3,il2) = dconst*ag(il3,il1)+hg(il3,il2)
                  ENDDO
                  !$acc end parallel
#else
                  zconst = vr_distr(ig1,il2)
                  CALL ZAXPY(npwq,zconst,ag(1,il1),1,hg(1,il2),1)
#endif
                  !
               ENDDO
            ENDDO
            !
            ! Cycle ag
            !
            CALL mp_circular_shift_left(ag,icycl,inter_image_comm)
            !
         ENDDO
         !
         !$acc update host(hg) wait
         !$acc exit data delete(ag,hg)
         !
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 > nselect) THEN
               ag(:,il2) = 0._DP
            ELSE
               ag(:,il2) = hg(:,il2)
            ENDIF
         ENDDO
         !
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
    SUBROUTINE refresh_with_vr_distr_complex(ag,nselect,n,lda,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,inter_bgrp_comm,&
                                     & my_bgrp_id,inter_pool_comm,my_pool_id
      USE mp,                   ONLY : mp_circular_shift_left,mp_bcast
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq,npwqx
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwqx,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect,n,lda
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: il1,il2,il3,ig1,ig2
      INTEGER :: icycl,idx,nloc,pert_nlocx
      COMPLEX(DP) :: zconst
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
#if defined(__CUDA)
      ATTRIBUTES(PINNED) :: hg
#endif
      !
#if defined(__CUDA)
      CALL start_clock_gpu('refresh_vr')
#else
      CALL start_clock('refresh_vr')
#endif
      !
      IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
         !
         pert_nlocx = pert%nlocx
         !
         ALLOCATE(hg(npwqx,pert%nlocx))
         !
         !$acc enter data create(ag,hg)
         !
         !$acc parallel loop collapse(2)
         DO il2 = 1,pert_nlocx
            DO il1 = 1,npwqx
               hg(il1,il2) = 0._DP
            ENDDO
         ENDDO
         !$acc end parallel
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
#if defined(__CUDA)
                  zconst = vr_distr(ig1,il2)
                  !
                  !$acc parallel loop async
                  DO il3 = 1,npwq
                     hg(il3,il2) = zconst*ag(il3,il1)+hg(il3,il2)
                  ENDDO
                  !$acc end parallel
#else
                  CALL ZAXPY(npwq,vr_distr(ig1,il2),ag(1,il1),1,hg(1,il2),1)
#endif
                  !
               ENDDO
            ENDDO
            !
            ! Cycle ag
            !
            CALL mp_circular_shift_left(ag,icycl,inter_image_comm)
            !
         ENDDO
         !
         !$acc update host(hg) wait
         !$acc exit data delete(ag,hg)
         !
         DO il2 = 1,pert%nloc
            ig2 = pert%l2g(il2)
            IF(ig2 > nselect) THEN
               ag(:,il2) = 0._DP
            ELSE
               ag(:,il2) = hg(:,il2)
            ENDIF
         ENDDO
         !
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
END MODULE
