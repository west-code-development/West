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
MODULE wstat_tools
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
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
  INTERFACE symm_hr_distr
     MODULE PROCEDURE symm_hr_distr_real, symm_hr_distr_complex
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
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE diagox_dsy( nbase, nvec, hr_distr, nvecx, ew, vr_distr  ) 
      !------------------------------------------------------------------------
      !
      USE io_global,             ONLY : stdout
      USE io_push,               ONLY : io_push_title,io_push_bar
      USE mp_global,             ONLY : nproc_bgrp,inter_image_comm,nimage 
      USE mp_world,              ONLY : nproc
      USE westcom,               ONLY : n_pdep_eigen
      USE io_push,               ONLY : io_push_title,io_push_bar
      USE distribution_center,   ONLY : pert
      USE mp,                    ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER  :: nbase, nvec, nvecx
      REAL(DP) :: hr_distr(nvecx,pert%nlocx)
      REAL(DP) :: ew(nvecx)
      REAL(DP) :: vr_distr(nvecx,pert%nlocx)
      !
      REAL(DP) :: time_spent(2)
      REAL(DP), EXTERNAL :: GET_CLOCK
      CHARACTER(20),EXTERNAL :: human_readable_time
      LOGICAL :: l_parallel 
      INTEGER :: npur, npuc
      CHARACTER(LEN=8) :: aux_label_npur
      CHARACTER(LEN=8) :: aux_label_npuc
      !
      !
#if defined __SCALAPACK
      l_parallel = nproc > 3 
#else
      l_parallel = .FALSE.
#endif
      !
      CALL start_clock( 'diagox' )
      !
      IF( l_parallel ) THEN 
         !
         time_spent(1)=get_clock( 'diagox' )
#if defined __SCALAPACK
         !
         npur = MAX( MIN(INT( SQRT(DBLE(nproc))),INT(DBLE(nbase)/32.)), 2)
         npuc = npur
         CALL parallel_distributed_diago_dsy ( nvec, nbase, nvecx, hr_distr, vr_distr, ew(1:nvec), npur, npuc, pert)
         !
#endif
         time_spent(2)=get_clock( 'diagox' )
         WRITE(aux_label_npur,'(i8)') npur
         WRITE(aux_label_npuc,'(i8)') npuc
         WRITE(stdout, "(5x, 'p-DIAGOX done in ',a,' with a SCALAPACK grid ',a)") &
          & TRIM(ADJUSTL(human_readable_time(time_spent(2)-time_spent(1)))), &
          & "("//TRIM(ADJUSTL(aux_label_npur))//"x"//TRIM(ADJUSTL(aux_label_npuc))//")"
         !
      ELSE
         !
         time_spent(1)=get_clock( 'diagox' )
         !
         CALL serial_diagox_dsy( nvec, nbase, nvecx, hr_distr, ew(1:nvec), vr_distr )
         !
         time_spent(2)=get_clock( 'diagox' )
         WRITE(stdout, "(5x, 's-DIAGOX done in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
         !
      ENDIF
      !
      CALL stop_clock( 'diagox' )
      !
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE diagox_zhe( nbase, nvec, hr_distr, nvecx, ew, vr_distr  ) 
      !------------------------------------------------------------------------
      !
      USE io_global,             ONLY : stdout
      USE io_push,               ONLY : io_push_title,io_push_bar
      USE mp_global,             ONLY : nproc_bgrp,inter_image_comm,nimage 
      USE mp_world,              ONLY : nproc
      USE westcom,               ONLY : n_pdep_eigen
      USE io_push,               ONLY : io_push_title,io_push_bar
      USE distribution_center,   ONLY : pert
      USE mp,                    ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER  :: nbase, nvec, nvecx
      COMPLEX(DP) :: hr_distr(nvecx,pert%nlocx)
      REAL(DP) :: ew(nvecx)
      COMPLEX(DP) :: vr_distr(nvecx,pert%nlocx)
      !
      REAL(DP) :: time_spent(2)
      REAL(DP), EXTERNAL :: GET_CLOCK
      CHARACTER(20),EXTERNAL :: human_readable_time
      LOGICAL :: l_parallel 
      INTEGER :: npur, npuc
      CHARACTER(LEN=8) :: aux_label_npur
      CHARACTER(LEN=8) :: aux_label_npuc
      !
      !
#if defined __SCALAPACK
      l_parallel = nproc > 3 
#else
      l_parallel = .FALSE.
#endif
      !
      CALL start_clock( 'diagox' )
      !
      IF( l_parallel ) THEN 
         !
         time_spent(1)=get_clock( 'diagox' )
#if defined __SCALAPACK
         !
         npur = MAX( MIN(INT( SQRT(DBLE(nproc))),INT(DBLE(nbase)/32.)), 2)
         npuc = npur
         CALL parallel_distributed_diago_zhe ( nvec, nbase, nvecx, hr_distr, vr_distr, ew(1:nvec), npur, npuc, pert)
         !
#endif
         time_spent(2)=get_clock( 'diagox' )
         WRITE(aux_label_npur,'(i8)') npur
         WRITE(aux_label_npuc,'(i8)') npuc
         WRITE(stdout, "(5x, 'p-DIAGOX done in ',a,' with a SCALAPACK grid ',a)") &
          & TRIM(ADJUSTL(human_readable_time(time_spent(2)-time_spent(1)))), &
          & "("//TRIM(ADJUSTL(aux_label_npur))//"x"//TRIM(ADJUSTL(aux_label_npuc))//")"
         !
      ELSE
         !
         time_spent(1)=get_clock( 'diagox' )
         !
         CALL serial_diagox_zhe( nvec, nbase, nvecx, hr_distr, ew(1:nvec), vr_distr )
         !
         time_spent(2)=get_clock( 'diagox' )
         WRITE(stdout, "(5x, 's-DIAGOX done in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
         !
      ENDIF
      !
      CALL stop_clock( 'diagox' )
      !
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE serial_diagox_dsy ( nselect, n, lda, hr_distr, e, vr_distr ) 
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
      USE mp_global,             ONLY : me_bgrp, root_bgrp, intra_bgrp_comm, inter_image_comm
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
      INTEGER :: i, j, il1, il2, ig1, ig2
      !
      ALLOCATE(zz(n,n))
      ALLOCATE(ee(n))
      !
      ! Create zz from hr_distr
      !
      zz = 0._DP
      DO il2 = 1, pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2>n) CYCLE
         DO ig1 = 1, pert%nglob
            IF(ig1>n) CYCLE
            zz( ig1,ig2 ) = hr_distr( ig1, il2)
         ENDDO
      ENDDO
      CALL mp_sum(zz,inter_image_comm)
      ee = 0._DP
      !
      IF(me_bgrp == root_bgrp ) THEN 
         !
         CALL matdiago_dsy(n,zz,ee,.FALSE.)
         !
      ENDIF
      !
      CALL mp_bcast( ee, root_bgrp, intra_bgrp_comm )
      CALL mp_bcast( zz, root_bgrp, intra_bgrp_comm )
      !
      DO il2 = 1, pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2>nselect) CYCLE
         DO ig1 = 1, pert%nglob
            IF(ig1>n) CYCLE
            vr_distr(ig1,il2) = zz(ig1,ig2)
         ENDDO
      ENDDO
      DO i = 1, nselect
         e(i) = ee(i)
      ENDDO
      !
      DEALLOCATE(zz)
      DEALLOCATE(ee)
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE serial_diagox_zhe ( nselect, n, lda, hr_distr, e, vr_distr ) 
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
      USE mp_global,             ONLY : me_bgrp, root_bgrp, intra_bgrp_comm, inter_image_comm
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
      INTEGER :: i, j, il1, il2, ig1, ig2
      !
      ALLOCATE(zz(n,n))
      ALLOCATE(ee(n))
      !
      ! Create zz from hr_distr
      !
      zz = 0._DP
      DO il2 = 1, pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2>n) CYCLE
         DO ig1 = 1, pert%nglob
            IF(ig1>n) CYCLE
            zz( ig1,ig2 ) = hr_distr( ig1, il2)
         ENDDO
      ENDDO
      CALL mp_sum(zz,inter_image_comm)
      ee = 0._DP
      !
      IF(me_bgrp == root_bgrp ) THEN 
         !
         CALL matdiago_zhe(n,zz,ee,.FALSE.)
         !
      ENDIF
      !
      CALL mp_bcast( ee, root_bgrp, intra_bgrp_comm )
      CALL mp_bcast( zz, root_bgrp, intra_bgrp_comm )
      !
      DO il2 = 1, pert%nloc
         ig2 = pert%l2g(il2)
         IF(ig2>nselect) CYCLE
         DO ig1 = 1, pert%nglob
            IF(ig1>n) CYCLE
            vr_distr(ig1,il2) = zz(ig1,ig2)
         ENDDO
      ENDDO
      DO i = 1, nselect
         e(i) = ee(i)
      ENDDO
      !
      DEALLOCATE(zz)
      DEALLOCATE(ee)
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE build_hr_real( ag, bg, l2_s, l2_e, c_distr, g_s, g_e )
      !------------------------------------------------------------------------
      !
      !  c_distr = < ag | bg >
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq0,npwq0x
      USE gvect,                ONLY : gstart
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwq0x,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: bg(npwq0x,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      REAL(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_s, g_e
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: DDOT
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock ('build_hr')
      !
      ! Initialize to zero
      !
      c_distr(:,l2_s:l2_e)=0.0_DP
      !
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > g_e ) CYCLE
            !
            DO il2=l2_s,l2_e
               !
               !ig2 = pert%l2g(il2)
               !IF( ig2 < n1 .OR. ig2 > n2 ) CYCLE
               !
               IF(gstart==1) THEN
                  c_distr(ig1,il2) = 2.0_DP * DDOT(2*npwq0,ag(1,il1),1,bg(1,il2),1)
               ELSE
                  c_distr(ig1,il2) = 2.0_DP * DDOT(2*npwq0-2,ag(2,il1),1,bg(2,il2),1)
               ENDIF
               !
            ENDDO
         ENDDO
         !
         ! Cycle the ag array 
         ! 
         CALL mp_circular_shift_left( ag,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      ! Syncronize c_distr
      !
      CALL mp_sum(c_distr(:,l2_s:l2_e), intra_bgrp_comm)
      !
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "build_hr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE build_hr_complex( ag, bg, l2_s, l2_e, c_distr, g_s, g_e )
      !------------------------------------------------------------------------
      !
      !  c_distr = < ag | bg >
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq0,npwq0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(INOUT) :: ag(npwq0x,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: bg(npwq0x,pert%nlocx)
      INTEGER,INTENT(IN) :: l2_s,l2_e
      COMPLEX(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      INTEGER,INTENT(IN) :: g_s, g_e
      !
      ! Workspace
      !
      COMPLEX(DP),EXTERNAL :: ZDOTC
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock ('build_hr')
      !
      ! Initialize to zero
      !
      c_distr(:,l2_s:l2_e)=0.0_DP
      !
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > g_e ) CYCLE
            !
            DO il2=l2_s,l2_e
               !
               !ig2 = pert%l2g(il2)
               !IF( ig2 < n1 .OR. ig2 > n2 ) CYCLE
               !
               c_distr(ig1,il2) = ZDOTC(npwq0,ag(1,il1),1,bg(1,il2),1)
               !
            ENDDO
         ENDDO
         !
         ! Cycle the ag array 
         ! 
         CALL mp_circular_shift_left( ag,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      ! Syncronize c_distr
      !
      CALL mp_sum(c_distr(:,l2_s:l2_e), intra_bgrp_comm)
      !
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "build_hr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE symm_hr_distr_real( hr_distr, n, lda )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: n, lda
      REAL(DP),INTENT(INOUT) :: hr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'symm_hr' )
      !
      ALLOCATE( tmp_distr(lda,pert%nlocx) )
      tmp_distr = hr_distr
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nloc
            ig1 = pert%l2g(il1)
            IF( ig1 > n ) CYCLE
            !
            DO il2=1,pert%nlocx
               ig2 = tmp_l2g(il2)
               IF( .NOT.ig2 > ig1 ) CYCLE
               !
               hr_distr(ig2,il1) = tmp_distr(ig1,il2)
               !
            ENDDO
         ENDDO
         !
         ! Cycle the tmp_distr array 
         ! 
         CALL mp_circular_shift_left( tmp_distr ,        icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g   , icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DEALLOCATE( tmp_distr )
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "symm_hr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE symm_hr_distr_complex( hr_distr, n, lda )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: n, lda
      COMPLEX(DP),INTENT(INOUT) :: hr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'symm_hr' )
      !
      ALLOCATE( tmp_distr(lda,pert%nlocx) )
      tmp_distr = hr_distr
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nloc
            ig1 = pert%l2g(il1)
            IF( ig1 > n ) CYCLE
            !
            DO il2=1,pert%nlocx
               ig2 = tmp_l2g(il2)
               IF( .NOT.ig2 > ig1 ) CYCLE
               !
               hr_distr(ig2,il1) = tmp_distr(ig1,il2)
               !
            ENDDO
         ENDDO
         !
         ! Cycle the tmp_distr array 
         ! 
         CALL mp_circular_shift_left( tmp_distr ,        icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g   , icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DEALLOCATE( tmp_distr )
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "symm_hr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE redistribute_vr_distr_real( nselect, n, lda, vr_distr, ishift)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nselect, n, lda
      REAL(DP),INTENT(INOUT) :: vr_distr(lda,pert%nlocx)
      INTEGER,INTENT(IN) :: ishift(lda)
      !
      ! Workspace
      !
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'redistr_vr' )
      !
      ALLOCATE( tmp_distr(lda,pert%nlocx) )
      tmp_distr = vr_distr
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il2=1,pert%nloc
            ig2 = pert%l2g(il2)
            IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
            !
            DO il1=1,pert%nlocx
               ig1 = tmp_l2g(il1)
               !IF( ig1 .NE. ig2-n ) CYCLE
               IF( ig1 == 0 ) CYCLE
               IF( ig1 .NE. ishift(ig2) ) CYCLE
               !
               vr_distr(:,il2) = tmp_distr(:,il1)
               !
            ENDDO
         ENDDO
         !
         ! Cycle the tmp_distr array 
         ! 
         CALL mp_circular_shift_left( tmp_distr ,        icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g   , icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DEALLOCATE( tmp_distr )
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "redistr_vr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE redistribute_vr_distr_complex( nselect, n, lda, vr_distr, ishift)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: nselect, n, lda
      COMPLEX(DP),INTENT(INOUT) :: vr_distr(lda,pert%nlocx)
      INTEGER,INTENT(IN) :: ishift(lda)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'redistr_vr' )
      !
      ALLOCATE( tmp_distr(lda,pert%nlocx) )
      tmp_distr = vr_distr
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il2=1,pert%nloc
            ig2 = pert%l2g(il2)
            IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
            !
            DO il1=1,pert%nlocx
               ig1 = tmp_l2g(il1)
               !IF( ig1 .NE. ig2-n ) CYCLE
               IF( ig1 == 0 ) CYCLE
               IF( ig1 .NE. ishift(ig2) ) CYCLE
               !
               vr_distr(:,il2) = tmp_distr(:,il1)
               !
            ENDDO
         ENDDO
         !
         ! Cycle the tmp_distr array 
         ! 
         CALL mp_circular_shift_left( tmp_distr ,        icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g   , icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DEALLOCATE( tmp_distr )
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "redistr_vr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_real( ag, bg, nselect, n, lda, vr_distr, ew )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq0,npwq0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwq0x,pert%nlocx)
      COMPLEX(DP) :: bg(npwq0x,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      COMPLEX(DP) :: zconst
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'update_vr' )
      !
      ALLOCATE( hg(npwq0x,pert%nlocx) )
      hg = 0._DP 
      !
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > n ) CYCLE
            !
            DO il2=1,pert%nloc
               !
               ig2 = pert%l2g(il2)
               IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE 
               !
               zconst = CMPLX( vr_distr(ig1,il2), 0_DP, KIND=DP )
               CALL ZAXPY(npwq0,zconst,ag(1,il1),1,hg(1,il2),1)
               !dhg(:,il2) = dhg(:,il2) + amat(:,il1) * z(ig1,ig2-n) 
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array
         ! 
         CALL mp_circular_shift_left( ag,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         ag(:,il2) = - ew(ig2) * hg(:,il2)
      ENDDO
      !
      hg = 0._DP
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > n ) CYCLE
            !
            DO il2=1,pert%nloc
               !
               ig2 = pert%l2g(il2)
               IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE 
               !
               zconst = CMPLX( vr_distr(ig1,il2), 0_DP, KIND=DP )
               CALL ZAXPY(npwq0,zconst,bg(1,il1),1,hg(1,il2),1)
               !dhg(:,il2) = dhg(:,il2) + bmat(:,il1) * z(ig1,ig2-n) 
               !
            ENDDO
         ENDDO
         !
         ! Cycle the bmat array 
         ! 
         CALL mp_circular_shift_left( bg,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         ag(:,il2) = ag(:,il2) + hg(:,il2)
      ENDDO
      !
      DEALLOCATE( tmp_l2g )
      DEALLOCATE( hg )
      !
      CALL stop_clock( "update_vr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_complex( ag, bg, nselect, n, lda, vr_distr, ew )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq0,npwq0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwq0x,pert%nlocx)
      COMPLEX(DP) :: bg(npwq0x,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'update_vr' )
      !
      ALLOCATE( hg(npwq0x,pert%nlocx) )
      hg = 0._DP 
      !
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > n ) CYCLE
            !
            DO il2=1,pert%nloc
               !
               ig2 = pert%l2g(il2)
               IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE 
               !
               CALL ZAXPY(npwq0,vr_distr(ig1,il2),ag(1,il1),1,hg(1,il2),1)
               !dhg(:,il2) = dhg(:,il2) + amat(:,il1) * z(ig1,ig2-n) 
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array
         ! 
         CALL mp_circular_shift_left( ag,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         ag(:,il2) = - ew(ig2) * hg(:,il2)
      ENDDO
      !
      hg = 0._DP
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > n ) CYCLE
            !
            DO il2=1,pert%nloc
               !
               ig2 = pert%l2g(il2)
               IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE 
               !
               CALL ZAXPY(npwq0,vr_distr(ig1,il2),bg(1,il1),1,hg(1,il2),1)
               !dhg(:,il2) = dhg(:,il2) + bmat(:,il1) * z(ig1,ig2-n) 
               !
            ENDDO
         ENDDO
         !
         ! Cycle the bmat array 
         ! 
         CALL mp_circular_shift_left( bg,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         ag(:,il2) = ag(:,il2) + hg(:,il2)
      ENDDO
      !
      DEALLOCATE( tmp_l2g )
      DEALLOCATE( hg )
      !
      CALL stop_clock( "update_vr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE refresh_with_vr_distr_real( ag, nselect, n, lda, vr_distr )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq0,npwq0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwq0x,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      COMPLEX(DP) :: zconst
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'refresh_vr' )
      !
      ALLOCATE( hg(npwq0x,pert%nlocx) )
      hg = 0._DP 
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > n ) CYCLE
            !
            DO il2=1,pert%nloc
               !
               ig2 = pert%l2g(il2)
               IF( ig2 > nselect ) CYCLE 
               !
               zconst = CMPLX( vr_distr(ig1,il2), 0_DP, KIND=DP )
               CALL ZAXPY(npwq0,zconst,ag(1,il1),1,hg(1,il2),1)
               !dhg(:,il2) = dhg(:,il2) + amat(:,il1) * z(ig1,ig2) 
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array 
         ! 
         CALL mp_circular_shift_left( ag  ,           icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 > nselect ) THEN
            ag(:,il2) = 0._DP
         ELSE
            ag(:,il2) = hg(:,il2)
         ENDIF
      ENDDO
      !
      DEALLOCATE( hg )
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "refresh_vr" )
      !
    END SUBROUTINE
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE refresh_with_vr_distr_complex( ag, nselect, n, lda, vr_distr )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE westcom,              ONLY : npwq0,npwq0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwq0x,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl
      !
      CALL mp_barrier( world_comm )
      !  
      CALL start_clock( 'refresh_vr' )
      !
      ALLOCATE( hg(npwq0x,pert%nlocx) )
      hg = 0._DP 
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc 
         tmp_l2g(il1) = pert%l2g(il1)
      ENDDO
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            IF( ig1 == 0 .OR. ig1 > n ) CYCLE
            !
            DO il2=1,pert%nloc
               !
               ig2 = pert%l2g(il2)
               IF( ig2 > nselect ) CYCLE 
               !
               CALL ZAXPY(npwq0,vr_distr(ig1,il2),ag(1,il1),1,hg(1,il2),1)
               !dhg(:,il2) = dhg(:,il2) + amat(:,il1) * z(ig1,ig2) 
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array 
         ! 
         CALL mp_circular_shift_left( ag  ,           icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 > nselect ) THEN
            ag(:,il2) = 0._DP
         ELSE
            ag(:,il2) = hg(:,il2)
         ENDIF
      ENDDO
      !
      DEALLOCATE( hg )
      DEALLOCATE( tmp_l2g )
      !
      CALL stop_clock( "refresh_vr" )
      !
    END SUBROUTINE
    !
END MODULE
