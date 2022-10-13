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
MODULE wbse_tools
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE parallel_include
  USE mp,        ONLY : mp_stop
  !
  IMPLICIT NONE
  !
  SAVE
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
     MODULE PROCEDURE preconditioner_complex_2nd
  END INTERFACE
  !
  CONTAINS

      SUBROUTINE mp_circular_shift_left_c4d( buf, itag, gid )
       IMPLICIT NONE
       COMPLEX(DP) :: buf( :, :, :, : )
       INTEGER, INTENT(IN) :: itag
       INTEGER, INTENT(IN) :: gid
       INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)

       INTEGER :: istatus( mpi_status_size )
       !
       group = gid
       !
       CALL mpi_comm_size( group, npe, ierr )
       IF (ierr/=0) CALL mp_stop( 8100 )
       CALL mpi_comm_rank( group, mype, ierr )
       IF (ierr/=0) CALL mp_stop( 8101 )
       !
       sour = mype + 1
       IF( sour == npe ) sour = 0
       dest = mype - 1
       IF( dest == -1 ) dest = npe - 1
       !
       CALL MPI_Sendrecv_replace( buf, SIZE(buf), MPI_DOUBLE_COMPLEX, &
            dest, itag, sour, itag, group, istatus, ierr)
       !
       IF (ierr/=0) CALL mp_stop( 8102 )
       !
#else
       ! do nothing
#endif
       RETURN
    END SUBROUTINE
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
      USE pwcom,                ONLY : nks,npwx,npw
      USE westcom,              ONLY : nbnd_occ, nbndval0x
      USE gvect,                ONLY : gstart
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN)       :: l2_s,l2_e
      INTEGER, INTENT(IN)       :: g_s, g_e
      COMPLEX(DP),INTENT(IN)    :: bg(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      REAL(DP),INTENT(INOUT)    :: c_distr(pert%nglob,pert%nlocx)
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: DDOT
      REAL(DP):: norm
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl, ibnd, iks, nbndval
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
               norm = 0.0_DP
               !
               DO iks=1, nks
                  !
                  nbndval = nbnd_occ(iks)
                  !
                  DO ibnd=1, nbndval
                     !
                     norm = norm + 2.0_DP * DDOT(2*npw,ag(1,ibnd,iks,il1),1,bg(1,ibnd,iks,il2),1)
                     !
                     IF(gstart==2) norm = norm -  REAL(ag(1,ibnd,iks,il1),KIND=DP)*REAL(bg(1,ibnd,iks,il2),KIND=DP)
                     !
                  ENDDO
                  !
               ENDDO
               !
               c_distr(ig1,il2) = norm
               !
            ENDDO
            !
         ENDDO
         !
         ! Cycle the ag array
         !
         CALL mp_circular_shift_left_c4d( ag,      icycl,        inter_image_comm)
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
      USE pwcom,                ONLY : nks,npwx,npw
      USE westcom,              ONLY : nbnd_occ, nbndval0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN)        :: l2_s,l2_e
      INTEGER,INTENT(IN)        :: g_s, g_e
      COMPLEX(DP),INTENT(IN)    :: bg(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),EXTERNAL :: ZDOTC
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl, ibnd, iks, nbndval
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
               DO iks  = 1, nks
                  !
                  nbndval = nbnd_occ(iks)
                  !
                  DO ibnd = 1, nbndval
                     !
                     c_distr(ig1,il2) = c_distr(ig1,il2) + ZDOTC(npw,ag(1,ibnd,iks,il1),1,bg(1,ibnd,iks,il2),1)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
         ENDDO
         !
         ! Cycle the ag array
         !
         CALL mp_circular_shift_left_c4d( ag,      icycl,        inter_image_comm)
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
    SUBROUTINE update_with_vr_distr_real( ag, bg, nselect, n, lda, vr_distr, ew )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx,npw
      USE westcom,              ONLY : nbnd_occ, nbndval0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,    INTENT(IN)    :: nselect, n, lda
      REAL(DP),   INTENT(IN)    :: vr_distr(lda,pert%nlocx)
      REAL(DP),   INTENT(IN)    :: ew(lda)
      COMPLEX(DP),INTENT(INOUT) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP),INTENT(INOUT) :: bg(npwx,nbndval0x,nks,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      INTEGER,    ALLOCATABLE :: tmp_l2g(:)
      INTEGER     :: il1, il2, ig1, ig2, icycl, ibnd, iks, nbndval
      COMPLEX(DP) :: zconst
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock( 'update_vr' )
      !
      ALLOCATE( hg(npwx,nbndval0x,nks,pert%nlocx) )
      hg = 0._DP
      !
      ALLOCATE( tmp_l2g(1:pert%nlocx) )
      !
      tmp_l2g = 0
      !
      DO il1 = 1, pert%nloc
         !
         tmp_l2g(il1) = pert%l2g(il1)
         !
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
               !
               DO iks  = 1,nks
                  !
                  nbndval = nbnd_occ(iks)
                  !
                  DO ibnd = 1, nbndval
                     !
                     CALL ZAXPY(npw,zconst,ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array
         !
         CALL mp_circular_shift_left_c4d( ag,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         !
         ig2 = pert%l2g(il2)
         !
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         !
         DO iks  = 1, nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd = 1, nbndval
               !
               ag(:,ibnd,iks,il2) = - ew(ig2) * hg(:,ibnd,iks,il2)
               !
            ENDDO
            !
         ENDDO
         !
      ENDDO
      !
      hg = 0._DP
      !
      DO icycl=0,nimage-1
         !
         DO il1=1,pert%nlocx
            !
            ig1 = tmp_l2g(il1)
            !
            IF( ig1 == 0 .OR. ig1 > n ) CYCLE
            !
            DO il2=1,pert%nloc
               !
               ig2 = pert%l2g(il2)
               !
               IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
               !
               zconst = CMPLX( vr_distr(ig1,il2), 0_DP, KIND=DP )
               !
               DO iks  = 1, nks
                  !
                  nbndval = nbnd_occ(iks)
                  !
                  DO ibnd = 1, nbndval
                     !
                     CALL ZAXPY(npw,zconst,bg(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
         ENDDO
         !
         ! Cycle the bmat array
         !
         CALL mp_circular_shift_left_c4d( bg,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         !
         DO iks  = 1, nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd = 1, nbndval
               !
               ag(:,ibnd,iks,il2) = ag(:,ibnd,iks,il2) + hg(:,ibnd,iks,il2)
               !
            ENDDO
            !
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE( tmp_l2g )
      DEALLOCATE( hg )
      !
      CALL stop_clock( "update_vr" )
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_complex( ag, bg, nselect, n, lda, vr_distr, ew )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx,npw
      USE westcom,              ONLY : nbnd_occ, nbndval0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP) :: bg(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      INTEGER,    ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl, ibnd, iks, nbndval
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock( 'update_vr' )
      !
      ALLOCATE( hg(npwx,nbndval0x,nks,pert%nlocx) )
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
               DO iks  = 1, nks
                  !
                  nbndval = nbnd_occ(iks)
                  !
                  DO ibnd = 1, nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array
         !
         CALL mp_circular_shift_left_c4d( ag,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         !
         DO iks =1, nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd=1, nbndval
               ag(:,ibnd,iks,il2) = - ew(ig2) * hg(:,ibnd,iks,il2)
            ENDDO
         ENDDO
         !
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
               DO iks  = 1, nks
                  nbndval = nbnd_occ(iks)
                  DO ibnd = 1, nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),bg(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle the bmat array
         !
         CALL mp_circular_shift_left_c4d( bg,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         !
         DO iks  = 1, nks
            nbndval = nbnd_occ(iks)
            DO ibnd = 1, nbndval
               ag(:,ibnd,iks,il2) = ag(:,ibnd,iks,il2) + hg(:,ibnd,iks,il2)
            ENDDO
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE( tmp_l2g )
      DEALLOCATE( hg )
      !
      CALL stop_clock( "update_vr" )
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_with_vr_distr_complex_2nd( ag, bg, nselect, n, lda, vr_distr, ew, eps_ref )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx,npw
      USE westcom,              ONLY : nbnd_occ,nbndval0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      COMPLEX(DP) :: bg(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      REAL(DP),INTENT(IN) :: ew(lda), eps_ref
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      INTEGER,    ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl, ibnd, iks, nbndval
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock( 'update_vr' )
      !
      ALLOCATE( hg(npwx,nbndval0x,nks,pert%nlocx) )
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
               DO iks  = 1, nks
                  !
                  nbndval = nbnd_occ(iks)
                  !
                  DO ibnd = 1, nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array
         !
         CALL mp_circular_shift_left_c4d( ag,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         !
         DO iks =1, nks
            !
            nbndval = nbnd_occ(iks)
            !
            DO ibnd=1, nbndval
               ag(:,ibnd,iks,il2) = - (ew(ig2)-eps_ref)**2 * hg(:,ibnd,iks,il2)
            ENDDO
         ENDDO
         !
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
               DO iks  = 1, nks
                  nbndval = nbnd_occ(iks)
                  DO ibnd = 1, nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),bg(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle the bmat array
         !
         CALL mp_circular_shift_left_c4d( bg,      icycl,        inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
         !
         DO iks  = 1, nks
            nbndval = nbnd_occ(iks)
            DO ibnd = 1, nbndval
               ag(:,ibnd,iks,il2) = ag(:,ibnd,iks,il2) + hg(:,ibnd,iks,il2)
            ENDDO
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE( tmp_l2g )
      DEALLOCATE( hg )
      !
      CALL stop_clock( "update_vr" )
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE refresh_with_vr_distr_real( ag, nselect, n, lda, vr_distr )
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
      USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
      USE distribution_center,  ONLY : pert
      USE pwcom,                ONLY : nks,npwx,npw
      USE westcom,              ONLY : nbnd_occ,nbndval0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl, ibnd, iks, nbndval
      COMPLEX(DP) :: zconst
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock( 'refresh_vr' )
      !
      ALLOCATE( hg(npwx,nbndval0x,nks,pert%nlocx) )
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
               !
               DO iks  = 1, nks
                  nbndval = nbnd_occ(iks)
                  DO ibnd = 1, nbndval
                     CALL ZAXPY(npw,zconst,ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array
         !
         CALL mp_circular_shift_left_c4d( ag  ,           icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 > nselect ) THEN
           !
           DO iks  = 1, nks
              nbndval = nbnd_occ(iks)
              DO ibnd = 1, nbndval
                 ag(:,ibnd,iks,il2) = 0._DP
              ENDDO
           ENDDO
           !
         ELSE
           !
           DO iks  = 1, nks
              nbndval = nbnd_occ(iks)
              DO ibnd = 1, nbndval
                 ag(:,ibnd,iks,il2) = hg(:,ibnd,iks,il2)
              ENDDO
           ENDDO
           !
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
      USE pwcom,                ONLY : nks,npwx,npw
      USE westcom,              ONLY : nbnd_occ,nbndval0x
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN) :: nselect, n, lda
      COMPLEX(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: hg(:,:,:,:)
      INTEGER,ALLOCATABLE :: tmp_l2g(:)
      INTEGER :: il1, il2, ig1, ig2, icycl, ibnd, iks, nbndval
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock( 'refresh_vr' )
      !
      ALLOCATE( hg(npwx,nbndval0x,nks,pert%nlocx) )
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
               DO iks  = 1, nks
                  nbndval = nbnd_occ(iks)
                  DO ibnd = 1, nbndval
                     CALL ZAXPY(npw,vr_distr(ig1,il2),ag(1,ibnd,iks,il1),1,hg(1,ibnd,iks,il2),1)
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
         !
         ! Cycle the amat array
         !
         CALL mp_circular_shift_left_c4d( ag  ,           icycl, inter_image_comm)
         CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
         !
      ENDDO
      !
      DO il2=1,pert%nloc
         ig2 = pert%l2g(il2)
         IF( ig2 > nselect ) THEN
           !
           DO iks  = 1, nks
              nbndval = nbnd_occ(iks)
              DO ibnd = 1, nbndval
                 ag(:,ibnd,iks,il2) = 0._DP
              ENDDO
           ENDDO
           !
         ELSE
           !
           DO iks  = 1, nks
              nbndval = nbnd_occ(iks)
              DO ibnd = 1, nbndval
                 ag(:,ibnd,iks,il2) = hg(:,ibnd,iks,il2)
              ENDDO
           ENDDO
           !
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
    SUBROUTINE preconditioner_complex( ag, nselect, n, lda, ew, turn_shift )
      !
      USE kinds,                ONLY : dp
      USE mp_global,            ONLY : my_image_id,inter_image_comm,world_comm
      USE mp,                   ONLY : mp_bcast,mp_barrier,mp_max
      USE wavefunctions,        ONLY : evc
      USE distribution_center,  ONLY : pert
      USE buffers,              ONLY : get_buffer
      USE pwcom,                ONLY : nks,npwx
      USE westcom,              ONLY : nbnd_occ,lrwfc,iuwfc, nbndval0x
      USE wvfct,                ONLY : g2kin,et
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN)  :: nselect, n, lda
      REAL(DP),INTENT(IN) :: ew(lda)
      LOGICAL, INTENT(IN) :: turn_shift
      !
      ! Workspace
      !
      INTEGER :: il1, ig1, ig, ibnd, nbndval
      INTEGER :: iks, current_k, current_spin
      INTEGER :: mloc, mstart, max_mloc
      REAL(DP):: temp, minimum
      REAL(DP),ALLOCATABLE :: eprec(:)
      !
      minimum=0.001d0
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock( 'precd_ag' )
      !
      mloc = 0
      mstart = 1
      DO il1 = 1, pert%nloc
         ig1 = pert%l2g(il1)
         IF( ig1 <= n .OR. ig1 > n+nselect ) CYCLE
         IF( mloc==0 ) mstart = il1
         mloc = mloc + 1
      ENDDO
      !
      ! Apply Liouville operator
      !
      max_mloc = mloc
      CALL mp_max (max_mloc, inter_image_comm)
      !
      DO il1 = mstart, mstart+max_mloc-1
         !
         ig1 = pert%l2g(il1)
         !
         DO iks  = 1, nks
            !
            nbndval = nbnd_occ(iks)
            !
            ! ... Set k-point, spin, kinetic energy, needed by Hpsi
            !
            current_k = iks
            !
            CALL g2_kin( iks )
            !
            ! ... read in GS wavefunctions iks
            !
            IF (nks>1) THEN
               !
               IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
               CALL mp_bcast(evc,0,inter_image_comm)
               !
            ENDIF
            !
            IF (.NOT.( ig1 <= n .OR. ig1 > n+nselect )) THEN
               !
               DO ibnd = 1, nbndval
                  !
                  DO ig = 1, npwx
                     !
                     IF (turn_shift) THEN
                        temp = (g2kin(ig) - et(ibnd,iks))
                     ELSE
                        temp =  g2kin(ig)
                     ENDIF
                     !
                     IF(ABS(temp) < minimum ) temp = SIGN(minimum,temp)
                     ag(ig,ibnd,iks,il1) = ag(ig,ibnd,iks,il1)/temp
                     !
                  ENDDO
                  !
               ENDDO
               !
               ! Pc amat
               !
!               CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,ag(:,:,iks,il1),(1.0_DP,0.0_DP))
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      !
      CALL stop_clock( "precd_ag" )
      !
    END SUBROUTINE
    !
    SUBROUTINE preconditioner_complex_2nd( ag, nselect, n, lda, ew, turn_shift, eps_ref)
      !
      USE kinds,                ONLY : dp
      USE mp_global,            ONLY : my_image_id,inter_image_comm,world_comm
      USE mp,                   ONLY : mp_bcast,mp_barrier,mp_max
      USE wavefunctions,        ONLY : evc
      USE distribution_center,  ONLY : pert
      USE buffers,              ONLY : get_buffer
      USE pwcom,                ONLY : nks,npwx
      USE westcom,              ONLY : nbnd_occ,lrwfc,iuwfc,nbndval0x
      USE wvfct,                ONLY : g2kin,et
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP) :: ag(npwx,nbndval0x,nks,pert%nlocx)
      INTEGER,INTENT(IN)  :: nselect, n, lda
      REAL(DP),INTENT(IN) :: ew(lda), eps_ref
      LOGICAL, INTENT(IN) :: turn_shift
      !
      ! Workspace
      !
      INTEGER :: il1, ig1, ig, ibnd, nbndval
      INTEGER :: iks, current_k, current_spin
      INTEGER :: mloc, mstart, max_mloc
      REAL(DP):: temp, minimum
      REAL(DP),ALLOCATABLE :: eprec(:)
      !
      minimum=0.001d0
      !
      CALL mp_barrier( world_comm )
      !
      CALL start_clock( 'precd_ag' )
      !
      mloc = 0
      mstart = 1
      DO il1 = 1, pert%nloc
         ig1 = pert%l2g(il1)
         IF( ig1 <= n .OR. ig1 > n+nselect ) CYCLE
         IF( mloc==0 ) mstart = il1
         mloc = mloc + 1
      ENDDO
      !
      ! Apply Liouville operator
      !
      max_mloc = mloc
      CALL mp_max (max_mloc, inter_image_comm)
      !
      DO il1 = mstart, mstart+max_mloc-1
         !
         ig1 = pert%l2g(il1)
         !
         DO iks  = 1, nks
            !
            nbndval = nbnd_occ(iks)
            !
            ! ... Set k-point, spin, kinetic energy, needed by Hpsi
            !
            current_k = iks
            !
            CALL g2_kin( iks )
            !
            ! ... read in GS wavefunctions iks
            !
            IF (nks>1) THEN
               !
               IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
               CALL mp_bcast(evc,0,inter_image_comm)
               !
            ENDIF
            !
            IF (.NOT.( ig1 <= n .OR. ig1 > n+nselect )) THEN
               !
               DO ibnd = 1, nbndval
                  !
                  DO ig = 1, npwx
                     !
                     IF (turn_shift) THEN
                        temp = (g2kin(ig)-eps_ref)**2 - et(ibnd,iks)
                     ELSE
                        temp = (g2kin(ig)-eps_ref)**2
                     ENDIF
                     !
                     IF(ABS(temp) < minimum**2 ) temp = SIGN(minimum,temp)
                     ag(ig,ibnd,iks,il1) = ag(ig,ibnd,iks,il1)/SQRT(temp)
                     !
                  ENDDO
                  !
               ENDDO
               !
               ! Pc amat
               !
               CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,ag(:,:,iks,il1),(1.0_DP,0.0_DP))
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      !
      CALL stop_clock( "precd_ag" )
      !
    END SUBROUTINE
END MODULE
