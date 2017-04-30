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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
#if defined __SCALAPACK
!----------------------------------------------------------------------------
SUBROUTINE parallel_distributed_diago_dsy ( glob_nselect,glob_ndim,glob_ndimx,a_distr,v_distr,e,NPROW,NPCOL,idistr ) 
  !----------------------------------------------------------------------------
  !
  ! Diagox -- parallel
  !   global_nselect : number of wanted ev
  !   global_ndim    : actual dimension of a
  !   global_ndimx   : leading dimension of a
  !   a_distr        : matrix to be diago
  !   v_distr        : unitary trans. 
  !   e              : eigenval 
  !
  USE kinds,                 ONLY : DP
  USE mp,                    ONLY : mp_bcast, mp_sum, mp_circular_shift_left
  USE mp_global,             ONLY : me_bgrp, root_bgrp, intra_bgrp_comm, nproc_bgrp, inter_image_comm, nimage
  USE mp_world,              ONLY : nproc, world_comm, mpime
  USE class_idistribute,     ONLY : idistribute
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute),INTENT(IN) :: idistr
  INTEGER,INTENT(IN) :: glob_nselect, glob_ndim, glob_ndimx 
  REAL(DP),INTENT(IN) :: a_distr(glob_ndimx,idistr%nlocx)
  REAL(DP),INTENT(OUT) :: v_distr(glob_ndimx,idistr%nlocx)
  REAL(DP),INTENT(OUT) :: e(glob_nselect)
  INTEGER,INTENT(IN) :: NPROW,NPCOL
  !
  ! Workspace
  !
  INTEGER :: il1, ig1
  INTEGER :: il2, ig2
  INTEGER :: i_loc, j_loc, i_glob, j_glob
  INTEGER :: icycl
  INTEGER,ALLOCATABLE :: tmp_igl(:)
  !
  ! Workspace for SCALAPACK
  !
  REAL(DP),ALLOCATABLE :: la(:,:),lv(:,:)
  REAL(DP) :: ge(glob_ndim)
  INTEGER :: NB, MYROW, MYCOL, LDROW, LDCOL
  INTEGER :: BLACS_CONTEXT, INFO
  INTEGER :: DESCA(9),DESCV(9)
  INTEGER :: ifail,LWORK
  REAL(DP),ALLOCATABLE :: work(:)
  INTEGER,EXTERNAL :: INDXL2G
  INTEGER,EXTERNAL :: NUMROC
  !
  ! ===============================
  ! 1) Definisco una griglia locale
  ! =============================== 
  !
  NB= MAX(glob_ndim / NPCOL,1)
  !
  CALL BLACS_GET( -1, 0, BLACS_CONTEXT )
  CALL BLACS_GRIDINIT( BLACS_CONTEXT, 'R', NPROW, NPCOL )
  CALL BLACS_GRIDINFO( BLACS_CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
  !
  IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
     !
     LDROW = NUMROC(glob_ndim, NB, MYROW, 0, NPROW)   ! number of elements in this element of the grid 
     LDCOL = NUMROC(glob_ndim, NB, MYCOL, 0, NPCOL)   ! number of elements in this element of the grid
     CALL DESCINIT( DESCA, glob_ndim, glob_ndim, NB, NB, 0, 0, BLACS_CONTEXT, LDROW, INFO )
     CALL DESCINIT( DESCV, glob_ndim, glob_ndim, NB, NB, 0, 0, BLACS_CONTEXT, LDROW, INFO )
     !
     ALLOCATE( la( LDROW, LDCOL), lv( LDROW, LDCOL) )
     !
  ENDIF
  !
  ! ========================
  ! 2) Ripartisco la matrice
  ! ========================
  !
  v_distr = a_distr
  ALLOCATE( tmp_igl(1:idistr%nlocx) ) 
  tmp_igl = 0
  DO il1 = 1, idistr%nloc 
     tmp_igl(il1) = idistr%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
        !
        ! Filling of the local array distributed in the processor grid
        !
        DO il2 = 1,idistr%nloc
           ig2 = tmp_igl(il2)
           IF(ig2<1) CYCLE 
           ! 
           DO j_loc=1,LDCOL
              j_glob=INDXL2G(j_loc, NB, MYCOL, 0, NPCOL)
              IF(j_glob.NE.ig2) CYCLE 
              DO i_loc=1,LDROW
                 i_glob=INDXL2G(i_loc, NB, MYROW, 0, NPROW)
                 la(i_loc,j_loc) = v_distr(i_glob,il2)
              ENDDO
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
     ! Cycle the v_distr array 
     ! 
     CALL mp_circular_shift_left( v_distr,        icycl, inter_image_comm)
     CALL mp_circular_shift_left( tmp_igl, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( tmp_igl )
  !
  ! ===============
  ! 3) DIAGONALIZZO
  ! ===============
  !
  ge = 0._DP
  IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
     !
     LWORK=-1 
     ALLOCATE( work(1) )
     CALL PDSYEV('V','U',glob_ndim,la,1,1,DESCA,ge,lv,1,1,DESCV,work,LWORK,IFAIL )
     ! 
     LWORK=INT(work(1))+1
     DEALLOCATE( work )
     !
     ALLOCATE( work(LWORK) )
     CALL PDSYEV('V','U',glob_ndim,la,1,1,DESCA,ge,lv,1,1,DESCV,work,LWORK,IFAIL )
     !
     DEALLOCATE( work, la )
     !
  ENDIF
  !
  ! =========================
  ! 4) DISTRIBUISCO RISULTATO
  ! =========================
  !
  ! Distribuisco eigenvalues
  !
  DO i_glob = 1, glob_nselect
     e(i_glob) = ge(i_glob)
  ENDDO
  CALL mp_bcast(e,0,world_comm)
  !
  ! Distribuisco eigenvectors
  !
  v_distr = 0._DP
  ALLOCATE( tmp_igl(1:idistr%nlocx) ) 
  tmp_igl = 0
  DO il1 = 1, idistr%nloc 
     tmp_igl(il1) = idistr%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
        !
        ! Filling of the local array distributed in the processor grid
        !
        DO il2 = 1,idistr%nloc
           ig2 = tmp_igl(il2) 
           IF(ig2<1) CYCLE 
           ! 
           DO j_loc=1,LDCOL
              j_glob=INDXL2G(j_loc, NB, MYCOL, 0, NPCOL)
              IF(j_glob.NE.ig2) CYCLE 
              DO i_loc=1,LDROW
                 i_glob=INDXL2G(i_loc, NB, MYROW, 0, NPROW)
                 v_distr(i_glob,il2) = lv(i_loc,j_loc)
              ENDDO
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
     ! Cycle the v_distr array 
     ! 
     CALL mp_circular_shift_left( v_distr,        icycl, inter_image_comm)
     CALL mp_circular_shift_left( tmp_igl, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( tmp_igl )
  !
  IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
     DEALLOCATE( lv )
     CALL BLACS_GRIDEXIT(BLACS_CONTEXT)
  ENDIF
  !
  CALL mp_sum(v_distr,intra_bgrp_comm)
  !
END SUBROUTINE
!
!----------------------------------------------------------------------------
SUBROUTINE parallel_distributed_diago_zhe ( glob_nselect,glob_ndim,glob_ndimx,a_distr,v_distr,e,NPROW,NPCOL,idistr ) 
  !----------------------------------------------------------------------------
  !
  ! Diagox -- parallel
  !   global_nselect : number of wanted ev
  !   global_ndim    : actual dimension of a
  !   global_ndimx   : leading dimension of a
  !   a_distr        : matrix to be diago
  !   v_distr        : unitary trans. 
  !   e              : eigenval 
  !
  USE kinds,                 ONLY : DP
  USE mp,                    ONLY : mp_bcast, mp_sum, mp_circular_shift_left
  USE mp_global,             ONLY : me_bgrp, root_bgrp, intra_bgrp_comm, nproc_bgrp, inter_image_comm, nimage
  USE mp_world,              ONLY : nproc, world_comm, mpime
  USE class_idistribute,     ONLY : idistribute
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute),INTENT(IN) :: idistr
  INTEGER,INTENT(IN) :: glob_nselect, glob_ndim, glob_ndimx
  COMPLEX(DP),INTENT(IN) :: a_distr(glob_ndimx,idistr%nlocx)
  COMPLEX(DP),INTENT(OUT) :: v_distr(glob_ndimx,idistr%nlocx)
  REAL(DP),INTENT(OUT) :: e(glob_nselect)
  INTEGER,INTENT(IN) :: NPROW,NPCOL
  !
  ! Workspace
  !
  INTEGER :: il1, ig1
  INTEGER :: il2, ig2
  INTEGER :: i_loc, j_loc, i_glob, j_glob
  INTEGER :: icycl
  INTEGER,ALLOCATABLE :: tmp_igl(:)
  !
  ! Workspace for SCALAPACK
  !
  COMPLEX(DP),ALLOCATABLE :: la(:,:),lv(:,:)
  REAL(DP) :: ge(glob_ndim)
  INTEGER :: NB, MYROW, MYCOL, LDROW, LDCOL
  INTEGER :: BLACS_CONTEXT, INFO
  INTEGER :: DESCA(9),DESCV(9)
  INTEGER :: ifail,LWORK,LRWORK
  COMPLEX(DP),ALLOCATABLE :: rwork(:)
  COMPLEX(DP),ALLOCATABLE :: work(:)
  INTEGER,EXTERNAL :: INDXL2G
  INTEGER,EXTERNAL :: NUMROC
  !
  ! ===============================
  ! 1) Definisco una griglia locale
  ! =============================== 
  !
  NB= MAX(glob_ndim / NPCOL,1)
  !
  CALL BLACS_GET( -1, 0, BLACS_CONTEXT )
  CALL BLACS_GRIDINIT( BLACS_CONTEXT, 'R', NPROW, NPCOL )
  CALL BLACS_GRIDINFO( BLACS_CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
  !
  IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
     !
     LDROW = NUMROC(glob_ndim, NB, MYROW, 0, NPROW)   ! number of elements in this element of the grid 
     LDCOL = NUMROC(glob_ndim, NB, MYCOL, 0, NPCOL)   ! number of elements in this element of the grid
     CALL DESCINIT( DESCA, glob_ndim, glob_ndim, NB, NB, 0, 0, BLACS_CONTEXT, LDROW, INFO )
     CALL DESCINIT( DESCV, glob_ndim, glob_ndim, NB, NB, 0, 0, BLACS_CONTEXT, LDROW, INFO )
     !
     ALLOCATE( la( LDROW, LDCOL), lv( LDROW, LDCOL) )
     !
  ENDIF
  !
  ! ========================
  ! 2) Ripartisco la matrice
  ! ========================
  !
  v_distr = a_distr
  ALLOCATE( tmp_igl(1:idistr%nlocx) ) 
  tmp_igl = 0
  DO il1 = 1, idistr%nloc 
     tmp_igl(il1) = idistr%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
        !
        ! Filling of the local array distributed in the processor grid
        !
        DO il2 = 1,idistr%nloc
           ig2 = tmp_igl(il2)
           IF(ig2<1) CYCLE 
           ! 
           DO j_loc=1,LDCOL
              j_glob=INDXL2G(j_loc, NB, MYCOL, 0, NPCOL)
              IF(j_glob.NE.ig2) CYCLE 
              DO i_loc=1,LDROW
                 i_glob=INDXL2G(i_loc, NB, MYROW, 0, NPROW)
                 la(i_loc,j_loc) = v_distr(i_glob,il2)
              ENDDO
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
     ! Cycle the v_distr array 
     ! 
     CALL mp_circular_shift_left( v_distr,        icycl, inter_image_comm)
     CALL mp_circular_shift_left( tmp_igl, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( tmp_igl )
  !
  ! ===============
  ! 3) DIAGONALIZZO
  ! ===============
  !
  ge = 0._DP
  IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
     !
     LWORK=-1 
     LRWORK=-1 
     ALLOCATE( work(1) )
     ALLOCATE( rwork(1) )
     CALL PZHEEV('V','U',glob_ndim,la,1,1,DESCA,ge,lv,1,1,DESCV,work,LWORK,rwork,LRWORK,IFAIL )
     ! 
     LWORK=INT(DBLE(work(1)))+1
     DEALLOCATE( work )
     LRWORK=INT(DBLE(rwork(1)))+1
     DEALLOCATE( rwork )
     !
     ALLOCATE( work(LWORK) )
     ALLOCATE( rwork(LRWORK) )
     CALL PZHEEV('V','U',glob_ndim,la,1,1,DESCA,ge,lv,1,1,DESCV,work,LWORK,rwork,LRWORK,IFAIL )
     !
     DEALLOCATE( work, rwork, la )
     !
  ENDIF
  !
  ! =========================
  ! 4) DISTRIBUISCO RISULTATO
  ! =========================
  !
  ! Distribuisco eigenvalues
  !
  DO i_glob = 1, glob_nselect
     e(i_glob) = ge(i_glob)
  ENDDO
  CALL mp_bcast(e,0,world_comm)
  !
  ! Distribuisco eigenvectors
  !
  v_distr = 0._DP
  ALLOCATE( tmp_igl(1:idistr%nlocx) ) 
  tmp_igl = 0
  DO il1 = 1, idistr%nloc 
     tmp_igl(il1) = idistr%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
        !
        ! Filling of the local array distributed in the processor grid
        !
        DO il2 = 1,idistr%nloc
           ig2 = tmp_igl(il2) 
           IF(ig2<1) CYCLE 
           ! 
           DO j_loc=1,LDCOL
              j_glob=INDXL2G(j_loc, NB, MYCOL, 0, NPCOL)
              IF(j_glob.NE.ig2) CYCLE 
              DO i_loc=1,LDROW
                 i_glob=INDXL2G(i_loc, NB, MYROW, 0, NPROW)
                 v_distr(i_glob,il2) = lv(i_loc,j_loc)
              ENDDO
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
     ! Cycle the v_distr array 
     ! 
     CALL mp_circular_shift_left( v_distr,        icycl, inter_image_comm)
     CALL mp_circular_shift_left( tmp_igl, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( tmp_igl )
  !
  IF( MYROW/=-1 .OR. MYCOL/=-1) THEN 
     DEALLOCATE( lv )
     CALL BLACS_GRIDEXIT(BLACS_CONTEXT)
  ENDIF
  !
  CALL mp_sum(v_distr,intra_bgrp_comm)
  !
END SUBROUTINE
#endif
