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
#define ZERO (0.0_DP,0.0_DP)
#define HALF (0.5_DP,0.0_DP)
#define ONE  (1.0_DP,0.0_DP)
!
#if defined(__SCALAPACK)
!----------------------------------------------------------------------------
SUBROUTINE parallel_distributed_diago_dsy(glob_nselect,glob_ndim,glob_ndimx,a_distr,v_distr,e,NPROW,NPCOL,idistr)
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
  USE mp,                    ONLY : mp_bcast
  USE mp_global,             ONLY : intra_bgrp_comm,inter_bgrp_comm,my_bgrp_id,me_bgrp,nbgrp,nproc_bgrp
  USE mp_world,              ONLY : world_comm,mpime,nproc
  USE class_idistribute,     ONLY : idistribute
  USE sort_tools,            ONLY : heapsort,i8b
  USE parallel_include
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute),INTENT(IN) :: idistr
  INTEGER,INTENT(IN) :: glob_nselect,glob_ndim,glob_ndimx
  REAL(DP),INTENT(IN) :: a_distr(glob_ndimx,idistr%nlocx)
  REAL(DP),INTENT(OUT) :: v_distr(glob_ndimx,idistr%nlocx)
  REAL(DP),INTENT(OUT) :: e(glob_nselect)
  INTEGER,INTENT(INOUT) :: NPROW,NPCOL
  !
  ! Workspace
  !
  INTEGER :: i_loc,j_loc,i_glob,j_glob
  INTEGER :: n_j_loc,i_val,this_dest,this_sour,dummy
  INTEGER :: i_proc,j_proc
  INTEGER,ALLOCATABLE :: dest(:)
  INTEGER,ALLOCATABLE :: swap(:)
  INTEGER,ALLOCATABLE :: send_count(:)
  INTEGER,ALLOCATABLE :: recv_count(:)
  INTEGER,ALLOCATABLE :: send_displ(:)
  INTEGER,ALLOCATABLE :: recv_displ(:)
  INTEGER(i8b) :: this_idx
  INTEGER(i8b),ALLOCATABLE :: tmp_i8(:)
  INTEGER(i8b),ALLOCATABLE :: idx_send(:)
  INTEGER(i8b),ALLOCATABLE :: idx_recv(:)
  REAL(DP),ALLOCATABLE :: tmp_r(:)
  REAL(DP),ALLOCATABLE :: val_send(:)
  REAL(DP),ALLOCATABLE :: val_recv(:)
  !
  ! Workspace for ScaLAPACK
  !
  REAL(DP),ALLOCATABLE :: la(:,:),lv(:,:)
  REAL(DP) :: ge(glob_ndim)
  INTEGER :: NB,MYROW,MYCOL,LDROW,LDCOL
  INTEGER :: BLACS_CONTEXT
  INTEGER :: DESC(9)
  INTEGER :: LWORK
  INTEGER :: ierr
  INTEGER :: color
  REAL(DP),ALLOCATABLE :: work(:)
  INTEGER,EXTERNAL :: INDXL2G
  INTEGER,EXTERNAL :: INDXG2P
  INTEGER,EXTERNAL :: NUMROC
  !
  ! ========================
  ! 1) Define processor grid
  ! ========================
  !
  ! Find suitable block size (between 1 and 64)
  !
  NB = 1
  DO WHILE(2*NB*NPCOL <= glob_ndim)
     NB = NB*2
     IF(NB == 64) EXIT
  ENDDO
  !
  i_proc = NPROW
  j_proc = NPCOL
  !
  CALL BLACS_GET(-1,0,BLACS_CONTEXT)
  CALL BLACS_GRIDINIT(BLACS_CONTEXT,'R',NPROW,NPCOL)
  CALL BLACS_GRIDINFO(BLACS_CONTEXT,NPROW,NPCOL,MYROW,MYCOL)
  !
  NPROW = i_proc
  NPCOL = j_proc
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     LDROW = NUMROC(glob_ndim,NB,MYROW,0,NPROW) ! number of elements in this element of the grid
     LDCOL = NUMROC(glob_ndim,NB,MYCOL,0,NPCOL) ! number of elements in this element of the grid
     CALL DESCINIT(DESC,glob_ndim,glob_ndim,NB,NB,0,0,BLACS_CONTEXT,LDROW,ierr)
  ELSE
     LDROW = 1
     LDCOL = 1
  ENDIF
  !
  ! ======================
  ! 2) Redistribute matrix
  ! ======================
  !
  CALL start_clock('diagox_redis')
  !
  n_j_loc = 0
  !
  DO j_loc = 1,idistr%nloc
     j_glob = idistr%l2g(j_loc)
     IF(j_glob > 0 .AND. j_glob <= glob_ndim) n_j_loc = n_j_loc+1
  ENDDO
  !
  IF(my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
     ALLOCATE(val_send(n_j_loc*glob_ndim))
     ALLOCATE(idx_send(n_j_loc*glob_ndim))
  ELSE
     ALLOCATE(val_send(1))
     ALLOCATE(idx_send(1))
  ENDIF
  !
  ALLOCATE(send_count(nproc))
  ALLOCATE(recv_count(nproc))
  !
  send_count = 0
  recv_count = 0
  !
  IF(my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
     ALLOCATE(dest(n_j_loc*glob_ndim))
     ALLOCATE(swap(n_j_loc*glob_ndim))
     !
     i_val = 0
     !
     DO j_loc = 1,n_j_loc
        !
        j_glob = idistr%l2g(j_loc)
        j_proc = INDXG2P(j_glob,NB,dummy,0,NPCOL)
        this_idx = INT(j_glob-1,KIND=i8b)*INT(glob_ndim,KIND=i8b)
        !
        DO i_glob = 1,glob_ndim
           i_val = i_val+1
           val_send(i_val) = a_distr(i_glob,j_loc)
           idx_send(i_val) = this_idx+INT(i_glob,KIND=i8b)
           i_proc = INDXG2P(i_glob,NB,dummy,0,NPROW)
           dest(i_val) = j_proc+i_proc*NPROW
           send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
        ENDDO
        !
     ENDDO
     !
     CALL heapsort(n_j_loc*glob_ndim,dest,swap)
     !
     DEALLOCATE(dest)
     ALLOCATE(tmp_r(n_j_loc*glob_ndim))
     ALLOCATE(tmp_i8(n_j_loc*glob_ndim))
     !
     DO i_val = 1,n_j_loc*glob_ndim
        tmp_r(i_val) = val_send(swap(i_val))
        tmp_i8(i_val) = idx_send(swap(i_val))
     ENDDO
     !
     val_send = tmp_r
     idx_send = tmp_i8
     !
     DEALLOCATE(swap)
     DEALLOCATE(tmp_r)
     DEALLOCATE(tmp_i8)
  ENDIF
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     DO j_loc = 1,LDCOL
        j_glob = INDXL2G(j_loc,NB,MYCOL,0,NPCOL)
        !
        CALL idistr%g2l(j_glob,dummy,this_sour)
        !
        this_sour = this_sour*nbgrp*nproc_bgrp
        recv_count(this_sour+1) = recv_count(this_sour+1)+LDROW
     ENDDO
     !
     ALLOCATE(idx_recv(LDROW*LDCOL))
     ALLOCATE(val_recv(LDROW*LDCOL))
  ELSE
     ALLOCATE(idx_recv(1))
     ALLOCATE(val_recv(1))
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
  CALL MPI_ALLTOALLV(idx_send,send_count,send_displ,MPI_INTEGER8,idx_recv,recv_count,&
  & recv_displ,MPI_INTEGER8,world_comm,ierr)
  !
  CALL MPI_ALLTOALLV(val_send,send_count,send_displ,MPI_DOUBLE_PRECISION,val_recv,recv_count,&
  & recv_displ,MPI_DOUBLE_PRECISION,world_comm,ierr)
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     ALLOCATE(swap(LDROW*LDCOL))
     ALLOCATE(la(LDROW,LDCOL))
     !
     CALL heapsort(LDROW*LDCOL,idx_recv,swap)
     !
     i_val = 0
     !
     DO j_loc = 1,LDCOL
        DO i_loc = 1,LDROW
           i_val = i_val+1
           la(i_loc,j_loc) = val_recv(swap(i_val))
        ENDDO
     ENDDO
     !
     DEALLOCATE(swap)
     ALLOCATE(lv(LDROW,LDCOL))
  ENDIF
  !
  CALL stop_clock('diagox_redis')
  !
  ! ==============
  ! 3) Diagonalize
  ! ==============
  !
  ge = 0.0_DP
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     !
     LWORK = -1
     ALLOCATE(work(1))
     CALL PDSYEV('V','U',glob_ndim,la,1,1,DESC,ge,lv,1,1,DESC,work,LWORK,ierr)
     !
     LWORK = INT(work(1))+1
     DEALLOCATE(work)
     !
     ALLOCATE(work(LWORK))
     CALL PDSYEV('V','U',glob_ndim,la,1,1,DESC,ge,lv,1,1,DESC,work,LWORK,ierr)
     !
     DEALLOCATE(work)
     DEALLOCATE(la)
     !
  ENDIF
  !
  ! =======================
  ! 4) Redistribute results
  ! =======================
  !
  CALL start_clock('diagox_redis')
  !
  ! Redistribute eigenvalues
  !
  e = ge(1:glob_nselect)
  CALL mp_bcast(e,0,world_comm)
  !
  ! Redistribute eigenvectors
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     ALLOCATE(dest(LDROW*LDCOL))
     ALLOCATE(swap(LDROW*LDCOL))
     !
     i_val = 0
     !
     DO j_loc = 1,LDCOL
        !
        j_glob = INDXL2G(j_loc,NB,MYCOL,0,NPCOL)
        this_idx = INT(j_glob-1,KIND=i8b)*INT(glob_ndim,KIND=i8b)
        !
        CALL idistr%g2l(j_glob,dummy,this_dest)
        !
        this_dest = this_dest*nbgrp*nproc_bgrp
        !
        DO i_loc = 1,LDROW
           i_val = i_val+1
           i_glob = INDXL2G(i_loc,NB,MYROW,0,NPROW)
           val_recv(i_val) = lv(i_loc,j_loc)
           idx_recv(i_val) = this_idx+INT(i_glob,KIND=i8b)
           dest(i_val) = this_dest
        ENDDO
        !
     ENDDO
     !
     CALL heapsort(LDROW*LDCOL,dest,swap)
     !
     DEALLOCATE(dest)
     ALLOCATE(tmp_r(LDROW*LDCOL))
     ALLOCATE(tmp_i8(LDROW*LDCOL))
     !
     DO i_val = 1,LDROW*LDCOL
        tmp_r(i_val) = val_recv(swap(i_val))
        tmp_i8(i_val) = idx_recv(swap(i_val))
     ENDDO
     !
     val_recv = tmp_r
     idx_recv = tmp_i8
     !
     DEALLOCATE(swap)
     DEALLOCATE(tmp_r)
     DEALLOCATE(tmp_i8)
  ENDIF
  !
  ! Compared to the forward redistribution, receive and send buffers are simply
  ! swapped, so no need to recalculate any of them
  !
  CALL MPI_ALLTOALLV(idx_recv,recv_count,recv_displ,MPI_INTEGER8,idx_send,send_count,&
  & send_displ,MPI_INTEGER8,world_comm,ierr)
  !
  CALL MPI_ALLTOALLV(val_recv,recv_count,recv_displ,MPI_DOUBLE_PRECISION,val_send,send_count,&
  & send_displ,MPI_DOUBLE_PRECISION,world_comm,ierr)
  !
  DEALLOCATE(send_count)
  DEALLOCATE(recv_count)
  DEALLOCATE(send_displ)
  DEALLOCATE(recv_displ)
  DEALLOCATE(val_recv)
  DEALLOCATE(idx_recv)
  !
  v_distr = 0.0_DP
  !
  IF(my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
     ALLOCATE(swap(n_j_loc*glob_ndim))
     !
     CALL heapsort(n_j_loc*glob_ndim,idx_send,swap)
     !
     i_val = 0
     !
     DO j_loc = 1,n_j_loc
        DO i_glob = 1,glob_ndim
           i_val = i_val+1
           v_distr(i_glob,j_loc) = val_send(swap(i_val))
        ENDDO
     ENDDO
     !
     DEALLOCATE(swap)
  ENDIF
  !
  DEALLOCATE(idx_send)
  DEALLOCATE(val_send)
  !
  CALL stop_clock('diagox_redis')
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     DEALLOCATE(lv)
     CALL BLACS_GRIDEXIT(BLACS_CONTEXT)
  ENDIF
  !
  ! v_distr only needed by band group 0 in the next step
  !
  IF(my_bgrp_id == 0) CALL mp_bcast(v_distr,0,intra_bgrp_comm)
  !
END SUBROUTINE
!
!----------------------------------------------------------------------------
SUBROUTINE parallel_distributed_diago_zhe(glob_nselect,glob_ndim,glob_ndimx,a_distr,v_distr,e,NPROW,NPCOL,idistr)
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
  USE kinds,                 ONLY : DP,i8b
  USE mp,                    ONLY : mp_bcast
  USE mp_global,             ONLY : intra_bgrp_comm,inter_bgrp_comm,my_bgrp_id,me_bgrp,nbgrp,nproc_bgrp
  USE mp_world,              ONLY : world_comm,mpime,nproc
  USE class_idistribute,     ONLY : idistribute
  USE sort_tools,            ONLY : heapsort
  USE parallel_include
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute),INTENT(IN) :: idistr
  INTEGER,INTENT(IN) :: glob_nselect,glob_ndim,glob_ndimx
  COMPLEX(DP),INTENT(IN) :: a_distr(glob_ndimx,idistr%nlocx)
  COMPLEX(DP),INTENT(OUT) :: v_distr(glob_ndimx,idistr%nlocx)
  REAL(DP),INTENT(OUT) :: e(glob_nselect)
  INTEGER,INTENT(INOUT) :: NPROW,NPCOL
  !
  ! Workspace
  !
  INTEGER :: i_loc,j_loc,i_glob,j_glob
  INTEGER :: n_j_loc,i_val,this_dest,this_sour,dummy
  INTEGER :: i_proc,j_proc
  INTEGER,ALLOCATABLE :: dest(:)
  INTEGER,ALLOCATABLE :: swap(:)
  INTEGER,ALLOCATABLE :: send_count(:)
  INTEGER,ALLOCATABLE :: recv_count(:)
  INTEGER,ALLOCATABLE :: send_displ(:)
  INTEGER,ALLOCATABLE :: recv_displ(:)
  INTEGER(i8b) :: this_idx
  INTEGER(i8b),ALLOCATABLE :: tmp_i8(:)
  INTEGER(i8b),ALLOCATABLE :: idx_send(:)
  INTEGER(i8b),ALLOCATABLE :: idx_recv(:)
  COMPLEX(DP),ALLOCATABLE :: tmp_c(:)
  COMPLEX(DP),ALLOCATABLE :: val_send(:)
  COMPLEX(DP),ALLOCATABLE :: val_recv(:)
  !
  ! Workspace for ScaLAPACK
  !
  COMPLEX(DP),ALLOCATABLE :: la(:,:),lv(:,:)
  REAL(DP) :: ge(glob_ndim)
  INTEGER :: NB,MYROW,MYCOL,LDROW,LDCOL
  INTEGER :: BLACS_CONTEXT
  INTEGER :: DESC(9)
  INTEGER :: LWORK,LRWORK
  INTEGER :: ierr
  INTEGER :: color
  REAL(DP),ALLOCATABLE :: rwork(:)
  COMPLEX(DP),ALLOCATABLE :: work(:)
  INTEGER,EXTERNAL :: INDXL2G
  INTEGER,EXTERNAL :: INDXG2P
  INTEGER,EXTERNAL :: NUMROC
  !
  ! ========================
  ! 1) Define processor grid
  ! ========================
  !
  ! Find suitable block size (between 1 and 64)
  !
  NB = 1
  DO WHILE(2*NB*NPCOL <= glob_ndim)
     NB = NB*2
     IF(NB == 64) EXIT
  ENDDO
  !
  i_proc = NPROW
  j_proc = NPCOL
  !
  CALL BLACS_GET(-1,0,BLACS_CONTEXT)
  CALL BLACS_GRIDINIT(BLACS_CONTEXT,'R',NPROW,NPCOL)
  CALL BLACS_GRIDINFO(BLACS_CONTEXT,NPROW,NPCOL,MYROW,MYCOL)
  !
  NPROW = i_proc
  NPCOL = j_proc
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     LDROW = NUMROC(glob_ndim,NB,MYROW,0,NPROW) ! number of elements in this element of the grid
     LDCOL = NUMROC(glob_ndim,NB,MYCOL,0,NPCOL) ! number of elements in this element of the grid
     CALL DESCINIT(DESC,glob_ndim,glob_ndim,NB,NB,0,0,BLACS_CONTEXT,LDROW,ierr)
  ELSE
     LDROW = 1
     LDCOL = 1
  ENDIF
  !
  ! ======================
  ! 2) Redistribute matrix
  ! ======================
  !
  CALL start_clock('diagox_redis')
  !
  n_j_loc = 0
  !
  DO j_loc = 1,idistr%nloc
     j_glob = idistr%l2g(j_loc)
     IF(j_glob > 0 .AND. j_glob <= glob_ndim) n_j_loc = n_j_loc+1
  ENDDO
  !
  IF(my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
     ALLOCATE(val_send(n_j_loc*glob_ndim))
     ALLOCATE(idx_send(n_j_loc*glob_ndim))
  ELSE
     ALLOCATE(val_send(1))
     ALLOCATE(idx_send(1))
  ENDIF
  !
  ALLOCATE(send_count(nproc))
  ALLOCATE(recv_count(nproc))
  !
  send_count = 0
  recv_count = 0
  !
  IF(my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
     ALLOCATE(dest(n_j_loc*glob_ndim))
     ALLOCATE(swap(n_j_loc*glob_ndim))
     !
     i_val = 0
     !
     DO j_loc = 1,n_j_loc
        !
        j_glob = idistr%l2g(j_loc)
        j_proc = INDXG2P(j_glob,NB,dummy,0,NPCOL)
        this_idx = INT(j_glob-1,KIND=i8b)*INT(glob_ndim,KIND=i8b)
        !
        DO i_glob = 1,glob_ndim
           i_val = i_val+1
           val_send(i_val) = a_distr(i_glob,j_loc)
           idx_send(i_val) = this_idx+INT(i_glob,KIND=i8b)
           i_proc = INDXG2P(i_glob,NB,dummy,0,NPROW)
           dest(i_val) = j_proc+i_proc*NPROW
           send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
        ENDDO
        !
     ENDDO
     !
     CALL heapsort(n_j_loc*glob_ndim,dest,swap)
     !
     DEALLOCATE(dest)
     ALLOCATE(tmp_c(n_j_loc*glob_ndim))
     ALLOCATE(tmp_i8(n_j_loc*glob_ndim))
     !
     DO i_val = 1,n_j_loc*glob_ndim
        tmp_c(i_val) = val_send(swap(i_val))
        tmp_i8(i_val) = idx_send(swap(i_val))
     ENDDO
     !
     val_send = tmp_c
     idx_send = tmp_i8
     !
     DEALLOCATE(swap)
     DEALLOCATE(tmp_c)
     DEALLOCATE(tmp_i8)
  ENDIF
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     DO j_loc = 1,LDCOL
        j_glob = INDXL2G(j_loc,NB,MYCOL,0,NPCOL)
        !
        CALL idistr%g2l(j_glob,dummy,this_sour)
        !
        this_sour = this_sour*nbgrp*nproc_bgrp
        recv_count(this_sour+1) = recv_count(this_sour+1)+LDROW
     ENDDO
     !
     ALLOCATE(idx_recv(LDROW*LDCOL))
     ALLOCATE(val_recv(LDROW*LDCOL))
  ELSE
     ALLOCATE(idx_recv(1))
     ALLOCATE(val_recv(1))
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
  CALL MPI_ALLTOALLV(idx_send,send_count,send_displ,MPI_INTEGER8,idx_recv,recv_count,&
  & recv_displ,MPI_INTEGER8,world_comm,ierr)
  !
  CALL MPI_ALLTOALLV(val_send,send_count,send_displ,MPI_DOUBLE_COMPLEX,val_recv,recv_count,&
  & recv_displ,MPI_DOUBLE_COMPLEX,world_comm,ierr)
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     ALLOCATE(swap(LDROW*LDCOL))
     ALLOCATE(la(LDROW,LDCOL))
     !
     CALL heapsort(LDROW*LDCOL,idx_recv,swap)
     !
     i_val = 0
     !
     DO j_loc = 1,LDCOL
        DO i_loc = 1,LDROW
           i_val = i_val+1
           la(i_loc,j_loc) = val_recv(swap(i_val))
        ENDDO
     ENDDO
     !
     DEALLOCATE(swap)
     ALLOCATE(lv(LDROW,LDCOL))
  ENDIF
  !
  CALL stop_clock('diagox_redis')
  !
  ! ==============
  ! 3) Diagonalize
  ! ==============
  !
  ge = 0.0_DP
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     !
     LWORK = -1
     LRWORK = -1
     ALLOCATE(work(1))
     ALLOCATE(rwork(1))
     CALL PZHEEV('V','U',glob_ndim,la,1,1,DESC,ge,lv,1,1,DESC,work,LWORK,rwork,LRWORK,ierr)
     !
     LWORK = CEILING(REAL(work(1),KIND=DP))+1
     DEALLOCATE(work)
     LRWORK = CEILING(rwork(1))+1
     DEALLOCATE(rwork)
     !
     ALLOCATE(work(LWORK))
     ALLOCATE(rwork(LRWORK))
     CALL PZHEEV('V','U',glob_ndim,la,1,1,DESC,ge,lv,1,1,DESC,work,LWORK,rwork,LRWORK,ierr)
     !
     DEALLOCATE(work)
     DEALLOCATE(rwork)
     DEALLOCATE(la)
     !
  ENDIF
  !
  ! =======================
  ! 4) Redistribute results
  ! =======================
  !
  CALL start_clock('diagox_redis')
  !
  ! Redistribute eigenvalues
  !
  e = ge(1:glob_nselect)
  CALL mp_bcast(e,0,world_comm)
  !
  ! Redistribute eigenvectors
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     ALLOCATE(dest(LDROW*LDCOL))
     ALLOCATE(swap(LDROW*LDCOL))
     !
     i_val = 0
     !
     DO j_loc = 1,LDCOL
        !
        j_glob = INDXL2G(j_loc,NB,MYCOL,0,NPCOL)
        this_idx = INT(j_glob-1,KIND=i8b)*INT(glob_ndim,KIND=i8b)
        !
        CALL idistr%g2l(j_glob,dummy,this_dest)
        !
        this_dest = this_dest*nbgrp*nproc_bgrp
        !
        DO i_loc = 1,LDROW
           i_val = i_val+1
           i_glob = INDXL2G(i_loc,NB,MYROW,0,NPROW)
           val_recv(i_val) = lv(i_loc,j_loc)
           idx_recv(i_val) = this_idx+INT(i_glob,KIND=i8b)
           dest(i_val) = this_dest
        ENDDO
        !
     ENDDO
     !
     CALL heapsort(LDROW*LDCOL,dest,swap)
     !
     DEALLOCATE(dest)
     ALLOCATE(tmp_c(LDROW*LDCOL))
     ALLOCATE(tmp_i8(LDROW*LDCOL))
     !
     DO i_val = 1,LDROW*LDCOL
        tmp_c(i_val) = val_recv(swap(i_val))
        tmp_i8(i_val) = idx_recv(swap(i_val))
     ENDDO
     !
     val_recv = tmp_c
     idx_recv = tmp_i8
     !
     DEALLOCATE(swap)
     DEALLOCATE(tmp_c)
     DEALLOCATE(tmp_i8)
  ENDIF
  !
  ! Compared to the forward redistribution, receive and send buffers are simply
  ! swapped, so no need to recalculate any of them
  !
  CALL MPI_ALLTOALLV(idx_recv,recv_count,recv_displ,MPI_INTEGER8,idx_send,send_count,&
  & send_displ,MPI_INTEGER8,world_comm,ierr)
  !
  CALL MPI_ALLTOALLV(val_recv,recv_count,recv_displ,MPI_DOUBLE_COMPLEX,val_send,send_count,&
  & send_displ,MPI_DOUBLE_COMPLEX,world_comm,ierr)
  !
  DEALLOCATE(send_count)
  DEALLOCATE(recv_count)
  DEALLOCATE(send_displ)
  DEALLOCATE(recv_displ)
  DEALLOCATE(val_recv)
  DEALLOCATE(idx_recv)
  !
  v_distr = 0.0_DP
  !
  IF(my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
     ALLOCATE(swap(n_j_loc*glob_ndim))
     !
     CALL heapsort(n_j_loc*glob_ndim,idx_send,swap)
     !
     i_val = 0
     !
     DO j_loc = 1,n_j_loc
        DO i_glob = 1,glob_ndim
           i_val = i_val+1
           v_distr(i_glob,j_loc) = val_send(swap(i_val))
        ENDDO
     ENDDO
     !
     DEALLOCATE(swap)
  ENDIF
  !
  DEALLOCATE(idx_send)
  DEALLOCATE(val_send)
  !
  CALL stop_clock('diagox_redis')
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     DEALLOCATE(lv)
     CALL BLACS_GRIDEXIT(BLACS_CONTEXT)
  ENDIF
  !
  ! v_distr only needed by band group 0 in the next step
  !
  IF(my_bgrp_id == 0) CALL mp_bcast(v_distr,0,intra_bgrp_comm)
  !
END SUBROUTINE
#endif
