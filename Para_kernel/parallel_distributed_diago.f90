!
! Copyright (C) 2015-2025 M. Govoni
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
  USE kinds,                 ONLY : DP,i8b
  USE mp,                    ONLY : mp_bcast,mp_comm_split,mp_comm_free
  USE mp_global,             ONLY : my_pool_id,npool,intra_bgrp_comm,my_bgrp_id,me_bgrp,nbgrp,nproc_bgrp
  USE mp_world,              ONLY : world_comm,mpime,nproc
  USE class_idistribute,     ONLY : idistribute
  USE sort_tools,            ONLY : heapsort
  USE west_mp,               ONLY : west_mp_alltoallv
#if defined(__ELPA)
  USE elpa
#endif
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
  REAL(DP),ALLOCATABLE :: ge(:),la(:,:),lv(:,:)
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
  ! Workspace for ELPA
  !
#if defined(__ELPA)
  INTEGER :: elpa_comm
  CLASS(elpa_t),POINTER :: elpa_h => NULL()
#endif
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
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
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
  send_count(:) = 0
  recv_count(:) = 0
  !
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
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
     val_send(:) = tmp_r
     idx_send(:) = tmp_i8
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
        this_sour = this_sour*npool*nbgrp*nproc_bgrp
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
  send_displ(:) = 0
  recv_displ(:) = 0
  !
  DO i_proc = 2,nproc
     send_displ(i_proc) = SUM(send_count(1:i_proc-1))
     recv_displ(i_proc) = SUM(recv_count(1:i_proc-1))
  ENDDO
  !
  CALL west_mp_alltoallv(idx_send,send_count,send_displ,idx_recv,recv_count,recv_displ,world_comm)
  CALL west_mp_alltoallv(val_send,send_count,send_displ,val_recv,recv_count,recv_displ,world_comm)
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
  ALLOCATE(ge(glob_ndim))
  ge(:) = 0.0_DP
  !
#if defined(__ELPA)
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     color = 1
  ELSE
     color = 0
  ENDIF
  !
  CALL mp_comm_split(world_comm,color,mpime,elpa_comm)
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     !
     ! Set up
     !
     ierr = elpa_init(20200417)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa_init',ierr)
     elpa_h => elpa_allocate(ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa_allocate',ierr)
     !
     CALL elpa_h%set('na',glob_ndim,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set na',ierr)
     CALL elpa_h%set('nev',glob_ndim,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set nev',ierr)
     CALL elpa_h%set('nblk',NB,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set nblk',ierr)
     CALL elpa_h%set('local_nrows',LDROW,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set local_nrows',ierr)
     CALL elpa_h%set('local_ncols',LDCOL,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set local_ncols',ierr)
     CALL elpa_h%set('mpi_comm_parent',elpa_comm,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set mpi_comm_parent',ierr)
     CALL elpa_h%set('process_row',MYROW,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set process_row',ierr)
     CALL elpa_h%set('process_col',MYCOL,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set process_col',ierr)
     !
     ierr = elpa_h%setup()
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa setup',ierr)
     !
     CALL elpa_h%set('solver',ELPA_SOLVER_2STAGE,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa set solver',ierr)
     !
     ! Fill in lower from upper
     !
     DO j_loc = 1,LDCOL
        j_glob = INDXL2G(j_loc,NB,MYCOL,0,NPCOL)
        DO i_loc = 1,LDROW
           i_glob = INDXL2G(i_loc,NB,MYROW,0,NPROW)
           IF(i_glob == j_glob) THEN
              la(i_loc,j_loc) = 0.5_DP*la(i_loc,j_loc)
           ELSEIF(i_glob > j_glob) THEN
              la(i_loc,j_loc) = 0.0_DP
           ENDIF
        ENDDO
     ENDDO
     !
     CALL PDTRAN(glob_ndim,glob_ndim,1.0_DP,la,1,1,DESC,0.0_DP,lv,1,1,DESC)
     !
     la = la+lv
     !
     ! Solve
     !
     CALL elpa_h%eigenvectors(la,ge,lv,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa eigenvectors',ierr)
     !
     ! Clean up
     !
     CALL elpa_deallocate(elpa_h,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_dsy','elpa_deallocate',ierr)
     !
     NULLIFY(elpa_h)
     DEALLOCATE(la)
     !
  ENDIF
  !
  CALL mp_comm_free(elpa_comm)
  !
#else
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
#endif
  !
  ! =======================
  ! 4) Redistribute results
  ! =======================
  !
  CALL start_clock('diagox_redis')
  !
  ! Redistribute eigenvalues
  !
  e(:) = ge(1:glob_nselect)
  CALL mp_bcast(e,0,world_comm)
  DEALLOCATE(ge)
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
        this_dest = this_dest*npool*nbgrp*nproc_bgrp
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
     val_recv(:) = tmp_r
     idx_recv(:) = tmp_i8
     !
     DEALLOCATE(swap)
     DEALLOCATE(tmp_r)
     DEALLOCATE(tmp_i8)
  ENDIF
  !
  ! Compared to the forward redistribution, receive and send buffers are simply
  ! swapped, so no need to recalculate any of them
  !
  CALL west_mp_alltoallv(idx_recv,recv_count,recv_displ,idx_send,send_count,send_displ,world_comm)
  CALL west_mp_alltoallv(val_recv,recv_count,recv_displ,val_send,send_count,send_displ,world_comm)
  !
  DEALLOCATE(send_count)
  DEALLOCATE(recv_count)
  DEALLOCATE(send_displ)
  DEALLOCATE(recv_displ)
  DEALLOCATE(val_recv)
  DEALLOCATE(idx_recv)
  !
  v_distr(:,:) = 0.0_DP
  !
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
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
  ! v_distr only needed by pool 0 and band group 0 in the next step
  !
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0) CALL mp_bcast(v_distr,0,intra_bgrp_comm)
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
  USE mp,                    ONLY : mp_bcast,mp_comm_split,mp_comm_free
  USE mp_global,             ONLY : my_pool_id,npool,intra_bgrp_comm,my_bgrp_id,me_bgrp,nbgrp,nproc_bgrp
  USE mp_world,              ONLY : world_comm,mpime,nproc
  USE class_idistribute,     ONLY : idistribute
  USE sort_tools,            ONLY : heapsort
  USE west_mp,               ONLY : west_mp_alltoallv
#if defined(__ELPA)
  USE elpa
#endif
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
  REAL(DP),ALLOCATABLE :: ge(:)
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
  ! Workspace for ELPA
  !
#if defined(__ELPA)
  INTEGER :: elpa_comm
  CLASS(elpa_t),POINTER :: elpa_h => NULL()
#endif
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
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
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
  send_count(:) = 0
  recv_count(:) = 0
  !
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
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
     val_send(:) = tmp_c
     idx_send(:) = tmp_i8
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
        this_sour = this_sour*npool*nbgrp*nproc_bgrp
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
  send_displ(:) = 0
  recv_displ(:) = 0
  !
  DO i_proc = 2,nproc
     send_displ(i_proc) = SUM(send_count(1:i_proc-1))
     recv_displ(i_proc) = SUM(recv_count(1:i_proc-1))
  ENDDO
  !
  CALL west_mp_alltoallv(idx_send,send_count,send_displ,idx_recv,recv_count,recv_displ,world_comm)
  CALL west_mp_alltoallv(val_send,send_count,send_displ,val_recv,recv_count,recv_displ,world_comm)
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
  ALLOCATE(ge(glob_ndim))
  ge(:) = 0.0_DP
  !
#if defined(__ELPA)
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     color = 1
  ELSE
     color = 0
  ENDIF
  !
  CALL mp_comm_split(world_comm,color,mpime,elpa_comm)
  !
  IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
     !
     ! Set up
     !
     ierr = elpa_init(20200417)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa_init',ierr)
     elpa_h => elpa_allocate(ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa_allocate',ierr)
     !
     CALL elpa_h%set('na',glob_ndim,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set na',ierr)
     CALL elpa_h%set('nev',glob_ndim,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set nev',ierr)
     CALL elpa_h%set('nblk',NB,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set nblk',ierr)
     CALL elpa_h%set('local_nrows',LDROW,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set local_nrows',ierr)
     CALL elpa_h%set('local_ncols',LDCOL,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set local_ncols',ierr)
     CALL elpa_h%set('mpi_comm_parent',elpa_comm,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set mpi_comm_parent',ierr)
     CALL elpa_h%set('process_row',MYROW,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set process_row',ierr)
     CALL elpa_h%set('process_col',MYCOL,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set process_col',ierr)
     !
     ierr = elpa_h%setup()
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa setup',ierr)
     !
     CALL elpa_h%set('solver',ELPA_SOLVER_2STAGE,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa set solver',ierr)
     !
     ! Fill in lower from upper
     !
     DO j_loc = 1,LDCOL
        j_glob = INDXL2G(j_loc,NB,MYCOL,0,NPCOL)
        DO i_loc = 1,LDROW
           i_glob = INDXL2G(i_loc,NB,MYROW,0,NPROW)
           IF(i_glob == j_glob) THEN
              la(i_loc,j_loc) = HALF*la(i_loc,j_loc)
           ELSEIF(i_glob > j_glob) THEN
              la(i_loc,j_loc) = ZERO
           ENDIF
        ENDDO
     ENDDO
     !
     CALL PZTRANC(glob_ndim,glob_ndim,ONE,la,1,1,DESC,ZERO,lv,1,1,DESC)
     !
     la = la+lv
     !
     DO j_loc = 1,LDCOL
        j_glob = INDXL2G(j_loc,NB,MYCOL,0,NPCOL)
        DO i_loc = 1,LDROW
           i_glob = INDXL2G(i_loc,NB,MYROW,0,NPROW)
           IF(i_glob == j_glob) THEN
              la(i_loc,j_loc) = REAL(la(i_loc,j_loc),KIND=DP)
           ENDIF
        ENDDO
     ENDDO
     !
     ! Solve
     !
     CALL elpa_h%eigenvectors(la,ge,lv,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa eigenvectors',ierr)
     !
     ! Clean up
     !
     CALL elpa_deallocate(elpa_h,ierr)
     if(ierr /= 0) CALL errore('parallel_distributed_diago_zhe','elpa_deallocate',ierr)
     !
     NULLIFY(elpa_h)
     DEALLOCATE(la)
     !
  ENDIF
  !
  CALL mp_comm_free(elpa_comm)
  !
#else
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
#endif
  !
  ! =======================
  ! 4) Redistribute results
  ! =======================
  !
  CALL start_clock('diagox_redis')
  !
  ! Redistribute eigenvalues
  !
  e(:) = ge(1:glob_nselect)
  CALL mp_bcast(e,0,world_comm)
  DEALLOCATE(ge)
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
        this_dest = this_dest*npool*nbgrp*nproc_bgrp
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
  CALL west_mp_alltoallv(idx_recv,recv_count,recv_displ,idx_send,send_count,send_displ,world_comm)
  CALL west_mp_alltoallv(val_recv,recv_count,recv_displ,val_send,send_count,send_displ,world_comm)
  !
  DEALLOCATE(send_count)
  DEALLOCATE(recv_count)
  DEALLOCATE(send_displ)
  DEALLOCATE(recv_displ)
  DEALLOCATE(val_recv)
  DEALLOCATE(idx_recv)
  !
  v_distr(:,:) = 0.0_DP
  !
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0 .AND. me_bgrp == 0) THEN
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
  ! v_distr only needed by pool 0 and band group 0 in the next step
  !
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0) CALL mp_bcast(v_distr,0,intra_bgrp_comm)
  !
END SUBROUTINE
#endif
