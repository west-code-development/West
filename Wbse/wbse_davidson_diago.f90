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
!----------------------------------------------------------------------------
SUBROUTINE wbse_davidson_diago ( )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  ! ... ( L - ev ) * dvg = 0
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE distribution_center,  ONLY : pert, aband
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE pwcom,                ONLY : nks,npwx
  USE westcom,              ONLY : n_pdep_eigen,trev_pdep,n_pdep_maxiter,n_pdep_basis,&
                                 & wstat_calculation,ev,conv,n_pdep_read_from_file,trev_pdep_rel,&
                                 & l_is_wstat_converged,dvg_exc,dng_exc,nbndval0x,&
                                 & l_preconditioning,nbnd_occ,lrwfc,iuwfc
  USE plep_db,              ONLY : plep_db_write,plep_db_read
  USE wbse_restart,         ONLY : wbse_restart_write, wbse_restart_clear, wbse_restart_read
  USE mp_global,            ONLY : inter_image_comm, my_image_id
  USE mp,                   ONLY : mp_sum,mp_max,mp_bcast
  USE wstat_tools,          ONLY : diagox,serial_diagox,redistribute_vr_distr
  USE wbse_tools,           ONLY : wbse_build_hr,wbse_update_with_vr_distr,&
                                 & wbse_refresh_with_vr_distr,apply_preconditioning_dvg
  USE bse_module,           ONLY : bse_calc,size_index_matrix_lz
  USE distribution_center,  ONLY : bseparal
  USE buffers,              ONLY : get_buffer
  USE wavefunctions,        ONLY : evc
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
  INTEGER :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  INTEGER :: kter, nbase, np, n, ip
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: ierr,mloc,mstart,max_mloc
  INTEGER, ALLOCATABLE  :: ishift(:)
  REAL(DP), ALLOCATABLE :: ew(:)
  REAL(DP), ALLOCATABLE :: hr_distr(:,:), vr_distr(:,:)
  COMPLEX(DP), ALLOCATABLE :: dng_exc_tmp(:,:,:), dvg_exc_tmp(:,:,:)
  !
  INTEGER :: il1,ig1
  INTEGER :: iks, nbndval
  INTEGER :: size_index_matrix
  REAL(DP) :: time_spent(2)
  CHARACTER(LEN=8) :: iter_label
  !
  REAL(DP), EXTERNAL :: GET_CLOCK
  !
  ! ... INITIALIZATION
  !
  l_is_wstat_converged = .FALSE.
  nvec  = n_pdep_eigen
  nvecx = n_pdep_basis
  !
  CALL start_clock( 'chidiago' )
  time_spent(1)=get_clock( 'chidiago' )
  !
  ! ... DISTRIBUTE nvecx
  !
  pert = idistribute()
  CALL pert%init(nvecx,'i','nvecx',.TRUE.)
  !
  ! ... DISTRIBUTE nband
  !
  aband = idistribute()
  CALL aband%init(nbndval0x,'b','band_paralel',.TRUE.)
  !
  ! ... DISTRIBUTE bse_kernel
  !
  size_index_matrix = MAXVAL(size_index_matrix_lz(:))
  IF (bse_calc) THEN
     !
     bseparal  = idistribute()
     CALL bseparal%init(size_index_matrix,'i','bse_kernel',.TRUE.)
     !
  ENDIF
  !
  CALL wbse_memory_report() ! Before allocating I report the memory required.
  !
  ! ... MEMORY ALLOCATION
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'chidiago', 'nvecx is too small', 1 )
  !
  ALLOCATE( dvg_exc( npwx, nbndval0x, nks, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dvg ', ABS(ierr) )
  !
  ALLOCATE( dvg_exc_tmp( npwx, nbndval0x, nks), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dvg ', ABS(ierr) )
  !
  ALLOCATE( dng_exc( npwx, nbndval0x, nks, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dng ', ABS(ierr) )
  !
  ALLOCATE( dng_exc_tmp( npwx, nbndval0x, nks ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dvg ', ABS(ierr) )
  !
  ALLOCATE( hr_distr( nvecx, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate hr_distr ', ABS(ierr) )
  !
  ALLOCATE( vr_distr( nvecx, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate vr_distr ', ABS(ierr) )
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate ew ', ABS(ierr) )
  !
  ALLOCATE( ev( nvec  ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate ev ', ABS(ierr) )
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate conv ', ABS(ierr) )
  !
  nbase  = nvec
  conv   = .FALSE.
  ev     = 0._DP
  ew     = 0._DP
  dng_exc= 0._DP
  dvg_exc= 0._DP
  hr_distr(:,:) = 0._DP
  vr_distr(:,:) = 0._DP
  notcnv  = nvec
  dav_iter = -2
  !
  ! KIND OF CALCULATION
  !
  SELECT CASE(wstat_calculation)
  CASE('r','R')
     !
     ! RESTART
     !
     CALL wbse_restart_read( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
     !
  CASE('s','S')
     !
     ! FROM SCRATCH
     !
     ! ... Eventually read from file
     !
     IF(n_pdep_read_from_file>0) CALL plep_db_read( n_pdep_read_from_file )
     !
     ! ... Eventually randomize
     !
     IF(n_pdep_read_from_file<nvec) CALL wbse_do_randomize ( dvg_exc, n_pdep_read_from_file+1, nvec  )
     !
     ! ... MGS
     !
     CALL wbse_do_mgs( dvg_exc, 1, nvec)
     !
     WRITE(stdout, "( /,5x,'                  *----------*              *----------*               *----------*') ")
     WRITE(stdout, &
         & "(   5x,'#     Iteration = | ', a8,' |','   ','WBSE_dim = | ', i8,' |', '   ','Diago_dim = | ', i8,' |  x 1/2')")&
         & 'starting', nbase, nbase
     WRITE(stdout, "(   5x,'                  *----------*              *----------*               *----------*') ")
     !
     ! Apply Liouville operator
     !
     mloc = 0
     mstart = 1
     DO il1 = 1, pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 < 1 .OR. ig1 > nvec ) CYCLE
        IF( mloc==0 ) mstart = il1
        mloc = mloc + 1
     ENDDO
     !
     ! Apply Liouville operator
     !
     max_mloc = mloc
     CALL mp_max (max_mloc, inter_image_comm)
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
           !
        ELSE
           !
           dvg_exc_tmp(:,:,:) = (0.0_DP, 0.0_DP)
           !
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp(:,:,:), dng_exc_tmp(:,:,:))
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:)
           !
        ENDIF
        !
     ENDDO
     !
     dav_iter = -1
     !
  CASE DEFAULT
     CALL errore('chidiago', 'Wrong wstat_calculation',1)
  END SELECT
  !
  IF( dav_iter == -2 ) CALL errore( 'chidiago','Cannot find the 1st starting loop',1)
  !
  IF( dav_iter == -1 ) THEN
     !
     ! < EXTRA STEP >
     !
     dvg_exc = dng_exc
     CALL wbse_do_mgs( dvg_exc, 1, nvec)
     !
     WRITE(stdout, "( /,5x,'                  *----------*              *----------*               *----------*') ")
     WRITE(stdout, &
         & "(   5x,'#     Iteration = | ', a8,' |','   ','WBSE_dim = | ', i8,' |', '   ','Diago_dim = | ', i8,' |  x 2/2')")&
         & 'starting', nbase, nbase
     WRITE(stdout, "(   5x,'                  *----------*              *----------*               *----------*') ")
     !
     ! Apply Liouville operator
     !
     mloc = 0
     mstart = 1
     DO il1 = 1, pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 < 1 .OR. ig1 > nvec ) CYCLE
        IF( mloc==0 ) mstart = il1
        mloc = mloc + 1
     ENDDO
     !
     ! Apply Liouville operator
     !
     max_mloc = mloc
     CALL mp_max (max_mloc, inter_image_comm)
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
           !
        ELSE
           !
           dvg_exc_tmp(:,:,:) = (0.0_DP, 0.0_DP)
           !
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp(:,:,:), dng_exc_tmp(:,:,:))
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:)
           !
        ENDIF
        !
     ENDDO
     !
     ! </ EXTRA STEP >
     !
     ! hr = <dvg|dng>
     !
     CALL wbse_build_hr( dvg_exc, dng_exc, mstart, mstart+mloc-1, hr_distr, 1, nvec )
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL diagox( nbase, nvec, hr_distr, nvecx, ew, vr_distr )
     time_spent(2)=get_clock( 'chidiago' )
     ev(1:nvec) = ew(1:nvec)
     !
     ! Write the eigenvalues & time spent
     !
     CALL wbse_output_ev_and_time(nvec,ev,time_spent)
     !
     dav_iter = 0
     !CALL wbse_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr)
     !
  ENDIF
  !
  ! --------------------
  !
  ! ... iterate
  !
  ! --------------------
  !
  iterate: DO kter = 1,  n_pdep_maxiter
     !
     time_spent(1) = get_clock( 'chidiago' )
     !
     dav_iter = dav_iter + 1
     !
     WRITE(stdout, "( /,5x,'                  *----------*              *----------*               *----------*') ")
     WRITE(stdout, "(   5x,'#     Iteration = | ', i8,' |','   ','WBSE_dim = | ', i8,' |', '   ','Diago_dim = | ', i8,' |')") &
         &dav_iter, notcnv, nbase+notcnv
     WRITE(stdout, "(   5x,'                  *----------*              *----------*               *----------*') ")
     !
     ALLOCATE( ishift( nvecx ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( 'chidiago',' cannot allocate ishift ', ABS(ierr) )
     ishift=0
     np = 0
     !
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ...
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix
           ! ... multiplications to set a new basis vector (see below)
           !
           ishift(nbase+np) = n
           !
           ew(nbase+np) = ev(n)
           !
        END IF
        !
     END DO
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL redistribute_vr_distr( notcnv, nbase, nvecx, vr_distr, ishift )
     DEALLOCATE(ishift)
     CALL wbse_update_with_vr_distr(dvg_exc, dng_exc, notcnv, nbase, nvecx, vr_distr, ew )
     !
     IF (l_preconditioning) THEN
        !
        IF (dav_iter < 4) THEN
           CALL apply_preconditioning_dvg( dvg_exc, notcnv, nbase, nvecx, ew, .FALSE. )
        ELSE
           CALL apply_preconditioning_dvg( dvg_exc, notcnv, nbase, nvecx, ew, .TRUE. )
        ENDIF
        !
     ENDIF
     !
     ! apply P_c to the new basis vectors
     !
     mloc = 0
     mstart = 1
     DO il1 = 1, pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 <= nbase .OR. ig1 > nbase+notcnv ) CYCLE
        IF( mloc==0 ) mstart = il1
        mloc = mloc + 1
     ENDDO
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
           ! ... read in GS wavefunctions iks
           !
           IF (nks>1) THEN
              !
              IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
              CALL mp_bcast(evc,0,inter_image_comm)
              !
           ENDIF
           !
           IF (.NOT.( ig1 <= nbase .OR. ig1 > nbase+notcnv )) THEN
              !
              ! Pc amat
              !
              CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,dvg_exc(:,:,iks,il1),(1.0_DP,0.0_DP))
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
     ! ... MGS
     !
     CALL wbse_do_mgs(dvg_exc,nbase+1,nbase+notcnv)
     !
     ! apply the response function to new vectors
     !
     ! determine image that actually compute dng first
     !
     mloc = 0
     mstart = 1
     DO il1 = 1, pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 < nbase+1 .OR. ig1 > nbase+notcnv ) CYCLE
        IF( mloc==0 ) mstart = il1
        mloc = mloc + 1
     ENDDO
     !
     ! Apply Liouville operator
     !
     max_mloc = mloc
     CALL mp_max (max_mloc, inter_image_comm)
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
           !
        ELSE
           !
           dvg_exc_tmp(:,:,:) = (0.0_DP, 0.0_DP)
           !
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp(:,:,:), dng_exc_tmp(:,:,:))
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:)
           !
        ENDIF
        !
     ENDDO
     !
     ! ... update the reduced Liouville hamiltonian
     !
     ! hr = <dvg|dng>
     !
     CALL wbse_build_hr( dvg_exc, dng_exc, mstart, mstart+mloc-1, hr_distr, nbase+1, nbase+notcnv )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced Liouville hamiltonian
     !
     CALL diagox( nbase, nvec, hr_distr, nvecx, ew, vr_distr )
     time_spent(2)=get_clock( 'chidiago' )
     !
     ! ... test for convergence
     !
     conv(1:nvec) = ( ( ABS( (ew(1:nvec) - ev(1:nvec))/ev(1:nvec) ) < trev_pdep_rel ) &
                 .AND. ( ABS( ew(1:nvec) - ev(1:nvec) ) < trev_pdep ) )
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     ev(1:nvec) = ew(1:nvec)
     !
     ! Write the eigenvalues & time spent
     !
     CALL wbse_output_ev_and_time(nvec,ev,time_spent)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. kter == n_pdep_maxiter ) THEN
        !
        CALL start_clock( 'chidiago:last' )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'chidiago:last' )
           !
           CALL wbse_refresh_with_vr_distr( dvg_exc, nvec, nbase, nvecx, vr_distr )
           !
           CALL plep_db_write( )
           CALL wbse_restart_clear()
           CALL wbse_output_a_report(-1)
           !
           WRITE(iter_label,'(i8)') kter
           CALL io_push_title("Convergence achieved !!! in "//TRIM(iter_label)//" steps")
           l_is_wstat_converged = .TRUE.
           !
           EXIT iterate
           !
        ELSE IF ( kter == n_pdep_maxiter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           CALL stop_clock( 'chidiago:last' )
           !
           CALL wbse_refresh_with_vr_distr( dvg_exc, nvec, nbase, nvecx, vr_distr )
           !
           CALL plep_db_write( )
           CALL wbse_restart_clear()
           !
           WRITE( stdout, '(5X,"WARNING: ",I5, &
                &   " eigenvalues not converged in chidiago")' ) notcnv
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        WRITE(stdout,'(/,7x,"Refresh the basis set")')
        !
        CALL wbse_refresh_with_vr_distr( dvg_exc, nvec, nbase, nvecx, vr_distr )
        CALL wbse_refresh_with_vr_distr( dng_exc, nvec, nbase, nvecx, vr_distr )
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hr_distr = 0._DP
        vr_distr = 0._DP
        !
        DO il1 = 1, pert%nloc
           ig1 = pert%l2g(il1)
           IF( ig1 > nbase ) CYCLE
           hr_distr(ig1,il1) = ev(ig1)
           vr_distr(ig1,il1) = 1._DP
        ENDDO
        !
        CALL stop_clock( 'chidiago:last' )
        !
     ENDIF
     !
     IF (MOD(dav_iter,10)==0) THEN
        CALL wbse_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr)
     ENDIF
     CALL wbse_output_a_report(dav_iter)
     !
  END DO iterate
  !
! xxx  CALL wstat_xml_dump( )
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( ev )
  DEALLOCATE( hr_distr )
  DEALLOCATE( vr_distr )
  !
  DEALLOCATE( dng_exc )
  DEALLOCATE( dvg_exc )
  !
  DEALLOCATE( dng_exc_tmp )
  DEALLOCATE( dvg_exc_tmp )
  !
  CALL stop_clock( 'chidiago' )
  !
END SUBROUTINE
!
SUBROUTINE wbse_do_mgs (amat,m_global_start,m_global_end)
  !
  ! MGS of the vectors beloging to the interval [ m_global_start, m_global_end ]
  !    also with respect to the vectors belonging to the interval [ 1, m_global_start -1 ]
  !
  USE kinds,                  ONLY : DP
  USE mp_global,              ONLY : intra_bgrp_comm,inter_image_comm,my_image_id,world_comm
  USE gvect,                  ONLY : gstart
  USE mp,                     ONLY : mp_sum,mp_barrier,mp_bcast
  USE pwcom,                  ONLY : nks,npw,npwx
  USE westcom,                ONLY : nbnd_occ,nbndval0x
  USE control_flags,          ONLY : gamma_only
  USE io_push,                ONLY : io_push_title
  USE distribution_center,    ONLY : pert
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m_global_start,m_global_end
  COMPLEX(DP) :: amat(npwx,nbndval0x,nks,pert%nlocx)
  !
  ! Workspace
  !
  INTEGER :: ig, ip, ibnd, iks
  INTEGER :: k_global, k_local, j_local, k_id, nbndval
  REAL(DP) :: anorm
  REAL(DP),ALLOCATABLE :: braket(:)
  COMPLEX(DP),ALLOCATABLE :: vec(:,:,:),zbraket(:)
  LOGICAL :: unfinished
  COMPLEX(DP) :: za, anormc
  INTEGER :: m_local_start,m_local_end
  !
  REAL(DP),EXTERNAL :: DDOT
  COMPLEX(DP),EXTERNAL :: ZDOTC
  !
  CALL mp_barrier(world_comm)
  !
  CALL start_clock ('paramgs')
  !
  ALLOCATE( vec(npwx,nbndval0x,nks), zbraket(pert%nloc), braket(pert%nloc) )
  !
  ! 1) Run some checks
  !
  IF( m_global_start < 1 .OR. m_global_start > m_global_end .OR. m_global_end > pert%nglob ) &
     & CALL errore( 'mgs', 'wbse_do_mgs problem', 1 )
  !
  ! 2) Localize m_global_start
  !
  m_local_start = 1
  DO ip = 1, pert%nloc
     ig = pert%l2g(ip)
     IF( ig < m_global_start ) CYCLE
     m_local_start = ip
     EXIT
  ENDDO
  !
  ! 3) Localize m_global_end
  !
  m_local_end = pert%nloc
  DO ip = pert%nloc, 1, -1
     ig = pert%l2g(ip)
     IF( ig > m_global_end ) CYCLE
     m_local_end = ip
     EXIT
  ENDDO
  !
  j_local=1
  unfinished=.TRUE.
  !
  DO k_global=1,m_global_end
     !
     CALL pert%g2l(k_global,k_local,k_id)
     !
     IF(my_image_id==k_id) THEN
        !
        ! 4) Eventually, normalize the current vector
        !
        IF( k_global >= m_global_start ) THEN
           !
           ! anorm = < k_l | k_l >
           !
           IF(gamma_only) THEN
             !
             anorm = 0.0_DP
             !
             DO iks  = 1, nks
                !
                nbndval = nbnd_occ(iks)
                !
                DO ibnd = 1, nbndval
                   !
                   anorm = anorm + 2.0_DP * DDOT(2*npw,amat(1,ibnd,iks,k_local),1,amat(1,ibnd,iks,k_local),1)
                   !
                   IF(gstart==2) THEN
                     !
                     anorm = anorm - REAL(amat(1,ibnd,iks,k_local),KIND=DP) * REAL(amat(1,ibnd,iks,k_local),KIND=DP)

                     !
                   ENDIF
                   !
                ENDDO
                !
             ENDDO
             !
           ELSE
             !
             anorm = 0.0_DP
             !
             DO iks  = 1, nks
                !
                nbndval = nbnd_occ(iks)
                !
                DO ibnd = 1, nbndval
                   !
                   anorm  = anorm + DDOT(2*npw,amat(1,ibnd,iks,k_local),1,amat(1,ibnd,iks,k_local),1)
                   !
                ENDDO
                !
             ENDDO
             !
           ENDIF
           !
           CALL mp_sum(anorm,intra_bgrp_comm)
           !
           ! normalize | k_l >
           !
           za = CMPLX( 1._DP/SQRT(anorm), 0._DP,KIND=DP)
           !
           DO iks  = 1, nks
              !
              nbndval = nbnd_occ(iks)
              !
              DO ibnd = 1, nbndval
                 !
                 CALL ZSCAL(npw,za,amat(1,ibnd,iks,k_local),1)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
        ! 5) Copy the current vector into V
        !
        DO iks  = 1, nks
           !
           nbndval = nbnd_occ(iks)
           !
           DO ibnd = 1, nbndval
              !
              CALL ZCOPY(npwx,amat(1,ibnd,iks,k_local),1,vec(1,ibnd,iks),1)
              !
           ENDDO
           !
        ENDDO
        !
        j_local=MAX(k_local+1,m_local_start)
        !
        IF(j_local>m_local_end) unfinished=.FALSE.
        !
     ENDIF
     !
     ! BCAST | vec >
     !
     CALL mp_bcast(vec,k_id,inter_image_comm)
     !
     ! Update when needed
     !
     IF(unfinished) THEN
        !
        ! IN the range ip=j_local:pert%nloc    = >    | ip > = | ip > - | vec > * < vec | ip >
        !
        IF (gamma_only) THEN
           !
           DO ip = j_local,m_local_end !pert%nloc
              !
              anorm = 0.0_DP
              !
              DO iks  = 1, nks
                 !
                 nbndval = nbnd_occ(iks)
                 !
                 DO ibnd = 1, nbndval
                    !
                    anorm = anorm + 2._DP * DDOT(2*npw,vec(1,ibnd,iks),1,amat(1,ibnd,iks,ip),1)
                    !
                    IF (gstart==2) THEN
                       !
                       anorm = anorm - REAL(vec(1,ibnd,iks),KIND=DP)*REAL(amat(1,ibnd,iks,ip),KIND=DP)
                       !
                    ENDIF
                    !
                 ENDDO
                 !
              ENDDO
              !
              braket(ip) = anorm
              !
           ENDDO
           !
           CALL mp_sum(braket(j_local:m_local_end),intra_bgrp_comm)
           !
           DO ip = j_local,m_local_end !pert%nloc
              !
              zbraket(ip) = CMPLX( braket(ip), 0._DP, KIND=DP)
              !
           ENDDO
           !
        ELSE
           !
           DO ip = j_local,m_local_end !pert%nloc
              !
              anormc = (0.0_DP,0.0_DP)
              !
              DO iks  = 1, nks
                 !
                 nbndval = nbnd_occ(iks)
                 !
                 DO ibnd = 1, nbndval
                    !
                    anormc = anormc +  ZDOTC(npw,vec(1,ibnd,iks),1,amat(1,ibnd,iks,ip),1)
                    !
                 ENDDO
                 !
              ENDDO
              !
              zbraket(ip) = anormc
              !
           ENDDO
           !
           CALL mp_sum(zbraket(j_local:m_local_end),intra_bgrp_comm)
           !
        ENDIF
        !
        DO ip = j_local, m_local_end
           !
           amat(:,:,:,ip) = amat(:,:,:,ip) - vec(:,:,:)*zbraket(ip)
           !
        ENDDO
        !
        !
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE( vec,zbraket,braket )
  !
  CALL mp_barrier(world_comm)
  !
  CALL stop_clock ('paramgs')
  !
END SUBROUTINE
!
SUBROUTINE wbse_do_randomize ( amat, mglobalstart, mglobalend  )
  !
  ! Randomize in dvg the vectors belonging to [ mglobalstart, mglobalend ]
  !
  USE kinds,                ONLY : DP
  USE random_numbers,       ONLY : randy
  USE gvect,                ONLY : g,gstart,ngm_g,ig_l2g
  USE pwcom,                ONLY : nks,npw,npwx
  USE westcom,              ONLY : lrwfc,iuwfc,nbnd_occ,nbndval0x
  USE constants,            ONLY : tpi
  USE distribution_center,  ONLY : pert
  USE mp_global,            ONLY : my_image_id,inter_image_comm,world_comm
  USE mp,                   ONLY : mp_bcast,mp_barrier,mp_max
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,    INTENT(IN)   :: mglobalstart, mglobalend
  COMPLEX(DP),INTENT(INOUT):: amat(npwx,nbndval0x,nks,pert%nlocx)
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: random_num_debug(:,:)
  INTEGER  :: il1,ig1,ig,ibnd,iks,nbndval
  INTEGER  :: mloc,mstart,max_mloc
  REAL(DP) :: aux_real
  REAL(DP) :: rr, arg
  !
  CALL mp_barrier(world_comm)
  !
  CALL start_clock ('randomize')
  !
  ! Random numbers are generated according to G-global
  !
  ALLOCATE(random_num_debug(2,ngm_g))
  !
  mloc = 0
  mstart = 1
  DO il1 = 1, pert%nloc
     ig1 = pert%l2g(il1)
     IF( ig1 < mglobalstart .OR. ig1 > mglobalend ) CYCLE
     IF( mloc==0 ) mstart = il1
     mloc = mloc + 1
  ENDDO
  !
  max_mloc = mloc
  CALL mp_max (max_mloc, inter_image_comm)
  !
  DO il1 = mstart, mstart+max_mloc-1
     !
     ig1=pert%l2g(il1)
     !
     ! Initialize the sequence
     !
     aux_real=randy(ig1)
     !
     DO iks  = 1, nks
        !
        nbndval = nbnd_occ(iks)
        !
        DO ibnd = 1, nbndval
           !
           DO ig=1,ngm_g
              random_num_debug(1:2,ig) = (/ randy(), randy() /)
           ENDDO
           !
           amat(:,ibnd,iks,il1) = 0._DP
!$OMP PARALLEL private(ig,rr,arg)
!$OMP DO
           DO ig=gstart,npw
              rr = random_num_debug(1,ig_l2g(ig))
              arg = tpi * random_num_debug(2,ig_l2g(ig))
              amat(ig,ibnd,iks,il1) = CMPLX( rr*COS( arg ), rr*SIN( arg ), KIND=DP) / &
                            ( g(1,ig)*g(1,ig) + &
                              g(2,ig)*g(2,ig) + &
                              g(3,ig)*g(3,ig) + 1.0_DP )
           ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
        ENDDO
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
        ! Pc amat
        !
        IF (.NOT.( ig1 < mglobalstart .OR. ig1 > mglobalend )) THEN
           !
           CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,amat(:,:,iks,il1),(1.0_DP,0.0_DP))
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(random_num_debug)
  !
  CALL stop_clock ('randomize')
  !
END SUBROUTINE
!
SUBROUTINE wbse_output_a_report(iteration)
   !
   USE kinds,        ONLY : DP
   USE westcom,      ONLY : n_pdep_eigen,ev,conv
   USE west_io ,     ONLY : serial_table_output
   USE mp_world,     ONLY : mpime,root
   !
   IMPLICIT NONE
   !
   INTEGER,INTENT(IN) :: iteration
   CHARACTER(LEN=9) :: pref
   INTEGER :: ip
   REAL(DP) :: out_tab(n_pdep_eigen,3)
   !
   DO ip=1,n_pdep_eigen
      out_tab(ip,1) = REAL(ip,DP)
      out_tab(ip,2) = ev(ip)
      out_tab(ip,3) = 0._DP
      IF(conv(ip)) out_tab(ip,3) = 1._DP
   ENDDO
   IF(iteration>=0) THEN
      WRITE(pref,"('itr_',i5.5)") iteration
   ELSE
      pref='converged'
   ENDIF
   CALL serial_table_output(mpime==root,'wbse.'//TRIM(ADJUSTL(pref)),out_tab,n_pdep_eigen,3,(/'   iprt','eigenv.','  conv.'/))
   !
END SUBROUTINE
!
SUBROUTINE wbse_output_ev_and_time(nvec,ev,time)
   !
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE io_push,              ONLY : io_push_title,io_push_bar
   !
   IMPLICIT NONE
   !
   INTEGER,INTENT(IN)  :: nvec
   REAL(DP),INTENT(IN) :: ev(nvec)
   REAL(DP),INTENT(IN) :: time(2)
   !
   INTEGER :: i,j
   CHARACTER(20),EXTERNAL :: human_readable_time
   !
   WRITE(stdout,'(7X,a)') ' '
   DO i = 1, INT( nvec / 9 )
      WRITE(stdout,'(6X, 9(f9.5,1x))') (ev(j), j=9*(i-1)+1,9*i)
   ENDDO
   IF( MOD(nvec,9) > 0 ) WRITE(stdout,'(6X, 9(f9.5,1x))') (ev(j), j=9*INT(nvec/9)+1,nvec)
   WRITE(stdout,'(7X,a)') ' '
   CALL io_push_bar()
   WRITE(stdout, "(5x,'Tot. elapsed time ',a,',  time spent in last iteration ',a) ") &
   TRIM(human_readable_time(time(2))), TRIM(human_readable_time(time(2)-time(1)))
   CALL io_push_bar()
   !
END SUBROUTINE
