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
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE mp,                   ONLY : mp_max,mp_bcast
  USE io_global,            ONLY : stdout
  USE pwcom,                ONLY : nks,npw,npwx,ngk
  USE distribution_center,  ONLY : pert,aband,bandpair
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title
  USE westcom,              ONLY : n_pdep_eigen,trev_pdep,n_pdep_maxiter,n_pdep_basis,&
                                 & wstat_calculation,ev,conv,n_pdep_read_from_file,&
                                 & n_steps_write_restart,trev_pdep_rel,l_is_wstat_converged,&
                                 & nbnd_occ,lrwfc,iuwfc,dvg_exc,dng_exc,nbndval0x,&
                                 & l_preconditioning,l_bse_calculation,n_bse_idx
  USE plep_db,              ONLY : plep_db_write,plep_db_read
  USE davidson_restart,     ONLY : davidson_restart_write,davidson_restart_clear,&
                                 & davidson_restart_read
  USE wstat_tools,          ONLY : diagox,redistribute_vr_distr
  USE wbse_tools,           ONLY : wbse_build_hr,wbse_update_with_vr_distr,&
                                 & wbse_refresh_with_vr_distr,apply_preconditioning_dvg
  USE buffers,              ONLY : get_buffer
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d
  USE wvfct_gpum,           ONLY : using_et,using_et_d
  USE west_gpu,             ONLY : caux1,allocate_gpu,deallocate_gpu,allocate_bse_gpu,&
                                 & deallocate_bse_gpu,reallocate_ps_gpu,memcpy_H2D,memcpy_D2H
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
  INTEGER :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  INTEGER :: kter, nbase, np, n, ip
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: ierr,mloc,mstart,max_mloc
  INTEGER, ALLOCATABLE :: ishift(:)
  REAL(DP), ALLOCATABLE :: ew(:)
  REAL(DP), ALLOCATABLE :: hr_distr(:,:), vr_distr(:,:)
  COMPLEX(DP), ALLOCATABLE :: dng_exc_tmp(:,:,:), dvg_exc_tmp(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: dng_exc_tmp
#endif
  !
  INTEGER :: il1,ig1
  INTEGER :: iks, nbndval
  INTEGER :: do_idx
  REAL(DP) :: time_spent(2)
  CHARACTER(LEN=8) :: iter_label
  !
  REAL(DP), EXTERNAL :: GET_CLOCK
  !
  ! ... INITIALIZATION
  !
  l_is_wstat_converged = .FALSE.
  nvec = n_pdep_eigen
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
  CALL aband%init(nbndval0x,'b','nbndval',.TRUE.)
  !
  ! ... DISTRIBUTE bse_kernel
  !
  do_idx = MAXVAL(n_bse_idx)
  IF (l_bse_calculation) THEN
     bandpair = idistribute()
     CALL bandpair%init(do_idx,'i','n_pairs',.TRUE.)
  ENDIF
  !
  CALL wbse_memory_report() ! Before allocating I report the memory required.
  !
#if defined(__CUDA)
  CALL allocate_gpu()
  !
  CALL using_et(2)
  CALL using_et_d(0)
  IF(nks == 1) THEN
     CALL using_evc(2)
     CALL using_evc_d(0)
  ENDIF
#endif
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
  !$acc enter data create(dvg_exc_tmp)
  !
  ALLOCATE( dng_exc( npwx, nbndval0x, nks, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dng ', ABS(ierr) )
  !
  ALLOCATE( dng_exc_tmp( npwx, nbndval0x, nks ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dng ', ABS(ierr) )
  !$acc enter data create(dng_exc_tmp)
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
  ALLOCATE( ev( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate ev ', ABS(ierr) )
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate conv ', ABS(ierr) )
  !
  nbase = nvec
  conv = .FALSE.
  ev = 0._DP
  ew = 0._DP
  dng_exc = 0._DP
  dvg_exc = 0._DP
  hr_distr = 0._DP
  vr_distr = 0._DP
  notcnv = nvec
  dav_iter = -2
  !
  ! KIND OF CALCULATION
  !
  SELECT CASE(wstat_calculation)
  CASE('r','R')
     !
     ! RESTART
     !
     CALL davidson_restart_read( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
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
#if defined(__CUDA)
     CALL allocate_bse_gpu(aband%nloc)
#endif
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_H2D(dvg_exc_tmp,dvg_exc(:,:,:,ip),npwx*nbndval0x*nks)
#else
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
#endif
        ELSE
           !$acc kernels present(dvg_exc_tmp)
           dvg_exc_tmp(:,:,:) = (0._DP, 0._DP)
           !$acc end kernels
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp, dng_exc_tmp)
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_D2H(dng_exc(:,:,:,ip),dng_exc_tmp,npwx*nbndval0x*nks)
#else
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:)
#endif
        ENDIF
        !
     ENDDO
     !
#if defined(__CUDA)
     CALL deallocate_bse_gpu()
#endif
     !
     dav_iter = -1
     IF(n_steps_write_restart == 1) CALL davidson_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
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
#if defined(__CUDA)
     CALL allocate_bse_gpu(aband%nloc)
#endif
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_H2D(dvg_exc_tmp,dvg_exc(:,:,:,ip),npwx*nbndval0x*nks)
#else
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
#endif
        ELSE
           !$acc kernels present(dvg_exc_tmp)
           dvg_exc_tmp(:,:,:) = (0._DP, 0._DP)
           !$acc end kernels
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp, dng_exc_tmp)
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_D2H(dng_exc(:,:,:,ip),dng_exc_tmp,npwx*nbndval0x*nks)
#else
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:)
#endif
        ENDIF
        !
     ENDDO
     !
#if defined(__CUDA)
     CALL deallocate_bse_gpu()
#endif
     !
     ! </ EXTRA STEP >
     !
     ! hr = <dvg|dng>
     !
     CALL wbse_build_hr( dvg_exc, dng_exc, mstart, mstart+mloc-1, hr_distr, nvec )
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
     IF(n_steps_write_restart == 1) CALL davidson_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
     !
  ENDIF
  !
  ! --------------------
  !
  ! ... iterate
  !
  ! --------------------
  !
  iterate: DO kter = 1, n_pdep_maxiter
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
        ENDIF
        !
     ENDDO
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
           CALL apply_preconditioning_dvg( dvg_exc, notcnv, nbase, .FALSE. )
        ELSE
           CALL apply_preconditioning_dvg( dvg_exc, notcnv, nbase, .TRUE. )
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
#if defined(__CUDA)
     ALLOCATE(caux1(npwx,nbndval0x))
#endif
     !
     DO il1 = mstart, mstart+max_mloc-1
        !
        ig1 = pert%l2g(il1)
        !
        DO iks  = 1, nks
           !
           nbndval = nbnd_occ(iks)
           npw = ngk(iks)
           !
           ! ... read in GS wavefunctions iks
           !
           IF(nks > 1) THEN
              IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
              CALL mp_bcast(evc,0,inter_image_comm)
              !
#if defined(__CUDA)
              CALL using_evc(2)
              CALL using_evc_d(0)
#endif
           ENDIF
           !
           ! Pc amat
           !
           IF (.NOT.( ig1 <= nbase .OR. ig1 > nbase+notcnv )) THEN
#if defined(__CUDA)
              CALL memcpy_H2D(caux1,dvg_exc(:,:,iks,il1),npwx*nbndval0x)
              !
              CALL reallocate_ps_gpu(nbndval,nbndval)
              CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,caux1,(1._DP,0._DP))
              !
              CALL memcpy_D2H(dvg_exc(:,:,iks,il1),caux1,npwx*nbndval0x)
#else
              CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,dvg_exc(:,:,iks,il1),(1._DP,0._DP))
#endif
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
#if defined(__CUDA)
     DEALLOCATE(caux1)
#endif
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
#if defined(__CUDA)
     CALL allocate_bse_gpu(aband%nloc)
#endif
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_H2D(dvg_exc_tmp,dvg_exc(:,:,:,ip),npwx*nbndval0x*nks)
#else
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
#endif
        ELSE
           !$acc kernels present(dvg_exc_tmp)
           dvg_exc_tmp(:,:,:) = (0._DP, 0._DP)
           !$acc end kernels
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp, dng_exc_tmp)
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_D2H(dng_exc(:,:,:,ip),dng_exc_tmp,npwx*nbndval0x*nks)
#else
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:)
#endif
        ENDIF
        !
     ENDDO
     !
#if defined(__CUDA)
     CALL deallocate_bse_gpu()
#endif
     !
     ! ... update the reduced Liouville hamiltonian
     !
     ! hr = <dvg|dng>
     !
     CALL wbse_build_hr( dvg_exc, dng_exc, mstart, mstart+mloc-1, hr_distr, nbase+notcnv )
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
     ! ... eigenvectors; set the basis dimension to nvec.
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
           CALL davidson_restart_clear()
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
           WRITE( stdout, '(5X,"WARNING: ",I5," eigenvalues not converged in chidiago")' ) notcnv
           !
           CALL stop_clock( 'chidiago:last' )
           !
           EXIT iterate
           !
        ENDIF
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
     IF(n_steps_write_restart > 0 .AND. MOD(dav_iter,n_steps_write_restart) == 0) &
        CALL davidson_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
     !
  ENDDO iterate
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
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
  !$acc exit data delete(dng_exc_tmp,dvg_exc_tmp)
  DEALLOCATE( dng_exc_tmp )
  DEALLOCATE( dvg_exc_tmp )
  !
  CALL stop_clock( 'chidiago' )
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------
SUBROUTINE wbse_do_mgs (amat,m_global_start,m_global_end)
  !----------------------------------------------------------------------------
  !
  ! MGS of the vectors beloging to the interval [ m_global_start, m_global_end ]
  !    also with respect to the vectors belonging to the interval [ 1, m_global_start -1 ]
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_pool_comm,my_pool_id,intra_bgrp_comm,inter_bgrp_comm,&
                                 & my_bgrp_id,inter_image_comm,my_image_id
  USE gvect,                ONLY : gstart
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE pwcom,                ONLY : nks,npw,npwx,ngk
  USE westcom,              ONLY : nbnd_occ,nbndval0x
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : pert
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m_global_start,m_global_end
  COMPLEX(DP),INTENT(INOUT) :: amat(npwx,nbndval0x,nks,pert%nlocx)
  !
  ! Workspace
  !
  LOGICAL :: unfinished
  INTEGER :: ig,ip,ncol,ibnd,iks,nbndval
  INTEGER :: k_global,k_local,j_local,k_id
  INTEGER :: m_local_start,m_local_end
  REAL(DP) :: anorm
  COMPLEX(DP) :: za,anormc
  COMPLEX(DP),ALLOCATABLE :: zbraket(:)
  COMPLEX(DP),ALLOCATABLE :: vec(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: vec
#endif
  COMPLEX(DP),PARAMETER :: mone = (-1._DP,0._DP)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('paramgs')
#else
  CALL start_clock('paramgs')
#endif
  !
  ! 1) Run some checks
  !
  IF( m_global_start < 1 .OR. m_global_start > m_global_end .OR. m_global_end > pert%nglob ) &
  & CALL errore( 'mgs', 'wbse_do_mgs problem', 1 )
  !
  IF(my_pool_id == 0 .AND. my_bgrp_id == 0) THEN
     !
     ALLOCATE(vec(npwx,nbndval0x,nks))
     ALLOCATE(zbraket(pert%nloc))
     !
     !$acc enter data create(vec,zbraket) copyin(amat)
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
              anorm = 0._DP
              DO iks = 1, nks
                 !
                 nbndval = nbnd_occ(iks)
                 npw = ngk(iks)
                 !
                 !$acc parallel loop collapse(2) reduction(+:anorm) present(amat) copy(anorm)
                 DO ibnd = 1, nbndval
                    DO ig = 1, npw
                       anorm = anorm+REAL(amat(ig,ibnd,iks,k_local),KIND=DP)**2 &
                       & +AIMAG(amat(ig,ibnd,iks,k_local))**2
                    ENDDO
                 ENDDO
                 !$acc end parallel
                 !
                 IF(gamma_only) THEN
                    anorm = 2._DP*anorm
                    IF(gstart == 2) THEN
                       !$acc parallel loop reduction(+:anorm) present(amat) copy(anorm)
                       DO ibnd = 1, nbndval
                          anorm = anorm-REAL(amat(1,ibnd,iks,k_local),KIND=DP)**2
                       ENDDO
                       !$acc end parallel
                    ENDIF
                 ENDIF
                 !
              ENDDO
              !
              CALL mp_sum(anorm,intra_bgrp_comm)
              !
              ! normalize | k_l >
              !
              za = CMPLX(1._DP/SQRT(anorm),KIND=DP)
              !
              !$acc host_data use_device(amat)
              CALL ZSCAL(npwx*nbndval0x*nks,za,amat(1,1,1,k_local),1)
              !$acc end host_data
              !
           ENDIF
           !
           ! 5) Copy the current vector into V
           !
           !$acc host_data use_device(amat,vec)
           CALL ZCOPY(npwx*nbndval0x*nks,amat(1,1,1,k_local),1,vec,1)
           !$acc end host_data
           !
           !$acc update host(vec)
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
           !$acc update device(vec)
           !
           ! IN the range ip=j_local:pert%nloc    = >    | ip > = | ip > - | vec > * < vec | ip >
           !
           DO ip = j_local,m_local_end
              !
              IF(gamma_only) THEN
                 anorm = 0._DP
              ELSE
                 anormc = (0._DP,0._DP)
              ENDIF
              !
              DO iks = 1, nks
                 !
                 nbndval = nbnd_occ(iks)
                 npw = ngk(iks)
                 !
                 IF(gamma_only) THEN
                    !$acc parallel loop collapse(2) reduction(+:anorm) present(vec,amat) copy(anorm)
                    DO ibnd = 1, nbndval
                       DO ig = 1, npw
                          anorm = anorm+REAL(vec(ig,ibnd,iks),KIND=DP)*REAL(amat(ig,ibnd,iks,ip),KIND=DP) &
                          & +AIMAG(vec(ig,ibnd,iks))*AIMAG(amat(ig,ibnd,iks,ip))
                       ENDDO
                    ENDDO
                    !$acc end parallel
                    !
                    anorm = 2._DP*anorm
                    IF(gstart == 2) THEN
                       !$acc parallel loop reduction(+:anorm) present(vec,amat) copy(anorm)
                       DO ibnd = 1, nbndval
                          anorm = anorm-REAL(vec(1,ibnd,iks),KIND=DP)*REAL(amat(1,ibnd,iks,ip),KIND=DP)
                       ENDDO
                       !$acc end parallel
                    ENDIF
                 ELSE
                    !$acc parallel loop collapse(2) reduction(+:anormc) present(vec,amat) copy(anormc)
                    DO ibnd = 1, nbndval
                       DO ig = 1, npw
                          anormc = anormc+CONJG(vec(ig,ibnd,iks))*amat(ig,ibnd,iks,ip)
                       ENDDO
                    ENDDO
                    !$acc end parallel
                 ENDIF
                 !
              ENDDO
              !
              IF(gamma_only) THEN
                 zbraket(ip) = CMPLX(anorm,KIND=DP)
              ELSE
                 zbraket(ip) = anormc
              ENDIF
              !
           ENDDO
           !
           CALL mp_sum(zbraket(j_local:m_local_end),intra_bgrp_comm)
           !
           !$acc update device(zbraket)
           !
           ncol=m_local_end-j_local+1
           !
           !$acc host_data use_device(vec,zbraket,amat)
           CALL ZGERU(npwx*nbndval0x*nks,ncol,mone,vec,1,zbraket(j_local),1,amat(1,1,1,j_local),npwx*nbndval0x*nks)
           !$acc end host_data
           !
        ENDIF
        !
     ENDDO
     !
     !$acc exit data delete(vec,zbraket) copyout(amat)
     !
     DEALLOCATE(vec)
     DEALLOCATE(zbraket)
     !
  ENDIF
  !
  CALL mp_bcast(amat,0,inter_bgrp_comm)
  CALL mp_bcast(amat,0,inter_pool_comm)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('paramgs')
#else
  CALL stop_clock('paramgs')
#endif
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------
SUBROUTINE wbse_do_randomize ( amat, mglobalstart, mglobalend )
  !----------------------------------------------------------------------------
  !
  ! Randomize in dvg the vectors belonging to [ mglobalstart, mglobalend ]
  !
  USE kinds,                ONLY : DP
  USE random_numbers,       ONLY : randy
  USE gvect,                ONLY : g,gstart,ngm_g,ig_l2g
  USE pwcom,                ONLY : nks,npw,npwx,ngk
  USE westcom,              ONLY : lrwfc,iuwfc,nbnd_occ,nbndval0x
  USE constants,            ONLY : tpi
  USE distribution_center,  ONLY : pert
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_bcast,mp_max
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d
  USE west_gpu,             ONLY : caux1,reallocate_ps_gpu,memcpy_H2D,memcpy_D2H
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: mglobalstart, mglobalend
  COMPLEX(DP),INTENT(INOUT) :: amat(npwx,nbndval0x,nks,pert%nlocx)
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: random_num_debug(:,:)
  INTEGER :: il1,ig1,ig,ibnd,iks,nbndval
  INTEGER :: mloc,mstart,max_mloc
  REAL(DP) :: aux_real
  REAL(DP) :: rr, arg
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
#if defined(__CUDA)
  ALLOCATE(caux1(npwx,nbndval0x))
#endif
  !
  DO il1 = mstart, mstart+max_mloc-1
     !
     ig1=pert%l2g(il1)
     !
     ! Initialize the sequence
     !
     aux_real=randy(ig1)
     !
     DO iks = 1, nks
        !
        nbndval = nbnd_occ(iks)
        npw = ngk(iks)
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
                              g(3,ig)*g(3,ig) + 1._DP )
           ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
        ENDDO
        !
        ! ... read in GS wavefunctions iks
        !
        IF(nks > 1) THEN
           IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc,0,inter_image_comm)
           !
#if defined(__CUDA)
           CALL using_evc(2)
           CALL using_evc_d(0)
#endif
        ENDIF
        !
        ! Pc amat
        !
        IF (.NOT.( ig1 < mglobalstart .OR. ig1 > mglobalend )) THEN
#if defined(__CUDA)
           CALL memcpy_H2D(caux1,amat(:,:,iks,il1),npwx*nbndval0x)
           !
           CALL reallocate_ps_gpu(nbndval,nbndval)
           CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,caux1,(1._DP,0._DP))
           !
           CALL memcpy_D2H(amat(:,:,iks,il1),caux1,npwx*nbndval0x)
#else
           CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,amat(:,:,iks,il1),(1._DP,0._DP))
#endif
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(random_num_debug)
  !
#if defined(__CUDA)
  DEALLOCATE(caux1)
#endif
  !
  CALL stop_clock ('randomize')
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------
SUBROUTINE wbse_output_ev_and_time(nvec,ev,time)
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_push,              ONLY : io_push_bar
  !
  ! I/O
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nvec
  REAL(DP),INTENT(IN) :: ev(nvec)
  REAL(DP),INTENT(IN) :: time(2)
  !
  ! Workspace
  !
  INTEGER :: i,j
  CHARACTER(20),EXTERNAL :: human_readable_time
  !
  WRITE(stdout,*)
  DO i = 1, INT( nvec / 9 )
     WRITE(stdout,'(6X, 9(f9.5,1x))') (ev(j), j=9*(i-1)+1,9*i)
  ENDDO
  IF( MOD(nvec,9) > 0 ) WRITE(stdout,'(6X, 9(f9.5,1x))') (ev(j), j=9*INT(nvec/9)+1,nvec)
  WRITE(stdout,*)
  CALL io_push_bar()
  WRITE(stdout, "(5x,'Tot. elapsed time ',a,',  time spent in last iteration ',a) ") &
  TRIM(human_readable_time(time(2))), TRIM(human_readable_time(time(2)-time(1)))
  CALL io_push_bar()
  !
END SUBROUTINE
