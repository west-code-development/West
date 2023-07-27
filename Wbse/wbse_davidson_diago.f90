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
!----------------------------------------------------------------------------
SUBROUTINE wbse_davidson_diago ( )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  ! ... ( L - ev ) * dvg = 0
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_image_comm,my_image_id,nimage,inter_pool_comm,&
                                 & inter_bgrp_comm,nbgrp
  USE mp,                   ONLY : mp_max,mp_bcast
  USE west_mp,              ONLY : west_mp_get
  USE io_global,            ONLY : stdout
  USE pwcom,                ONLY : npw,npwx,ngk
  USE distribution_center,  ONLY : pert,kpt_pool,band_group
  USE class_idistribute,    ONLY : idistribute,IDIST_BLK
  USE io_push,              ONLY : io_push_title
  USE westcom,              ONLY : n_pdep_eigen,trev_pdep,n_pdep_maxiter,n_pdep_basis,ev,conv,&
                                 & wstat_calculation,n_pdep_read_from_file,n_steps_write_restart,&
                                 & trev_pdep_rel,l_is_wstat_converged,nbnd_occ,lrwfc,iuwfc,dvg_exc,&
                                 & dng_exc,nbndval0x,n_trunc_bands,l_preconditioning,l_pre_shift,&
                                 & l_spin_flip,l_forces,forces_state
  USE plep_db,              ONLY : plep_db_write,plep_db_read
  USE davidson_restart,     ONLY : davidson_restart_write,davidson_restart_clear,&
                                 & davidson_restart_read
  USE wstat_tools,          ONLY : diagox,redistribute_vr_distr
  USE wbse_tools,           ONLY : wbse_build_hr,wbse_update_with_vr_distr,&
                                 & wbse_refresh_with_vr_distr,wbse_precondition_dvg
  USE buffers,              ONLY : get_buffer
  USE wavefunctions,        ONLY : evc
  USE wbse_bgrp,            ONLY : init_gather_bands
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu,allocate_bse_gpu,deallocate_bse_gpu,&
                                 & reallocate_ps_gpu,memcpy_H2D,memcpy_D2H
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
  ATTRIBUTES(PINNED) :: dng_exc_tmp, dvg_exc_tmp
#endif
  !
  INTEGER :: iks,il1,ig1,lbnd,ibnd,iks_do
  INTEGER :: nbndval,nbnd_do,flnbndval
  INTEGER :: owner
  REAL(DP) :: time_spent(2)
  CHARACTER(LEN=8) :: iter_label
  !
  REAL(DP), EXTERNAL :: GET_CLOCK
  !
  INTEGER, PARAMETER :: flks(2) = [2,1]
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
  IF(nimage > nvecx) CALL errore('chidiago','nimage>nvecx',1)
  !
  pert = idistribute()
  CALL pert%init(nvecx,'i','nvecx',.TRUE.)
  !
  ! ... DISTRIBUTE nbndval
  !
  IF(nbgrp > nbndval0x-n_trunc_bands) CALL errore('chidiago','nbgrp>nbndval',1)
  !
  band_group = idistribute()
  CALL band_group%init(nbndval0x-n_trunc_bands,'b','nbndval',.TRUE.,IDIST_BLK)
  !
  CALL init_gather_bands()
  !
  CALL wbse_memory_report() ! Before allocating I report the memory required.
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  ! ... MEMORY ALLOCATION
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'chidiago', 'nvecx is too small', 1 )
  !
  ALLOCATE( dvg_exc( npwx, band_group%nlocx, kpt_pool%nloc, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dvg ', ABS(ierr) )
  !
  ALLOCATE( dvg_exc_tmp( npwx, band_group%nlocx, kpt_pool%nloc), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dvg ', ABS(ierr) )
  !$acc enter data create(dvg_exc_tmp)
  !
  ALLOCATE( dng_exc( npwx, band_group%nlocx, kpt_pool%nloc, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dng ', ABS(ierr) )
  !
  ALLOCATE( dng_exc_tmp( npwx, band_group%nlocx, kpt_pool%nloc ), STAT=ierr )
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
     !$acc enter data copyin(dvg_exc,dng_exc)
     !
  CASE('s','S')
     !
     ! FROM SCRATCH
     !
     ! ... Eventually read from file
     !
     IF(n_pdep_read_from_file>0) CALL plep_db_read( n_pdep_read_from_file )
     !
     ! ... Eventually initialize or randomize
     !
     IF(n_pdep_read_from_file<nvec) THEN
        CALL wbse_vc_initialize ( dvg_exc, n_pdep_read_from_file+1, nvec, l_spin_flip )
     ENDIF
     !
     dav_iter = -1
     !
  CASE DEFAULT
     CALL errore('chidiago', 'Wrong wstat_calculation',1)
  END SELECT
  !
  IF( dav_iter == -1 ) THEN
     !
     ! < EXTRA STEP >
     !
     !$acc enter data copyin(dvg_exc)
     CALL wbse_do_mgs( dvg_exc, 1, nvec, l_spin_flip )
     !$acc exit data copyout(dvg_exc)
     !
     WRITE(stdout, "(/,5x,'                  *----------*              *----------*               *----------*')")
     WRITE(stdout, "(  5x,'#     Iteration = | ',a8,' |   ','WBSE_dim = | ',i8,' |   ','Diago_dim = | ',i8,' |')") &
          & 'starting', nbase, nbase
     WRITE(stdout, "(  5x,'                  *----------*              *----------*               *----------*')")
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
     CALL allocate_bse_gpu(band_group%nlocx)
#endif
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_H2D(dvg_exc_tmp,dvg_exc(:,:,:,ip),npwx*band_group%nlocx*kpt_pool%nloc)
#else
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
#endif
        ELSE
           !$acc kernels present(dvg_exc_tmp)
           dvg_exc_tmp(:,:,:) = (0._DP, 0._DP)
           !$acc end kernels
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp, dng_exc_tmp, l_spin_flip)
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_D2H(dng_exc(:,:,:,ip),dng_exc_tmp,npwx*band_group%nlocx*kpt_pool%nloc)
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
     !$acc enter data copyin(dvg_exc,dng_exc)
     CALL wbse_build_hr( dvg_exc, dng_exc, mstart, mstart+mloc-1, hr_distr, nvec, l_spin_flip )
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL diagox( nbase, nvec, hr_distr, nvecx, ew, vr_distr )
     time_spent(2)=get_clock( 'chidiago' )
     ev(1:nvec) = ew(1:nvec)
     !
     ! Write the eigenvalues & time spent
     !
     CALL wbse_output_ev_and_time(nvec,ev,conv,time_spent,dav_iter,notcnv)
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
     WRITE(stdout, "(/,5x,'                  *----------*              *----------*               *----------*')")
     WRITE(stdout, "(  5x,'#     Iteration = | ',i8,' |   ','WBSE_dim = | ',i8,' |   ','Diago_dim = | ',i8,' |')") &
         & dav_iter, notcnv, nbase+notcnv
     WRITE(stdout, "(  5x,'                  *----------*              *----------*               *----------*')")
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
     CALL mp_bcast(vr_distr,0,inter_bgrp_comm)
     CALL mp_bcast(vr_distr,0,inter_pool_comm)
     DEALLOCATE(ishift)
     CALL wbse_update_with_vr_distr(dvg_exc, dng_exc, notcnv, nbase, nvecx, vr_distr, ew, l_spin_flip )
     !
     IF (l_preconditioning) THEN
        !
        IF (dav_iter < 4) THEN
           CALL wbse_precondition_dvg( dvg_exc, notcnv, nbase, .FALSE., l_spin_flip )
        ELSE
           CALL wbse_precondition_dvg( dvg_exc, notcnv, nbase, l_pre_shift, l_spin_flip )
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
        DO iks = 1, kpt_pool%nloc
           !
           IF(l_spin_flip) THEN
              iks_do = flks(iks)
           ELSE
              iks_do = iks
           ENDIF
           !
           flnbndval = nbnd_occ(iks_do)
           nbndval = nbnd_occ(iks)
           npw = ngk(iks)
           !
           nbnd_do = 0
           DO lbnd = 1, band_group%nloc
              ibnd = band_group%l2g(lbnd)+n_trunc_bands
              IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
           ENDDO
           !
           ! ... read in GS wavefunctions iks
           !
           IF(kpt_pool%nloc > 1) THEN
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
              CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
              CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,dvg_exc(:,:,iks,il1),(1._DP,0._DP))
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
     ! ... MGS
     !
     CALL wbse_do_mgs(dvg_exc,nbase+1,nbase+notcnv,l_spin_flip)
     !$acc exit data delete(dng_exc) copyout(dvg_exc)
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
     CALL allocate_bse_gpu(band_group%nlocx)
#endif
     !
     DO ip = mstart, mstart+max_mloc-1
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_H2D(dvg_exc_tmp,dvg_exc(:,:,:,ip),npwx*band_group%nlocx*kpt_pool%nloc)
#else
           dvg_exc_tmp(:,:,:) = dvg_exc(:,:,:,ip)
#endif
        ELSE
           !$acc kernels present(dvg_exc_tmp)
           dvg_exc_tmp(:,:,:) = (0._DP, 0._DP)
           !$acc end kernels
        ENDIF
        !
        CALL west_apply_liouvillian (dvg_exc_tmp, dng_exc_tmp, l_spin_flip)
        !
        IF (mstart <= ip .AND. ip <= mstart+mloc-1) THEN
#if defined(__CUDA)
           CALL memcpy_D2H(dng_exc(:,:,:,ip),dng_exc_tmp,npwx*band_group%nlocx*kpt_pool%nloc)
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
     !$acc enter data copyin(dvg_exc,dng_exc)
     CALL wbse_build_hr( dvg_exc, dng_exc, mstart, mstart+mloc-1, hr_distr, nbase+notcnv, l_spin_flip )
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
     ! Print max difference
     !
     WRITE(stdout,'(5X,A,E10.3)') 'Max diff = ', MAXVAL(ABS(ew(1:nvec) - ev(1:nvec)))
     !
     ev(1:nvec) = ew(1:nvec)
     !
     ! Write the eigenvalues & time spent
     !
     CALL wbse_output_ev_and_time(nvec,ev,conv,time_spent,dav_iter,notcnv)
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
           CALL wbse_refresh_with_vr_distr( dvg_exc, nvec, nbase, nvecx, vr_distr, l_spin_flip )
           !$acc update host(dvg_exc)
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
        ELSEIF ( kter == n_pdep_maxiter ) THEN
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
        CALL wbse_refresh_with_vr_distr( dvg_exc, nvec, nbase, nvecx, vr_distr, l_spin_flip )
        !$acc update host(dvg_exc)
        CALL wbse_refresh_with_vr_distr( dng_exc, nvec, nbase, nvecx, vr_distr, l_spin_flip )
        !$acc update host(dng_exc)
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
  !$acc exit data delete(dvg_exc,dng_exc)
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( ev )
  DEALLOCATE( hr_distr )
  DEALLOCATE( vr_distr )
  !
  DEALLOCATE( dng_exc )
  !
  CALL stop_clock( 'chidiago' )
  !
  IF(l_forces) THEN
     !
#if defined(__CUDA)
     CALL allocate_bse_gpu(band_group%nlocx)
#endif
     !
     ! send forces_state to root image
     !
     CALL pert%g2l(forces_state,il1,owner)
     !
     CALL west_mp_get(dvg_exc_tmp,dvg_exc(:,:,:,il1),my_image_id,0,owner,owner,inter_image_comm)
     !
     !$acc update device(dvg_exc_tmp)
     !
     DEALLOCATE( dvg_exc )
     !
     ! root image computes forces
     !
     CALL wbse_calc_forces( dvg_exc_tmp )
     !
     !$acc exit data delete(dvg_exc_tmp)
     DEALLOCATE( dvg_exc_tmp )
     !
#if defined(__CUDA)
     CALL deallocate_bse_gpu()
#endif
     !
  ELSE
     !
     DEALLOCATE( dvg_exc )
     !$acc exit data delete(dvg_exc_tmp)
     DEALLOCATE( dvg_exc_tmp )
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------
SUBROUTINE wbse_do_mgs (amat,m_global_start,m_global_end,sf)
  !----------------------------------------------------------------------------
  !
  ! MGS of the vectors beloging to the interval [ m_global_start, m_global_end ]
  !    also with respect to the vectors belonging to the interval [ 1, m_global_start -1 ]
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,&
                                 & intra_bgrp_comm
  USE gvect,                ONLY : gstart
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE pwcom,                ONLY : npw,npwx,ngk
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands
  USE distribution_center,  ONLY : pert,kpt_pool,band_group
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m_global_start,m_global_end
  COMPLEX(DP),INTENT(INOUT) :: amat(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
  LOGICAL,INTENT(IN) :: sf
  !
  ! Workspace
  !
  LOGICAL :: unfinished
  INTEGER :: ig,ip,ncol,lbnd,ibnd,iks,nbndval,nbnd_do,iks_do
  INTEGER :: k_global,k_local,j_local,k_id
  INTEGER :: m_local_start,m_local_end
  REAL(DP) :: anorm
  COMPLEX(DP) :: za
  COMPLEX(DP),ALLOCATABLE :: zbraket(:)
  COMPLEX(DP),ALLOCATABLE :: vec(:,:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: vec
#endif
  COMPLEX(DP),PARAMETER :: mone = (-1._DP,0._DP)
  INTEGER,PARAMETER :: flks(2) = [2,1]
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
  ALLOCATE(vec(npwx,band_group%nlocx,kpt_pool%nloc))
  ALLOCATE(zbraket(pert%nloc))
  !
  !$acc enter data create(vec,zbraket)
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
           DO iks = 1, kpt_pool%nloc
              !
              IF(sf) THEN
                 iks_do = flks(iks)
              ELSE
                 iks_do = iks
              ENDIF
              !
              nbndval = nbnd_occ(iks_do)
              npw = ngk(iks)
              !
              nbnd_do = 0
              DO lbnd = 1, band_group%nloc
                 ibnd = band_group%l2g(lbnd)+n_trunc_bands
                 IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
              ENDDO
              !
              !$acc parallel loop collapse(2) reduction(+:anorm) present(amat) copy(anorm)
              DO lbnd = 1, nbnd_do
                 DO ig = 1, npw
                    anorm = anorm+2._DP*REAL(amat(ig,lbnd,iks,k_local),KIND=DP)**2 &
                    & +2._DP*AIMAG(amat(ig,lbnd,iks,k_local))**2
                 ENDDO
              ENDDO
              !$acc end parallel
              !
              IF(gstart == 2) THEN
                 !$acc parallel loop reduction(+:anorm) present(amat) copy(anorm)
                 DO lbnd = 1, nbnd_do
                    anorm = anorm-REAL(amat(1,lbnd,iks,k_local),KIND=DP)**2
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ENDDO
           !
           CALL mp_sum(anorm,intra_bgrp_comm)
           CALL mp_sum(anorm,inter_bgrp_comm)
           CALL mp_sum(anorm,inter_pool_comm)
           !
           ! normalize | k_l >
           !
           za = CMPLX(1._DP/SQRT(anorm),KIND=DP)
           !
           !$acc host_data use_device(amat)
           CALL ZSCAL(npwx*band_group%nlocx*kpt_pool%nloc,za,amat(1,1,1,k_local),1)
           !$acc end host_data
           !
        ENDIF
        !
        ! 5) Copy the current vector into V
        !
        !$acc host_data use_device(amat,vec)
        CALL ZCOPY(npwx*band_group%nlocx*kpt_pool%nloc,amat(1,1,1,k_local),1,vec,1)
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
           anorm = 0._DP
           !
           DO iks = 1, kpt_pool%nloc
              !
              IF(sf) THEN
                 iks_do = flks(iks)
              ELSE
                 iks_do = iks
              ENDIF
              !
              nbndval = nbnd_occ(iks_do)
              npw = ngk(iks)
              !
              nbnd_do = 0
              DO lbnd = 1, band_group%nloc
                 ibnd = band_group%l2g(lbnd)+n_trunc_bands
                 IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
              ENDDO
              !
              !$acc parallel loop collapse(2) reduction(+:anorm) present(vec,amat) copy(anorm)
              DO lbnd = 1, nbnd_do
                 DO ig = 1, npw
                    anorm = anorm+2._DP*REAL(vec(ig,lbnd,iks),KIND=DP)*REAL(amat(ig,lbnd,iks,ip),KIND=DP) &
                    & +2._DP*AIMAG(vec(ig,lbnd,iks))*AIMAG(amat(ig,lbnd,iks,ip))
                 ENDDO
              ENDDO
              !$acc end parallel
              !
              IF(gstart == 2) THEN
                 !$acc parallel loop reduction(+:anorm) present(vec,amat) copy(anorm)
                 DO lbnd = 1, nbnd_do
                    anorm = anorm-REAL(vec(1,lbnd,iks),KIND=DP)*REAL(amat(1,lbnd,iks,ip),KIND=DP)
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ENDDO
           !
           zbraket(ip) = CMPLX(anorm,KIND=DP)
           !
        ENDDO
        !
        CALL mp_sum(zbraket(j_local:m_local_end),intra_bgrp_comm)
        CALL mp_sum(zbraket(j_local:m_local_end),inter_bgrp_comm)
        CALL mp_sum(zbraket(j_local:m_local_end),inter_pool_comm)
        !
        !$acc update device(zbraket)
        !
        ncol=m_local_end-j_local+1
        !
        !$acc host_data use_device(vec,zbraket,amat)
        CALL ZGERU(npwx*band_group%nlocx*kpt_pool%nloc,ncol,mone,vec,1,zbraket(j_local),1,&
        & amat(1,1,1,j_local),npwx*band_group%nlocx*kpt_pool%nloc)
        !$acc end host_data
        !
     ENDIF
     !
  ENDDO
  !
  !$acc exit data delete(vec,zbraket)
  !
  DEALLOCATE(vec)
  DEALLOCATE(zbraket)
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
SUBROUTINE wbse_output_ev_and_time(nvec,ev_,conv_,time,dav_iter,notcnv)
  !----------------------------------------------------------------------------
  !
  USE json_module,          ONLY : json_file
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_push,              ONLY : io_push_bar
  USE mp_world,             ONLY : mpime,root
  USE westcom,              ONLY : logfile
  !
  ! I/O
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nvec
  REAL(DP),INTENT(IN) :: ev_(nvec)
  LOGICAL,INTENT(IN) :: conv_(nvec)
  REAL(DP),INTENT(IN) :: time(2)
  INTEGER,INTENT(IN) :: dav_iter,notcnv
  !
  ! Workspace
  !
  TYPE(json_file) :: json
  INTEGER :: i,j
  CHARACTER(20),EXTERNAL :: human_readable_time
  CHARACTER(LEN=6) :: cdav
  INTEGER :: ndav,iunit
  LOGICAL :: found
  !
  WRITE(stdout,*)
  DO i = 1, INT( nvec / 9 )
     WRITE(stdout,'(6X, 9(f9.5,1x))') (ev_(j), j=9*(i-1)+1,9*i)
  ENDDO
  IF( MOD(nvec,9) > 0 ) WRITE(stdout,'(6X, 9(f9.5,1x))') (ev_(j), j=9*INT(nvec/9)+1,nvec)
  WRITE(stdout,*)
  CALL io_push_bar()
  WRITE(stdout, "(5x,'Tot. elapsed time ',a,',  time spent in last iteration ',a) ") &
  TRIM(human_readable_time(time(2))), TRIM(human_readable_time(time(2)-time(1)))
  CALL io_push_bar()
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     !
     CALL json%load(filename=TRIM(logfile))
     !
     CALL json%get('exec.ndav', ndav, found)
     !
     IF( found ) THEN
        ndav = ndav+1
     ELSE
        ndav = 1
     ENDIF
     !
     WRITE(cdav,'(i6)') ndav
     !
     CALL json%update('exec.ndav', ndav, found)
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').dav_iter', dav_iter)
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').ev', ev_(1:nvec))
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').conv', conv_(1:nvec))
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').notcnv', notcnv)
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').time_elap:sec', time(2))
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').time_elap:hum', TRIM(human_readable_time(time(2))))
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').time_iter:sec', (time(2)-time(1)))
     CALL json%add('exec.davitr('//TRIM(ADJUSTL(cdav))//').time_iter:hum', TRIM(human_readable_time(time(2)-time(1))))
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------
SUBROUTINE wbse_vc_initialize(amat,mglobalstart,mglobalend,sf)
  !----------------------------------------------------------------------------
  !
  ! Adapted from lr_dav_set_init in Turbo-TDDFPT
  !
  ! Use couples of occ/vir states as initial guess
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer
  USE pwcom,                ONLY : npwx,nbnd,et,nspin
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,my_pool_id,&
                                 & my_bgrp_id
  USE mp,                   ONLY : mp_sum,mp_bcast,mp_max
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,n_trunc_bands
  USE distribution_center,  ONLY : pert,kpt_pool,band_group
  USE sort_tools,           ONLY : heapsort
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: mglobalstart,mglobalend
  COMPLEX(DP),INTENT(INOUT) :: amat(npwx,band_group%nlocx,kpt_pool%nloc,pert%nlocx)
  LOGICAL,INTENT(IN) :: sf
  !
  ! Workspace
  !
  INTEGER :: iv,ic,ib,lv,nbndval,flnbndval,iks,iks_do,is,is_g
  INTEGER :: nbnd_v_window,nbnd_c_window,npair
  INTEGER :: il1,ig1,nvec
  INTEGER :: itmp1,itmp2
  INTEGER :: mloc,mstart,max_mloc
  INTEGER :: owner
  INTEGER,ALLOCATABLE :: occ(:)
  INTEGER,ALLOCATABLE :: e_diff_order(:)
  REAL(DP),ALLOCATABLE :: e_diff(:)
  INTEGER,PARAMETER :: flks(2) = [2,1]
  !
  CALL start_clock ('vc_init')
  !
  WRITE(stdout,'(5x,"Using lowest energy electron-hole pairs as initial guess")')
  !
  ! get global copy of occupation
  !
  ALLOCATE(occ(nspin))
  !
  occ(:) = 0._DP
  DO is = 1, kpt_pool%nloc
     is_g = kpt_pool%l2g(is)
     occ(is_g) = nbnd_occ(is)
  ENDDO
  !
  CALL mp_sum(occ,inter_pool_comm)
  !
  nvec = mglobalend-mglobalstart+1
  nbnd_v_window = MIN(CEILING(SQRT(REAL(4*nvec,KIND=DP))), MINVAL(occ))
  nbnd_c_window = MIN(CEILING(SQRT(REAL(4*nvec,KIND=DP))), nbnd-MAXVAL(occ))
  npair = nbnd_v_window*nbnd_c_window
  !
  ALLOCATE(e_diff(nspin*npair))
  ALLOCATE(e_diff_order(nspin*npair))
  !
  IF(my_pool_id == 0) THEN
     !
     ib = 0
     DO iks = 1,nspin
        !
        IF(sf) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        !
        nbndval = occ(iks)
        flnbndval = occ(iks_do)
        !
        DO iv = flnbndval-nbnd_v_window+1, flnbndval
           DO ic = nbndval+1, nbndval+nbnd_c_window
              ib = ib + 1
              e_diff(ib) = et(ic,iks) - et(iv,iks_do)
           ENDDO
        ENDDO
        !
     ENDDO
     !
     CALL heapsort(nspin*npair,e_diff,e_diff_order)
     !
  ENDIF
  !
  CALL mp_bcast(e_diff,0,inter_pool_comm)
  CALL mp_bcast(e_diff_order,0,inter_pool_comm)
  !
  mloc = 0
  mstart = 1
  DO il1 = 1,pert%nloc
     ig1 = pert%l2g(il1)
     IF(ig1 < mglobalstart .OR. ig1 > mglobalend) CYCLE
     IF(mloc == 0) mstart = il1
     mloc = mloc + 1
  ENDDO
  !
  max_mloc = mloc
  CALL mp_max(max_mloc,inter_image_comm)
  !
  amat(:,:,:,mstart:mstart+max_mloc-1) = (0._DP,0._DP)
  !
  DO is = 1,kpt_pool%nloc
     !
     is_g = kpt_pool%l2g(is)
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,is)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     DO il1 = mstart,mstart+max_mloc-1
        !
        ig1 = pert%l2g(il1)
        IF(ig1 < mglobalstart .OR. ig1 > mglobalend) CYCLE
        !
        itmp1 = e_diff_order(ig1)
        itmp2 = MOD(itmp1-1, npair) + 1
        !
        iks = (itmp1-1)/npair + 1
        IF(iks /= is_g) CYCLE
        !
        IF(sf) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        !
        nbndval = occ(iks)
        flnbndval = occ(iks_do)
        !
        iv = (itmp2-1)/nbnd_c_window + 1 + flnbndval - nbnd_v_window - n_trunc_bands
        ic = MOD(itmp2-1, nbnd_c_window) + 1 + nbndval
        !
        CALL band_group%g2l(iv,lv,owner)
        !
        IF(owner == my_bgrp_id) amat(:,lv,is,il1) = evc(:,ic)
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(e_diff)
  DEALLOCATE(e_diff_order)
  DEALLOCATE(occ)
  !
  CALL stop_clock ('vc_init')
  !
END SUBROUTINE
