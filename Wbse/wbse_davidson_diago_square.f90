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
!----------------------------------------------------------------------------
SUBROUTINE wbse_davidson_diago_square ( )
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
  USE pwcom,                ONLY : nks,npw,npwx
  USE lsda_mod,             ONLY : nspin
  USE westcom,              ONLY : n_pdep_eigen,trev_pdep,n_pdep_maxiter,n_pdep_basis,wstat_calculation,ev,conv,&
                                 & n_pdep_restart_from_itr,n_pdep_read_from_file,n_steps_write_restart,n_pdep_times,&
                                 & trev_pdep_rel,tr2_dfpt,l_is_wstat_converged,dvg_exc,dng_exc,nbndval0x,&
                                 & l_preconditioning
  USE plep_db,              ONLY : plep_db_write,plep_db_read
  USE wbse_restart,         ONLY : wbse_restart_write, wbse_restart_clear, wbse_restart_read
  USE mp_world,             ONLY : mpime
  USE mp_global,            ONLY : inter_image_comm
  USE mp,                   ONLY : mp_sum,mp_max
  USE gvect,                ONLY : gstart
  USE wstat_tools,          ONLY : diagox,serial_diagox,symm_hr_distr,redistribute_vr_distr
  USE wbse_tools,           ONLY : wbse_build_hr,wbse_update_with_vr_distr,&
                                   wbse_refresh_with_vr_distr,apply_preconditioning_dvg
  USE bse_module,           ONLY : bse_calc,size_index_matrix_lz
  USE distribution_center,  ONLY : bseparal
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
  INTEGER :: kter, nbase, np, n, m, ip
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: ierr,mloc,mstart,mend,max_mloc
  INTEGER, ALLOCATABLE  :: ishift(:)
  REAL(DP), ALLOCATABLE :: ew(:)
  REAL(DP), ALLOCATABLE :: hr_distr(:,:), vr_distr(:,:)
  COMPLEX(DP), ALLOCATABLE :: dng_exc_tmp(:,:,:), dvg_exc_tmp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dng_exc_tmp_tmp(:,:,:)
  !
  INTEGER :: il1,il2,ig1,ig2,i
  INTEGER :: size_index_matrix
  REAL(DP) :: time_spent(2)
  REAL(DP) :: epsilon_ref, norm_tmp(nspin)
  CHARACTER(LEN=8) :: iter_label
  !
  REAL(DP), EXTERNAL :: GET_CLOCK
  !
  ! ... INITIALIZATION
  epsilon_ref = 0.38403
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
     CALL errore( 'chidiago',' cannot allocate dng ', ABS(ierr) )
  !
  ALLOCATE( dng_exc_tmp_tmp( npwx, nbndval0x, nks ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dng ', ABS(ierr) )
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
        !dng_exc_tmp_tmp(:,:,:) = dng_exc_tmp(:,:,:) - epsilon_ref * dvg_exc_tmp(:,:,:)
        !
        !CALL west_apply_liouvillian (dng_exc_tmp_tmp(:,:,:), dng_exc_tmp(:,:,:))
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:)! - epsilon_ref * dng_exc_tmp_tmp(:,:,:)
           !
        ENDIF
        !
     ENDDO
     !
     dav_iter = -1
     !
  CASE DEFAULT
     !
     CALL errore('chidiago', 'Wrong wstat_calculation',1)
     !
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
        dng_exc_tmp_tmp(:,:,:) = dng_exc_tmp(:,:,:) - epsilon_ref * dvg_exc_tmp(:,:,:)
        !
        CALL west_apply_liouvillian (dng_exc_tmp_tmp(:,:,:), dng_exc_tmp(:,:,:))
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:) - epsilon_ref * dng_exc_tmp_tmp(:,:,:)
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
     !new version call differ old version of west
     !CALL diagox( pert, nbase, nvec, hr_distr, nvecx, ew, vr_distr )
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
           !IF ( np /= n ) vr(:,np) = vr(:,n)
           ishift(nbase+np) = n
           !
           ew(nbase+np) = ev(n)
           !
        ENDIF
        !
     ENDDO
     !
     ! ... expand the basis set with new basis vectors ( H**2 - e*S )|psi> ...
     !
     CALL redistribute_vr_distr_real(notcnv, nbase, nvecx, vr_distr, ishift)
     !old version need pert
     !CALL redistribute_vr_distr(pert, notcnv, nbase, nvecx, vr_distr, ishift )
     DEALLOCATE(ishift)
     CALL wbse_update_with_vr_distr(dvg_exc, dng_exc, notcnv, nbase, nvecx, vr_distr, ew )
     !
     IF (l_preconditioning) THEN
        !
        IF (dav_iter < 4) THEN
           CALL apply_preconditioning_dvg( dvg_exc, notcnv, nbase, nvecx, ew, .FALSE., epsilon_ref)
        ELSE
           CALL apply_preconditioning_dvg( dvg_exc, notcnv, nbase, nvecx, ew, .TRUE., epsilon_ref )
        ENDIF
        !
     ENDIF
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
        dng_exc_tmp_tmp(:,:,:) = dng_exc_tmp(:,:,:) - epsilon_ref * dvg_exc_tmp(:,:,:)
        !
        CALL west_apply_liouvillian (dng_exc_tmp_tmp(:,:,:), dng_exc_tmp(:,:,:))
        !
        IF ((mstart <= ip).AND.(ip <= mstart+mloc-1)) THEN
           !
           dng_exc(:,:,:,ip) = dng_exc_tmp(:,:,:) - epsilon_ref * dng_exc_tmp_tmp(:,:,:)
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
     !old version of west need pert
     CALL symm_hr_distr(hr_distr,nbase,nvecx)
     !CALL symm_hr_distr(pert,hr_distr,nbase,nvecx)
     !
     ! ... diagonalize the reduced Liouville hamiltonian
     !
     !old version west need pert
     CALL diagox(nbase, nvec, hr_distr, nvecx, ew, vr_distr )
     !CALL diagox( pert, nbase, nvec, hr_distr, nvecx, ew, vr_distr )
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
           mloc = 0
           mstart = 1
           DO il1 = 1, pert%nloc
              ig1 = pert%l2g(il1)
              IF( ig1 < 1 .OR. ig1 > nvec ) CYCLE
              IF( mloc==0 ) mstart = il1
              mloc = mloc + 1
           ENDDO
           !
           max_mloc = mloc
           CALL mp_max (max_mloc, inter_image_comm)
           !
           ev(:) = 0.0_DP
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
                 ig1 = pert%l2g(ip)
                 !
                 CALL wbse_dot (dvg_exc_tmp,dng_exc_tmp,npwx,nbndval0x,nks,norm_tmp(nspin))
                 !
                 ev(ig1) = sum(norm_tmp(1:nspin))
                 !
              ENDIF
              !
           ENDDO
           !
           CALL mp_sum(ev,inter_image_comm)
           !
           CALL io_push_title("Original eigenvalues at the last step")
           !
           DO ip = 1, nvec
              !
              WRITE(stdout,*) ip, ev(ip)
              !
           ENDDO
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
  RETURN
  !
END SUBROUTINE
