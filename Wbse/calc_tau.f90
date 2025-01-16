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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE calc_tau()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : isk,npw,ngk
  USE wavefunctions,        ONLY : evc
  USE westcom,              ONLY : lrwfc,iuwfc,ev,dvg,n_pdep_eigen_to_use,npwqx,nbnd_occ,l_pdep,&
                                 & spin_channel,l_bse
  USE lsda_mod,             ONLY : nspin
  USE pdep_db,              ONLY : pdep_db_read
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE buffers,              ONLY : get_buffer
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : pert,kpt_pool
  USE qbox_interface,       ONLY : init_qbox,finalize_qbox
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iks,current_spin
  INTEGER :: nbndval
  LOGICAL :: spin_resolve
  !
#if defined(__CUDA)
  IF(.NOT. l_pdep) CALL errore('calc_tau','GPU not implemented for FF_Qbox',1)
  !
  CALL allocate_gpu()
#endif
  !
  IF(l_bse) THEN
     IF(l_pdep) THEN
        pert = idistribute()
        CALL pert%init(n_pdep_eigen_to_use,'b','nvecx',.TRUE.)
        !
        ALLOCATE(dvg(npwqx,pert%nlocx))
        ALLOCATE(ev(n_pdep_eigen_to_use))
        !
        CALL pdep_db_read(n_pdep_eigen_to_use)
     ELSE
        CALL init_qbox()
     ENDIF
  ENDIF
  !
  spin_resolve = spin_channel > 0 .AND. nspin > 1
  !
  DO iks = 1,kpt_pool%nloc
     !
     current_spin = isk(iks)
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     IF((.NOT. spin_resolve) .OR. (spin_resolve .AND. current_spin == spin_channel)) THEN
        CALL calc_tau_single_q(current_spin,nbndval)
     ENDIF
     !
  ENDDO
  !
  IF(l_bse) THEN
     IF(l_pdep) THEN
        DEALLOCATE(dvg)
        DEALLOCATE(ev)
     ELSE
        CALL finalize_qbox()
     ENDIF
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE calc_tau_single_q(current_spin,nbndval)
  !-----------------------------------------------------------------------
  !
  ! Compute and store screened exchange integrals tau in case of BSE, or unscreened exchange
  ! integrals tau_u in case of hybrid TDDFT
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE io_push,              ONLY : io_push_title
  USE types_coulomb,        ONLY : pot3D
  USE westcom,              ONLY : ev,dvg,wbse_init_calculation,wbse_init_save_dir,l_bse,l_pdep,&
                                 & chi_kernel,l_local_repr,overlap_thr,n_trunc_bands
  USE fft_base,             ONLY : dffts
  USE noncollin_module,     ONLY : npol
  USE pwcom,                ONLY : npw,npwx,lsda
  USE pdep_io,              ONLY : pdep_merge_and_write_G
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : bandpair
  USE function3d,           ONLY : write_function3d,read_function3d
  USE mp,                   ONLY : mp_barrier,mp_sum
  USE lsda_mod,             ONLY : nspin
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_invfft_gamma
  USE mp_global,            ONLY : intra_image_comm,my_image_id,inter_bgrp_comm,nbgrp,&
                                 & intra_bgrp_comm,me_bgrp
  USE conversions,          ONLY : itoa
  USE qbox_interface,       ONLY : sleep_and_wait_for_lock_to_be_removed
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wbse_dv,              ONLY : wbse_dv_setup,wbse_dv_of_drho
  USE distribution_center,  ONLY : pert
  USE wavefunctions,        ONLY : evc,psic
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin,nbndval
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,ibnd_g,jbnd_g,tmp_size
  INTEGER :: il1,ig1,ir,do_idx
  INTEGER :: iu,ig,ierr,ip,ip_g
  INTEGER :: nbnd_do
  INTEGER :: dffts_nnr
  REAL(DP) :: factor
  REAL(DP) :: ovl_value
  REAL(DP), ALLOCATABLE :: ovl_matrix(:,:)
  INTEGER, ALLOCATABLE :: idx_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_loc(:,:)
  COMPLEX(DP), ALLOCATABLE :: tau(:)
  !
  REAL(DP), ALLOCATABLE :: dotp(:)
  REAL(DP), ALLOCATABLE :: aux_rr(:)
  COMPLEX(DP), ALLOCATABLE :: aux_r(:),aux1_r(:,:),aux1_g(:)
  REAL(DP), ALLOCATABLE :: frspin(:,:)
  !
  CHARACTER(LEN=:), ALLOCATABLE :: lockfile
  CHARACTER(LEN=:), ALLOCATABLE :: fname
  CHARACTER :: slabel
  CHARACTER(LEN=6) :: ilabel,jlabel
  !
  LOGICAL :: l_xcchi,l_skip,l_restart
  !
  TYPE(bar_type) :: barra
  !
  SELECT CASE(wbse_init_calculation)
  CASE('r','R')
     l_restart = .TRUE.
  CASE('s','S')
     l_restart = .FALSE.
  CASE DEFAULT
     CALL errore('calc_tau','invalid wbse_init_calculation',1)
  END SELECT
  !
  IF(.NOT. l_pdep) THEN
     !
     SELECT CASE(chi_kernel)
     CASE('XC_CHI','XC_CHI_RPA')
        l_xcchi = .TRUE.
        l_skip = .FALSE.
     CASE DEFAULT
        l_xcchi = .FALSE.
        l_skip = .TRUE.
     END SELECT
     !
     CALL wbse_dv_setup(l_skip)
     !
  ENDIF
  !
  nbnd_do = nbndval-n_trunc_bands
  dffts_nnr = dffts%nnr
  !
  tmp_size = nbnd_do*nbnd_do
  ALLOCATE(idx_matrix(tmp_size,2))
  ALLOCATE(ovl_matrix(nbnd_do,nbnd_do))
  !
  ovl_matrix(:,:) = 0._DP
  idx_matrix(:,:) = 0
  !
  IF(l_local_repr) THEN
     !
     ALLOCATE(evc_loc(npwx,nbnd_do))
     !$acc enter data create(evc_loc)
     !
     CALL wbse_localization(current_spin,n_trunc_bands+1,nbndval,evc_loc,ovl_matrix,l_restart)
     !
  ENDIF
  !
  ! compute idx_matrix
  !
  do_idx = 0
  !
  DO jbnd = 1,nbnd_do
     DO ibnd = 1,nbnd_do
        !
        ovl_value = ovl_matrix(ibnd,jbnd)
        !
        IF(.NOT. l_local_repr .OR. ovl_value >= overlap_thr) THEN
           IF(jbnd >= ibnd) THEN
              do_idx = do_idx+1
              !
              idx_matrix(do_idx,1) = ibnd+n_trunc_bands
              idx_matrix(do_idx,2) = jbnd+n_trunc_bands
           ENDIF
        ENDIF
        !
     ENDDO
  ENDDO
  !
  DEALLOCATE(ovl_matrix)
  !
  ! initialize the paralellization
  !
  bandpair = idistribute()
  CALL bandpair%init(do_idx,'i','n_pairs',.TRUE.)
  !
  ALLOCATE(tau(npwx))
  ALLOCATE(aux_r(dffts%nnr))
  ALLOCATE(aux1_g(npwx))
  !$acc enter data create(tau,aux_r,aux1_g) copyin(dvg)
  IF(l_pdep) THEN
     ALLOCATE(dotp(pert%nloc))
     !$acc enter data create(dotp)
  ENDIF
  !
  IF(l_pdep) THEN
     CALL io_push_title('Applying CHI kernel with PDEP')
  ELSE
     CALL io_push_title('Applying '//TRIM(chi_kernel)//' kernel with FF_Qbox')
  ENDIF
  !
  CALL start_bar_type(barra,'tau',bandpair%nlocx)
  !
  ! parallel loop
  !
  DO il1 = 1,bandpair%nlocx
     !
     ig1  = bandpair%l2g(il1) ! global index of n_total
     !
     ibnd_g = idx_matrix(ig1,1)
     jbnd_g = idx_matrix(ig1,2)
     ibnd = ibnd_g-n_trunc_bands
     jbnd = jbnd_g-n_trunc_bands
     !
     l_skip = .FALSE.
     !
     IF(ig1 < 1 .OR. ig1 > do_idx) l_skip = .TRUE.
     !
     IF(.NOT. l_skip) THEN
        !
        IF(l_local_repr) THEN
           CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
        ELSE
           CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd_g),evc(:,jbnd_g),psic,'Wave')
        ENDIF
        !
        !$acc parallel loop present(aux_r)
        DO ir = 1,dffts_nnr
           aux_r(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
        ENDDO
        !$acc end parallel
        !
        ! aux_r -> aux1_g
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,aux_r,aux1_g,'Wave')
        !
        ! vc in fock like term
        !
        !$acc kernels present(tau)
        tau(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        !$acc parallel loop present(tau,aux1_g,pot3D,pot3D%sqvc)
        DO ig = 1,npw
           tau(ig) = aux1_g(ig)*(pot3D%sqvc(ig)**2)
        ENDDO
        !$acc end parallel
        !
        IF(l_bse) THEN
           !
           IF(l_pdep) THEN
              !
              !$acc parallel loop present(aux1_g,pot3D,pot3D%sqvc)
              DO ig = 1,npw
                 aux1_g(ig) = aux1_g(ig)*pot3D%sqvc(ig)
              ENDDO
              !$acc end parallel
              !
              CALL glbrak_gamma(aux1_g,dvg,dotp,npw,npwx,1,pert%nloc,1,npol)
              !
              !$acc update host(dotp)
              !
              CALL mp_sum(dotp,intra_bgrp_comm)
              !
              !$acc kernels present(aux1_g)
              aux1_g(:) = (0._DP,0._DP)
              !$acc end kernels
              !
              DO ip = 1,pert%nloc
                 !
                 ip_g = pert%l2g(ip)
                 factor = dotp(ip)*ev(ip_g)/(1._DP-ev(ip_g))
                 !
                 !$acc parallel loop present(aux1_g,dvg)
                 DO ig = 1,npw
                    aux1_g(ig) = aux1_g(ig)+dvg(ig,ip)*factor
                 ENDDO
                 !$acc end parallel
                 !
              ENDDO
              !
              IF(nbgrp > 1) THEN
                 !$acc update host(aux1_g)
                 CALL mp_sum(aux1_g,inter_bgrp_comm)
                 !$acc update device(aux1_g)
              ENDIF
              !
              !$acc parallel loop present(tau,aux1_g,pot3D,pot3D%sqvc)
              DO ig = 1,npw
                 tau(ig) = tau(ig)+aux1_g(ig)*pot3D%sqvc(ig)
              ENDDO
              !$acc end parallel
              !
           ELSE
              !
              ALLOCATE(aux1_r(dffts%nnr,nspin))
              ALLOCATE(aux_rr(dffts%nnr))
              !
              ! vc in correlation like term G->0 = 0
              !
              ! aux1_g -> aux_r
              !
#if !defined(__CUDA)
              CALL single_invfft_gamma(dffts,npw,npwx,aux1_g,aux_r,'Wave')
#endif
              !
              aux1_r(:,:) = (0._DP,0._DP)
              aux1_r(:,current_spin) = aux_r
              !
              ! aux1_r = vc*aux1_r()
              !
              CALL wbse_dv_of_drho(aux1_r,.TRUE.,.FALSE.)
              !
              aux_r(:) = aux1_r(:,current_spin)
              !
              ! Send data to Qbox to compute X|vc rho>
              !
              aux_rr(:) = REAL(aux_r,KIND=DP)/SQRT(omega) ! scale down pert.
              aux_rr(:) = 0.5_DP*aux_rr ! change from rydberg to hartree
              !
              ! Write aux_rr --> fname.xml
              !
              fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml'
              !
              CALL write_function3d(fname,aux_rr,dffts)
              !
              ! DUMP A LOCK FILE
              !
              IF(me_bgrp == 0) THEN
                 lockfile = 'I.'//itoa(my_image_id)//'.lock'
                 OPEN(NEWUNIT=iu,FILE=lockfile)
                 fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml'
                 WRITE(iu,'(A)') fname
                 CLOSE(iu)
                 !
                 CALL sleep_and_wait_for_lock_to_be_removed(lockfile,'["response"]')
              ENDIF
              !
              CALL mp_barrier(intra_image_comm)
              !
              ! READ RESPONSES
              !
              IF(lsda) THEN
                 ALLOCATE(frspin(dffts%nnr,2))
                 !
                 fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response.spin0'
                 CALL read_function3d(fname,frspin(:,1),dffts)
                 fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response.spin1'
                 CALL read_function3d(fname,frspin(:,2),dffts)
                 !
                 aux_rr(:) = frspin(:,1)+frspin(:,2)
                 !
                 DEALLOCATE(frspin)
              ELSE
                 fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response'
                 CALL read_function3d(fname,aux_rr,dffts)
              ENDIF
              !
              DO ir = 1,dffts%nnr
                 aux_r(ir) = CMPLX(aux_rr(ir)*SQRT(omega),KIND=DP) ! rescale response
              ENDDO
              !
              aux1_r(:,:) = (0._DP,0._DP)
              aux1_r(:,current_spin) = aux_r
              !
              ! aux1_r = vc*aux1_r()
              !
              IF(l_xcchi) THEN
                 CALL wbse_dv_of_drho(aux1_r,.FALSE.,.FALSE.)
              ELSE
                 CALL wbse_dv_of_drho(aux1_r,.TRUE.,.FALSE.)
              ENDIF
              !
              aux_r(:) = aux1_r(:,current_spin)
              !
              ! aux_r -> aux_g
              !
#if !defined(__CUDA)
              CALL single_fwfft_gamma(dffts,npw,npwx,aux_r,aux1_g,'Wave')
#endif
              !
              ! vc + vc/fxc X vc
              !
              tau(:) = tau+aux1_g
              !
              DEALLOCATE(aux1_r)
              DEALLOCATE(aux_rr)
              !
           ENDIF
           !
        ENDIF
        !
        ! write dvg vc_rho + vc_rho X vc_rho to disk
        !
        WRITE(ilabel,'(i6.6)') ibnd_g
        WRITE(jlabel,'(i6.6)') jbnd_g
        WRITE(slabel,'(i1)') current_spin
        !
        IF(l_bse) THEN
           fname = TRIM(wbse_init_save_dir)//'/int_W'//ilabel//'_'//jlabel//'_'//slabel//'.dat'
        ELSE
           fname = TRIM(wbse_init_save_dir)//'/int_v'//ilabel//'_'//jlabel//'_'//slabel//'.dat'
        ENDIF
        !
        !$acc update host(tau)
        !
        CALL pdep_merge_and_write_G(fname,tau)
        !
     ENDIF
     !
     ! clean up
     !
     IF(.NOT. l_pdep) THEN
        IF(me_bgrp == 0) THEN
           fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml'
           OPEN(NEWUNIT=iu,IOSTAT=ierr,FILE=fname,STATUS='OLD')
           IF(ierr == 0) CLOSE(iu,STATUS='DELETE')
           IF(lsda) THEN
              fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response.spin0'
              OPEN(NEWUNIT=iu,IOSTAT=ierr,FILE=fname,STATUS='OLD')
              IF(ierr == 0) CLOSE(iu,STATUS='DELETE')
              fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response.spin1'
              OPEN(NEWUNIT=iu,IOSTAT=ierr,FILE=fname,STATUS='OLD')
              IF(ierr == 0) CLOSE(iu,STATUS='DELETE')
           ELSE
              fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response'
              OPEN(NEWUNIT=iu,IOSTAT=ierr,FILE=fname,STATUS='OLD')
              IF(ierr == 0) CLOSE(iu,STATUS='DELETE')
           ENDIF
        ENDIF
     ENDIF
     !
     CALL update_bar_type(barra,'tau',1)
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'tau')
  !
  !$acc exit data delete(dvg,tau,aux_r,aux1_g)
  DEALLOCATE(tau)
  DEALLOCATE(aux_r)
  DEALLOCATE(aux1_g)
  IF(l_pdep) THEN
     !$acc exit data delete(dotp)
     DEALLOCATE(dotp)
  ENDIF
  !
  DEALLOCATE(idx_matrix)
  !
  IF(l_local_repr) THEN
     !$acc exit data delete(evc_loc)
     DEALLOCATE(evc_loc)
  ENDIF
  !
END SUBROUTINE
