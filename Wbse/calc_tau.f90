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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE calc_tau()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : isk,nks,npw,ngk
  USE wavefunctions,        ONLY : evc
  USE westcom,              ONLY : lrwfc,iuwfc,ev,dvg,n_pdep_eigen_to_use,npwqx,nbnd_occ,&
                                 & wbse_init_calculation,bse_method,spin_channel
  USE lsda_mod,             ONLY : nspin
  USE pdep_db,              ONLY : pdep_db_read
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE buffers,              ONLY : get_buffer
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : pert
  USE qbox_interface,       ONLY : init_qbox,finalize_qbox
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iks,current_spin
  INTEGER :: iq,nkq,ikq
  LOGICAL :: l_restart_calc,l_pdep,spin_resolve
  !
  SELECT CASE(wbse_init_calculation)
  CASE('r','R')
     l_restart_calc = .TRUE.
  CASE('s','S')
     l_restart_calc = .FALSE.
  CASE DEFAULT
     CALL errore('wbse_init','invalid wbse_init_calculation',1)
  END SELECT
  !
  SELECT CASE(TRIM(bse_method))
  CASE('PDEP','pdep')
     l_pdep = .TRUE.
  CASE('FF_QBOX','FF_Qbox','ff_qbox')
     l_pdep = .FALSE.
  END SELECT
  !
  IF(l_pdep) THEN
     pert = idistribute()
     CALL pert%init(n_pdep_eigen_to_use,'b','nvecx',.TRUE.)
     !
     ALLOCATE(dvg(npwqx,pert%nlocx))
     ALLOCATE(ev(n_pdep_eigen_to_use))
     !
     CALL pdep_db_read(n_pdep_eigen_to_use)
     CALL start_clock('tau_pdep')
  ELSE
     CALL init_qbox()
     CALL start_clock('tau_qbox')
  ENDIF
  !
  spin_resolve = spin_channel > 0 .AND. nspin > 1
  !
  nkq = 1
  !
  DO iq = 1,nkq
     DO iks = 1,nks
        !
        current_spin = isk(iks)
        ikq = 1
        npw = ngk(iks)
        !
        IF(nks > 1) THEN
           IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc,0,inter_image_comm)
        ENDIF
        !
        IF((.NOT. spin_resolve) .OR. (spin_resolve .AND. current_spin == spin_channel)) THEN
           CALL calc_tau_single_q(iks,iq,current_spin,nbnd_occ(iks),l_restart_calc)
        ENDIF
        !
     ENDDO
  ENDDO
  !
  IF(l_pdep) THEN
     CALL stop_clock('tau_pdep')
     !
     DEALLOCATE(dvg)
     DEALLOCATE(ev)
  ELSE
     CALL stop_clock('tau_qbox')
     CALL finalize_qbox()
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE calc_tau_single_q(iks,ikq,current_spin,nbndval,l_restart_calc)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE io_push,              ONLY : io_push_title
  USE types_coulomb,        ONLY : pot3D
  USE westcom,              ONLY : ev,dvg,wbse_init_save_dir,bse_method,chi_kernel,l_local_repr,&
                                 & overlap_thr
  USE wavefunctions,        ONLY : evc,psic
  USE fft_base,             ONLY : dffts
  USE noncollin_module,     ONLY : npol
  USE pwcom,                ONLY : npw,npwx,lsda
  USE pdep_io,              ONLY : pdep_merge_and_write_G
  USE wbse_init_restart,    ONLY : wbse_status_restart_read,wbse_status_restart_write,&
                                 & wbse_index_matrix_read,wbse_index_matrix_write
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : bandpair
  USE function3d,           ONLY : write_function3d,read_function3d
  USE mp,                   ONLY : mp_barrier,mp_sum
  USE lsda_mod,             ONLY : nspin
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_invfft_gamma
  USE mp_global,            ONLY : inter_image_comm,intra_image_comm,my_image_id,inter_bgrp_comm,&
                                 & intra_bgrp_comm,me_bgrp
  USE conversions,          ONLY : itoa
  USE qbox_interface,       ONLY : sleep_and_wait_for_lock_to_be_removed
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wbse_dv,              ONLY : wbse_dv_setup,wbse_dv_of_drho
  USE distribution_center,  ONLY : pert
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: iks,ikq,current_spin,nbndval
  LOGICAL, INTENT(IN) :: l_restart_calc
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,tmp_size
  INTEGER :: il1,ig1,ir,do_idx
  INTEGER :: iu,ig,ierr,ip,ip_g
  REAL(DP) :: factor
  REAL(DP) :: ovl_value
  REAL(DP), ALLOCATABLE :: ovl_matrix(:,:)
  INTEGER, ALLOCATABLE :: restart_matrix(:)
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
  CHARACTER(LEN=40) :: flabel
  !
  LOGICAL :: l_xcchi
  LOGICAL :: calc_is_done
  LOGICAL :: l_skip
  LOGICAL :: l_pdep
  !
  TYPE(bar_type) :: barra
  !
  SELECT CASE(TRIM(bse_method))
  CASE('PDEP','pdep')
     l_pdep = .TRUE.
  CASE('FF_QBOX','FF_Qbox','ff_qbox')
     l_pdep = .FALSE.
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
  WRITE(flabel,'(A,I6.6,A,I6.6,A,I1,A)') '_iq',ikq,'_ik',iks,'_spin',current_spin,'.dat'
  !
  tmp_size = nbndval*nbndval
  ALLOCATE(idx_matrix(tmp_size,2))
  ALLOCATE(ovl_matrix(nbndval,nbndval))
  !
  ovl_matrix(:,:) = 0._DP
  idx_matrix(:,:) = 0
  !
  IF(l_local_repr) THEN
     ALLOCATE(evc_loc(npwx,nbndval))
     CALL wbse_localization(current_spin,nbndval,evc_loc,ovl_matrix,l_restart_calc)
  ENDIF
  !
  ! compute idx_matrix
  !
  do_idx = 0
  !
  IF(.NOT. l_restart_calc) THEN
     !
     DO jbnd = 1,nbndval
        DO ibnd = 1,nbndval
           !
           ovl_value = ovl_matrix(ibnd,jbnd)
           !
           IF(ovl_value >= overlap_thr) THEN
              IF(jbnd >= ibnd) THEN
                 do_idx = do_idx + 1
                 !
                 idx_matrix(do_idx,1) = ibnd
                 idx_matrix(do_idx,2) = jbnd
              ENDIF
           ENDIF
           !
        ENDDO
     ENDDO
     !
     fname = TRIM(wbse_init_save_dir)//'/index_matrix'//flabel
     CALL wbse_index_matrix_write(fname,do_idx,2,idx_matrix(1:do_idx,:))
     !
  ELSE
     !
     fname = TRIM(wbse_init_save_dir)//'/index_matrix'//flabel
     CALL wbse_index_matrix_read(fname,tmp_size,do_idx,2,idx_matrix)
     !
  ENDIF
  !
  DEALLOCATE(ovl_matrix)
  ALLOCATE(restart_matrix(do_idx))
  !
  restart_matrix(:) = 0
  !
  calc_is_done = .FALSE.
  IF(l_restart_calc) THEN
     fname = TRIM(wbse_init_save_dir)//'/restart_matrix'//flabel
     CALL wbse_status_restart_read(fname,do_idx,restart_matrix,calc_is_done)
  ENDIF
  !
  IF(.NOT. calc_is_done) THEN
     !
     ! initialize the paralellization
     !
     bandpair = idistribute()
     CALL bandpair%init(do_idx,'i','n_pairs',.TRUE.)
     !
     IF(l_pdep) THEN
        CALL io_push_title('Applying CHI kernel with PDEP')
     ELSE
        CALL io_push_title('Applying '//TRIM(chi_kernel)//' kernel with FF_Qbox')
     ENDIF
     !
     CALL start_bar_type(barra,'CHI',bandpair%nlocx)
     !
     ! parallel loop
     !
     DO il1 = 1,bandpair%nlocx
        !
        ig1  = bandpair%l2g(il1) ! global index of n_total
        !
        ibnd = idx_matrix(ig1,1)
        jbnd = idx_matrix(ig1,2)
        !
        l_skip = .FALSE.
        IF(l_restart_calc) THEN
           IF(restart_matrix(ig1) > 0) l_skip = .TRUE.
        ENDIF
        !
        IF(ig1 < 1 .OR. ig1 > do_idx) l_skip = .TRUE.
        !
        IF(.NOT. l_skip) THEN
           !
           ALLOCATE(tau(npwx))
           ALLOCATE(aux_r(dffts%nnr))
           ALLOCATE(aux1_g(npwx))
           !
           IF(l_local_repr) THEN
              CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
           ELSE
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
           ENDIF
           !
           aux_r(:) = CMPLX(REAL(psic,KIND=DP)*AIMAG(psic)/omega,KIND=DP)
           !
           ! aux_r -> aux1_g
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,aux_r,aux1_g,'Wave')
           !
           ! vc in fock like term
           !
           tau(:) = (0._DP,0._DP)
           DO ig = 1,npw
              tau(ig) = aux1_g(ig) * pot3D%sqvc(ig)**2
           ENDDO
           !
           IF(l_pdep) THEN
              !
              ALLOCATE(dotp(pert%nloc))
              !
              DO ig = 1,npw
                 aux1_g(ig) = aux1_g(ig) * pot3D%sqvc(ig)
              ENDDO
              !
              CALL glbrak_gamma(aux1_g,dvg,dotp,npw,npwx,1,pert%nloc,1,npol)
              CALL mp_sum(dotp,intra_bgrp_comm)
              !
              aux1_g(:) = (0._DP,0._DP)
              !
              DO ip = 1,pert%nloc
                 !
                 ip_g = pert%l2g(ip)
                 factor = dotp(ip)*ev(ip_g)/(1._DP-ev(ip_g))
                 !
                 DO ig = 1,npw
                    aux1_g(ig) = aux1_g(ig)+dvg(ig,ip)*factor
                 ENDDO
                 !
              ENDDO
              !
              CALL mp_sum(aux1_g,inter_bgrp_comm)
              !
              DO ig = 1,npw
                 tau(ig) = tau(ig) + aux1_g(ig)*pot3D%sqvc(ig)
              ENDDO
              !
              DEALLOCATE(dotp)
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
              CALL single_invfft_gamma(dffts,npw,npwx,aux1_g,aux_r,'Wave')
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
              aux_rr(:) = 0.5_DP * aux_rr ! change from rydberg to hartree
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
                 aux_rr(:) = frspin(:,1) + frspin(:,2)
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
              CALL single_fwfft_gamma(dffts,npw,npwx,aux_r,aux1_g,'Wave')
              !
              ! vc + vc/fxc X vc
              !
              tau(:) = tau + aux1_g
              !
              DEALLOCATE(aux1_r)
              DEALLOCATE(aux_rr)
              !
           ENDIF
           !
           ! write dvg vc_rho + vc_rho X vc_rho to disk
           !
           WRITE(ilabel,'(i6.6)') ibnd
           WRITE(jlabel,'(i6.6)') jbnd
           WRITE(slabel,'(i1)') current_spin
           !
           fname = TRIM(wbse_init_save_dir)//'/E'//ilabel//'_'//jlabel//'_'//slabel//'.dat'
           CALL pdep_merge_and_write_G(fname,tau)
           !
           DEALLOCATE(tau)
           DEALLOCATE(aux_r)
           DEALLOCATE(aux1_g)
           !
           restart_matrix(ig1) = 1
           !
        ENDIF
        !
        ! for restarting, update status of restart_matrix
        !
        CALL mp_sum(restart_matrix(1:do_idx),inter_image_comm)
        !
        DO ir = 1,do_idx
           IF(restart_matrix(ir) > 0) restart_matrix(ir) = 1
        ENDDO
        !
        fname = TRIM(wbse_init_save_dir)//'/restart_matrix'//flabel
        CALL wbse_status_restart_write(fname,do_idx,restart_matrix)
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
        CALL update_bar_type(barra,'CHI',1)
        !
     ENDDO
     !
     fname = TRIM(wbse_init_save_dir)//'/restart_matrix'//flabel
     CALL wbse_status_restart_write(fname,do_idx,restart_matrix)
     !
     CALL stop_bar_type(barra,'CHI')
     !
  ENDIF
  !
  DEALLOCATE(idx_matrix)
  DEALLOCATE(restart_matrix)
  IF(l_local_repr) DEALLOCATE(evc_loc)
  !
END SUBROUTINE
