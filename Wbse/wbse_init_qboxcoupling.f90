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
SUBROUTINE wbse_init_qboxcoupling_single_q(iks,ikq,current_spin,nbndval,l_restart_calc)
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE io_push,              ONLY : io_push_title
  USE types_coulomb,        ONLY : pot3D
  USE westcom,              ONLY : wbse_init_save_dir,chi_kernel,l_xcchi,l_use_localise_repr,&
                                 & overlap_thr
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc,psic
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : igk_k,npw,npwx,lsda
  USE pdep_io,              ONLY : pdep_merge_and_write_G
  USE wbse_init_restart,    ONLY : wbse_status_restart_read,wbse_status_restart_write,&
                                 & wbse_index_matrix_read,wbse_index_matrix_write
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : aband,bseparal
  USE function3d,           ONLY : write_function3d,read_function3d
  USE mp,                   ONLY : mp_barrier,mp_sum
  USE lsda_mod,             ONLY : nspin
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE mp_global,            ONLY : inter_image_comm,intra_image_comm,my_image_id,me_bgrp
  USE conversions,          ONLY : itoa
  USE qbox_interface,       ONLY : sleep_and_wait_for_lock_to_be_removed
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iks,ikq,current_spin,nbndval
  LOGICAL, INTENT(IN) :: l_restart_calc
  !
  INTEGER :: ibnd, jbnd, tmp_size
  INTEGER :: il1, ig1, ir, do_index
  REAL(DP):: ovl_value
  REAL(DP), ALLOCATABLE :: rho_aux(:)
  REAL(DP), ALLOCATABLE :: ovl_matrix(:,:)
  INTEGER, ALLOCATABLE :: restart_matrix(:)
  INTEGER, ALLOCATABLE :: index_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_loc(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvg(:), psic_aux(:)
  !
  INTEGER :: iu, ig, stat
  CHARACTER(LEN=:), ALLOCATABLE :: lockfile
  CHARACTER(LEN=:), ALLOCATABLE :: fname
  REAL(DP), ALLOCATABLE :: aux_rr(:)
  COMPLEX(DP), ALLOCATABLE :: aux_r(:),aux1_r(:,:),aux1_g(:)
  REAL(DP), ALLOCATABLE :: frspin(:,:)
  !
  CHARACTER :: my_spin
  CHARACTER(LEN=6) :: my_labeliq, my_labelik
  CHARACTER(LEN=6) :: my_label1, my_label2
  CHARACTER(LEN=256) :: kernel, driver
  !
  REAL(DP), EXTERNAL :: get_clock
  CHARACTER(20), EXTERNAL :: human_readable_time
  !
  LOGICAL :: calc_is_done
  LOGICAL :: l_skip
  REAL(DP), EXTERNAL :: DDOT
  !
  TYPE(bar_type) :: barra
  !
  CALL wbse_init_memory_report()
  !
  CALL start_clock('wbse_qbox_coupling')
  !
  IF(chi_kernel == 'XC_CHI') THEN
     kernel  = 'CHI'
     l_xcchi = .TRUE.
  ELSEIF(chi_kernel == 'XC_CHI_RPA') THEN
     kernel  = 'CHI_RPA'
     l_xcchi = .TRUE.
  ELSE
     kernel  = chi_kernel
     l_xcchi = .FALSE.
  ENDIF
  !
  driver = 'FF_Qbox'
  !
  WRITE(my_labeliq,'(i6.6)') ikq
  WRITE(my_labelik,'(i6.6)') iks
  WRITE(my_spin,'(i1)') current_spin
  !
  CALL io_push_title('Wbse_init for '//TRIM(kernel)//' with '//TRIM(driver)//' driver'//&
       & ' ik'//my_labelik//' iq'//my_labeliq//' ispin'//my_spin)
  !
  aband = idistribute()
  CALL aband%init(nbndval,'i','bse_nbndval',.TRUE.)
  !
  IF(.NOT. ALLOCATED(psic)) ALLOCATE(psic(dffts%nnr))
  IF(.NOT. gamma_only) ALLOCATE(psic_aux(dffts%nnr))
  !
  tmp_size = nbndval*nbndval
  ALLOCATE(index_matrix(tmp_size,2))
  ALLOCATE(ovl_matrix(nbndval,nbndval))
  !
  index_matrix(:,:) = 0
  !
  IF(l_use_localise_repr) THEN
     ALLOCATE(evc_loc(npwx,nbndval))
     CALL wbse_localization(current_spin, nbndval, evc_loc, ovl_matrix, l_restart_calc)
  ENDIF
  !
  ! compute index_matrix
  !
  do_index = 0
  !
  IF(.NOT. l_restart_calc) THEN
     DO ibnd = 1, nbndval
        DO jbnd = 1, nbndval
           !
           ovl_value = ovl_matrix(ibnd,jbnd)
           !
           IF(ovl_value >= overlap_thr) THEN
              IF(gamma_only) THEN
                 IF(jbnd >= ibnd) THEN
                    do_index = do_index + 1
                    !
                    index_matrix(do_index, 1) = ibnd
                    index_matrix(do_index, 2) = jbnd
                 ENDIF
              ELSE
                 do_index = do_index + 1
                 !
                 index_matrix(do_index, 1) = ibnd
                 index_matrix(do_index, 2) = jbnd
              ENDIF
           ENDIF
           !
        ENDDO
     ENDDO
     !
     fname = TRIM(wbse_init_save_dir)//'/index_matrix_iq'//my_labeliq//'_ik'//&
             & my_labelik//'_spin'//my_spin//'.dat'
     CALL wbse_index_matrix_write(fname,do_index,2,index_matrix(1:do_index,:))
  ELSE
     fname = TRIM(wbse_init_save_dir)//'/index_matrix_iq'//my_labeliq//'_ik'//&
             & my_labelik//'_spin'//my_spin//'.dat'
     CALL wbse_index_matrix_read(fname,tmp_size,do_index,2,index_matrix)
  ENDIF
  !
  ALLOCATE(restart_matrix(do_index))
  !
  restart_matrix(:) = 0
  !
  calc_is_done = .FALSE.
  IF(l_restart_calc) THEN
     fname = TRIM(wbse_init_save_dir)//'/restart_matrix_iq'//my_labeliq//'_ik'//&
             & my_labelik//'_spin'//my_spin//'.dat'
     CALL wbse_status_restart_read(fname,do_index,restart_matrix,calc_is_done)
  ENDIF
  !
  IF(.NOT. calc_is_done) THEN
     !
     ! initialize the paralellization
     !
     bseparal = idistribute()
     CALL bseparal%init(do_index,'i','number_pairs', .TRUE.)
     !
     CALL io_push_title('Applying ' // TRIM(kernel) // ' kernel with FF_Qbox ...')
     !
     CALL start_bar_type(barra,'FF_Qbox',bseparal%nlocx)
     !
     ! parallel loop
     !
     DO il1 = 1, bseparal%nlocx
        !
        ig1  = bseparal%l2g(il1) ! global index of n_total
        !
        ibnd = index_matrix(ig1,1)
        jbnd = index_matrix(ig1,2)
        !
        l_skip = .FALSE.
        IF(l_restart_calc) THEN
           IF(restart_matrix(ig1) > 0) l_skip = .TRUE.
        ENDIF
        !
        IF(ig1 < 1 .OR. ig1 > do_index) l_skip = .TRUE.
        !
        IF(.NOT. l_skip) THEN
           !
           ALLOCATE(rho_aux(dffts%nnr))
           ALLOCATE(dvg(npwx))
           !
           IF(gamma_only) THEN
              IF(l_use_localise_repr) THEN
                 CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
              ELSE
                 CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
              ENDIF
              !
              rho_aux(:) = REAL(psic,KIND=DP) * AIMAG(psic)
           ELSE
              IF(l_use_localise_repr) THEN
                 CALL single_invfft_k(dffts,npw,npwx,evc_loc(:,ibnd),psic,'Wave',igk_k(:,1)) ! only 1 kpoint
                 CALL single_invfft_k(dffts,npw,npwx,evc_loc(:,jbnd),psic_aux,'Wave',igk_k(:,1)) ! only 1 kpoint
              ELSE
                 CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,1)) ! only 1 kpoint
                 CALL single_invfft_k(dffts,npw,npwx,evc(:,jbnd),psic_aux,'Wave',igk_k(:,1)) ! only 1 kpoint
              ENDIF
              !
              rho_aux(:) = REAL(CONJG(psic)*psic_aux,KIND=DP)
           ENDIF
           !
           rho_aux(:) = rho_aux/omega
           !
           ALLOCATE(aux1_g(npwx))
           ALLOCATE(aux_r(dffts%nnr))
           ALLOCATE(aux1_r(dffts%nnr,nspin))
           ALLOCATE(aux_rr(dffts%nnr))
           !
           aux_r(:) = CMPLX(rho_aux,KIND=DP)
           !
           ! aux_r -> aux1_g
           !
           IF(gamma_only) THEN
              CALL single_fwfft_gamma(dffts,npw,npwx,aux_r,aux1_g,'Wave')
           ELSE
              CALL single_fwfft_k(dffts,npw,npwx,aux_r,aux1_g,'Wave')
           ENDIF
           !
           ! vc in fock like term
           !
           dvg(:) = (0._DP,0._DP)
           DO ig = 1, npw
              dvg(ig) = aux1_g(ig) * pot3D%sqvc(ig) * pot3D%sqvc(ig)
           ENDDO
           !
           ! vc in correlation like term G->0 = 0
           !
           ! aux1_g -> aux_r
           !
           IF(gamma_only) THEN
              CALL single_invfft_gamma(dffts,npw,npwx,aux1_g,aux_r,'Wave')
           ELSE
              CALL single_invfft_k(dffts,npw,npwx,aux1_g,aux_r,'Wave')
           ENDIF
           !
           aux1_r(:,:) = (0._DP,0._DP)
           aux1_r(:,current_spin) = aux_r
           !
           ! aux1_r = vc*aux1_r()
           !
           CALL wbse_dv_of_drho(aux1_r, .TRUE., .FALSE.)
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
              CALL sleep_and_wait_for_lock_to_be_removed(lockfile, '["response"]')
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
           DO ir = 1, dffts%nnr
              aux_r(ir) = CMPLX(aux_rr(ir)*SQRT(omega), KIND=DP) ! rescale response
           ENDDO
           !
           aux1_r(:,:) = (0._DP,0._DP)
           aux1_r(:,current_spin) = aux_r
           !
           ! aux1_r = vc*aux1_r()
           !
           IF(l_xcchi) THEN
              CALL wbse_dv_of_drho(aux1_r, .FALSE., .FALSE.)
           ELSE
              CALL wbse_dv_of_drho(aux1_r, .TRUE., .FALSE.)
           ENDIF
           !
           aux_r(:) = aux1_r(:,current_spin)
           !
           ! aux_r -> aux_g
           !
           aux1_g(:) = (0._DP, 0._DP)
           !
           IF(gamma_only) THEN
              CALL single_fwfft_gamma(dffts,npw,npwx,aux_r,aux1_g,'Wave')
           ELSE
              CALL single_fwfft_k(dffts,npw,npwx,aux_r,aux1_g,'Wave')
           ENDIF
           !
           ! vc + vc/fxc X vc
           !
           dvg(:) = dvg + aux1_g
           !
           ! write dvg vc_rho + vc_rho X vc_rho to disk
           !
           WRITE(my_label1,'(i6.6)') ibnd
           WRITE(my_label2,'(i6.6)') jbnd
           WRITE(my_spin,'(i1)') current_spin
           !
           fname = TRIM(wbse_init_save_dir)//'/E'//my_label1//'_'//my_label2//'_'//my_spin//'.dat'
           CALL pdep_merge_and_write_G(fname,dvg)
           !
           DEALLOCATE(rho_aux, dvg)
           DEALLOCATE(aux1_g)
           DEALLOCATE(aux_r, aux1_r, aux_rr)
           !
           restart_matrix(ig1) = 1
           !
        ENDIF
        !
        ! for restarting, update status of restart_matrix
        !
        CALL mp_sum(restart_matrix(1:do_index), inter_image_comm)
        !
        DO ir = 1, do_index
           IF(restart_matrix(ir) > 0) restart_matrix(ir) = 1
        ENDDO
        !
        fname = TRIM(wbse_init_save_dir)//'/restart_matrix_iq'//my_labeliq//'_ik'//&
                & my_labelik//'_spin'//my_spin//'.dat'
        CALL wbse_status_restart_write(fname,do_index,restart_matrix)
        !
        ! clean up
        !
        IF(me_bgrp == 0) THEN
           fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml'
           OPEN(NEWUNIT=iu, IOSTAT=stat, FILE=fname, STATUS='OLD')
           IF(stat == 0) CLOSE(iu, STATUS='DELETE')
           IF(lsda) THEN
              fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response.spin0'
              OPEN(NEWUNIT=iu, IOSTAT=stat, FILE=fname, STATUS='OLD')
              IF(stat == 0) CLOSE(iu, STATUS='DELETE')
              fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response.spin1'
              OPEN(NEWUNIT=iu, IOSTAT=stat, FILE=fname, STATUS='OLD')
              IF(stat == 0) CLOSE(iu, STATUS='DELETE')
           ELSE
              fname = 'I.'//itoa(my_image_id)//'_P.'//itoa(il1)//'.xml.response'
              OPEN(NEWUNIT=iu, IOSTAT=stat, FILE=fname, STATUS='OLD')
              IF(stat == 0) CLOSE(iu, STATUS='DELETE')
           ENDIF
        ENDIF
        !
        CALL update_bar_type(barra,'FF_Qbox',1)
        !
     ENDDO
     !
     fname = TRIM(wbse_init_save_dir)//'/restart_matrix_iq'//my_labeliq//'_ik'//&
             & my_labelik//'_spin'//my_spin//'.dat'
     CALL wbse_status_restart_write(fname,do_index,restart_matrix)
     !
     CALL stop_bar_type(barra,'FF_Qbox')
     !
  ENDIF
  !
  DEALLOCATE(index_matrix)
  IF(ALLOCATED(restart_matrix)) DEALLOCATE(restart_matrix)
  DEALLOCATE(ovl_matrix)
  IF(ALLOCATED(psic)) DEALLOCATE(psic)
  IF(ALLOCATED(psic_aux)) DEALLOCATE(psic_aux)
  IF(ALLOCATED(evc_loc)) DEALLOCATE(evc_loc)
  !
  CALL stop_clock('wbse_qbox_coupling')
  !
END SUBROUTINE
