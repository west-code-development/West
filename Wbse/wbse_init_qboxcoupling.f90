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
SUBROUTINE wbse_init_qboxcoupling_single_q (iks,ikq,xq,current_spin,nbndval,l_restart_calc)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : e2, fpi
  USE cell_base,            ONLY : alat, tpiba2, omega
  USE gvect,                ONLY : nl,ngm,g,nlm,gstart
  USE io_push,              ONLY : io_push_title
   !sqvc not in westcom pot3D%sqvc  TODO: pot3d init
  USE types_coulomb,         ONLY : pot3D
  !USE westcom,              ONLY : wstat_save_dir,sqvc,fftdriver,npwq,npwqx
  USE westcom,              ONLY : wstat_save_dir,fftdriver,chi_driver, chi_kernel,l_xcchi
  !USE westcom,              ONLY : wstat_save_dir,sqvc,fftdriver,npwq0,npwq0x,chi_driver
  !wbsecom combined with westcom
  !USE wbsecom,              ONLY : chi_kernel,l_xcchi
  USE pwcom,                ONLY : omega
  USE gvect,                ONLY : gstart,ngm,ngmx
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc,psic
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE pwcom,                ONLY : wk,nks,nelup,neldw,isk,g,igk_k,ngm,tpiba2,xk,omega,npw,npwx,lsda,nkstot,&
                                 & current_k,ngk,nbnd,wg
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_global,            ONLY : inter_image_comm
  USE mp,                   ONLY : mp_sum
  USE pdep_io,              ONLY : pdep_merge_and_write_G
  USE bse_module,           ONLY : ovl_thr, l_wannier_repr
  USE wbse_init_restart,    ONLY : wbse_stat_restart_read, wbse_stat_restart_write
  USE wbse_init_restart,    ONLY : wbse_index_matrix_read, wbse_index_matrix_write
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : aband, bseparal
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN) :: iks,ikq,current_spin,nbndval
  REAL(DP),    INTENT(IN) :: xq(3)
  LOGICAL,     INTENT(IN) :: l_restart_calc
  !
  INTEGER :: ibnd, jbnd, tmp_size
  INTEGER :: il1, ig1, ir, do_index
  REAL(DP):: ovl_value
  REAL(DP),    ALLOCATABLE  :: rho_aux(:)
  REAL(DP),    ALLOCATABLE  :: ovl_matrix(:,:)
  REAL(DP),    ALLOCATABLE  :: restart_matrix(:)
  REAL(DP),    ALLOCATABLE  :: index_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE  :: evc_loc(:,:)
  COMPLEX(DP), ALLOCATABLE  :: dvg(:), psic_aux(:)
  !
  CHARACTER(LEN=1)          :: my_spin
  CHARACTER(LEN=6)          :: my_labeliq, my_labelik
  CHARACTER(LEN=6)          :: my_label1, my_label2
  CHARACTER(LEN=256)        :: kernel, driver
  CHARACTER(LEN=256)        :: filename
  !
  REAL(DP), EXTERNAL        :: get_clock
  CHARACTER(20),EXTERNAL    :: human_readable_time
  REAL(DP)                  :: wtime(2)
  !
  LOGICAL                   :: calc_is_done
  REAL(kind=dp), EXTERNAL   :: DDOT
  !
  CALL wbse_init_memory_report()
  !
  CALL start_clock( 'wbse_qbox_coupling' )
  !
  IF (chi_kernel == 'XC_CHI') THEN
     !
     kernel  = "CHI"
     l_xcchi = .true.
     !
  ELSEIF (chi_kernel == 'XC_CHI_RPA') THEN
     !
     kernel  = "CHI_RPA"
     l_xcchi = .true.
     !
  ELSE
     !
     kernel  = chi_kernel
     l_xcchi = .false.
     !
  ENDIF
  !
  driver = "FF_QBOX"
  !
  WRITE(my_labeliq,'(i6.6)') ikq
  WRITE(my_labelik,'(i6.6)') iks
  WRITE(my_spin,'(i1)') current_spin
  !
  CALL io_push_title( "Wbse_init for " // TRIM(kernel) // &
               & " with " // TRIM(driver) // " driver" // &
               & " ik"// TRIM(my_labelik) // &
               & " iq"// TRIM(my_labeliq) // &
               & " ispin"// TRIM(my_spin) )
  !
  aband = idistribute()
  CALL aband%init(nbndval,'i','bse_nbndval',.TRUE.)
  !
  IF (.NOT.ALLOCATED(psic)) ALLOCATE (psic(dffts%nnr))
  IF (.NOT.gamma_only) ALLOCATE (psic_aux(dffts%nnr))
  !
  tmp_size = nbndval*nbndval
  ALLOCATE (index_matrix(tmp_size,2))
  ALLOCATE (ovl_matrix(nbndval,nbndval))
  !
  ovl_matrix(:,:)   = 0.0_DP
  index_matrix(:,:) = 0.0_DP
  !
  IF (l_wannier_repr) THEN
     !
     ALLOCATE (evc_loc(npwx,nbndval))
     !
     CALL bse_do_localization (current_spin, nbndval, evc_loc, ovl_matrix, l_restart_calc)
     !
  ENDIF
  !
  ! compute index_matrix
  !
  do_index = 0
  !
  IF (.not. l_restart_calc) THEN
     !
     DO ibnd = 1, nbndval
        !
        DO jbnd = 1, nbndval
           !
           ovl_value = ovl_matrix(ibnd,jbnd)
           !
           IF (ovl_value >= ovl_thr ) THEN
              !
              IF (gamma_only) THEN
                 !
                 IF (jbnd >= ibnd)  THEN
                    !
                    do_index = do_index + 1
                    !
                    index_matrix(do_index, 1) = ibnd
                    index_matrix(do_index, 2) = jbnd
                    !
                 ENDIF
                 !
              ELSE
                 !
                 do_index = do_index + 1
                 !
                 index_matrix(do_index, 1) = ibnd
                 index_matrix(do_index, 2) = jbnd
                 !
              ENDIF
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
     filename = TRIM( wstat_save_dir )//"/index_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
     CALL wbse_index_matrix_write(filename,do_index,2,index_matrix(1:do_index,:))
     !
  ELSE
     !
     filename = TRIM( wstat_save_dir )//"/index_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
     CALL wbse_index_matrix_read (filename,tmp_size,do_index,2,index_matrix)
     !
  ENDIF
  !
  ALLOCATE (restart_matrix(do_index))
  !
  restart_matrix(:) = 0.0_DP
  !
  calc_is_done = .FALSE.
  IF (l_restart_calc) THEN
     !
     filename = TRIM( wstat_save_dir )//"/restart_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
     CALL wbse_stat_restart_read (filename,do_index,restart_matrix,calc_is_done)
     !
  ENDIF
  !
  IF (calc_is_done) GOTO 2222
  !
  ! initialize the paralellization
  !
  bseparal = idistribute()
  CALL bseparal%init( do_index,'i','number_pairs', .TRUE.)
  !
  ! parallel loop
  !
  DO il1 = 1, bseparal%nlocx
     !
     ig1  = bseparal%l2g(il1) ! global index of n_total
     !
     ibnd = INT(index_matrix(ig1,1))
     jbnd = INT(index_matrix(ig1,2))
     !
     IF (l_restart_calc) THEN
        !
        IF (INT(restart_matrix(ig1)) > 0) GOTO 1111
        !
     ENDIF
     !
     IF ((ig1 < 1).or.(ig1 > do_index)) GOTO 1111
     !
     ! Here is the main part of this code
     !
     ALLOCATE (rho_aux(dffts%nnr))
     ALLOCATE (dvg(ngmx))
     !
     IF (gamma_only) THEN
        !
        IF (l_wannier_repr) THEN
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(1,ibnd),evc_loc(1,jbnd), psic,'Wave')
           !
        ELSE
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,jbnd), psic,'Wave')
           !
        ENDIF
        !
        rho_aux(:) = DBLE(psic(:)) * AIMAG(psic(:))
        !
     ELSE
        !
        IF (l_wannier_repr) THEN
           !
           CALL single_invfft_k(dffts,npw,npwx,evc_loc(1,ibnd),psic,'Wave',igk_k(1,1)) !only 1 kpoint
           CALL single_invfft_k(dffts,npw,npwx,evc_loc(1,jbnd),psic_aux,'Wave',igk_k(1,1)) !only 1 kpoint
           !
        ELSE
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,1)) !only 1 kpoint
           CALL single_invfft_k(dffts,npw,npwx,evc(1,jbnd),psic_aux,'Wave',igk_k(1,1)) !only 1 kpoint
           !
        ENDIF
        !
        rho_aux(:) = DBLE(CONJG(psic(:)) * psic_aux(:))
        !
     ENDIF
     !
     rho_aux(:) = rho_aux(:)/omega
     !
     CALL couple_with_qbox_routine (kernel, current_spin, rho_aux, dvg)
     !
     ! write dvg vc_rho + vc_rho X vc_rho to disk
     !
     WRITE(my_label1,'(i6.6)') ibnd
     WRITE(my_label2,'(i6.6)') jbnd
     WRITE(my_spin,'(i1)') current_spin
     !
     filename = TRIM( wstat_save_dir )//"/E"//TRIM(ADJUSTL(my_label1))//"_"//&
             TRIM(ADJUSTL(my_label2))//"_"//TRIM(ADJUSTL(my_spin))//".dat"
     CALL pdep_merge_and_write_G(filename,dvg(:))
     !
     DEALLOCATE(rho_aux, dvg)
     !
     restart_matrix(ig1)  = 1.0
     !
1111 CONTINUE
     !
     ! for restarting, update status of restart_matrix
     !
     CALL mp_sum(restart_matrix(1:do_index), inter_image_comm)
     !
     DO ir = 1, do_index
        !
        IF (restart_matrix(ir) > 0) restart_matrix(ir) = 1.0
        !
     ENDDO
     !
     calc_is_done = .FALSE.
     filename = TRIM( wstat_save_dir )//"/restart_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
     CALL wbse_stat_restart_write (filename,do_index,restart_matrix,calc_is_done)
     !
  ENDDO
  !
  calc_is_done = .TRUE.
  filename = TRIM( wstat_save_dir )//"/restart_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
             TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
  CALL wbse_stat_restart_write (filename,do_index,restart_matrix,calc_is_done)
  !
2222 CONTINUE
  !
  DEALLOCATE (index_matrix)
  IF (ALLOCATED(restart_matrix))   DEALLOCATE (restart_matrix)
  DEALLOCATE (ovl_matrix)
  !
  IF (ALLOCATED(psic))    DEALLOCATE(psic)
  IF (ALLOCATED(psic_aux))DEALLOCATE(psic_aux)
  IF (ALLOCATED(evc_loc)) DEALLOCATE(evc_loc)
  !
  CALL stop_clock( 'wbse_qbox_coupling' )
  !
  RETURN
  !
ENDSUBROUTINE
!
!
SUBROUTINE couple_with_qbox_routine (kernel, current_spin, dnr, dvg)
  !
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE lsda_mod,              ONLY : nspin
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE westcom,               ONLY : fftdriver, l_xcchi
  !old west has sqvc, new version use pot3D
  !USE westcom,               ONLY : fftdriver, sqvc
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,              ONLY : single_fwfft_k,single_invfft_k
  USE gvect,                 ONLY : gstart,g,ngm,ngmx
  USE constants,             ONLY : fpi, e2
  USE cell_base,             ONLY : tpiba2, omega
  USE mp,                    ONLY : mp_barrier, mp_sum
  USE mp_global,             ONLY : intra_image_comm
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE pwcom,                 ONLY : omega
  USE control_flags,         ONLY : gamma_only
  USE io_push,               ONLY : io_push_title
  USE martyna_tuckerman,     ONLY : wg_corr_h, do_comp_mt
  USE qbox_interface,        ONLY : apply_kernel_by_qbox
  !wbsecom combined with westcom
  !USE wbsecom,               ONLY : l_xcchi
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  CHARACTER(LEN=256), INTENT(IN)  :: kernel
  INTEGER,            INTENT(IN)  :: current_spin
  REAL(DP),           INTENT(IN)  :: dnr(dffts%nnr)
  COMPLEX(DP),        INTENT(INOUT) :: dvg(ngmx)
  !
  ! Workspace
  !
  INTEGER :: ig, ir
  !
  REAL(DP),   ALLOCATABLE :: aux_rr(:)
  COMPLEX(DP),ALLOCATABLE :: aux_r(:),aux1_r(:,:),aux1_g(:)
  COMPLEX(DP),ALLOCATABLE :: dvaux_mt(:)
  !
  REAL(DP), EXTERNAL      :: get_clock
  CHARACTER(20),EXTERNAL  :: human_readable_time
  REAL(DP)                :: rtime(2)
  REAL(DP)                :: itime(2), otime(2), wtime(2)
  REAL(DP)                :: gnorm2, eh_corr
  !
  TYPE(bar_type) :: barra
  !
  CALL io_push_title("Applying " // TRIM(kernel) // " kernel with FF_QBOX ... ")
  !
  CALL mp_barrier( intra_image_comm )
  !
  CALL start_clock( 'kernel_r' )
  !
  rtime(1) = MAX(0.0_DP, get_clock('kernel_r'))
  itime(1) = MAX(0.0_DP, get_clock('read_resp'))
  otime(1) = MAX(0.0_DP, get_clock('write_vext'))
  wtime(1) = MAX(0.0_DP, get_clock('wait_qbox'))
  !
  CALL start_bar_type( barra, 'apply_kernel_r', 1 )
  !
  ALLOCATE( aux1_g(ngmx))
  ALLOCATE( aux_r(dffts%nnr) )
  ALLOCATE( aux1_r(dffts%nnr,nspin) )
  ALLOCATE( aux_rr(dffts%nnr) )
  !
  CALL report_dynamical_memory()
  !
  aux_r(:) = CMPLX(dnr(:), 0.0_DP)
  !
  ! aux_r -> aux1_g
  !
  IF (gamma_only) THEN
     !
     CALL single_fwfft_gamma(dffts,ngm,ngmx,aux_r,aux1_g,'Dense')
     !
  ELSE
     !
     CALL single_fwfft_k(dffts,ngm,ngmx,aux_r,aux1_g,'Dense')
     !
  ENDIF
  !
  ! vc in fock like term
  !
  dvg(:) = (0.0_DP, 0.0_DP)
  DO ig = 1, ngm
     !
     dvg(ig) = aux1_g(ig) * pot3D%sqvc(ig) * pot3D%sqvc(ig)
     !
  ENDDO
  !
  ! vc in correlation like term G->0 = 0
  !
  ! aux1_g -> aux_r
  !
  IF (gamma_only) THEN
     !
     CALL single_invfft_gamma(dffts,ngm,ngmx,aux1_g,aux_r,'Dense')
     !
  ELSE
     !
     CALL single_invfft_k(dffts,ngm,ngmx,aux1_g,aux_r,'Dense')
     !
  ENDIF
  !
  aux1_r(:,:) = (0.0_DP,0.0_DP)
  aux1_r(:,current_spin) = aux_r(:)
  !
  ! aux1_r = vc*aux1_r()
  !
  CALL west_dv_of_drho(aux1_r, .true., .false.)
  !
  aux_r(:) = aux1_r(:,current_spin)
  !
  ! Send data to QBOX to compute X|vc rho>
  !
  aux_rr(:) = DBLE(aux_r(:))/DSQRT(omega) ! scale down pert.
  !
  CALL apply_kernel_by_qbox(kernel, aux_rr)
  !
  DO ir = 1, dffts%nnr
     !
     aux_r(ir) = CMPLX( aux_rr(ir)*DSQRT(omega), 0._DP, KIND=DP) ! rescale response
     !
  ENDDO
  !
  aux1_r(:,:) = (0.0_DP,0.0_DP)
  aux1_r(:,current_spin) = aux_r(:)
  !
  ! aux1_r = vc*aux1_r()
  !
  IF (l_xcchi) THEN
     !
     CALL west_dv_of_drho(aux1_r, .false., .false.)
     !
  ELSE
     !
     CALL west_dv_of_drho(aux1_r, .true., .false.)
     !
  ENDIF
  !
  aux_r(:) = aux1_r(:,current_spin)
  !
  ! aux_r -> aux_g
  !
  aux1_g(:) = (0.0_DP, 0.0_DP)
  !
  IF (gamma_only) THEN
     !
     CALL single_fwfft_gamma(dffts,ngm,ngmx,aux_r,aux1_g,'Dense')
     !
  ELSE
     !
     CALL single_fwfft_k(dffts,ngm,ngmx,aux_r,aux1_g,'Dense')
     !
  ENDIF
  !
  ! vc + vc/fxc X vc
  !
  dvg(:) = dvg(:) + aux1_g(:)
  !
  DEALLOCATE (aux1_g)
  DEALLOCATE (aux_r, aux1_r, aux_rr)
  !
  CALL mp_barrier( intra_image_comm )
  !
  CALL stop_bar_type( barra, 'apply_kernel_r' )
  !
  CALL stop_clock( 'kernel_r' )
  !
  rtime(2) = get_clock( 'kernel_r' )
  !
  WRITE(stdout, "(5x,'apply_kernel_to_a_vectors_r run time ',a)") &
                & human_readable_time(rtime(2) - rtime(1) )
  !
  itime(2) = get_clock( 'read_resp' )
  otime(2) = get_clock( 'write_vext' )
  wtime(2) = get_clock( 'wait_qbox' )
  !
  WRITE(stdout, "(5x,'qbox interface timing:')")
  WRITE(stdout, "(7x,'read ',a,'write ',a,'wait ',a)") &
                   & human_readable_time(itime(2) - itime(1)), &
                   & human_readable_time(otime(2) - otime(1)), &
                   & human_readable_time(wtime(2) - wtime(1))
  !
  RETURN
  !
END SUBROUTINE couple_with_qbox_routine
