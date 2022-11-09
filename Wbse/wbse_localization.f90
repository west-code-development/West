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
SUBROUTINE wbse_localization(current_spin, nbndval, evc_loc, ovl_matrix, l_restart)
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : aband
  USE westcom,              ONLY : l_use_bisection_thr,wfc_from_qbox,wbse_init_save_dir
  USE wavefunctions,        ONLY : evc,psic
  USE plep_io,              ONLY : plep_merge_and_write_G,plep_read_G_and_distribute
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k
  USE pwcom,                ONLY : igk_k,npw,npwx
  USE gvect,                ONLY : gstart
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE qbox_interface,       ONLY : load_qbox_wfc
  USE check_ovl_wfc,        ONLY : check_ovl_wannier,read_bisection_loc,check_ovl_bisection
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: current_spin, nbndval
  REAL(DP),INTENT(OUT) :: ovl_matrix(nbndval,nbndval)
  COMPLEX(DP), INTENT(OUT) :: evc_loc(npwx,nbndval)
  LOGICAL, INTENT(IN) :: l_restart
  !
  INTEGER :: ibnd, jbnd
  INTEGER :: il1, bisec_i, bisec_j
  INTEGER, ALLOCATABLE :: bisec_loc(:)
  REAL(DP):: ovl_value
  COMPLEX(DP), ALLOCATABLE :: psic1(:)
  COMPLEX(DP), ALLOCATABLE :: u_matrix(:,:)
  CHARACTER(LEN=256) :: fname
  CHARACTER :: labels
  LOGICAL :: l_load_west_loc_wfc
  LOGICAL :: l_load_qbox_loc_wfc
  REAL(KIND=DP), EXTERNAL :: DDOT
  !
  CALL start_clock('local')
  !
  IF(.NOT. gamma_only) ALLOCATE(psic1(dffts%nnr))
  !
  l_load_west_loc_wfc = .FALSE.
  l_load_qbox_loc_wfc = .TRUE.
  !
  IF(.NOT. l_restart) THEN
     !
     ALLOCATE(u_matrix(nbndval,nbndval))
     !
     ! Use wannier as representation
     !
     IF(l_load_qbox_loc_wfc) CALL load_qbox_wfc(evc_loc, wfc_from_qbox, current_spin, nbndval)
!     IF(l_load_west_loc_wfc) CALL read_pwscf_wannier_orbs(nbndval, npwx, evc_loc, wfc_from_qbox)
     !
     WRITE(labels,'(i1)') current_spin
     fname = TRIM(wbse_init_save_dir)//'/evc_loc.'//labels//'.dat'
     CALL plep_merge_and_write_G(fname,evc_loc,nbndval)
     !
     ! Compute unitary rotation matrix
     !
     u_matrix(:,:) = (0._DP,0._DP)
     !
     DO ibnd = 1, nbndval
        DO jbnd = 1, nbndval
           IF(gamma_only) THEN
              u_matrix(ibnd,jbnd) = 2.0_DP*DDOT(2*npwx,evc(:,ibnd),1,evc_loc(:,jbnd),1)
              !
              IF(gstart == 2) THEN
                 u_matrix(ibnd,jbnd) = u_matrix(ibnd,jbnd) &
                 & - CMPLX(REAL(evc(1,ibnd),KIND=DP)*REAL(evc_loc(1,jbnd),KIND=DP),KIND=DP)
              ENDIF
           ELSE
              u_matrix(ibnd,jbnd) = DDOT(2*npwx,evc(:,ibnd),1,evc_loc(:,jbnd),1)
           ENDIF
        ENDDO
     ENDDO
     !
     CALL mp_sum(u_matrix,intra_bgrp_comm)
     !
     ! Compute overlap matrix
     !
     ovl_matrix(:,:) = 0._DP
     !
     IF(l_use_bisection_thr) THEN
        !
        ALLOCATE(bisec_loc(nbndval))
        !
        ! read bisection localization from file
        !
        CALL read_bisection_loc(current_spin, nbndval, bisec_loc)
        !
        ! compute orbital overlap matrix, using besection indication
        !
        DO ibnd = 1, nbndval
           DO jbnd = 1, nbndval
              bisec_i = bisec_loc(ibnd)
              bisec_j = bisec_loc(jbnd)
              !
              CALL check_ovl_bisection(bisec_i, bisec_j, ovl_value)
              !
              ovl_matrix(ibnd,jbnd) = ovl_value
           ENDDO
        ENDDO
        !
        DEALLOCATE (bisec_loc)
        !
     ELSE
        !
        ! compute orbital overlap matrix, using our method (see paper)
        !
        DO il1 = 1, aband%nloc
           !
           ibnd = aband%l2g(il1) ! global index of n_total
           !
           DO jbnd = 1, nbndval
              !
              IF(gamma_only) THEN
                 CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
                 CALL check_ovl_wannier(REAL(psic,KIND=DP),AIMAG(psic),ovl_value)
              ELSE
                 CALL single_invfft_k(dffts,npw,npwx,evc_loc(:,ibnd),psic,'Wave',igk_k(:,1)) ! only 1 kpoint
                 CALL single_invfft_k(dffts,npw,npwx,evc_loc(:,jbnd),psic1,'Wave',igk_k(:,1)) ! only 1 kpoint
                 CALL check_ovl_wannier(psic,psic1,ovl_value)
              ENDIF
              !
              ovl_matrix(ibnd,jbnd) = ovl_value
              !
           ENDDO
           !
        ENDDO
        !
        CALL mp_sum(ovl_matrix,inter_image_comm)
        !
     ENDIF
     !
     ! Save u_matrix and ovl_matrix to file
     !
     CALL write_umatrix_and_omatrix(nbndval, current_spin, u_matrix, ovl_matrix)
     !
     DEALLOCATE(u_matrix)
     !
  ELSE
     !
     WRITE(labels,'(i1)') current_spin
     fname = TRIM(wbse_init_save_dir)//'/evc_loc.'//labels//'.dat'
     CALL plep_read_G_and_distribute(fname,evc_loc,nbndval)
     !
  ENDIF
  !
  IF(ALLOCATED(psic1)) DEALLOCATE(psic1)
  !
  CALL stop_clock('local')
  !
END SUBROUTINE
