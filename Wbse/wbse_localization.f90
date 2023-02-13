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
!----------------------------------------------------------------------------
SUBROUTINE wbse_localization(current_spin,nbndval,evc_loc,ovl_matrix,l_restart)
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : localization,wbse_init_save_dir
  USE wavefunctions,        ONLY : evc,psic
  USE plep_io,              ONLY : plep_merge_and_write_G,plep_read_G_and_distribute
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma
  USE pwcom,                ONLY : npw,npwx
  USE gvect,                ONLY : gstart
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE qbox_interface,       ONLY : load_qbox_wfc
  USE check_ovl_wfc,        ONLY : check_ovl_wannier,read_bisection_loc,check_ovl_bisection
  USE wbse_io,              ONLY : write_umatrix_and_omatrix
  USE wann_loc_wfc,         ONLY : wann_calc_proj,wann_joint_d
  USE distribution_center,  ONLY : aband
  USE class_idistribute,    ONLY : idistribute
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: current_spin,nbndval
  REAL(DP),INTENT(OUT) :: ovl_matrix(nbndval,nbndval)
  COMPLEX(DP),INTENT(OUT) :: evc_loc(npwx,nbndval)
  LOGICAL,INTENT(IN) :: l_restart
  !
  LOGICAL :: l_wann
  INTEGER :: ibnd_l,ibnd,jbnd,ir,il
  INTEGER :: bisec_i,bisec_j
  INTEGER,ALLOCATABLE :: bisec_loc(:)
  REAL(DP) :: val(6)
  REAL(DP) :: ovl_val
  REAL(DP),ALLOCATABLE :: proj(:,:)
  REAL(DP),ALLOCATABLE :: a_matrix(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: u_matrix(:,:)
  CHARACTER(LEN=256) :: fname
  CHARACTER :: labels
  REAL(DP),EXTERNAL :: DDOT
  !
  CALL start_clock('local')
  !
  SELECT CASE(TRIM(localization))
  CASE('W','w')
     l_wann = .TRUE.
  CASE('B','b')
     l_wann = .FALSE.
  END SELECT
  !
  IF(.NOT. l_restart) THEN
     !
     ovl_matrix(:,:) = 0._DP
     !
     ALLOCATE(u_matrix(nbndval,nbndval))
     !
     IF(l_wann) THEN
        !
        aband = idistribute()
        CALL aband%init(nbndval,'i','wann_local',.TRUE.)
        !
        ! compute unitary rotation matrix
        !
        ALLOCATE(a_matrix(nbndval,nbndval,6))
        ALLOCATE(proj(dffts%nnr,6))
        !
        CALL wann_calc_proj(proj)
        !
        a_matrix(:,:,:) = 0._DP
        !
        DO ibnd_l = 1,aband%nloc
           !
           ibnd = aband%l2g(ibnd_l)
           !
           DO jbnd = ibnd,nbndval
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
              !
              val = 0._DP
              DO il = 1,6
                 DO ir = 1,dffts%nnr
                    val(il) = val(il) + REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))*proj(ir,il)
                 ENDDO
              ENDDO
              !
              a_matrix(ibnd,jbnd,1:6) = val(1:6)
              IF(ibnd /= jbnd) a_matrix(jbnd,ibnd,1:6) = val(1:6)
              !
           ENDDO
           !
        ENDDO
        !
        CALL mp_sum(a_matrix,intra_bgrp_comm)
        CALL mp_sum(a_matrix,inter_image_comm)
        !
        CALL wann_joint_d(nbndval,a_matrix,6,u_matrix)
        !
        DEALLOCATE(a_matrix)
        DEALLOCATE(proj)
        !
        ! compute localized wfc
        !
        CALL ZGEMM('N','N',npw,nbndval,nbndval,(1._DP,0._DP),evc,npwx,u_matrix,nbndval,&
        & (0._DP,0._DP),evc_loc,npwx)
        !
        ! compute overlap
        !
        DO ibnd_l = 1,aband%nloc
           !
           ibnd = aband%l2g(ibnd_l)
           !
           DO jbnd = 1,nbndval
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
              CALL check_ovl_wannier(REAL(psic,KIND=DP),AIMAG(psic),ovl_val)
              !
              ovl_matrix(ibnd,jbnd) = ovl_val
              !
           ENDDO
           !
        ENDDO
        !
        CALL mp_sum(ovl_matrix,inter_image_comm)
        !
     ELSE
        !
        ! load localized wfc
        !
        CALL load_qbox_wfc(current_spin,nbndval,evc_loc)
        !
        ! compute unitary rotation matrix
        !
        u_matrix(:,:) = (0._DP,0._DP)
        !
        DO ibnd = 1,nbndval
           DO jbnd = 1,nbndval
              !
              u_matrix(ibnd,jbnd) = 2._DP*DDOT(2*npwx,evc(:,ibnd),1,evc_loc(:,jbnd),1)
              !
              IF(gstart == 2) THEN
                 u_matrix(ibnd,jbnd) = u_matrix(ibnd,jbnd) &
                 & - CMPLX(REAL(evc(1,ibnd),KIND=DP)*REAL(evc_loc(1,jbnd),KIND=DP),KIND=DP)
              ENDIF
              !
           ENDDO
        ENDDO
        !
        CALL mp_sum(u_matrix,intra_bgrp_comm)
        !
        ALLOCATE(bisec_loc(nbndval))
        !
        ! read bisection localization from file
        !
        CALL read_bisection_loc(current_spin,nbndval,bisec_loc)
        !
        ! compute overlap
        !
        DO ibnd = 1,nbndval
           DO jbnd = 1,nbndval
              !
              bisec_i = bisec_loc(ibnd)
              bisec_j = bisec_loc(jbnd)
              !
              CALL check_ovl_bisection(bisec_i,bisec_j,ovl_val)
              !
              ovl_matrix(ibnd,jbnd) = ovl_val
              !
           ENDDO
        ENDDO
        !
        DEALLOCATE (bisec_loc)
        !
     ENDIF
     !
     ! save to file
     !
     WRITE(labels,'(i1)') current_spin
     fname = TRIM(wbse_init_save_dir)//'/evc_loc.'//labels//'.dat'
     CALL plep_merge_and_write_G(fname,evc_loc,nbndval)
     !
     CALL write_umatrix_and_omatrix(nbndval,current_spin,u_matrix,ovl_matrix)
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
  CALL stop_clock('local')
  !
END SUBROUTINE
