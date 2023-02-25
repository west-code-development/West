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
SUBROUTINE wbse_localization(current_spin,nbnd_s,nbnd_e,evc_loc,ovl_matrix,l_restart)
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : localization,wbse_init_save_dir
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
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : evc=>evc_d,psic=>psic_d
  USE cublas
#else
  USE wavefunctions,        ONLY : evc,psic
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: current_spin,nbnd_s,nbnd_e
  REAL(DP),INTENT(OUT) :: ovl_matrix(nbnd_e-nbnd_s+1,nbnd_e-nbnd_s+1)
  COMPLEX(DP),INTENT(OUT) :: evc_loc(npwx,nbnd_e-nbnd_s+1)
  LOGICAL,INTENT(IN) :: l_restart
  !
  LOGICAL :: l_wann
  INTEGER :: ibnd,jbnd,ibnd_l,ibnd_g,jbnd_g,ir,il
  INTEGER :: bisec_i,bisec_j
  INTEGER :: dffts_nnr
  INTEGER :: nbnd_do
  INTEGER,ALLOCATABLE :: bisec_loc(:)
  REAL(DP) :: reduce
  REAL(DP) :: val(6)
  REAL(DP) :: ovl_val
  REAL(DP),ALLOCATABLE :: proj(:,:)
  REAL(DP),ALLOCATABLE :: a_matrix(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: u_matrix(:,:)
  COMPLEX(DP),ALLOCATABLE :: evc_tmp(:,:)
  CHARACTER(LEN=256) :: fname
  CHARACTER :: labels
#if !defined(__CUDA)
  REAL(DP),EXTERNAL :: DDOT
#endif
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
  nbnd_do = nbnd_e-nbnd_s+1
  dffts_nnr = dffts%nnr
  !
  IF(.NOT. l_restart) THEN
     !
     ovl_matrix(:,:) = 0._DP
     !
     ALLOCATE(u_matrix(nbnd_do,nbnd_do))
     !
     IF(l_wann) THEN
        !
        aband = idistribute()
        CALL aband%init(nbnd_do,'i','wann_local',.TRUE.)
        !
        ! compute unitary rotation matrix
        !
        ALLOCATE(a_matrix(nbnd_do,nbnd_do,6))
        ALLOCATE(proj(dffts%nnr,6))
        !
        CALL wann_calc_proj(proj)
        !
        !$acc enter data copyin(proj)
        !
        a_matrix(:,:,:) = 0._DP
        !
        DO ibnd_l = 1,aband%nloc
           !
           ibnd = aband%l2g(ibnd_l)
           ibnd_g = ibnd+nbnd_s-1
           !
           DO jbnd = ibnd,nbnd_do
              !
              jbnd_g = jbnd+nbnd_s-1
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd_g),evc(:,jbnd_g),psic,'Wave')
              !
              DO il = 1,6
                 !
                 reduce = 0._DP
                 !
                 !$acc parallel loop reduction(+:reduce) present(proj)
                 DO ir = 1,dffts_nnr
                    reduce = reduce+REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))*proj(ir,il)
                 ENDDO
                 !$acc end parallel
                 !
                 val(il) = reduce
                 !
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
        CALL wann_joint_d(nbnd_do,a_matrix,6,u_matrix)
        !
        DEALLOCATE(a_matrix)
        !$acc exit data delete(proj)
        DEALLOCATE(proj)
        !
        ! compute localized wfc
        !
        !$acc enter data copyin(u_matrix)
        !
        !$acc host_data use_device(u_matrix,evc_loc)
        CALL ZGEMM('N','N',npw,nbnd_do,nbnd_do,(1._DP,0._DP),evc(:,nbnd_s:nbnd_e),npwx,&
        & u_matrix,nbnd_do,(0._DP,0._DP),evc_loc,npwx)
        !$acc end host_data
        !
        !$acc update host(evc_loc)
        !$acc exit data delete(u_matrix)
        !
        ! compute overlap
        !
        DO ibnd_l = 1,aband%nloc
           !
           ibnd = aband%l2g(ibnd_l)
           !
           DO jbnd = 1,nbnd_do
              !
              !$acc host_data use_device(evc_loc)
              CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
              !$acc end host_data
              !
              CALL check_ovl_wannier(psic,ovl_val)
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
        ALLOCATE(evc_tmp(npwx,nbnd_e))
        !
        CALL load_qbox_wfc(current_spin,nbnd_e,evc_tmp)
        !
        evc_loc(:,:) = evc_tmp(:,nbnd_s:nbnd_e)
        !
        DEALLOCATE(evc_tmp)
        !
        ! compute unitary rotation matrix
        !
        u_matrix(:,:) = (0._DP,0._DP)
        !
        DO ibnd = 1,nbnd_do
           !
           ibnd_g = ibnd+nbnd_s-1
           !
           DO jbnd = 1,nbnd_do
              !
#if !defined(__CUDA)
              u_matrix(ibnd,jbnd) = 2._DP*DDOT(2*npwx,evc(:,ibnd_g),1,evc_loc(:,jbnd),1)
#endif
              !
              IF(gstart == 2) THEN
                 u_matrix(ibnd,jbnd) = u_matrix(ibnd,jbnd) &
                 & - CMPLX(REAL(evc(1,ibnd_g),KIND=DP)*REAL(evc_loc(1,jbnd),KIND=DP),KIND=DP)
              ENDIF
              !
           ENDDO
        ENDDO
        !
        CALL mp_sum(u_matrix,intra_bgrp_comm)
        !
        ALLOCATE(bisec_loc(nbnd_e))
        !
        ! read bisection localization from file
        !
        CALL read_bisection_loc(current_spin,nbnd_e,bisec_loc)
        !
        ! compute overlap
        !
        DO ibnd = 1,nbnd_do
           ibnd_g = ibnd+nbnd_s-1
           !
           DO jbnd = 1,nbnd_do
              !
              jbnd_g = jbnd+nbnd_s-1
              !
              bisec_i = bisec_loc(ibnd_g)
              bisec_j = bisec_loc(jbnd_g)
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
     CALL plep_merge_and_write_G(fname,evc_loc,nbnd_do)
     !
     CALL write_umatrix_and_omatrix(nbnd_do,current_spin,u_matrix,ovl_matrix)
     !
     DEALLOCATE(u_matrix)
     !
  ELSE
     !
     WRITE(labels,'(i1)') current_spin
     fname = TRIM(wbse_init_save_dir)//'/evc_loc.'//labels//'.dat'
     CALL plep_read_G_and_distribute(fname,evc_loc,nbnd_do)
     !
  ENDIF
  !
  CALL stop_clock('local')
  !
END SUBROUTINE
