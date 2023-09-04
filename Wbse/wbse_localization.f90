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
SUBROUTINE wbse_localization(current_spin,nbnd_s,nbnd_e,evc_loc,ovl_matrix,l_restart)
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : localization
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma,single_invfft_gamma
  USE pwcom,                ONLY : npw,npwx
  USE gvect,                ONLY : gstart
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE qbox_interface,       ONLY : load_qbox_wfc
  USE check_ovl_wfc,        ONLY : check_ovl_wannier,read_bisection_loc,check_ovl_bisection
  USE wbse_io,              ONLY : write_umatrix_and_omatrix,read_umatrix_and_omatrix
  USE wann_loc_wfc,         ONLY : wann_calc_proj,wann_jade
  USE distribution_center,  ONLY : band_group
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
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
  INTEGER :: barra_load
  INTEGER,ALLOCATABLE :: bisec_loc(:)
  REAL(DP) :: reduce
  REAL(DP) :: val(6)
  REAL(DP) :: ovl_val
  REAL(DP),ALLOCATABLE :: proj(:,:)
  REAL(DP),ALLOCATABLE :: u_real(:,:)
  REAL(DP),ALLOCATABLE :: a_matrix(:,:,:)
  REAL(DP),ALLOCATABLE :: aux(:)
  REAL(DP),ALLOCATABLE :: aux2(:)
  !$acc declare device_resident(aux,aux2)
  COMPLEX(DP),ALLOCATABLE :: u_matrix(:,:)
  COMPLEX(DP),ALLOCATABLE :: evc_tmp(:,:)
  TYPE(bar_type) :: barra
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
        band_group = idistribute()
        CALL band_group%init(nbnd_do,'i','wann_local',.TRUE.)
        !
        ! compute unitary rotation matrix
        !
        ALLOCATE(proj(dffts_nnr,6))
        ALLOCATE(a_matrix(nbnd_do,nbnd_do,6))
        ALLOCATE(aux(dffts_nnr))
        ALLOCATE(aux2(dffts_nnr))
        ALLOCATE(u_real(nbnd_do,nbnd_do))
        !
        CALL wann_calc_proj(proj)
        !
        !$acc enter data copyin(proj)
        !
        a_matrix(:,:,:) = 0._DP
        !
        barra_load = 0
        DO ibnd_l = 1,band_group%nloc
           ibnd = band_group%l2g(ibnd_l)
           DO jbnd = ibnd,nbnd_do
              barra_load = barra_load+1
           ENDDO
        ENDDO
        !
        CALL io_push_title('Wannier (A matrices)')
        !
        CALL start_bar_type(barra,'wann',barra_load)
        !
        DO ibnd_l = 1,band_group%nloc
           !
           ibnd = band_group%l2g(ibnd_l)
           ibnd_g = ibnd+nbnd_s-1
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd_g),psic,'Wave')
           !
           !$acc kernels present(aux)
           aux(:) = REAL(psic,KIND=DP)
           !$acc end kernels
           !
           DO jbnd = ibnd,nbnd_do,2
              !
              jbnd_g = jbnd+nbnd_s-1
              !
              IF(jbnd < nbnd_do) THEN
                 !
                 CALL double_invfft_gamma(dffts,npw,npwx,evc(:,jbnd_g),evc(:,jbnd_g+1),psic,'Wave')
                 !
                 DO il = 1,6
                    !
                    reduce = 0._DP
                    !
                    !$acc parallel loop reduction(+:reduce) present(proj) copy(reduce)
                    DO ir = 1,dffts_nnr
                       reduce = reduce+aux(ir)*REAL(psic(ir),KIND=DP)*proj(ir,il)
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
                 DO il = 1,6
                    !
                    reduce = 0._DP
                    !
                    !$acc parallel loop reduction(+:reduce) present(proj) copy(reduce)
                    DO ir = 1,dffts_nnr
                       reduce = reduce+aux(ir)*AIMAG(psic(ir))*proj(ir,il)
                    ENDDO
                    !$acc end parallel
                    !
                    val(il) = reduce
                    !
                 ENDDO
                 !
                 a_matrix(ibnd,jbnd+1,1:6) = val(1:6)
                 IF(ibnd /= jbnd+1) a_matrix(jbnd+1,ibnd,1:6) = val(1:6)
                 !
              ELSE
                 !
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jbnd_g),psic,'Wave')
                 !
                 DO il = 1,6
                    !
                    reduce = 0._DP
                    !
                    !$acc parallel loop reduction(+:reduce) present(proj) copy(reduce)
                    DO ir = 1,dffts_nnr
                       reduce = reduce+aux(ir)*REAL(psic(ir),KIND=DP)*proj(ir,il)
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
              ENDIF
              !
           ENDDO
           !
           CALL update_bar_type(barra,'wann',nbnd_do-ibnd+1)
           !
        ENDDO
        !
        CALL mp_sum(a_matrix,intra_bgrp_comm)
        CALL mp_sum(a_matrix,inter_image_comm)
        !
        CALL stop_bar_type(barra,'wann')
        !
        CALL wann_jade(nbnd_do,a_matrix,6,u_real)
        !
        u_matrix(:,:) = CMPLX(u_real,KIND=DP)
        !
        !$acc exit data delete(proj)
        DEALLOCATE(proj)
        DEALLOCATE(a_matrix)
        DEALLOCATE(u_real)
        !
        ! compute localized wfc
        !
        !$acc enter data copyin(u_matrix)
        !
        !$acc host_data use_device(u_matrix,evc_loc)
        CALL ZGEMM('N','N',npw,nbnd_do,nbnd_do,(1._DP,0._DP),evc(1,nbnd_s),npwx,u_matrix,nbnd_do,&
        & (0._DP,0._DP),evc_loc,npwx)
        !$acc end host_data
        !
        !$acc exit data delete(u_matrix)
        !
        ! compute overlap
        !
        CALL io_push_title('Wannier (overlap)')
        !
        CALL start_bar_type(barra,'wann',barra_load)
        !
        DO ibnd_l = 1,band_group%nloc
           !
           ibnd = band_group%l2g(ibnd_l)
           !
           !$acc host_data use_device(evc_loc)
           CALL single_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),psic,'Wave')
           !$acc end host_data
           !
           !$acc kernels present(aux)
           aux(:) = REAL(psic,KIND=DP)
           !$acc end kernels
           !
           DO jbnd = ibnd,nbnd_do,2
              !
              IF(jbnd < nbnd_do) THEN
                 !
                 !$acc host_data use_device(evc_loc)
                 CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,jbnd),evc_loc(:,jbnd+1),psic,'Wave')
                 !$acc end host_data
                 !
                 !$acc kernels present(aux2)
                 aux2(:) = REAL(psic,KIND=DP)
                 !$acc end kernels
                 !
                 CALL check_ovl_wannier(aux,aux2,ovl_val)
                 !
                 ovl_matrix(ibnd,jbnd) = ovl_val
                 ovl_matrix(jbnd,ibnd) = ovl_val
                 !
                 !$acc kernels present(aux2)
                 aux2(:) = AIMAG(psic)
                 !$acc end kernels
                 !
                 CALL check_ovl_wannier(aux,aux2,ovl_val)
                 !
                 ovl_matrix(ibnd,jbnd+1) = ovl_val
                 ovl_matrix(jbnd+1,ibnd) = ovl_val
                 !
              ELSE
                 !
                 !$acc host_data use_device(evc_loc)
                 CALL single_invfft_gamma(dffts,npw,npwx,evc_loc(:,jbnd),psic,'Wave')
                 !$acc end host_data
                 !
                 !$acc kernels present(aux2)
                 aux2(:) = REAL(psic,KIND=DP)
                 !$acc end kernels
                 !
                 CALL check_ovl_wannier(aux,aux2,ovl_val)
                 !
                 ovl_matrix(ibnd,jbnd) = ovl_val
                 ovl_matrix(jbnd,ibnd) = ovl_val
                 !
              ENDIF
              !
           ENDDO
           !
           CALL update_bar_type(barra,'wann',nbnd_do-ibnd+1)
           !
        ENDDO
        !
        CALL mp_sum(ovl_matrix,inter_image_comm)
        !
        CALL stop_bar_type(barra,'wann')
        !
        DEALLOCATE(aux)
        DEALLOCATE(aux2)
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
     ! save overlap matrix and unitary transformation matrix
     !
     CALL write_umatrix_and_omatrix(nbnd_do,current_spin,u_matrix,ovl_matrix)
     !
     DEALLOCATE(u_matrix)
     !
  ELSE
     !
     ALLOCATE(u_matrix(nbnd_do,nbnd_do))
     !
     ! read overlap matrix and unitary transformation matrix
     !
     CALL read_umatrix_and_omatrix(nbnd_do,current_spin,u_matrix,ovl_matrix)
     !
     !$acc enter data copyin(u_matrix)
     !
     !$acc host_data use_device(u_matrix,evc_loc)
     CALL ZGEMM('N','N',npw,nbnd_do,nbnd_do,(1._DP,0._DP),evc(1,nbnd_s),npwx,u_matrix,nbnd_do,&
     & (0._DP,0._DP),evc_loc,npwx)
     !$acc end host_data
     !
     !$acc exit data delete(u_matrix)
     DEALLOCATE(u_matrix)
     !
  ENDIF
  !
  CALL stop_clock('local')
  !
END SUBROUTINE
