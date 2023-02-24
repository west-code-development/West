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
! Marco Govoni
!
!------------------------------------------------------------------------
SUBROUTINE do_wann()
  !------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE cell_base,            ONLY : at,alat
  USE wavefunctions,        ONLY : evc,psic
  USE westcom,              ONLY : iuwfc,lrwfc,westpp_range,logfile
  USE mp_world,             ONLY : mpime,root
  USE mp,                   ONLY : mp_bcast,mp_sum
  USE mp_global,            ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma
  USE pwcom,                ONLY : npw,npwx,current_k,ngk
  USE buffers,              ONLY : get_buffer
  USE types_bz_grid,        ONLY : k_grid
  USE wann_loc_wfc,         ONLY : wann_calc_proj,wann_joint_d
  USE distribution_center,  ONLY : aband
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE json_module,          ONLY : json_file,json_core,json_value
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: nstate,local_ib,global_ib,global_jb,ib,jb
  INTEGER :: iks,ir,il
  INTEGER :: iunit
  REAL(DP) :: val(6)
  REAL(DP) :: tmp(3)
  REAL(DP) :: wan_center(3)
  REAL(DP), ALLOCATABLE :: proj(:,:)
  REAL(DP), ALLOCATABLE :: amat(:,:,:)
  REAL(DP), ALLOCATABLE :: umat_r(:)
  COMPLEX(DP), ALLOCATABLE :: umat(:,:)
  CHARACTER(LEN=5) :: label_k
  CHARACTER(LEN=9) :: label_b
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  TYPE(json_file) :: json
  TYPE(json_core) :: jcor
  TYPE(json_value), POINTER :: jval
  !
  IF(mpime == root) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
  !
  nstate = westpp_range(2)-westpp_range(1)+1
  !
  ALLOCATE(proj(dffts%nnr,6))
  ALLOCATE(amat(nstate,nstate,6))
  ALLOCATE(umat(nstate,nstate))
  !
  aband = idistribute()
  CALL aband%init(nstate,'i','westpp_range',.TRUE.)
  !
  barra_load = 0
  DO iks = 1,k_grid%nps
     DO local_ib = 1,aband%nloc
        global_ib = aband%l2g(local_ib)
        DO global_jb = global_ib,nstate
           barra_load = barra_load+1
        ENDDO
     ENDDO
  ENDDO
  !
  CALL io_push_title('(B)oys/Wannier localization')
  !
  CALL start_bar_type(barra,'westpp',barra_load)
  !
  DO iks = 1,k_grid%nps ! KPOINT-SPIN LOOP
     !
     ! ... set k-point etc
     !
     current_k = iks
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(k_grid%nps > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     ! compute unitary transformation matrix
     !
     CALL wann_calc_proj(proj)
     !
     amat(:,:,:) = 0._DP
     !
     DO local_ib = 1,aband%nloc
        !
        global_ib = aband%l2g(local_ib)
        ib = global_ib+westpp_range(1)-1
        !
        DO global_jb = global_ib,nstate
           !
           jb = global_jb+westpp_range(1)-1
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ib),evc(:,jb),psic,'Wave')
           !
           val = 0._DP
           DO il = 1,6
              DO ir = 1,dffts%nnr
                 val(il) = val(il) + REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))*proj(ir,il)
              ENDDO
           ENDDO
           !
           amat(global_ib,global_jb,1:6) = val(1:6)
           IF(ib /= jb) amat(global_jb,global_ib,1:6) = val(1:6)
           !
           CALL update_bar_type(barra,'westpp',1)
           !
        ENDDO
        !
     ENDDO
     !
     CALL mp_sum(amat,intra_bgrp_comm)
     CALL mp_sum(amat,inter_image_comm)
     !
     CALL wann_joint_d(nstate,amat,6,umat)
     !
     IF(mpime == root) THEN
        !
        WRITE(label_k,'(I5.5)') iks
        !
        ! output Wannier centers
        !
        CALL jcor%create_array(jval,'center')
        CALL json%add('output.B.K'//label_k//'.wan_center',jval)
        !
        DO ib = 1,nstate
           !
           tmp(1) = AIMAG(LOG(CMPLX(amat(ib,ib,1),amat(ib,ib,2),KIND=DP))) / tpi
           tmp(2) = AIMAG(LOG(CMPLX(amat(ib,ib,3),amat(ib,ib,4),KIND=DP))) / tpi
           tmp(3) = AIMAG(LOG(CMPLX(amat(ib,ib,5),amat(ib,ib,6),KIND=DP))) / tpi
           !
           tmp(1) = MODULO(tmp(1),1._DP)
           tmp(2) = MODULO(tmp(2),1._DP)
           tmp(3) = MODULO(tmp(3),1._DP)
           !
           wan_center(:) = tmp(1)*at(:,1)*alat + tmp(2)*at(:,2)*alat + tmp(3)*at(:,3)*alat
           !
           WRITE(label_b,'(I9)') ib
           !
           CALL json%add('output.B.K'//label_k//'.wan_center('//label_b//').x',wan_center(1))
           CALL json%add('output.B.K'//label_k//'.wan_center('//label_b//').y',wan_center(2))
           CALL json%add('output.B.K'//label_k//'.wan_center('//label_b//').z',wan_center(3))
           !
        ENDDO
        !
        ! output transformation matrix
        !
        ALLOCATE(umat_r(nstate*nstate))
        !
        umat_r = RESHAPE(REAL(umat,KIND=DP),(/nstate*nstate/))
        !
        CALL json%add('output.B.K'//label_k//'.trans_matrix',umat_r)
        !
        DEALLOCATE(umat_r)
        !
     ENDIF
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'westpp')
  !
  DEALLOCATE(umat)
  DEALLOCATE(amat)
  DEALLOCATE(proj)
  !
  IF(mpime == root) THEN
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
  ENDIF
  !
END SUBROUTINE
