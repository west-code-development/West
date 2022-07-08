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
!-----------------------------------------------------------------------
SUBROUTINE do_dip()
  !-----------------------------------------------------------------------
  !
  USE control_flags,        ONLY : gamma_only
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : iuwfc,lrwfc,l_skip_nl_part_of_hcomr,westpp_range,logfile
  USE mp_world,             ONLY : mpime,root
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_bcast,mp_sum
  USE pwcom,                ONLY : npw,npwx,current_spin,isk,xk,lsda,igk_k,current_k,ngk,nspin,et
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE uspp_init,            ONLY : init_us_2
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : npol
  USE cell_base,            ONLY : bg
  USE buffers,              ONLY : get_buffer
  USE types_bz_grid,        ONLY : k_grid
  USE json_module,          ONLY : json_file,json_core,json_value
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE uspp,                 ONLY : vkb,nkb,deeq,qq_at
  USE becmod_subs_gpum,     ONLY : using_becp_auto,using_becp_d_auto
  USE wvfct_gpum,           ONLY : using_et,using_et_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu,allocate_macropol_gpu,&
                                 & deallocate_macropol_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc
  USE uspp,                 ONLY : vkb,nkb
#endif
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iks
  INTEGER :: ipol
  INTEGER :: icart
  INTEGER :: istate
  INTEGER :: jstate
  INTEGER :: nstate
  INTEGER :: iaux
  INTEGER :: iunit
  INTEGER :: trans(2)
  REAL(DP) :: aux_r(3)
  REAL(DP), ALLOCATABLE :: dip_cryst_r(:,:,:)
  REAL(DP), ALLOCATABLE :: dip_cart_r(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dip_cryst_c(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dip_cart_c(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: Hx_psi(:,:)
  CHARACTER(LEN=5) :: label_k
  CHARACTER(LEN=9) :: label_d
  TYPE(bar_type) :: barra
  TYPE(json_file) :: json
  TYPE(json_core) :: jcor
  TYPE(json_value), POINTER :: jval
  !
  IF(nspin == 4) CALL errore('do_dip','nspin 4 not yet implemented',1)
  !
  nstate = westpp_range(2)-westpp_range(1)+1
  IF(gamma_only) THEN
     ALLOCATE(dip_cryst_r(nstate,nstate,3))
     !$acc enter data create(dip_cryst_r)
     ALLOCATE(dip_cart_r(nstate,nstate,3))
  ELSE
     ALLOCATE(dip_cryst_c(nstate,nstate,3))
     !$acc enter data create(dip_cryst_c)
     ALLOCATE(dip_cart_c(nstate,nstate,3))
  ENDIF
  ALLOCATE(Hx_psi(npwx*npol,nstate))
  !$acc enter data create(Hx_psi)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
  CALL allocate_macropol_gpu(nstate)
#endif
  !
  IF(mpime == root) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
  !
  CALL io_push_title('(D)ipole matrix')
  !
  CALL start_bar_type(barra,'westpp',k_grid%nps*3)
  !
  ! LOOP
  !
  DO iks = 1,k_grid%nps ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     npw = ngk(iks)
     IF(lsda) current_spin = isk(iks)
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
#if defined(__CUDA)
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... read in wavefunctions
     !
     IF(k_grid%nps > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_host,0,inter_image_comm)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
#if defined(__CUDA)
     CALL using_becp_auto(2)
     CALL using_becp_d_auto(0)
     CALL using_evc(2)
     CALL using_evc_d(0)
     CALL using_et(2)
     CALL using_et_d(0)
     !
     !$acc update device(deeq,qq_at)
#endif
     !
     DO ipol = 1,3
        !
        !$acc host_data use_device(Hx_psi)
        CALL commut_Hx_psi(iks,nstate,ipol,evc_work(:,westpp_range(1):westpp_range(2)),Hx_psi,&
        & l_skip_nl_part_of_hcomr)
        !$acc end host_data
        !
        IF(gamma_only) THEN
           !$acc host_data use_device(Hx_psi,dip_cryst_r)
           CALL glbrak_gamma(evc_work(:,westpp_range(1):westpp_range(2)),Hx_psi,&
           & dip_cryst_r(:,:,ipol),npw,npwx,nstate,nstate,nstate,npol)
           !$acc end host_data
        ELSE
           !$acc host_data use_device(Hx_psi,dip_cryst_c)
           CALL glbrak_k(evc_work(:,westpp_range(1):westpp_range(2)),Hx_psi,dip_cryst_c(:,:,ipol),&
           & npw,npwx,nstate,nstate,nstate,npol)
           !$acc end host_data
        ENDIF
        !
        CALL update_bar_type(barra,'westpp',1)
        !
     ENDDO
     !
     IF(gamma_only) THEN
        !$acc update host(dip_cryst_r)
        !
        CALL mp_sum(dip_cryst_r,intra_bgrp_comm)
        !
        dip_cart_r = 0._DP
        DO icart = 1,3
           DO ipol = 1,3
              dip_cart_r(:,:,icart) = dip_cart_r(:,:,icart)+bg(icart,ipol)*dip_cryst_r(:,:,ipol)
           ENDDO
        ENDDO
     ELSE
        !$acc update host(dip_cryst_c)
        !
        CALL mp_sum(dip_cryst_c,intra_bgrp_comm)
        !
        dip_cart_c = (0._DP,0._DP)
        DO icart = 1,3
           DO ipol = 1,3
              dip_cart_c(:,:,icart) = dip_cart_c(:,:,icart)+bg(icart,ipol)*dip_cryst_c(:,:,ipol)
           ENDDO
        ENDDO
     ENDIF
     !
     IF(mpime == root) THEN
        !
        WRITE(label_k,'(I5.5)') iks
        !
        CALL json%add('output.D.K'//label_k//'.weight',k_grid%weight(iks))
        CALL json%add('output.D.K'//label_k//'.energies',et(:,iks))
        !
        CALL jcor%create_array(jval,'dipole')
        CALL json%add('output.D.K'//label_k//'.dipole',jval)
        !
        iaux = 0
        trans = 0
        DO jstate = 1,nstate
           trans(2) = jstate+westpp_range(1)-1
           DO istate = 1,nstate
              trans(1) = istate+westpp_range(1)-1
              iaux = iaux+1
              WRITE(label_d,'(I9)') iaux
              !
              CALL json%add('output.D.K'//label_k//'.dipole('//label_d//').trans',trans)
              IF(gamma_only) THEN
                 aux_r = dip_cart_r(istate,jstate,:)
                 CALL json%add('output.D.K'//label_k//'.dipole('//label_d//').re',aux_r)
                 aux_r = 0._DP
                 CALL json%add('output.D.K'//label_k//'.dipole('//label_d//').im',aux_r)
              ELSE
                 aux_r = REAL(dip_cart_c(istate,jstate,:),KIND=DP)
                 CALL json%add('output.D.K'//label_k//'.dipole('//label_d//').re',aux_r)
                 aux_r = AIMAG(dip_cart_c(istate,jstate,:))
                 CALL json%add('output.D.K'//label_k//'.dipole('//label_d//').im',aux_r)
              ENDIF
           ENDDO
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'westpp')
  !
  IF(gamma_only) THEN
     !$acc exit data delete(dip_cryst_r)
     DEALLOCATE(dip_cryst_r)
     DEALLOCATE(dip_cart_r)
  ELSE
     !$acc exit data delete(dip_cryst_c)
     DEALLOCATE(dip_cryst_c)
     DEALLOCATE(dip_cart_c)
  ENDIF
  !$acc exit data delete(Hx_psi)
  DEALLOCATE(Hx_psi)
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
  CALL deallocate_macropol_gpu()
#endif
  !
  IF(mpime == root) THEN
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
  ENDIF
  !
END SUBROUTINE
