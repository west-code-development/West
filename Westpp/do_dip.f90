!
! Copyright (C) 2015-2021 M. Govoni
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
  USE wavefunctions,        ONLY : evc
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE uspp,                 ONLY : vkb,nkb
  USE uspp_init,            ONLY : init_us_2
  USE io_push,              ONLY : io_push_title
  USE noncollin_module,     ONLY : npol
  USE cell_base,            ONLY : bg
  USE buffers,              ONLY : get_buffer
  USE types_bz_grid,        ONLY : k_grid
  USE json_module,          ONLY : json_file,json_core,json_value
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
  CHARACTER(5) :: label_k
  CHARACTER(9) :: label_d
  TYPE(bar_type) :: barra
  TYPE(json_file) :: json
  TYPE(json_core) :: jcor
  TYPE(json_value), POINTER :: jval
  !
  IF(nspin == 4) CALL errore('do_dip','nspin 4 not yet implemented',1)
  !
  CALL io_push_title('(D)ipole matrix')
  !
  nstate = westpp_range(2)-westpp_range(1)+1
  IF(gamma_only) THEN
     ALLOCATE(dip_cryst_r(nstate,nstate,3))
     ALLOCATE(dip_cart_r(nstate,nstate,3))
  ELSE
     ALLOCATE(dip_cryst_c(nstate,nstate,3))
     ALLOCATE(dip_cart_c(nstate,nstate,3))
  ENDIF
  ALLOCATE(Hx_psi(npwx*npol,nstate))
  !
  IF(mpime == root) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
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
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb)
     !
     ! ... read in wavefunctions
     !
     IF(k_grid%nps > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     DO ipol = 1,3
        !
        CALL commut_Hx_psi(iks,nstate,ipol,evc(:,westpp_range(1):westpp_range(2)),Hx_psi,&
        & l_skip_nl_part_of_hcomr)
        !
        IF(gamma_only) THEN
           CALL glbrak_gamma(evc(:,westpp_range(1):westpp_range(2)),Hx_psi,dip_cryst_r(:,:,ipol),&
           & npw,npwx,nstate,nstate,nstate,npol)
        ELSE
           CALL glbrak_k(evc(:,westpp_range(1):westpp_range(2)),Hx_psi,dip_cryst_c(:,:,ipol),npw,&
           & npwx,nstate,nstate,nstate,npol)
        ENDIF
        !
        CALL update_bar_type(barra,'westpp',1)
        !
     END DO
     !
     IF(gamma_only) THEN
        CALL mp_sum(dip_cryst_r,intra_bgrp_comm)
        !
        dip_cart_r = 0._DP
        DO icart = 1,3
           DO ipol = 1,3
              dip_cart_r(:,:,icart) = dip_cart_r(:,:,icart)+bg(icart,ipol)*dip_cryst_r(:,:,ipol)
           ENDDO
        ENDDO
     ELSE
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
              CALL json%add('output.D.K'//label_k//'.dipole('//TRIM(ADJUSTL(label_d))//').trans',trans)
              IF(gamma_only) THEN
                 aux_r = dip_cart_r(istate,jstate,:)
                 CALL json%add('output.D.K'//label_k//'.dipole('//TRIM(ADJUSTL(label_d))//').re',aux_r)
                 aux_r = 0._DP
                 CALL json%add('output.D.K'//label_k//'.dipole('//TRIM(ADJUSTL(label_d))//').im',aux_r)
              ELSE
                 aux_r = REAL(dip_cart_c(istate,jstate,:),KIND=DP)
                 CALL json%add('output.D.K'//label_k//'.dipole('//TRIM(ADJUSTL(label_d))//').re',aux_r)
                 aux_r = AIMAG(dip_cart_c(istate,jstate,:))
                 CALL json%add('output.D.K'//label_k//'.dipole('//TRIM(ADJUSTL(label_d))//').im',aux_r)
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
     DEALLOCATE(dip_cryst_r)
     DEALLOCATE(dip_cart_r)
  ELSE
     DEALLOCATE(dip_cryst_c)
     DEALLOCATE(dip_cart_c)
  ENDIF
  DEALLOCATE(Hx_psi)
  !
  IF(mpime == root) THEN
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
  ENDIF
  !
END SUBROUTINE
