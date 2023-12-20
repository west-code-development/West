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
! Yu Jin
!
SUBROUTINE do_exc_comp()
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE io_push,               ONLY : io_push_title
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE pwcom,                 ONLY : npw,npwx,nks,current_k,ngk
  USE control_flags,         ONLY : gamma_only
  USE gvect,                 ONLY : gstart
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE mp_global,             ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE cell_base,             ONLY : bg
  USE buffers,               ONLY : get_buffer
  USE westcom,               ONLY : logfile,iuwfc,lrwfc,nbndval0x,nbnd_occ,dvg_exc,ev,westpp_range,&
                                  & westpp_n_liouville_to_use,westpp_l_spin_flip,westpp_l_compute_tdm,&
                                  & westpp_l_dipole_realspace,l_dipole_realspace,d0psi,alphapv_dfpt
  USE mp_world,              ONLY : mpime,root
  USE plep_db,               ONLY : plep_db_read
  USE distribution_center,   ONLY : pert,kpt_pool,band_group
  USE class_idistribute,     ONLY : idistribute
  USE types_bz_grid,         ONLY : k_grid
  USE json_module,           ONLY : json_file,json_core,json_value
  USE wvfct,                 ONLY : nbnd
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE wavefunctions,         ONLY : evc_host=>evc
  USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu
#else
  USE wavefunctions,         ONLY : evc_work=>evc
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: iks,iks_do,iexc,lexc,ig,iocc,iemp,iaux,iunit,ipol,icart
  INTEGER :: trans(4)
  INTEGER :: nbndx_occ,nbndx_emp
  INTEGER :: nbndval,flnbndval
  CHARACTER(5) :: label_exc
  CHARACTER(5) :: label_k
  CHARACTER(9) :: label_d
  REAL(DP) :: reduce
  REAL(DP), ALLOCATABLE :: projection_matrix(:,:,:,:)
  REAL(DP), ALLOCATABLE :: transition_dipole_cry(:,:), transition_dipole_cart(:,:)
  TYPE(bar_type) :: barra
  TYPE(json_file) :: json
  TYPE(json_core) :: jcor
  TYPE(json_value), POINTER :: jval
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  COMPLEX(DP), EXTERNAL :: get_alpha_pv
  !
  nbndx_occ = MAXVAL(nbnd_occ)
  nbndx_emp = nbnd - MINVAL(nbnd_occ)
  !
  IF(nbndx_emp < 1) &
     CALL errore('do_exc_comp', 'decomposition of eigenvectors needs empty states, rerun pwscf with nbnd > nbnd_occ', 1)
  !
  IF(westpp_range(1) > nbndx_occ) &
     CALL errore('do_exc_comp', 'westpp_range(1) should not exceed the number of occupied bands', 1)
  !
  IF(westpp_range(2) <= nbndx_occ) &
     CALL errore('do_exc_comp', 'westpp_range(2) should be greater than the number of occupied bands', 1)
  !
  IF(westpp_range(2) > nbnd) &
     CALL errore('do_exc_comp', 'westpp_range(2) should not exceed the number of bands', 1)
  !
  IF(.NOT. gamma_only) &
     CALL errore('do_exc_comp', 'decomposition of eigenvectors requires gamma_only', 1)
  !
  ! ... DISTRIBUTE
  !
  pert = idistribute()
  CALL pert%init(westpp_n_liouville_to_use,'i','nvec',.TRUE.)
  kpt_pool = idistribute()
  CALL kpt_pool%init(k_grid%nps,'p','kpt',.FALSE.)
  band_group = idistribute()
  CALL band_group%init(nbndval0x,'b','nbndval',.FALSE.)
  !
  ! READ EIGENVALUES AND VECTORS FROM OUTPUT
  !
  CALL plep_db_read(westpp_n_liouville_to_use)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  IF(mpime == root) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
  !
  IF((.NOT. westpp_l_spin_flip) .AND. westpp_l_compute_tdm) THEN
     !
     ALLOCATE(d0psi(npwx,nbndx_occ,nks,3))
     !$acc enter data create(d0psi)
     !
     l_dipole_realspace = westpp_l_dipole_realspace
     !
     ! Calculate ALPHA_PV
     !
     alphapv_dfpt = get_alpha_pv()
     !
     CALL solve_e_psi()
     !
     ALLOCATE(transition_dipole_cry(3,westpp_n_liouville_to_use))
     ALLOCATE(transition_dipole_cart(3,westpp_n_liouville_to_use))
     transition_dipole_cry(:,:) = 0._DP
     transition_dipole_cart(:,:) = 0._DP
     !
  ENDIF
  !
  ALLOCATE(projection_matrix(nks,nbndx_occ,nbndx_emp,westpp_n_liouville_to_use))
  !$acc enter data create(projection_matrix) copyin(dvg_exc)
  !
  !$acc kernels present(projection_matrix)
  projection_matrix(:,:,:,:) = 0._DP
  !$acc end kernels
  !
  CALL io_push_title('BSE/TDDFT Excited State De(C)omposition')
  !
  CALL start_bar_type(barra,'westpp',pert%nloc*k_grid%nps)
  !
  DO lexc = 1,pert%nloc
     !
     ! local -> global
     !
     iexc = pert%l2g(lexc)
     !
     DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
        !
        ! ... Set k-point, spin, kinetic energy, needed by Hpsi
        !
        current_k = iks
        npw = ngk(iks)
        !
        ! ... read in wavefunctions from the previous iteration
        !
        IF(k_grid%nps > 1) THEN
#if defined(__CUDA)
           IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc_host,0,inter_image_comm)
           !
           CALL using_evc(2)
           CALL using_evc_d(0)
#else
           IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
        ENDIF
        !
        IF(westpp_l_spin_flip) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        !
        nbndval = nbnd_occ(iks)
        flnbndval = nbnd_occ(iks_do)
        !
        !$acc parallel vector_length(1024) present(dvg_exc,projection_matrix)
        !$acc loop collapse(2)
        DO iocc = 1, flnbndval
           DO iemp = 1, nbnd - nbndval
              !
              reduce = 0._DP
              !$acc loop reduction(+:reduce)
              DO ig = 1, npw
                 reduce = reduce + 2._DP*REAL(dvg_exc(ig,iocc,iks,lexc),KIND=DP)*REAL(evc_work(ig,nbndval+iemp),KIND=DP) &
                 &               + 2._DP*AIMAG(dvg_exc(ig,iocc,iks,lexc))*AIMAG(evc_work(ig,nbndval+iemp))
              ENDDO
              !
              IF(gstart == 2) THEN
                 reduce = reduce - REAL(dvg_exc(1,iocc,iks,lexc),KIND=DP)*REAL(evc_work(1,nbndval+iemp),KIND=DP)
              ENDIF
              !
              projection_matrix(iks,iocc,iemp,iexc) = reduce
              !
           ENDDO
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'westpp',1)
        !
     ENDDO
     !
     ! Compute the transition dipole moment
     !
     IF((.NOT. westpp_l_spin_flip) .AND. westpp_l_compute_tdm) THEN
        !
        DO ipol = 1, 3
           !
           reduce = 0._DP
           !
           DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
              !
              ! ... Set k-point, spin, kinetic energy, needed by Hpsi
              !
              current_k = iks
              npw = ngk(iks)
              !
              IF(westpp_l_spin_flip) THEN
                 iks_do = flks(iks)
              ELSE
                 iks_do = iks
              ENDIF
              !
              nbndval = nbnd_occ(iks)
              flnbndval = nbnd_occ(iks_do)
              !
              !$acc parallel loop collapse(2) reduction(+:reduce) present(dvg_exc,d0psi) copy(reduce)
              DO iocc = 1, flnbndval
                 DO ig = 1, npw
                    reduce = reduce + 2._DP*REAL(dvg_exc(ig,iocc,iks,lexc),KIND=DP)*REAL(d0psi(ig,iocc,iks_do,ipol),KIND=DP) &
                    &               + 2._DP*AIMAG(dvg_exc(ig,iocc,iks,lexc))*AIMAG(d0psi(ig,iocc,iks_do,ipol))
                 ENDDO
              ENDDO
              !$acc end parallel
              !
              IF(gstart == 2) THEN
                 !$acc parallel loop reduction(+:reduce) present(dvg_exc,d0psi) copy(reduce)
                 DO iocc = 1, flnbndval
                    reduce = reduce - REAL(dvg_exc(1,iocc,iks,lexc),KIND=DP)*REAL(d0psi(1,iocc,iks_do,ipol),KIND=DP)
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ENDDO
           !
           transition_dipole_cry(ipol,iexc) = reduce
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  !$acc update host(projection_matrix)
  !
  CALL mp_sum(projection_matrix,intra_bgrp_comm)
  CALL mp_sum(projection_matrix,inter_image_comm)
  !
  IF((.NOT. westpp_l_spin_flip) .AND. westpp_l_compute_tdm) THEN
     !
     CALL mp_sum(transition_dipole_cry,intra_bgrp_comm)
     CALL mp_sum(transition_dipole_cry,inter_image_comm)
     !
     IF(nks == 1) transition_dipole_cry(:,:) = SQRT(2._DP) * transition_dipole_cry
     !
     IF(westpp_l_dipole_realspace) THEN
        !
        transition_dipole_cart(:,:) = transition_dipole_cry(:,:)
        !
     ELSE
        !
        transition_dipole_cart(:,:) = 0._DP
        DO icart = 1,3
           DO ipol = 1,3
              transition_dipole_cart(icart,:) = transition_dipole_cart(icart,:) &
                      & +bg(icart,ipol)*transition_dipole_cry(ipol,:)
           ENDDO
        ENDDO
        !
     ENDIF
     !
  ENDIF
  !
  CALL stop_bar_type(barra,'westpp')
  !
  ! ... Print out results
  !
  WRITE(stdout, "(/,5x,'*-------------* THE PRINCIPLE PROJECTED COMPONENTS *-------------*')")
  !
  DO iexc = 1,westpp_n_liouville_to_use
     !
     WRITE(stdout, "(/, 5x, '#     Exciton :   ', i8,' |','   ','Excitation energy :   ', f12.6)") iexc, ev(iexc)
     WRITE(stdout, "(   5x, '#     Transition_from      |   Transition_to       |    Coeffcient')")
     !
     DO iks = 1,nks
        !
        IF(westpp_l_spin_flip) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        !
        DO iocc = 1,nbnd_occ(iks_do)
           DO iemp = 1,(nbnd - nbnd_occ(iks))
              IF(ABS(projection_matrix(iks,iocc,iemp,iexc)) >= 0.1_DP) THEN
                 WRITE(stdout, "(4x, i8, 4x, i8, 8x, '|', i4, 4x, i8, 7x, '|', f13.6)") &
                 & iks_do,iocc,iks,iemp+nbnd_occ(iks),projection_matrix(iks,iocc,iemp,iexc)
              ENDIF
           ENDDO
        ENDDO
        !
     ENDDO
     !
     IF((.NOT. westpp_l_spin_flip) .AND. westpp_l_compute_tdm) THEN
        WRITE(stdout, "(5x, '#     TDM_x                |   TDM_y               |    TDM_z')")
        WRITE(stdout, "(9x, f18.9, 5x, '|', f17.9, 6x, '|', f16.9)") &
        & transition_dipole_cart(1,iexc),transition_dipole_cart(2,iexc),transition_dipole_cart(3,iexc)
     ENDIF
     !
  ENDDO
  !
  ! ... Write the results in the json file
  !
  IF(mpime == root) THEN
     !
     DO iexc = 1,westpp_n_liouville_to_use
        !
        WRITE(label_exc,'(I5.5)') iexc
        !
        CALL json%add('output.E'//label_exc//'.excitation_energy',ev(iexc))
        !
        IF((.NOT. westpp_l_spin_flip) .AND. westpp_l_compute_tdm) &
        & CALL json%add('output.E'//label_exc//'.transition_dipole_moment',transition_dipole_cart(:,iexc))
        !
        DO iks = 1,nks
           !
           IF(westpp_l_spin_flip) THEN
              iks_do = flks(iks)
           ELSE
              iks_do = iks
           ENDIF
           !
           WRITE(label_k,'(I5.5)') iks
           !
           CALL jcor%create_array(jval,'projection')
           CALL json%add('output.E'//label_exc//'.K'//label_k//'.projection',jval)
           !
           iaux = 0
           trans = 0
           !
           DO iocc = westpp_range(1),nbnd_occ(iks_do)
              !
              trans(1) = iks_do
              trans(2) = iocc
              !
              DO iemp = 1,(westpp_range(2) - nbnd_occ(iks))
                 !
                 trans(3) = iks
                 trans(4) = iemp + nbnd_occ(iks)
                 iaux = iaux+1
                 WRITE(label_d,'(I9)') iaux
                 !
                 CALL json%add('output.E'//label_exc//'.K'//label_k//'.projection('//label_d//').trans',trans)
                 !
                 reduce = projection_matrix(iks,iocc,iemp,iexc)
                 CALL json%add('output.E'//label_exc//'.K'//label_k//'.projection('//label_d//').value',reduce)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
  !$acc exit data delete(projection_matrix,dvg_exc)
  DEALLOCATE(projection_matrix)
  IF((.NOT. westpp_l_spin_flip) .AND. westpp_l_compute_tdm) THEN
     !$acc exit data delete(d0psi)
     DEALLOCATE(d0psi)
     DEALLOCATE(transition_dipole_cry)
     DEALLOCATE(transition_dipole_cart)
  ENDIF
  !
  IF(mpime == root) THEN
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
  ENDIF
  !
END SUBROUTINE
