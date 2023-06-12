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
SUBROUTINE do_exc_comp()
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE io_push,               ONLY : io_push_title
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE pwcom,                 ONLY : npw,npwx,nks,igk_k,current_k,ngk,wg
  USE control_flags,         ONLY : gamma_only
  USE gvect,                 ONLY : gstart
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE mp_global,             ONLY : my_image_id,inter_image_comm,intra_bgrp_comm
  USE buffers,               ONLY : get_buffer
  USE westcom,               ONLY : iuwfc,lrwfc,nbndval0x,nbnd_occ,dvg_exc,ev,westpp_range,&
                                  & westpp_n_liouville_to_use,westpp_save_dir,westpp_l_spin_flip,&
                                  & logfile
  USE mp_world,              ONLY : mpime,root
  USE plep_db,               ONLY : plep_db_read
  USE distribution_center,   ONLY : pert,kpt_pool,band_group
  USE class_idistribute,     ONLY : idistribute
  USE types_bz_grid,         ONLY : k_grid
  USE json_module,           ONLY : json_file,json_core,json_value
  USE wvfct,                 ONLY : nbnd
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
  USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: iks,iks_do,iexc,lexc,ig,iocc,iemp,iaux,iunit
  INTEGER :: trans(2)
  INTEGER :: barra_load
  INTEGER :: nbndx_occ,nbndx_emp
  CHARACTER(5) :: label_exc
  CHARACTER(5) :: label_k
  CHARACTER(9) :: label_d
  REAL(DP) :: reduce
  REAL(DP), ALLOCATABLE :: projection_matrix(:,:,:,:)
  TYPE(bar_type) :: barra
  TYPE(json_file) :: json
  TYPE(json_core) :: jcor
  TYPE(json_value), POINTER :: jval
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  nbndx_occ = MAXVAL(nbnd_occ)
  nbndx_emp = nbnd - MINVAL(nbnd_occ)
  !
  IF(nbndx_emp < 1) THEN
     CALL errore('do_exc_comp', 'Eigenvectors decomposition need KS empty stats, rerun pwscf with nbnd > nbnd_occ', 1)
  ENDIF
  !
  IF(westpp_range(1) >= nbndx_occ) THEN
     CALL errore('do_exc_comp', 'westpp_range(1) should be smaller than the number of occupied bands', 1)
  ENDIF
  !
  IF(westpp_range(2) <= nbndx_occ) THEN
     CALL errore('do_exc_comp', 'westpp_range(2) should be greater than the number of occupied bands', 1)
  ENDIF
  !
  IF(westpp_range(2) > nbnd) THEN
     CALL errore('do_exc_comp', 'westpp_range(2) should not exceed the number of bands', 1)
  ENDIF
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
  ALLOCATE(projection_matrix(nks,nbndx_occ,nbndx_emp,westpp_n_liouville_to_use))
  !
  CALL io_push_title('BSE/TDDFT Excited State De(C)omposition')
  !
  !barra_load = 0
  !DO lexc = 1,pert%nloc
  !   iexc = pert%l2g(lexc)
  !   IF(iexc < westpp_range(1) .OR. iexc > westpp_range(2)) CYCLE
  !   barra_load = barra_load+1
  !ENDDO
  !
  !CALL start_bar_type(barra,'westpp',barra_load*k_grid%nps*SUM(nbnd_occ))
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
#else
           IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks)
           CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
        ENDIF
        !
#if defined(__CUDA)
        CALL using_evc(2)
        CALL using_evc_d(0)
#endif
        !
        IF(westpp_l_spin_flip) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        !
        DO iocc = 1, nbnd_occ(iks_do)
           !
           DO iemp = 1, (nbnd - nbnd_occ(iks))
              !
              IF(gamma_only) THEN
                 !
                 reduce = 0._DP
                 !$acc loop reduction(+:reduce)
                 DO ig = 1, npw
                    reduce = reduce + 2._DP*REAL(dvg_exc(ig,iocc,iks,lexc),KIND=DP)*REAL(evc_work(ig,nbnd_occ(iks)+iemp),KIND=DP) &
                    &               + 2._DP*AIMAG(dvg_exc(ig,iocc,iks,lexc))*AIMAG(evc_work(ig,nbnd_occ(iks)+iemp))
                 ENDDO
                 !
                 IF (gstart==2) THEN
                    reduce = reduce - REAL(dvg_exc(1,iocc,iks,lexc),KIND=DP)*REAL(evc_work(1,nbnd_occ(iks)+iemp),KIND=DP)
                 ENDIF
                 !
              ELSE
                 !
                 CALL errore('do_exc_comp', 'gamma_only is required', 1)
                 !
              ENDIF
              !
              projection_matrix(iks,iocc,iemp,iexc) = reduce
              !
              !CALL update_bar_type(barra,'westpp',1)
              !
           ENDDO
           !
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(projection_matrix(:,:,:,:),intra_bgrp_comm)
  CALL mp_sum(projection_matrix(:,:,:,:),inter_image_comm)
  !
  ! ... Print out results  
  !
  WRITE(stdout, "( /,5x,' *-------------* THE PRINCIPLE PROJECTED COMPONENTS *-------------*')")
  !
  DO iexc = 1,westpp_n_liouville_to_use
     !
     WRITE(stdout, &
         & "(   /, 5x, ' #     Exciton : | ', i8,' |','   ','Excitation energy : | ', f12.6)") &
         iexc, ev(iexc)
     !
     DO iks = 1,nks
        IF(westpp_l_spin_flip) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        DO iocc = 1,nbnd_occ(iks_do)
           DO iemp = 1,(nbnd - nbnd_occ(iks))
              reduce = ABS(projection_matrix(iks,iocc,iemp,iexc))
              IF (reduce >= 0.1) THEN
                 WRITE(stdout, "(6x, i8, 4x, i8, 16x, i8, 7x, f12.6, 2x)" ) iks,iocc,iemp,reduce
              ENDIF
           ENDDO
        ENDDO
     ENDDO
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
              trans(1) = iocc
              DO iemp = 1,(westpp_range(2) - nbnd_occ(iks))
                 trans(2) = iemp
                 iaux = iaux+1
                 WRITE(label_d,'(I9)') iaux
                 !
                 CALL json%add('output.E'//label_exc//'.K'//label_k//'.projection('//label_d//').trans',trans)
                 IF(gamma_only) THEN
                    reduce = projection_matrix(iks,iocc,iemp,iexc)
                    CALL json%add('output.E'//label_exc//'.K'//label_k//'.projection('//label_d//').value',reduce)
                 ENDIF
                 !
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  !CALL stop_bar_type(barra,'westpp')
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
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
