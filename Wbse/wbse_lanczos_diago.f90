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
!-----------------------------------------------------------------------
SUBROUTINE wbse_lanczos_diago()
  !---------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE lsda_mod,             ONLY : nspin
  USE pwcom,                ONLY : npw,npwx,ngk,nks,isk,current_spin
  USE westcom,              ONLY : nbnd_occ,lrwfc,iuwfc,nbnd_occ,wbse_calculation,d0psi,wbse_ipol,&
                                 & n_lanczos,beta_store,zeta_store,nbndval0x,n_trunc_bands,&
                                 & n_steps_write_restart
  USE lanczos_restart,      ONLY : lanczos_restart_write,lanczos_restart_read,&
                                 & lanczos_restart_clear,lanczos_log
  USE mp_global,            ONLY : inter_image_comm,nimage,my_image_id,me_image,nbgrp
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer
  USE distribution_center,  ONLY : pert,aband
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d
  USE wvfct_gpum,           ONLY : using_et,using_et_d
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu,allocate_bse_gpu,&
                                 & deallocate_bse_gpu,reallocate_ps_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  LOGICAL :: l_from_scratch
  INTEGER :: lp,ip,iip,pol_index,nipol_input
  INTEGER :: iter
  INTEGER :: iks,is,nbndval,ig,lbnd,ibnd
  INTEGER :: ilan_restart,ilan_stopped,ipol_restart,ipol_stopped
  INTEGER :: nbnd_do
  INTEGER, PARAMETER :: n_ipol = 3
  INTEGER, ALLOCATABLE :: pol_index_input(:)
  CHARACTER(LEN=3), ALLOCATABLE :: pol_label_input(:)
  REAL(DP) :: factor
  REAL(DP) :: beta(nspin)
  COMPLEX(DP) :: dotp(nspin)
  COMPLEX(DP), ALLOCATABLE :: evc1(:,:,:),evc1_old(:,:,:),evc1_new(:,:,:)
  TYPE(bar_type) :: barra
  !
  SELECT CASE(wbse_ipol)
  CASE('XX','xx')
     nipol_input = 1
     ALLOCATE(pol_index_input(1))
     ALLOCATE(pol_label_input(1))
     pol_index_input(1) = 1
     pol_label_input(1) = 'XX'
  CASE('YY','yy')
     nipol_input = 1
     ALLOCATE(pol_index_input(1))
     ALLOCATE(pol_label_input(1))
     pol_index_input(1) = 2
     pol_label_input(1) = 'YY'
  CASE('ZZ','zz')
     nipol_input = 1
     ALLOCATE(pol_index_input(1))
     ALLOCATE(pol_label_input(1))
     pol_index_input(1) = 3
     pol_label_input(1) = 'ZZ'
  CASE('XYZ','xyz')
     nipol_input = 3
     ALLOCATE(pol_index_input(3))
     ALLOCATE(pol_label_input(3))
     pol_index_input(1) = 1
     pol_label_input(1) = 'XX'
     pol_index_input(2) = 2
     pol_label_input(2) = 'YY'
     pol_index_input(3) = 3
     pol_label_input(3) = 'ZZ'
  CASE DEFAULT
     CALL errore('wbse_lanczos_diago','wrong wbse_ipol',1)
  END SELECT
  !
  ! ... DISTRIBUTE nipol_input
  !
  IF(MOD(nipol_input,nimage) /= 0) CALL errore('wbse_lanczos_diago','bad nimage',1)
  !
  pert = idistribute()
  CALL pert%init(nipol_input,'i','nipol_input',.TRUE.)
  !
  ! ... DISTRIBUTE nbndval
  !
  IF(nbgrp > nbndval0x-n_trunc_bands) CALL errore('wbse_lanczos_diago','nbgrp>nbndval',1)
  !
  aband = idistribute()
  CALL aband%init(nbndval0x-n_trunc_bands,'b','nbndval',.TRUE.)
  !
  CALL wbse_memory_report()
  !
#if defined(__CUDA)
  CALL allocate_gpu()
  CALL allocate_bse_gpu(aband%nloc)
  !
  CALL using_et(2)
  CALL using_et_d(0)
  IF(nks == 1) THEN
     CALL using_evc(2)
     CALL using_evc_d(0)
  ENDIF
#endif
  !
  ! Main Lanczos code
  !
  ALLOCATE(beta_store(n_lanczos,nipol_input,nspin))
  ALLOCATE(zeta_store(n_lanczos,n_ipol,nipol_input,nspin))
  !
  beta_store(:,:,:) = 0._DP
  zeta_store(:,:,:,:) = 0._DP
  !
  ALLOCATE(d0psi(npwx,aband%nlocx,nks,n_ipol))
  ALLOCATE(evc1(npwx,aband%nlocx,nks))
  ALLOCATE(evc1_old(npwx,aband%nlocx,nks))
  ALLOCATE(evc1_new(npwx,aband%nlocx,nks))
  !$acc enter data create(d0psi,evc1,evc1_old,evc1_new)
  !
  SELECT CASE(wbse_calculation)
  CASE('l')
     !
     IF(nimage > 1) CALL errore('wbse_lanczos_diago','restart with nimage>1 not implemented',1)
     !
     ! RESTART
     !
     CALL lanczos_restart_read(nipol_input,ipol_stopped,ilan_stopped,evc1,evc1_old)
     !
     !$acc update device(evc1,evc1_old)
     !
     ipol_restart = ipol_stopped
     ilan_restart = ilan_stopped+1
     !
     l_from_scratch = .FALSE.
     !
  CASE('L')
     !
     ! FROM SCRATCH
     !
     ipol_restart = 1
     ilan_restart = 1
     !
     l_from_scratch = .TRUE.
     !
  CASE DEFAULT
     CALL errore('wbse_lanczos_diago','wrong wlzcos_calculation',1)
  END SELECT
  !
  CALL solve_e_psi()
  !
  !$acc update host(d0psi)
  !
  CALL io_push_title('Lanczos linear-response absorption spectrum calculation')
  WRITE(stdout,'(5x,"Using Tamm-Dancoff Liouvillian operator")')
  !
  polarization_loop : DO lp = ipol_restart,pert%nloc
     !
     ip = pert%l2g(lp)
     pol_index = pol_index_input(ip)
     !
     IF(l_from_scratch) THEN
        IF(nimage == 3) THEN
           CALL io_push_title('Starting new Lanczos loops (XX,YY,ZZ)')
        ELSE
           CALL io_push_title('Starting new Lanczos loop ('//TRIM(pol_label_input(ip))//')')
        ENDIF
        !
        !$acc kernels present(evc1_old,evc1_new,evc1,d0psi)
        evc1_old(:,:,:) = (0._DP,0._DP)
        evc1_new(:,:,:) = (0._DP,0._DP)
        evc1(:,:,:) = d0psi(:,:,:,pol_index)
        !$acc end kernels
     ELSE
        CALL io_push_title('Retarting Lanczos loop ('//TRIM(pol_label_input(ip))//')')
     ENDIF
     !
     CALL start_bar_type(barra,'lan_diago',n_lanczos-ilan_restart+1)
     !
     ! Loop on the Lanczos iterations
     !
     lancz_loop : DO iter = ilan_restart,n_lanczos
        !
        ! Application of the Liouvillian superoperator
        !
        CALL west_apply_liouvillian(evc1,evc1_new)
        !
        ! By construction <p|Lq>=0 should be 0, forcing this both conserves
        ! resources and increases stability.
        ! ( i.e., alpha = 0 )
        !
        ! Orthogonality requirement: <v|\bar{L}|v> = 1
        !
        CALL wbse_dot(evc1,evc1_new,aband%nlocx,nks,dotp)
        !
        beta(:) = REAL(dotp,KIND=DP)
        !
        ! beta<0 is a serious error for the pseudo-Hermitian algorithm
        !
        DO is = 1,nspin
           IF(beta(is) < 0._DP) THEN
              CALL errore('wbse_lanczos_diago','negative beta',1)
           ELSE
              beta(is) = SQRT(beta(is))
           ENDIF
        ENDDO
        !
        beta_store(iter,ip,:) = beta
        !
        ! Renormalize q(i) and Lq(i)
        !
        DO iks = 1,nks
           !
           npw = ngk(iks)
           current_spin = isk(iks)
           factor = 1._DP/beta(current_spin)
           nbnd_do = aband%nlocx
           !
           !$acc parallel loop collapse(2) present(evc1)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 evc1(ig,lbnd,iks) = factor*evc1(ig,lbnd,iks)
              ENDDO
           ENDDO
           !$acc end parallel
           !
           !$acc parallel loop collapse(2) present(evc1_new)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 evc1_new(ig,lbnd,iks) = factor*evc1_new(ig,lbnd,iks)
              ENDDO
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        ! Calculation of zeta coefficients.
        ! See Eq.(35) in Malcioglu et al., Comput. Phys. Commun. 182, 1744 (2011).
        !
        IF(MOD(iter,2) == 0) THEN
           DO iip = 1,n_ipol
              CALL wbse_dot(d0psi(:,:,:,iip),evc1,aband%nlocx,nks,dotp)
              !
              zeta_store(iter,iip,ip,:) = dotp
           ENDDO
        ELSE
           DO iip = 1,n_ipol
              zeta_store(iter,iip,ip,:) = 0._DP
           ENDDO
        ENDIF
        !
        DO iks = 1,nks
           !
           npw = ngk(iks)
           current_spin = isk(iks)
           factor = beta(current_spin)
           nbnd_do = aband%nlocx
           !
           !$acc parallel loop collapse(2) present(evc1_new)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 evc1_new(ig,lbnd,iks) = evc1_new(ig,lbnd,iks)-factor*evc1_old(ig,lbnd,iks)
              ENDDO
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        ! Apply P_c|evc1_new>
        !
        DO iks = 1,nks
           !
           nbndval = nbnd_occ(iks)
           npw = ngk(iks)
           !
           nbnd_do = 0
           DO lbnd = 1,aband%nloc
              ibnd = aband%l2g(lbnd)+n_trunc_bands
              IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
           ENDDO
           !
           ! ... Read GS wavefunctions
           !
           IF(nks > 1) THEN
              IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
              CALL mp_bcast(evc,0,inter_image_comm)
              !
#if defined(__CUDA)
              CALL using_evc(2)
              CALL using_evc_d(0)
#endif
           ENDIF
           !
#if defined(__CUDA)
           CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
           CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,evc1_new(:,:,iks),(1._DP,0._DP))
           !
        ENDDO
        !
        ! Throw away q(i-1),and make q(i+1) to be the current vector,
        ! be ready for the next iteration. evc1_new will be free again after this step
        !
        !$acc kernels present(evc1_old,evc1,evc1_new)
        evc1_old(:,:,:) = evc1
        evc1(:,:,:) = evc1_new
        !$acc end kernels
        !
        IF(n_steps_write_restart > 0 .AND. MOD(iter,n_steps_write_restart) == 0 &
           & .AND. nimage == 1) THEN
           !$acc update host(evc1,evc1_old)
           !
           CALL lanczos_restart_write(nipol_input,ip,iter,evc1,evc1_old)
        ENDIF
        !
        CALL update_bar_type(barra,'lan_diago',1)
        !
     ENDDO lancz_loop
     !
     IF(nimage == 1) CALL lanczos_log(ip,pol_label_input(ip))
     !
     ilan_restart = 1
     l_from_scratch = .TRUE.
     !
     CALL stop_bar_type(barra,'lan_diago')
     !
  ENDDO polarization_loop
  !
  IF(nimage > 1) THEN
     IF(me_image == 0) THEN
        CALL mp_sum(beta_store,inter_image_comm)
        CALL mp_sum(zeta_store,inter_image_comm)
     ENDIF
     !
     DO ip = 1,nipol_input
        CALL lanczos_log(ip,pol_label_input(ip))
     ENDDO
  ENDIF
  !
  CALL lanczos_restart_clear()
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
  CALL deallocate_bse_gpu()
#endif
  !
  DEALLOCATE(pol_index_input)
  DEALLOCATE(pol_label_input)
  !$acc exit data delete(d0psi,evc1,evc1_old,evc1_new)
  DEALLOCATE(d0psi)
  DEALLOCATE(evc1)
  DEALLOCATE(evc1_new)
  DEALLOCATE(evc1_old)
  DEALLOCATE(beta_store)
  DEALLOCATE(zeta_store)
  !
END SUBROUTINE
