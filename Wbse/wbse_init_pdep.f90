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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
SUBROUTINE wbse_init_methods ()
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : isk,nks
  USE wavefunctions_module, ONLY : evc
  USE westcom,              ONLY : lrwfc,iuwfc,ev,dvg,n_pdep_eigen,npwqx,&
                                   nbnd_occ,&
                                   wbse_init_calculation,which_spin_channel
  USE lsda_mod,             ONLY : nspin
  USE pdep_db,              ONLY : pdep_db_read
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE buffers,              ONLY : get_buffer
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : pert
  !wbsecom combined into westcom
  !USE wbsecom,              ONLY : wbse_init_calculation
  !USE wbsecom,              ONLY : which_spin_channel
  !
  IMPLICIT NONE
  !
  !LOGICAL,  INTENT(IN) :: ff_activate
  !
  INTEGER :: iks, current_spin
  INTEGER :: iq, nkq, ikq
  REAL(DP):: xq(3)
  LOGICAL :: l_restart_calc, spin_resolve
  !
  SELECT CASE(wbse_init_calculation)
  CASE('r','R')
     !
     ! RESTART
     !
     l_restart_calc = .True.
     !
  CASE('s','S')
     !
     ! FROM SCRATCH
     !
     l_restart_calc = .False.
     !
  CASE DEFAULT
     CALL errore('Wbse_init', 'Wrong wstat_init_calculation',1)
  END SELECT
  !
  !IF (.NOT.ff_activate) THEN
     !
     ! Activate band group parallel PDEP
     !
  !   pert=idistribute()
  !   CALL pert%init(n_pdep_eigen, 'B','nvecx',.TRUE.)
     !
  !ENDIF
  !
  ALLOCATE( dvg( npwqx, pert%nlocx ) )
  ALLOCATE( ev( n_pdep_eigen ) )
  !
  spin_resolve = ((which_spin_channel > 0).AND.(nspin>1))
  !
  nkq = 1
  !
  DO iq = 1, nkq
     !
     xq(:)  = 0.0_DP
     !
     !IF (.NOT.ff_activate) THEN
        !
     !   CALL pdep_db_read( n_pdep_eigen )
        !
     !ENDIF
     !
     DO iks = 1, nks
        !
        current_spin = isk(iks)
        !
        ikq    = 1 !grid_ikq (iq,iks)
        !
        ! ... read in GS wavefunctions iks
        !
        IF (nkq > 1) THEN
           !
       !    IF (my_image_id==0) CALL get_buffer (evc, lrwfc, iuwfc, iks)
       !    IF (my_image_id==0) CALL get_buffer (evq, lrwfc, iuwfc, ikq)
       !    CALL mp_bcast(evc,0,inter_image_comm)
       !    CALL mp_bcast(evq,0,inter_image_comm)
           !
        ELSE
           !
           IF (nks > 1) THEN
              !
              IF (my_image_id==0)  CALL get_buffer (evc, lrwfc, iuwfc, iks)
              CALL mp_bcast(evc,0,inter_image_comm)
              !
           ENDIF
           !
           IF (spin_resolve) THEN
              !
              IF (current_spin == which_spin_channel) THEN
                 !
                 !IF (ff_activate) THEN
                    !
                    CALL wbse_init_qboxcoupling_single_q(iks,iq,xq,current_spin,&
                                                   nbnd_occ(iks),l_restart_calc)
                    !
                 !ELSE
                    !
                 !   CALL wbse_init_pdep_single_q(iks,iq,xq,current_spin,nbnd_occ(iks),&
                 !                                 n_pdep_eigen,dvg,ev,l_restart_calc)
                    !
                 !ENDIF
                 !
              ENDIF
              !
           ELSE
              !
              !IF (ff_activate) THEN
                 !
                 CALL wbse_init_qboxcoupling_single_q(iks,iq,xq,current_spin, &
                                                nbnd_occ(iks),l_restart_calc)
                 !
              !ELSE
                 !
              !   CALL wbse_init_pdep_single_q(iks,iq,xq,current_spin,nbnd_occ(iks),&
              !                                 n_pdep_eigen,dvg,ev,l_restart_calc)
                 !
              !ENDIF
              !
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
ENDSUBROUTINE
!
!!----------------------------------------------------------------------------
!SUBROUTINE wbse_init_pdep_single_q (iks,ikq,xq,current_spin,nbndval,n_pdep_eigen,dvg,ev,l_restart_calc)
!  !----------------------------------------------------------------------------
!  !
!  USE kinds,                ONLY : DP
!  USE io_global,            ONLY : stdout
!  USE constants,            ONLY : e2, fpi
!  USE cell_base,            ONLY : alat, tpiba2, omega
!  USE gvect,                ONLY : nl,ngm,g,nlm,gstart
!  USE io_push,              ONLY : io_push_title
!  !sqvc not in westcom pot3D%sqvc  TODO: pot3d init
!  USE types_coulomb,         ONLY : pot3D
!  !USE westcom,              ONLY : wstat_save_dir,sqvc,fftdriver,npwq,npwqx
!  USE westcom,              ONLY : wbse_init_save_dir,fftdriver,npwq,npwqx
!  USE pwcom,                ONLY : omega
!  USE control_flags,        ONLY : gamma_only
!  USE wavefunctions_module, ONLY : evc,psic
!  USE fft_base,             ONLY : dfftp,dffts
!  USE fft_at_gamma,         ONLY : double_fwfft_gamma,double_invfft_gamma
!  USE pwcom,                ONLY : wk,nks,nelup,neldw,isk,g,igk_k,ngm,tpiba2,xk,omega,npw,npwx,lsda,nkstot,&
!                                 & current_k,ngk,nbnd,wg
!  USE mp_bands,             ONLY : intra_bgrp_comm
!  USE mp_global,            ONLY : inter_image_comm
!  USE mp,                   ONLY : mp_sum
!  USE bse_module,           ONLY : ovl_thr, l_wannier_repr
!  USE wbse_init_restart,    ONLY : wbse_stat_restart_read, wbse_stat_restart_write
!  USE wbse_init_restart,    ONLY : wbse_index_matrix_read, wbse_index_matrix_write
!  USE wbse_init_restart,    ONLY : wbse_pdep_coeffie_read, wbse_pdep_coeffie_write
!  USE wbse_init_restart,    ONLY : wbse_init_restart_para_read, wbse_init_restart_para_write
!  USE class_idistribute,    ONLY : idistribute
!  USE distribution_center,  ONLY : aband, bseparal
!  USE distribution_center,  ONLY : pert
!  !
!  IMPLICIT NONE
!  !
!  INTEGER,     INTENT(IN) :: iks,ikq,current_spin,nbndval,n_pdep_eigen
!  REAL(DP),    INTENT(IN) :: xq(3),ev(n_pdep_eigen)
!  COMPLEX(DP), INTENT(IN) :: dvg(npwqx, pert%nlocx)
!  LOGICAL,     INTENT(IN) :: l_restart_calc
!  !
!  INTEGER :: ibnd, jbnd, tmp_size, alnd
!  INTEGER :: il1, ig1, ir, ig, do_index
!  REAL(DP):: ovl_value, qg2,alpha_ija_vx, alpha_ija_vc
!  REAL(DP),    ALLOCATABLE  :: rho_aux(:)
!  COMPLEX(DP), ALLOCATABLE  :: phi_c(:), phi_x(:), psic_aux(:)
!  COMPLEX(DP), ALLOCATABLE  :: evc_loc(:,:)
!  REAL(DP),    ALLOCATABLE  :: ovl_matrix(:,:)
!  REAL(DP),    ALLOCATABLE  :: restart_matrix(:)
!  REAL(DP),    ALLOCATABLE  :: index_matrix(:,:)
!  REAL(DP),    ALLOCATABLE  :: alpha_list_vx(:)
!  REAL(DP),    ALLOCATABLE  :: alpha_list_vc(:)
!  !
!  CHARACTER(LEN=256)        :: filename
!  CHARACTER(LEN=6)          :: my_labeliq, my_labelik
!  CHARACTER(LEN=1)          :: my_spin
!  !
!  REAL(DP), EXTERNAL        :: get_clock
!  CHARACTER(20),EXTERNAL    :: human_readable_time
!  REAL(DP)                  :: wtime(2)
!  !
!  LOGICAL                   :: calc_is_done
!  REAL(kind=dp), EXTERNAL   :: DDOT
!  !
!  CALL wbse_init_memory_report()
!  !
!  CALL start_clock( "wbse_stat_pdep" )
!  !
!  CALL io_push_title( "Wbse_init for CHI represented with PDEP" )
!  !
!  IF (.NOT. ALLOCATED(psic)) ALLOCATE (psic(dffts%nnr))
!  WRITE(my_labeliq,'(i6.6)') ikq
!  WRITE(my_labelik,'(i6.6)') iks
!  WRITE(my_spin,'(i1)') current_spin
!  !
!  ! Set up a image parallel over the bands
!  !
!  aband = idistribute()
!  CALL aband%init(nbndval,'i','bse_nbndval',.TRUE.)
!  !
!  tmp_size = nbndval*nbndval*n_pdep_eigen
!  ALLOCATE (index_matrix(tmp_size,3))
!  ALLOCATE (ovl_matrix(nbndval,nbndval))
!  !
!  ovl_matrix(:,:) = 0.0_DP
!  index_matrix(:,:) = 0.0_DP
!  !
!  IF (l_wannier_repr) THEN
!     !
!     ALLOCATE (evc_loc(npwx,nbndval))
!     !
!     CALL bse_do_localization (current_spin, nbndval, evc_loc, ovl_matrix, l_restart_calc)
!     !
!  ENDIF
!  !
!  IF (.NOT. l_restart_calc) THEN
!     !
!     do_index = 0
!     !
!     DO alnd = 1, n_pdep_eigen
!        !
!        DO ibnd = 1, nbndval
!           !
!           DO jbnd = 1, nbndval
!              !
!              ovl_value = ovl_matrix(ibnd,jbnd)
!              !
!              IF (ovl_value >= ovl_thr ) THEN
!                 !
!                 IF (gamma_only) THEN
!                    !
!                    IF (jbnd >= ibnd)  THEN
!                       !
!                       do_index = do_index + 1
!                       !
!                       index_matrix(do_index,1) = ibnd
!                       index_matrix(do_index,2) = jbnd
!                       index_matrix(do_index,3) = alnd
!                       !
!                    ENDIF
!                    !
!                 ELSE
!                    !
!                    do_index = do_index + 1
!                    !
!                    index_matrix(do_index,1) = ibnd
!                    index_matrix(do_index,2) = jbnd
!                    index_matrix(do_index,3) = alnd
!                    !
!                 ENDIF
!                 !
!              ENDIF
!              !
!           ENDDO
!           !
!        ENDDO
!        !
!     ENDDO
!     !
!     filename = TRIM( wbse_init_save_dir )//"/index_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
!                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
!     CALL wbse_index_matrix_write(filename,do_index,3,index_matrix(1:do_index,:))
!     !
!  ELSE
!     !
!     filename = TRIM( wbse_init_save_dir )//"/index_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
!                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
!     CALL wbse_index_matrix_read (filename,tmp_size,do_index,3,index_matrix)
!     !
!  ENDIF
!  !
!  ALLOCATE (restart_matrix(do_index))
!  ALLOCATE ( alpha_list_vx(do_index))
!  ALLOCATE ( alpha_list_vc(do_index))
!  !
!  restart_matrix(:) = 0.0_DP
!  alpha_list_vc(:)  = 0.0_DP
!  alpha_list_vx(:)  = 0.0_DP
!  !
!  calc_is_done = .FALSE.
!  IF (l_restart_calc) THEN
!     !
!     filename = TRIM( wbse_init_save_dir )//"/restart_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
!                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
!     CALL wbse_stat_restart_read (filename,do_index,restart_matrix,calc_is_done)
!     !
!  ENDIF
!  !
!  IF (calc_is_done) GOTO 2222
!  !
!  ! initialize the paralellization
!  !
!  bseparal = idistribute()
!  CALL bseparal%init(do_index,'i','number_triples',.TRUE.)
!  !
!  IF (l_restart_calc) THEN
!     !
!     CALL wbse_init_restart_para_read (do_index,alpha_list_vx,alpha_list_vc)
!     !
!  ENDIF
!  !
!  ! parallel loop
!  !
!  DO il1 = 1, bseparal%nlocx
!     !
!     ig1  = bseparal%l2g(il1) ! global index of n_total
!     !
!     ibnd = INT(index_matrix(ig1,1))
!     jbnd = INT(index_matrix(ig1,2))
!     alnd = INT(index_matrix(ig1,3))
!     !
!     IF ((ig1 < 1).OR.(ig1 > do_index)) GOTO 1111
!     !
!     IF (l_restart_calc) THEN
!        !
!        IF (INT(restart_matrix(ig1)) > 0) GOTO 1111
!        !
!     ENDIF
!     !
!     IF ((ig1 < 1).OR.(ig1 > do_index)) GOTO 1111
!     !
!     ! Here is the main part of this code
!     !
!     ALLOCATE(phi_c(npwqx),phi_x(npwqx),psic_aux(dffts%nnr))
!     !
!     psic_aux(:) = (0.d0,0.d0)
!     psic(:) = (0.d0,0.d0)
!     !
!     ! 1) compute phi_c = phi_i * sqrt(v_c*)
!     !    G->0 : v_c(G) = 0.0
!     !
!     phi_c(:) = (0.d0,0.d0)
!     DO ig = 1, npwqx
!        !
!        qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
!        !
!        IF (qg2 > 1.d-8) then
!           !
!           ! 1) compute phi_c = phi_i * sqrt(v_c*)
!           !    G->0 : v_c(G) = 0.0
!           !
!           phi_c(ig) = dvg(ig,alnd) * DSQRT(e2*fpi/(tpiba2*qg2))
!           !
!        ENDIF
!        !
!     ENDDO
!     !
!     ! 2) compute phi_x = phi_i * sqrt(v_c)
!     ! using F-B correction
!     !
!     phi_x(:) = (0.d0,0.d0)
!     phi_x(:) = dvg(:,alnd) * pot3D%sqvc(:)
!     !
!     CALL double_invfft_gamma(dffts,npwq,npwqx,phi_x,phi_c,psic_aux,TRIM(fftdriver))
!     !
!     IF (l_wannier_repr) THEN
!        !
!        CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
!        !
!     ELSE
!        !
!        CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
!        !
!     ENDIF
!     !
!     alpha_ija_vx  = 0.0_DP
!     alpha_ija_vc  = 0.0_DP
!     alpha_ija_vx  = SUM(DBLE(psic(:)) * AIMAG(psic(:)) *  DBLE(psic_aux(:)))
!     alpha_ija_vc  = SUM(DBLE(psic(:)) * AIMAG(psic(:)) * AIMAG(psic_aux(:))) &
!                                       * ev(alnd)/(1.0_DP-ev(alnd))
!     !
!     alpha_ija_vx  = alpha_ija_vx * wg(ibnd,iks)/(dffts%nr1*dffts%nr2*dffts%nr3)
!     alpha_ija_vc  = alpha_ija_vc * wg(ibnd,iks)/(dffts%nr1*dffts%nr2*dffts%nr3)
!     !
!     CALL mp_sum(alpha_ija_vx,intra_bgrp_comm)
!     CALL mp_sum(alpha_ija_vc,intra_bgrp_comm)
!     !
!     alpha_list_vx(ig1) = alpha_ija_vx
!     alpha_list_vc(ig1) = alpha_ija_vc
!     !
!     restart_matrix(ig1)  = 1.0
!     !
!     DEALLOCATE (phi_c, phi_x, psic_aux)
!     !
!1111 CONTINUE
!     !
!     CALL wbse_init_restart_para_write(do_index,alpha_list_vx,alpha_list_vc)
!     !
!     ! for restarting, update status of restart_matrix
!     !
!     CALL mp_sum(restart_matrix(1:do_index), inter_image_comm)
!     !
!     DO ir = 1, do_index
!        !
!        IF (restart_matrix(ir) > 0) restart_matrix(ir) = 1.0
!        !
!     ENDDO
!     !
!     calc_is_done = .FALSE.
!     filename = TRIM( wbse_init_save_dir )//"/restart_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
!                TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
!     CALL wbse_stat_restart_write (filename,do_index,restart_matrix,calc_is_done)
!     !
!  ENDDO
!  !
!  CALL mp_sum(alpha_list_vx(1:do_index), inter_image_comm)
!  CALL mp_sum(alpha_list_vc(1:do_index), inter_image_comm)
!  !
!  filename = TRIM( wbse_init_save_dir )//"/PDEP_vc_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
!          TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
!  CALL wbse_pdep_coeffie_write (filename,do_index,alpha_list_vx,alpha_list_vc)
!  !
!  calc_is_done = .TRUE.
!  filename = TRIM( wbse_init_save_dir )//"/restart_matrix_iq"//TRIM(ADJUSTL(my_labeliq))//"_ik"//&
!             TRIM(ADJUSTL(my_labelik))//"_spin"//TRIM(ADJUSTL(my_spin))//".dat"
!  CALL wbse_stat_restart_write (filename,do_index,restart_matrix,calc_is_done)
!  !
!2222 CONTINUE
!  !
!  DEALLOCATE (index_matrix)
!  DEALLOCATE (restart_matrix)
!  DEALLOCATE (ovl_matrix)
!  DEALLOCATE (alpha_list_vx)
!  DEALLOCATE (alpha_list_vc)
!  !
!  IF (ALLOCATED(psic))    DEALLOCATE(psic)
!  IF (ALLOCATED(evc_loc)) DEALLOCATE(evc_loc)
!  !
!  CALL stop_clock( 'wbse_stat_pdep' )
!  !
!  RETURN
!  !
!ENDSUBROUTINE
