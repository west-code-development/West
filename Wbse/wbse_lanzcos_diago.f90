!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_lanzcos_diago ()
  !---------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : tmp_dir
  USE control_flags,        ONLY : gamma_only
  USE wvfct,                ONLY : npwx
  USE lsda_mod,             ONLY : nspin
  USE pwcom,                ONLY : nks,isk,current_spin
  USE westcom,              ONLY : nbnd_occ,west_prefix,lrwfc,iuwfc,nbnd_occ
  USE wbsecom,              ONLY : wlz_calculation, &
                                   d0psi,ipol_input,n_lzstep,&
                                   alpha_store,beta_store,& 
                                   gamma_store,zeta_store
  USE lanzcos_db,           ONLY : lanzcos_d0psi_read,&
                                   lanzcos_d0psi_write,&
                                   lanzcos_evcs_write,&
                                   lanzcos_evcs_read 
  USE lanzcos_restart,      ONLY : lanzcos_restart_write,&
                                   lanzcos_restart_read, &
                                   lanzcos_postpro_write 
  !
  USE bse_module,           ONLY : bse_calc,size_index_matrix_lz,&
                                   bseparal
  USE io_files,             ONLY : tmp_dir
  USE wbsecom,              ONLY : nbndval0x
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_bcast,mp_barrier
  USE wavefunctions_module, ONLY : evc
  USE buffers,              ONLY : get_buffer
  USE distribution_center,  ONLY : aband
  USE class_idistribute,    ONLY : idistribute
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  LOGICAL     :: l_from_scratch
  INTEGER     :: ip,iip,pol_index,nipol_input
  INTEGER     :: iteration, lz_iteration
  INTEGER     :: iks,is, nbndval
  INTEGER     :: iter_restart, ipol_restart
  INTEGER     :: lriter_restart, pliter_restart 
  INTEGER     :: pliter_stop, lriter_stop 
  INTEGER     :: size_index_matrix
  INTEGER, PARAMETER   :: n_ipol = 3 
  INTEGER, ALLOCATABLE :: pol_index_input(:)
  CHARACTER(LEN=3), ALLOCATABLE :: pol_label_input(:) 
  REAL(DP)    :: alpha(nspin),beta(nspin),gamma(nspin)
  COMPLEX(DP) :: zeta(nspin), wbse_dot_out(nspin)
  COMPLEX(DP), ALLOCATABLE :: evc1(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: evc1_old(:,:,:), evc1_new(:,:,:)
  !
  CHARACTER(LEN=256)       :: tmp_lz
  !
  tmp_lz = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.tmp_lz' 
  !
  ! ... DISTRIBUTE lanzcos
  !
  aband = idistribute()
  !
  CALL aband%init(nbndval0x,'i','band_paralel',.TRUE.)
  ! 
  ! ... DISTRIBUTE bse_kernel
  ! 
  size_index_matrix = MAXVAL(size_index_matrix_lz(:))
  !
  IF (bse_calc) THEN
     !
     ! initialize the paralellization
     !
     bseparal  = idistribute()
     CALL bseparal%init(size_index_matrix,'i','bse_kernel',.TRUE.)
     !  
  ENDIF
  !
  ! Main Lanzcos program 
  !
  IF (n_lzstep == 0) THEN
     !
     CALL errore('lanzcosdiago', ' n_lzstep shoudl be > 0 ', 1)
     ! 
  ENDIF
  !
  SELECT CASE(ipol_input)
  ! 
  CASE('XX', 'xx')
      !   
      nipol_input = 1
      ALLOCATE(pol_index_input(1))
      ALLOCATE(pol_label_input(1))
      pol_index_input(1) = 1
      pol_label_input(1) = 'XX'
      !
  CASE('YY', 'yy')
      !
      nipol_input = 1
      ALLOCATE(pol_index_input(1))
      ALLOCATE(pol_label_input(1))
      pol_index_input(1) = 2
      pol_label_input(1) = 'YY' 
      ! 
  CASE('ZZ', 'zz')
      !
      nipol_input = 1
      ALLOCATE(pol_index_input(1))
      ALLOCATE(pol_label_input(1))
      pol_index_input(1) = 3
      pol_label_input(1) = 'ZZ'
      !
  CASE('XYZ', 'xyz')
      !
      nipol_input = 3
      ALLOCATE(pol_index_input(3))
      ALLOCATE(pol_label_input(3))
      pol_index_input(1) = 1
      pol_label_input(1) = 'XX'
      pol_index_input(2) = 2
      pol_label_input(2) = 'YY'
      pol_index_input(3) = 3
      pol_label_input(3) = 'ZZ'
      !
  CASE DEFAULT
      !
      CALL errore('lanzcosdiago', 'Wrong ipol_input', 1)
      ! 
  END SELECT
  !
  ALLOCATE (d0psi(npwx,nbndval0x,nks,n_ipol))
  ALLOCATE (alpha_store(nipol_input,n_lzstep,nspin),beta_store(nipol_input,n_lzstep,nspin))
  ALLOCATE (gamma_store(nipol_input,n_lzstep,nspin),zeta_store(nipol_input,n_ipol,n_lzstep,nspin))
  !
  ALLOCATE(evc1(npwx,nbndval0x,nks),evc1_old(npwx,nbndval0x,nks),evc1_new(npwx,nbndval0x,nks))
  ! 
  SELECT CASE(wlz_calculation)
  !
  CASE('r','R')
      !
      ! RESTART
      !
      alpha_store(:,:,:) = 0.0_DP
      beta_store(:,:,:)  = 0.0_DP
      gamma_store(:,:,:) = 0.0_DP
      zeta_store(:,:,:,:)  = (0.0_DP,0.0_DP)
      !
      CALL lanzcos_restart_read (nipol_input, pliter_stop, lriter_stop)
      !
      ! 1) read pliter_stopped
      ! 2) read lriter_stopped
      ! 3) read store_alpha, store_beta, store_gamma
      !
      pliter_restart = pliter_stop 
      lriter_restart = lriter_stop+1 
      !
      ! 4) read d0psi evc1, evc_old, saved on files
      ! 
      CALL lanzcos_d0psi_read ()
      CALL lanzcos_evcs_read (evc1, evc1_old)
      !
      l_from_scratch = .false.
      ! 
  CASE('s','S')
      !
      ! FROM SCRATCH
      ! 
      CALL solve_e_psi()
      !
      CALL lanzcos_d0psi_write ()
      !
      alpha_store(:,:,:) = 0.0_DP
      beta_store(:,:,:)  = 0.0_DP
      gamma_store(:,:,:) = 0.0_DP
      zeta_store(:,:,:,:)  = (0.0_DP,0.0_DP)
      !
      lriter_restart = 1 
      pliter_restart = 1 
      !
      l_from_scratch = .true.
      !  
  CASE DEFAULT
      !
      CALL errore('lanzcosdiago', 'Wrong wlzcos_calculation',1)
      !
  END SELECT
  ! 
  WRITE(stdout,'(/,5X,"LANCZOS LINEAR-RESPONSE ADSORPTION SPECTRUM CALCULATION")')
  WRITE(stdout,'(/,10X,"USING TAMM-DANCOFF LIOUVILLIAN OPERATOR ")')
  WRITE(stdout,'(/,5x,"Number of Lanczos iterations = ",i6)') n_lzstep
  !
  polarization_loop : DO ip = pliter_restart, nipol_input 
     !
     pol_index = pol_index_input(ip)
     !
     IF (l_from_scratch) THEN
        ! 
        WRITE(stdout,'(/5x,"***Starting new Lanczos loop at ipol: ",1x,a8)') pol_label_input(ip) 
        !
        evc1_old(:,:,:) = (0.0d0,0.0d0)
        evc1_new(:,:,:) = (0.0d0,0.0d0)
        evc1(:,:,:) = d0psi(:,:,:,pol_index)
        !
     ELSE
        ! 
        WRITE(stdout,'(/5x,"***Restarting a Lanczos loop at ipol: ",1x,a8)') pol_label_input(ip)
        !
     ENDIF 
     !
     ! Loop on the Lanczos iterations
     ! 
     lancz_loop : DO iteration = lriter_restart, n_lzstep 
        !
        lz_iteration = iteration
        !
        WRITE(stdout,'(/5x,"**Lanczos iteration:",1x,i6,3x,"at Polar:",i5,a8)') lz_iteration, pol_index 
        !
        ! Application of the Liouvillian superoperator
        !
        CALL west_apply_liouvillian (evc1(:,:,:), evc1_new(:,:,:))
        !
        ! By construction <p|Lq>=0 should be 0, forcing this both conserves 
        ! resources and increases stability.
        !
        alpha = 0.0d0
        !
        alpha_store(ip,lz_iteration,:) = alpha
        ! 
        WRITE(stdout,'(5X,"^-^alpha(",i8.8,")=",f10.6)') lz_iteration, alpha
        !
        ! Orthogonality requirement: <v|\bar{L}|v> = 1
        !
        CALL wbse_dot(evc1(:,:,:),evc1_new(:,:,:),npwx,nbndval0x,nks,wbse_dot_out)
        !
        beta(:) = dble(wbse_dot_out(:))
        !
        ! beta<0 is a serious error for the pseudo-Hermitian algorithm
        !
        DO is = 1, nspin
           !
           IF ( beta(is) < 0.0d0) THEN
              !
              ! CALL error_
              !
           ELSE ! ( beta>0.0d0 ) 
              !
              beta(is)  = sqrt(beta(is))
              gamma(is) = beta(is)
              !
           ENDIF
           !
        ENDDO
        !
        beta_store (ip,lz_iteration,:) = beta(:)
        !
        gamma_store(ip,lz_iteration,:) = gamma(:)
        !
        DO is = 1, nspin
           ! 
           WRITE(stdout,'(5X,"ispin:",i2,5X,"beta (",i8.8,")=", f12.6)') is,lz_iteration, beta(is)
           WRITE(stdout,'(5X,"ispin:",i2,5X,"gamma(",i8.8,")=", f12.6)') is,lz_iteration, gamma(is)
           !
        ENDDO
        !
        ! Renormalize q(i) and Lq(i)
        !
        DO iks = 1, nks
           !
           current_spin = isk(iks)
           !
           evc1(:,:,iks) = (1.0_DP/beta(current_spin))*evc1(:,:,iks)
           evc1_new(:,:,iks) = (1.0_DP/beta(current_spin))*evc1_new(:,:,iks)
           !
        ENDDO
        !
        ! Calculation of zeta coefficients.
        ! See Eq.(35) in Malcioglu et al., Comput. Phys. Commun. 182, 1744 (2011).
        !
        IF (mod(lz_iteration,2)==0) THEN
           !
           DO iip = 1, n_ipol
              !
              CALL wbse_dot(d0psi(:,:,:,iip), evc1(:,:,:),npwx,nbndval0x,nks, wbse_dot_out)   
              !
              zeta(:) = wbse_dot_out(:)
              !
              zeta_store (ip,iip,lz_iteration,:) = zeta(:)
              ! 
              DO is = 1, nspin
                 ! 
                 WRITE(stdout,'(5X,"ispin:",i2,5X,"zeta = ",1x,i3,i3,2(1x,f18.13))') is, ip, iip, real(zeta(is)),aimag(zeta(is))
                 ! 
              ENDDO
              !
           ENDDO
           !
        ELSE
           !
           DO iip = 1, n_ipol
              !
              zeta(:) = (0.0d0,0.0d0)
              !
              zeta_store (ip,iip,lz_iteration,:) = zeta(:)
              !
              DO is = 1, nspin
                 ! 
                 WRITE(stdout,'(5X, "ispin:",i2,5X,"zeta = ",1x,i3,i3,2(1x,f18.13))') is,ip, iip, real(zeta(is)),aimag(zeta(is))
                 !
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
        DO iks = 1, nks
           !
           current_spin = isk(iks)
           !
           evc1_new(:,:,iks) = evc1_new(:,:,iks) - gamma(current_spin)*evc1_old(:,:,iks)
           !
        ENDDO
        !
        ! Apply P_c|evc1_new>  
        ! 
        DO iks=1, nks
           !
           nbndval = nbnd_occ(iks)
           !
           ! ... read in GS wavefunctions iks 
           !
           IF (nks>1) THEN
              !
              IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
              CALL mp_bcast(evc,0,inter_image_comm)
              !
           ENDIF
           !
           CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,evc1_new(1,1,iks),(1._DP,0._DP))
           ! 
        ENDDO
        !
        ! Throw away q(i-1), and make q(i+1) to be the current vector,
        ! be ready for the next iteration. evc1_new will be free again after this step
        !
        evc1_old(:,:,:) = (0.0_DP, 0.0_DP)
        evc1_old(:,:,:) = evc1(:,:,:)
        !
        evc1(:,:,:) = (0.0_DP, 0.0_DP)
        evc1(:,:,:) = evc1_new(:,:,:)
        !
        IF (mod(lz_iteration,100)==0) THEN
           !
           IF (lz_iteration>5) THEN
              !
              CALL my_mkdir(tmp_lz)
              CALL my_copy_lz(tmp_lz)
              !
           ENDIF
           !
           CALL lanzcos_restart_write (nipol_input, ip, lz_iteration)
           CALL lanzcos_evcs_write (evc1, evc1_old)
           !
        ENDIF
        !  
     ENDDO lancz_loop
     !
     CALL lanzcos_postpro_write (nipol_input, ip, pol_label_input(ip)) 
     !
     lriter_restart = 1
     l_from_scratch = .true.
     !
  ENDDO polarization_loop 
  ! 
  WRITE(stdout,'(5x,"End of Lanczos diagonalization")')
  !
  DEALLOCATE(pol_index_input)
  DEALLOCATE(pol_label_input)
  DEALLOCATE(d0psi)
  DEALLOCATE(evc1,evc1_new,evc1_old)
  DEALLOCATE(alpha_store,beta_store)
  DEALLOCATE(gamma_store,zeta_store)
  !
  RETURN
  !
ENDSUBROUTINE
!
!
SUBROUTINE my_copy_lz(tmp_lz)
  !
  USE io_global,            ONLY : stdout
  USE mp_world,             ONLY : root,mpime,world_comm
  USE mp,                   ONLY : mp_barrier,mp_bcast
  USE wrappers,             ONLY : f_copy
  USE westcom,              ONLY : wstat_save_dir 
  !
  IMPLICIT NONE
  ! 
  ! I/O
  ! 
  CHARACTER(LEN=256), INTENT(IN) :: tmp_lz
  !
  ! Workspace
  !
  CHARACTER(LEN=320) :: cp_source, cp_dest
  INTEGER            :: cp_status
  !
  ! BARRIER
  !
  CALL mp_barrier( world_comm )
  !
  ! ... clear the directory
  !
  IF(mpime==root) THEN
    !
    cp_source = TRIM( wstat_save_dir ) // "/EVC1.dat"
    cp_dest   = TRIM( tmp_lz ) // "/EVC1.dat"
    cp_status = f_copy(cp_source, cp_dest)
    !
  ENDIF
  !
  CALL mp_bcast( cp_status, root, world_comm) 
  CALL errore( 'my_copy_lz', 'cannot copy evc1', cp_status )
  !
  IF(mpime==root) THEN
    ! 
    cp_source = TRIM( wstat_save_dir ) // "/EVC1_OLD.dat"
    cp_dest   = TRIM( tmp_lz ) // "/EVC1_OLD.dat"
    cp_status = f_copy(cp_source, cp_dest)
    !
  ENDIF
  !
  CALL mp_bcast( cp_status, root, world_comm) 
  CALL errore( 'my_copy_lz', 'cannot copy evc1_old', cp_status )
  !
  IF(mpime==root) THEN
    !
    cp_source = TRIM( wstat_save_dir ) // '/' // TRIM( 'summary.xml' ) 
    cp_dest   = TRIM( tmp_lz ) // '/' // TRIM( 'summary.xml' ) 
    cp_status = f_copy(cp_source, cp_dest)
    ! 
  ENDIF
  !
  CALL mp_bcast( cp_status, root, world_comm)
  CALL errore( 'my_copy_lz', 'cannot copy summary', cp_status )
  !
  WRITE(stdout,*) " "
  WRITE(stdout, "(\5x, 'Done tmp save in location : ',a)") TRIM( tmp_lz )
  WRITE(stdout,*) " "
  !
ENDSUBROUTINE
