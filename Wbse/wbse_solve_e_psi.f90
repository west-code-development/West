!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE solve_e_psi ()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : dp
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_barrier
  USE uspp,                 ONLY : okvan
  USE gvect,                ONLY : gstart
  USE control_flags,        ONLY : gamma_only
  USE wbsecom,              ONLY : d0psi,macropol_dfpt
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !
  INTEGER, PARAMETER :: n_ipol = 3 
  !
  IF (okvan) THEN
     !
     WRITE(stdout,'(10x,"Real space dipole + USPP is not supported",/)')
#if defined(__MPI)
     CALL mp_barrier(intra_bgrp_comm)
#endif
     STOP
     !
  ENDIF
  !
  CALL start_clock ('solve_e_psi')
  !
  ! Compute dipole in the R space. This option can be used 
  ! only for finite systems (e.g. molecules).
  !
  IF (macropol_dfpt) THEN
     ! 
     CALL compute_d0psi_dfpt()
     ! 
  ELSE
     ! 
     CALL compute_d0psi_rs()
     !
  ENDIF
  !
  IF (gstart==2 .and. gamma_only) d0psi(1,:,:,:) = &
                         & CMPLX(DBLE(d0psi(1,:,:,:)),0.0d0,dp)
  !
  CALL stop_clock ('solve_e_psi')
  !
  RETURN
  !
ENDSUBROUTINE 

SUBROUTINE compute_d0psi_rs( )
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg, alat, omega
  USE fft_base,             ONLY : dfftp,dffts
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : me_bgrp,my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_sum,mp_barrier,mp_bcast 
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE buffers,              ONLY : get_buffer
  USE wvfct,                ONLY : nbnd,npwx
  USE klist,                ONLY : nks
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc,psic
  USE noncollin_module,     ONLY : noncolin,npol
  USE pwcom,                ONLY : isk,igk_k,ngk,lsda
  USE wbsecom,              ONLY : d0psi 
  USE westcom,              ONLY : nbnd_occ,iuwfc,lrwfc
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER  :: i, j, k, ip, ir, ir_end, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  INTEGER  :: ibnd, nbndval
  INTEGER  :: npw, iks, current_k, current_spin
  INTEGER, PARAMETER :: n_ipol = 3
  REAL(DP),    ALLOCATABLE :: r(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux_r(:)
  !
  WRITE(stdout,'(/,5X,"Calculation of the dipole in real space")')
  !
  IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  ALLOCATE(aux_r(dfftp%nnr))
  ALLOCATE(r(dfftp%nnr,n_ipol))
  ! 
  r(:,:) = 0.0d0 
  !
  ! Calculate r 
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
#if defined (__MPI)
  index0 = dfftp%nr1x*dfftp%nr2x*SUM(dfftp%npp(1:me_bgrp))
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
  index0 = 0
  ir_end = dfftp%nnr
#endif
  !
  DO ir = 1, ir_end 
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     !
     DO ip = 1, n_ipol
        !
        r(ir,ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                   DBLE( j )*inv_nr2*at(ip,2) + &
                   DBLE( k )*inv_nr3*at(ip,3)
        ! 
     ENDDO
     !
  ENDDO
  !
  IF (.true.) CALL shift_d0psi(r,n_ipol)
  !
  ! Calculate the product r * psi(r)
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, number pw
     !
     current_k = iks
     ! 
     IF ( lsda ) current_spin = isk(iks)
     !
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF (nks>1) THEN
        !  
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        CALL mp_bcast(evc,0,inter_image_comm)
        ! 
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     !
     IF (gamma_only) THEN
        !
        DO ibnd=1,nbndval,2
           !
           IF (ibnd<nbndval) THEN
              !
              ! double bands @ gamma
              ! 
              ! G -> R
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,ibnd+1),psic,'Wave')
              ! 
              ! R -> G
              !
              aux_r(:) = (0.0_DP, 0.0_DP)
              aux_r(:) = psic(:)
              !
              DO ip = 1, n_ipol
                 !
                 psic(:)  = (0.0_DP, 0.0_DP)
                 psic(:) =  aux_r(:) * r(:,ip) * alat 
                 !
                 CALL double_fwfft_gamma(dffts,npw,npwx,psic,d0psi(1,ibnd,iks,ip),d0psi(1,ibnd+1,iks,ip),'Wave')
                 !
              ENDDO
              !
           ELSE
              !
              IF (nbndval == nbnd) THEN  
                 !
                 ! single band @ gamma
                 !
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),psic,'Wave')
                 ! 
                 aux_r(:) = (0.0_DP, 0.0_DP)
                 aux_r(:) = psic(:)
                 ! 
                 DO ip = 1, n_ipol
                    !
                    psic(:)  = (0.0_DP, 0.0_DP)
                    psic(:) =  CMPLX(DBLE(aux_r(:)) * r(:,ip) * alat, 0.0_DP)
                    !
                    CALL single_fwfft_gamma(dffts,npw,npwx,psic,d0psi(1,ibnd,iks,ip),'Wave')
                    !
                 ENDDO
                 !
              ENDIF
              !
              IF (nbndval < nbnd) THEN
                 !
                 ! single band @ gamma
                 !
                 CALL double_invfft_gamma(dffts,npw,npwx,evc(1,ibnd),evc(1,ibnd+1),psic,'Wave')
                 ! 
                 aux_r(:) = (0.0_DP, 0.0_DP)
                 aux_r(:) = psic(:)
                 ! 
                 DO ip = 1, n_ipol
                    !
                    psic(:)  = (0.0_DP, 0.0_DP)
                    psic(:) =  CMPLX(DBLE(aux_r(:)) * r(:,ip) * alat, 0.0_DP)
                    !
                    CALL single_fwfft_gamma(dffts,npw,npwx,psic,d0psi(1,ibnd,iks,ip),'Wave')
                    !
                 ENDDO
                 !
              ENDIF
              !
           ENDIF
           !
        ENDDO
        !
     ELSE
        !
        ! only single bands
        !
        DO ibnd=1,nbndval
           !
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ibnd),psic,'Wave',igk_k(1,current_k))
           !
           aux_r(:) = (0.0_DP, 0.0_DP)
           aux_r(:) = psic(:)
           !
           DO ip = 1, n_ipol
              !
              psic(:)  = (0.0_DP, 0.0_DP)
              psic(:) =  aux_r(:) * r(:,ip) * alat      
              !
              CALL single_fwfft_k(dffts,npw,npwx,psic,d0psi(1,ibnd,iks,ip),'Wave',igk_k(1,current_k))
              !
           ENDDO
           !
        ENDDO
        !  
        IF (npol==2) THEN
           !
           DO ibnd=1,nbndval
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ibnd),psic,'Wave',igk_k(1,current_k))
              ! 
              aux_r(:) = (0.0_DP, 0.0_DP)
              aux_r(:) = psic(:)
              !
              DO ip = 1, n_ipol
                 !
                 psic(:)  = (0.0_DP, 0.0_DP)
                 psic(:) = aux_r(:) * r(:,ip) * alat
                 !  
                 CALL single_fwfft_k(dffts,npw,npwx,psic,d0psi(npwx+1,ibnd,iks,ip),'Wave',igk_k(1,current_k))
                 !
              ENDDO
              !
           ENDDO
           !    
        ENDIF
        !
     ENDIF
     !
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,d0psi(1,1,iks,ip),(1._DP,0._DP))
        !
     ENDDO
     !
  ENDDO
  !
  IF (ALLOCATED(psic)) DEALLOCATE(psic)
  DEALLOCATE(aux_r)
  DEALLOCATE(r)
  !
  RETURN
  !
END SUBROUTINE compute_d0psi_rs

SUBROUTINE shift_d0psi( r, n_ipol )
  !
  ! Shift a position operator r to the center of the molecule 
  ! for a proper calculation of d0psi in R-space.
  !
  USE fft_base,         ONLY : dfftp
  use kinds,            only : dp
  use ions_base,        only : nat,tau
  USE io_global,        ONLY : stdout
  use cell_base,        only : alat,at
  !
  IMPLICIT NONE
  !
  integer,intent(in) :: n_ipol
  real(dp),intent(inout) :: r(dfftp%nnr,n_ipol)
  !
  ! local vars
  !
  real(dp) :: mmin(3), mmax(3), center(3), origin(3), check_cell
  integer  :: ip, iatm, ir, ip1, ip2
  !
  WRITE(stdout,'(/,5X,"Dipole is shifted to the center of cell", &
                     & " for the calculation of d0psi. ")')
  !
  check_cell = 0.d0
  !
  do ip1 = 1,3
     ! 
     do ip2 = 1,3
        ! 
        if(.not. ip1 .eq. ip2) check_cell = check_cell + at(ip1,ip2)**2
        !
    enddo
    !
  enddo
  !
  ! XG: I am not sure that this type of super cell is supported now.
  !
  if (check_cell .gt. 1.d-5) call errore('shift_d0psi', &
        & "This type of the supercell is not supported",1)
  !
  mmin(:) = 2000.d0
  mmax(:)= -2000.d0
  !
  do ip = 1, n_ipol
     !
     do iatm = 1, nat
        ! 
        mmin(ip) = min(mmin(ip), tau(ip,iatm))
        !  
        mmax(ip) = max(mmax(ip), tau(ip,iatm))
        ! 
     enddo
     !  
  enddo
  !
  center(:)= 0.5d0*(mmin(:)+mmax(:))
  !
  do ip = 1, n_ipol
     !
     origin(ip)= center(ip)-0.5d0*at(ip,ip) 
     !
  enddo
  !
  do ir = 1, dfftp%nnr
     !
     do ip = 1, n_ipol 
        !
        r(ir,ip)= r(ir,ip) - origin(ip)
        !
     enddo
     ! 
     do ip = 1, n_ipol 
        !
        if(r(ir,ip) .lt. 0) r(ir,ip)=r(ir,ip)+at(ip,ip)
        !
        if(r(ir,ip) .gt. at(ip,ip)) r(ir,ip)=r(ir,ip)-at(ip,ip)
        !
     enddo
     ! 
  enddo
  !
  return
  !
ENDSUBROUTINE shift_d0psi
!
!
!
SUBROUTINE compute_d0psi_dfpt( )
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg, alat, omega
  USE fft_base,             ONLY : dfftp,dffts
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : me_bgrp,my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_sum,mp_barrier,mp_bcast
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE fft_at_k,             ONLY : single_fwfft_k,single_invfft_k
  USE buffers,              ONLY : get_buffer
  USE wvfct,                ONLY : nbnd,npwx
  USE klist,                ONLY : nks
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc,psic
  USE noncollin_module,     ONLY : noncolin,npol
  USE uspp,                 ONLY : vkb,nkb
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,current_k,wk,ngk 
  USE wbsecom,              ONLY : d0psi, l_lanzcos
  USE westcom,              ONLY : nbnd_occ,iuwfc,lrwfc,tr2_dfpt,&
                                   n_dfpt_maxiter, l_kinetic_only                              
  USE distribution_center,  ONLY : aband
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  INTEGER  :: il1, i, j, k, iv, ip, ir, ir_end, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  INTEGER  :: ibnd, nbndval
  INTEGER  :: iks, ierr
  INTEGER, PARAMETER :: n_ipol = 3
  LOGICAL, PARAMETER :: l_skip_nl_part_of_hcomr = .False.
  REAL(DP),    ALLOCATABLE :: eprec(:), e(:)
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), phi_tmp(:,:)
  !
  IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  !
  WRITE(stdout,'(/,5X,"Calculation of the dipole using DFPT method")')
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     !
     IF ( lsda ) current_spin = isk(iks)
     !  
     CALL g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF (nks>1) THEN
        !         
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        CALL mp_bcast(evc,0,inter_image_comm)
        !
     ENDIF
     !
     nbndval = nbnd_occ(iks)
     !
     ! LOOP over band states 
     !
     d0psi(:,:,iks,:) = (0.0_DP,0.0_DP)
     ! 
     DO il1 = 1, aband%nloc
        !
        iv = aband%l2g(il1)
        !
        ALLOCATE(phi(npwx*npol,3))
        ALLOCATE(phi_tmp(npwx*npol,3))
        !  
        CALL commutator_Hx_psi (iks, 1, 1, evc(1,iv), phi_tmp(1,1), l_skip_nl_part_of_hcomr)
        CALL commutator_Hx_psi (iks, 1, 2, evc(1,iv), phi_tmp(1,2), l_skip_nl_part_of_hcomr)
        CALL commutator_Hx_psi (iks, 1, 3, evc(1,iv), phi_tmp(1,3), l_skip_nl_part_of_hcomr)
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,3,phi_tmp,(1._DP,0._DP))
        !
        ALLOCATE( eprec(3) )
        ALLOCATE( e(3) )
        !
        CALL set_eprec( 1, evc(1,iv), eprec(1))
        eprec(2) = eprec(1)
        eprec(3) = eprec(1)
        !
        e(1) = et(iv,iks)
        e(2) = et(iv,iks)
        e(3) = et(iv,iks)
        !
        CALL precondition_m_wfcts( 3, phi_tmp, phi, eprec ) 
        !
        tr2_dfpt = 1.d-12
        n_dfpt_maxiter = 250
        l_kinetic_only = .false.
        !
        CALL linsolve_sternheimer_m_wfcts (nbndval, 3, phi_tmp, phi, e, eprec, tr2_dfpt, ierr )
        !
        IF (ierr/=0) THEN 
           !
           WRITE(stdout, '(7X,"** WARNING : MACROPOL not converged, ierr = ",i8)') ierr
           !
        ENDIF
        !
        d0psi(:,iv,iks,:) = phi(:,:)
        !
        DEALLOCATE( eprec )
        DEALLOCATE( e )
        DEALLOCATE( phi_tmp )
        DEALLOCATE( phi )
        !
     ENDDO
     !
     IF (l_lanzcos) THEN
        ! 
        CALL mp_sum (d0psi(:,:,iks,:),inter_image_comm)     
        !
     ENDIF
     ! 
     ! P_c|d0psi>
     !
     DO ip = 1, n_ipol
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,d0psi(1,1,iks,ip),(1._DP,0._DP))
        !
     ENDDO
     !
  ENDDO     
  !
  IF (ALLOCATED(psic)) DEALLOCATE(psic)
  !
  RETURN
  !
END SUBROUTINE
