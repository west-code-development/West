!
! Copyright (C) 2015-2016 M. Govoni 
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
SUBROUTINE calc_vxc( sigma_vxcl, sigma_vxcnl )
  !-----------------------------------------------------------------------
  !
  ! store in sigma_vxc(n,iks) = < n,iks | V_xc  | n,iks >     n = qp_bandrange(1):qp_bandrange(2)
  !
  USE kinds,                ONLY : DP 
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm,my_image_id,nimage
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE gvect,                ONLY : g,nl,gstart,ngm_g,ig_l2g,ngm
  USE gvecw,                ONLY : gcutw
  USE cell_base,            ONLY : tpiba2
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE constants,            ONLY : tpi,fpi,rytoev
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,nspin,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE wavefunctions_module, ONLY : evc,psic
  USE westcom,              ONLY : qp_bandrange,iuwfc,lrwfc
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol 
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : vkb,nkb
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE funct,                ONLY : dft_is_hybrid,get_exx_fraction
  USE class_idistribute,    ONLY : idistribute
  USE exx,                  ONLY : vexx,exxalfa
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP),INTENT(OUT) :: sigma_vxcl( qp_bandrange(1):qp_bandrange(2), nks )
  REAL(DP),INTENT(OUT) :: sigma_vxcnl( qp_bandrange(1):qp_bandrange(2), nks )
  !
  ! Workspace
  !
  REAL(DP) :: etxc_
  REAL(DP) :: vtxc_
  REAL(DP), ALLOCATABLE :: vxc(:,:)
  INTEGER :: ib,iv,i1,ir,iks
  COMPLEX(DP) :: braket
  REAL(DP) :: nnr 
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP), ALLOCATABLE :: xpsi(:,:),vxpsi(:,:)
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  INTEGER :: numbandegw
  TYPE(idistribute) :: gwbnd 
  !
  ALLOCATE( vxc(dfftp%nnr,nspin) )
  !
  WRITE(stdout,'(5x,a)') ' '
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') 'Vxc'
  CALL io_push_bar()
  !
  numbandegw = qp_bandrange(2)-qp_bandrange(1)+1
  gwbnd = idistribute()
  CALL gwbnd%init(numbandegw,'i','numbandegw',.FALSE.) 
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc_, vtxc_, vxc ) 
  !
  sigma_vxcl = 0._DP
  sigma_vxcnl = 0._DP
  !
  nnr = REAL( dfftp%nr1*dfftp%nr2*dfftp%nr3, KIND=DP )
  !
  barra_load = nks 
  CALL start_bar_type( barra, 'sigmavxc', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
     !npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks>1) THEN
        !iuwfc = 20
        !lrwfc = nbnd * npwx * npol 
        !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        !CALL mp_bcast(evc,0,inter_image_comm)
        !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
!     !
!     ! ... Needed for LDA+U
!     !
!     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
!          CALL get_buffer ( wfcU, nwordwfcU, iunhub, iks )
!     !
!     current_k = iks
!     current_spin = isk(iks)
!     !
!     CALL gk_sort(xk(1,iks),ngm,g,gcutw,npw,igk,g2kin)
!     g2kin=g2kin*tpiba2
!     !
!     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
!     !
!     !
!     CALL init_us_2 (npw, igk, xk (1, iks), vkb)
     !
     !nbndval = nbnd_occ(iks)
     !
     ! NON-HYBRID CONTRIBUTION TO VXC
     !
     IF( gwbnd%nloc>0 ) THEN
        !
        IF(gamma_only) THEN 
           !
           DO ib = 1, gwbnd%nloc
              CALL single_invfft_gamma(dfftp,npw,npwx,evc(1,qp_bandrange(1)+gwbnd%l2g(ib)-1),psic,'Wave')
              braket = 0._DP
              DO ir = 1, dfftp%nnr
                 braket = braket + psic(ir)*DCONJG(psic(ir)) * vxc(ir,current_spin) 
              ENDDO
              sigma_vxcl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) = REAL(braket,KIND=DP) / nnr
           ENDDO
           !
        ELSE
           !
           DO ib = 1, gwbnd%nloc
              CALL single_invfft_k(dffts,npw,npwx,evc(1,qp_bandrange(1)+gwbnd%l2g(ib)-1),psic,'Wave',igk_k(1,current_k))
              braket = 0._DP
              DO ir = 1, dfftp%nnr
                 braket = braket + psic(ir)*DCONJG(psic(ir)) * vxc(ir,current_spin) 
              ENDDO
              sigma_vxcl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) = REAL(braket,KIND=DP) / nnr
           ENDDO
           !
           IF(noncolin) THEN
              !
              DO ib = 1,gwbnd%nloc
                 CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx,qp_bandrange(1)+gwbnd%l2g(ib)-1),psic,'Wave',igk_k(1,current_k))
                 braket = 0._DP
                 DO ir = 1, dfftp%nnr
                    braket = braket + psic(ir)*DCONJG(psic(ir)) * vxc(ir,current_spin) 
                 ENDDO
                 sigma_vxcl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) = &
                &sigma_vxcl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) + REAL(braket,KIND=DP) / nnr
              ENDDO
              !
           ENDIF
           !
        ENDIF
        !
     ENDIF
     !
     ! HYBRID CONTRIBUTION TO VXC
     !
     IF(  dft_is_hybrid() ) THEN
        !
        IF( gwbnd%nloc>0 ) THEN
           !
           ALLOCATE( vxpsi(npwx*npol,gwbnd%nloc ) )
           ALLOCATE(  xpsi(npwx*npol,gwbnd%nloc ) )
           !
           xpsi = 0._DP
           DO ib=1,gwbnd%nloc
              xpsi(:,ib) = evc(:,qp_bandrange(1)+gwbnd%l2g(ib)-1)
           ENDDO
           vxpsi = 0._DP
           CALL vexx( npwx, npw, gwbnd%nloc, xpsi, vxpsi )
           !
           IF( gamma_only ) THEN
              !
              DO ib = 1,gwbnd%nloc
                 braket = 2._DP * DDOT( 2*npw, xpsi(1,ib), 1, vxpsi(1,ib), 1)
                 IF(gstart==2) braket = braket - REAL( xpsi(1,ib), KIND=DP) * REAL( vxpsi(1,ib), KIND=DP)
                 sigma_vxcnl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) = REAL( braket, KIND=DP )
              ENDDO
              !
           ELSE
              !
              DO ib = 1,gwbnd%nloc
                 braket = ZDOTC( npw, xpsi(1,ib),1,vxpsi(1,ib),1)
                 sigma_vxcnl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) = REAL( braket, KIND=DP )
              ENDDO
              !
              IF(noncolin) THEN
                 !
                 DO ib = 1, gwbnd%nloc
                    braket = ZDOTC( npw, xpsi(1+npwx,ib),1,vxpsi(1+npwx,ib),1)
                    sigma_vxcnl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) = &
                   &sigma_vxcnl(qp_bandrange(1)+gwbnd%l2g(ib)-1,iks) + REAL( braket, KIND=DP )
                 ENDDO
              ENDIF
              !
           ENDIF
           ! 
           DEALLOCATE( vxpsi )
           DEALLOCATE(  xpsi )
           !
        ENDIF
        !
     ENDIF
     !
     CALL update_bar_type( barra, 'sigmavxc', 1 )
     !
  ENDDO 
  !
  CALL stop_bar_type( barra, 'sigmavxc' )
  !
  CALL mp_sum( sigma_vxcl, intra_bgrp_comm )
  CALL mp_sum( sigma_vxcnl, intra_bgrp_comm )
  CALL mp_sum( sigma_vxcl, inter_image_comm )
  CALL mp_sum( sigma_vxcnl, inter_image_comm )
  !
  DEALLOCATE( vxc )
  !
END SUBROUTINE 
