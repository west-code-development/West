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
SUBROUTINE calc_exx2( sigma_exx, nb1, nb2 )
  !-----------------------------------------------------------------------
  !
  ! store in sigma_exx(n,iks) = < n,iks | V_exx | n,iks >     n = nb1, nb2
  !
  USE kinds,                ONLY : DP 
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm,my_image_id
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE gvect,                ONLY : g,nl,gstart,ngm_g,ig_l2g,ngm
  USE gvecs,                ONLY : ngms
  USE gvecw,                ONLY : gcutw
  USE cell_base,            ONLY : tpiba2,omega,tpiba,at,alat
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE constants,            ONLY : tpi,fpi,rytoev,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,nspin,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE wavefunctions_module, ONLY : evc,psic,psic_nc
  USE westcom,              ONLY : iuwfc,lrwfc,npwq0,nbnd_occ,div_kind_hf
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol 
  USE buffers,              ONLY : get_buffer
  USE funct,                ONLY : dft_is_hybrid,init_dft_exxrpa,stop_exx
  USE exx,                  ONLY : x_gamma_extrapolation,exxdiv_treatment,exx_grid_init,exx_div_check,&
                                   &deallocate_exx,exxinit,vexx,exx_grid_initialized
  USE uspp,                 ONLY : vkb,nkb
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE class_idistribute,    ONLY : idistribute
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get, vcut_destroy
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nb1, nb2
  REAL(DP),INTENT(OUT) :: sigma_exx( nb1:nb2, nks) 
  !
  ! Workspace
  !
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  REAL(DP), EXTERNAL :: DDOT
  INTEGER :: ib,iv,i1,ir,iks,ig,iv_glob
  INTEGER :: nbndval
  TYPE(idistribute) :: vband
  REAL(DP) :: peso
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  INTEGER :: nq1, nq2, nq3         ! integers defining the X integration mesh
  REAL(DP),ALLOCATABLE :: mysqvc(:)
  REAL(DP) :: q(3)
  REAL(DP) :: ecutvcut
  TYPE(vcut_type)   :: vcut
  REAL(DP) :: mydiv
  !
  WRITE(stdout,'(5x,a)') ' '
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') '(X)-Sigma'
  CALL io_push_bar()
  !
  ALLOCATE( pertg( ngms ) )
  ALLOCATE( mysqvc(ngms) )
  IF(noncolin) THEN 
     ALLOCATE( pertr_nc( dffts%nnr, npol ) )
  ELSE
     ALLOCATE( pertr( dffts%nnr ) )
  ENDIF
  !
  CALL store_sqvc(mysqvc,ngms,div_kind_hf,mydiv)
  !
  ! Set to zero
  !
  sigma_exx = 0._DP
  !
  barra_load = nks * ( nb2-nb1 + 1 )
  CALL start_bar_type( barra, 'sigmax', barra_load )
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
     nbndval = nbnd_occ(iks)
     !
     IF( gamma_only ) THEN 
        peso = 2._DP  
     ELSE
        peso = 1._DP
     ENDIF
     !
     vband = idistribute()
     CALL vband%init(nbndval,'i','nbndval',.FALSE.)
     !
     DO ib = nb1, nb2
        !
        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ib),psic,'Wave') 
        ELSEIF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1     ,ib),psic_nc(1,1),'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx,ib),psic_nc(1,2),'Wave',igk_k(1,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic,'Wave',igk_k(1,current_k))
        ENDIF
        !
        DO iv = 1, vband%nloc
           !
           iv_glob = vband%l2g(iv)
           !
           ! Bring it to R-space
           IF(gamma_only) THEN
              CALL single_invfft_gamma(dffts,npw,npwx,evc(1,iv_glob),pertr,'Wave')
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic(ir)*pertr(ir)
              ENDDO
              CALL single_fwfft_gamma(dffts,ngms,ngms,pertr,pertg,'Smooth')
           ELSEIF(noncolin) THEN
              CALL single_invfft_k(dffts,npw,npwx,evc(1     ,iv_glob),pertr_nc(1,1),'Wave',igk_k(1,current_k))
              CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx,iv_glob),pertr_nc(1,2),'Wave',igk_k(1,current_k))
              DO ir=1,dffts%nnr 
                 pertr_nc(ir,1)=DCONJG(psic_nc(ir,1))*pertr_nc(ir,1)+DCONJG(psic_nc(ir,2))*pertr_nc(ir,2)
              ENDDO
              CALL single_fwfft_k(dffts,ngms,ngms,pertr_nc(1,1),pertg,'Smooth') ! no igk
           ELSE
              CALL single_invfft_k(dffts,npw,npwx,evc(1,iv_glob),pertr,'Wave',igk_k(1,current_k))
              DO ir=1,dffts%nnr 
                 pertr(ir)=DCONJG(psic(ir))*pertr(ir)
              ENDDO
              CALL single_fwfft_k(dffts,ngms,ngms,pertr,pertg,'Smooth') ! no igk
           ENDIF 
           !
           DO ig = 1,ngms
              pertg(ig) = pertg(ig) * mysqvc(ig) 
           ENDDO
           sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - peso * DDOT( 2*ngms, pertg(1), 1, pertg(1), 1) / omega
           !IF(gstart==2) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) + REAL( pertg(1), KIND = DP )**2 / omega
           IF( ib == iv_glob .AND. gstart == 2 ) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - mydiv
           !
        ENDDO
        !
        CALL update_bar_type( barra, 'sigmax', 1 )
        !
     ENDDO ! ib
     !
  ENDDO ! iks
  !
  CALL stop_bar_type( barra, 'sigmax' )
  !
  CALL mp_sum( sigma_exx, intra_bgrp_comm )
  CALL mp_sum( sigma_exx, inter_image_comm ) 
  !
  DEALLOCATE( pertg ) 
  DEALLOCATE( mysqvc )
  IF( noncolin ) THEN 
    DEALLOCATE( pertr_nc ) 
  ELSE
    DEALLOCATE( pertr ) 
  ENDIF
  !
END SUBROUTINE 
