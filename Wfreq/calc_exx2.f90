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
  USE gvect,                ONLY : g,nl,gstart,ngm_g,ngm
  USE gvecw,                ONLY : gcutw
  USE cell_base,            ONLY : tpiba2,omega,tpiba,at,alat
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE constants,            ONLY : tpi,fpi,rytoev,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk_k,g2kin,nspin,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE wavefunctions_module, ONLY : evc,psic,psic_nc
  USE westcom,              ONLY : iuwfc,lrwfc,npwq,nbnd_occ
  USE westcom,              ONLY : l_enable_off_diagonal,sigma_exx_full,ijpmap,npair
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
  USE types_bz_grid,        ONLY : k_grid, q_grid, compute_phase
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nb1, nb2
  REAL(DP),INTENT(OUT) :: sigma_exx( nb1:nb2, k_grid%nps) 
  !
  ! Workspace
  !
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  COMPLEX(DP),ALLOCATABLE :: psic1(:), pertr1(:), pertg1(:)
  COMPLEX(DP), ALLOCATABLE :: evckmq(:,:), phase(:)
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  INTEGER :: ib,iv,i1,ir,iks,ik,is,ig,iv_glob,iq,ikqs,ikq,jb,index
  INTEGER :: nbndval
  INTEGER :: npwkq
  TYPE(idistribute) :: vband
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  LOGICAL :: l_gammaq
  REAL(DP) :: g0(3), peso
  COMPLEX(DP) :: braket
  !
  WRITE(stdout,'(5x,a)') ' '
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') '(X)-Sigma'
  CALL io_push_bar()
  !
  ALLOCATE( pertg( ngm ) )
  IF (gamma_only) THEN
     peso = 2._DP
     IF (l_enable_off_diagonal) THEN
        ALLOCATE( psic1(dffts%nnr) )
        ALLOCATE( pertr1(dffts%nnr) )
        ALLOCATE( pertg1(ngm) )
     ENDIF
  ELSE
     peso = 1._DP
     ALLOCATE( phase(dffts%nnr) )
     ALLOCATE( evckmq(npwx*npol,nbnd) )
  ENDIF
  IF(noncolin) THEN 
     ALLOCATE( pertr_nc( dffts%nnr, npol ) )
  ELSE
     ALLOCATE( pertr( dffts%nnr ) )
  ENDIF
  !
  !
  ! Set to zero
  !
  sigma_exx = 0._DP
  !
  IF (l_enable_off_diagonal) THEN
      barra_load = k_grid%nps * npair
  ELSE
      barra_load = k_grid%nps * ( nb2-nb1 + 1 )
  ENDIF
  CALL start_bar_type( barra, 'sigmax', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, k_grid%nps   ! KPOINT-SPIN
     !
     ik = k_grid%ip(iks) 
     is = k_grid%is(iks)
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), k_grid%p_cart(1,ik), vkb )
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(k_grid%nps>1) THEN
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
     !
     DO ib = nb1, nb2
        !
        !sigma_exx(ib,iks) = 0._DP
        !
        IF (gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ib),psic,'Wave')
        ELSEIF (noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1     ,ib),psic_nc(1,1),'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ib),psic_nc(1,2),'Wave',igk_k(1,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic,'Wave',igk_k(1,current_k))
        ENDIF
        !
        DO jb = nb1, nb2
           !
           IF ( l_enable_off_diagonal ) THEN
              IF ( jb > ib ) CYCLE
              index = ijpmap(jb,ib)
           ELSE
              IF ( jb /= ib ) CYCLE
           ENDIF
           !
           IF (gamma_only) THEN
              CALL single_invfft_gamma(dffts,npw,npwx,evc(1,jb),psic1,'Wave')
           ENDIF
           !
           DO iq = 1, q_grid%np
              !
              IF (gamma_only) THEN
                 l_gammaq = .TRUE.
                 CALL pot3D%init('Dense',.FALSE.,'gb')
                 nbndval = nbnd_occ(iks)
              ELSE
                 l_gammaq = q_grid%l_pIsGamma(iq)
                 CALL pot3D%init('Dense',.FALSE.,'gb',iq)
                 !
                 !CALL k_grid%find( k_grid%p_cart(:,ik) - q_grid%p_cart(:,iq), is, 'cart', ikqs, g0 )  !M
                 CALL k_grid%find( k_grid%p_cart(:,ik) - q_grid%p_cart(:,iq), 'cart', ikq, g0 )        
                 ikqs = k_grid%ipis2ips(ikq,is)                                                        
                 CALL compute_phase( g0, 'cart', phase )
                 !
                 nbndval = nbnd_occ(ikqs)
                 npwkq = ngk(ikqs)
                 IF ( my_image_id == 0 ) CALL get_buffer( evckmq, lrwfc, iuwfc, ikqs )
                 CALL mp_bcast( evckmq, 0, inter_image_comm )
              ENDIF
              !
              vband = idistribute()
              CALL vband%init(nbndval,'i','nbndval',.FALSE.)
              !
              DO iv = 1, vband%nloc
                 !
                 iv_glob = vband%l2g(iv)
                 !
                 ! Bring it to R-space
                 IF (gamma_only) THEN
                    CALL single_invfft_gamma(dffts,npw,npwx,evc(1,iv_glob),pertr,'Wave')
                    DO ir=1,dffts%nnr
                       pertr(ir)=psic(ir)*pertr(ir)
                       IF ( l_enable_off_diagonal ) pertr1(ir)=psic1(ir)*pertr(ir)
                    ENDDO
                    CALL single_fwfft_gamma(dffts,ngm,ngm,pertr,pertg,'Dense')
                    CALL single_fwfft_gamma(dffts,ngm,ngm,pertr1,pertg1,'Dense')
                 ELSEIF(noncolin) THEN
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1     ,iv_glob),pertr_nc(1,1),'Wave',igk_k(1,ikqs))
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1+npwx,iv_glob),pertr_nc(1,2),'Wave',igk_k(1,ikqs))
                    DO ir=1,dffts%nnr 
                       pertr_nc(ir,1)=DCONJG(pertr_nc(ir,1)*phase(ir))*psic_nc(ir,1)+DCONJG(pertr_nc(ir,2)*phase(ir))*psic_nc(ir,2)
                    ENDDO
                    CALL single_fwfft_k(dffts,ngm,ngm,pertr_nc(1,1),pertg,'Dense') ! no igk
                 ELSE
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1,iv_glob),pertr,'Wave',igk_k(1,ikqs))
                    DO ir=1,dffts%nnr 
                       pertr(ir)=DCONJG(pertr(ir)*phase(ir)) * psic(ir)
                    ENDDO
                    CALL single_fwfft_k(dffts,ngm,ngm,pertr,pertg,'Dense') ! no igk
                 ENDIF 
                 !
                 DO ig = 1,ngm
                    pertg(ig) = pertg(ig) * pot3D%sqvc(ig) 
                    IF ( l_enable_off_diagonal ) pertg1(ig) = pertg1(ig) * pot3D%sqvc(ig)
                 ENDDO
                 !
                 IF ( l_enable_off_diagonal .AND. jb < ib ) THEN
                    braket = peso*REAL( ZDOTC( ngm, pertg(1), 1, pertg1(1), 1 ) )/omega*q_grid%weight(iq)
                    sigma_exx_full( index, iks ) = sigma_exx_full( index, iks ) - braket
                 ELSEIF ( jb == ib ) THEN
                    braket = peso*DDOT( 2*ngm, pertg(1), 1, pertg(1), 1)/omega*q_grid%weight(iq)
                    sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - braket
                    IF ( l_enable_off_diagonal ) sigma_exx_full( index, iks ) = sigma_exx_full( index, iks ) - braket
                    !IF(gstart==2) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) + REAL( pertg(1), KIND = DP )**2 / omega
                    IF( ib == iv_glob .AND. gstart == 2 .AND. l_gammaq ) THEN
                       sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - pot3D%div
                       IF ( l_enable_off_diagonal ) sigma_exx_full( index, iks ) = sigma_exx_full( index, iks ) - pot3D%div
                    ENDIF 
                 ENDIF
                 !
              ENDDO ! iv
              !
            ENDDO ! iq
            !
            CALL update_bar_type( barra, 'sigmax', 1 )
            !
        ENDDO ! jb
        !
     ENDDO ! ib
     !
  ENDDO ! iks
  !
  CALL stop_bar_type( barra, 'sigmax' )
  !
  CALL mp_sum( sigma_exx, intra_bgrp_comm )
  CALL mp_sum( sigma_exx, inter_image_comm ) 
  IF (l_enable_off_diagonal) THEN
     CALL mp_sum( sigma_exx_full, intra_bgrp_comm)
     CALL mp_sum( sigma_exx_full, inter_image_comm)
  ENDIF
  !
  DEALLOCATE( pertg ) 
  IF (l_enable_off_diagonal) DEALLOCATE( pertg1 )
  IF( noncolin ) THEN 
    DEALLOCATE( pertr_nc ) 
  ELSE
    DEALLOCATE( pertr ) 
    IF (l_enable_off_diagonal) DEALLOCATE( pertr1, psic1 )
  ENDIF
  IF (.NOT.gamma_only) THEN
     DEALLOCATE( phase )
     DEALLOCATE( evckmq )
  ENDIF
  !
END SUBROUTINE 
