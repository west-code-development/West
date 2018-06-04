!
! Copyright (C) 2015-2017 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file: 
! Marco Govoni, Nicholas Brawand 
!
!----------------------------------------------------------------------------
SUBROUTINE do_sxx ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE uspp,                  ONLY : vkb,nkb
  USE io_global,             ONLY : stdout
  USE pwcom,                 ONLY : current_spin,wk,nks,nelup,neldw,isk,igk_k,xk,npw,npwx,lsda,nkstot,current_k,ngk,et
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE westcom,               ONLY : iuwfc,lrwfc,westpp_range,westpp_save_dir,nbnd_occ,iks_l2g,westpp_epsinfty,dvg,ev,&
                                  & npwq,npwqx,fftdriver,logfile,ngq,westpp_n_pdep_eigen_to_use
  USE mp_global,             ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_base,              ONLY : dfftp,dffts
  USE wvfct,                 ONLY : nbnd
  USE buffers,               ONLY : get_buffer
  USE wavefunctions_module,  ONLY : evc,psic,psic_nc
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k,single_fwfft_k
  USE distribution_center,   ONLY : pert
  USE control_flags,         ONLY : gamma_only 
  USE gvect,                 ONLY : g,nl,gstart,ngm_g,ngm
  USE cell_base,             ONLY : omega,at,alat
  USE noncollin_module,      ONLY : noncolin,npol 
  USE mp_world,              ONLY : mpime,root
  USE constants,             ONLY : rytoev
  USE json_module,           ONLY : json_file
  USE pdep_db,               ONLY : pdep_db_read
  USE types_bz_grid,         ONLY : k_grid, q_grid, compute_phase
  USE types_coulomb,         ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ir, ip, ig, iks, ib, iv, ip_glob, ik, is, ikqs, ikq, iq, nbndval, npwkq
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  COMPLEX(DP), ALLOCATABLE :: evckmq(:,:), phase(:)
  LOGICAL :: l_gammaq
  REAL(DP) :: g0(3)
  TYPE(bar_type) :: barra
  REAL(DP),ALLOCATABLE :: sigma_exx( :, : ) 
  REAL(DP),ALLOCATABLE :: sigma_sxx( :, : ) 
  REAL(DP) :: peso
  REAL(DP), EXTERNAL :: DDOT
  CHARACTER(LEN=6) :: myglobk
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  COMPLEX(DP),ALLOCATABLE :: zproj(:,:)
  REAL(DP),ALLOCATABLE :: dproj(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  LOGICAL :: l_print_pdep_read
  !
  CALL io_push_title("(S)creened eXact eXchange")
  !
  ALLOCATE( sigma_exx( westpp_range(1):westpp_range(2), k_grid%nps) )
  ALLOCATE( sigma_sxx( westpp_range(1):westpp_range(2), k_grid%nps) )
  !
  sigma_exx = 0._DP
  sigma_sxx = 0._DP
  !
  !ALLOCATE( pertg(npwqx) )
  ALLOCATE( pertg(ngm) )
  IF(noncolin) THEN 
     ALLOCATE( pertr_nc( dffts%nnr, npol ) )
  ELSE
     ALLOCATE( pertr( dffts%nnr ) )
  ENDIF
  !
  IF( gamma_only ) THEN 
     peso = 2._DP  
     ALLOCATE( dproj( 1, pert%nloc ) )
  ELSE
     peso = 1._DP
     ALLOCATE( zproj( 1, pert%nloc ) )
     ALLOCATE( phase(dffts%nnr) )
     ALLOCATE( evckmq(npwx*npol,nbnd) )
  ENDIF
  !
  CALL start_bar_type( barra, 'westpp', k_grid%nps * (westpp_range(2)-westpp_range(1)+1)  ) 
  !
  DO iks = 1, k_grid%nps  ! KPOINT-SPIN LOOP
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
     !IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
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
     !
     !nbndval = nbnd_occ(iks)
     !
     DO ib = 1, nbnd
        !
        IF( ib < westpp_range(1) .OR. ib > westpp_range(2) ) CYCLE 
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
        DO iq = 1, q_grid%np
           !
           IF ( gamma_only ) THEN
              !
              l_gammaq = .TRUE.
              nbndval = nbnd_occ(iks)
              CALL pot3D%init('Dense',.FALSE.,'gb')
              !
           ELSE
              !
              l_gammaq = q_grid%l_pIsGamma(iq)
              !
              !CALL k_grid%find( k_grid%p_cart(:,ik) - q_grid%p_cart(:,iq), is, 'cart', ikqs, g0 )  !MATTEO
              CALL k_grid%find( k_grid%p_cart(:,ik) - q_grid%p_cart(:,iq), 'cart', ikq, g0 )        !MARCO
              ikqs = k_grid%ipis2ips(ikq,is)                                                        !MARCO
              CALL compute_phase( g0, 'cart', phase )
              !
              nbndval = nbnd_occ(ikqs)
              npwq = ngq(iq)
              npwkq = ngk(ikqs)
              !
              CALL pot3D%init('Dense',.FALSE.,'gb',iq)
              !
              IF ( my_image_id == 0 ) CALL get_buffer( evckmq, lrwfc, iuwfc, ikqs )
              CALL mp_bcast( evckmq, 0, inter_image_comm )
              !
              IF (iks==1 .AND. ib==1 .AND. iq==1) THEN 
                 l_print_pdep_read = .TRUE.
              ELSE
                 l_print_pdep_read = .FALSE.
              ENDIF
              !
           ENDIF
           !
           DO iv = 1, nbndval
              !
              ! Bring it to R-space
              IF(gamma_only) THEN
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(1,iv),pertr,'Wave')
                 DO ir=1,dffts%nnr 
                    pertr(ir)=psic(ir)*pertr(ir)
                 ENDDO
                 CALL single_fwfft_gamma(dffts,npwq,npwqx,pertr,pertg,TRIM(fftdriver))
              ELSEIF(noncolin) THEN
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1     ,iv),pertr_nc(1,1),'Wave',igk_k(1,ikqs))
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1+npwx,iv),pertr_nc(1,2),'Wave',igk_k(1,ikqs))
                 DO ir=1,dffts%nnr 
                    pertr_nc(ir,1)=DCONJG(pertr_nc(ir,1)*phase(ir))*psic_nc(ir,1)+DCONJG(pertr_nc(ir,2)*phase(ir))*psic_nc(ir,2)
                 ENDDO
                 !CALL single_fwfft_k(dffts,npwq,npwqx,pertr_nc(1,1),pertg,'Dense') ! no igk
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr_nc(1,1),pertg,'Dense') ! no igk
              ELSE
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1,iv),pertr,'Wave',igk_k(1,ikqs))
                 DO ir=1,dffts%nnr 
                    pertr(ir)=DCONJG(pertr(ir)*phase(ir)) * psic(ir)
                 ENDDO
                 !CALL single_fwfft_k(dffts,npwq,npwqx,pertr,pertg,'Dense') ! no igk
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr,pertg,'Dense') ! no igk
              ENDIF 
              !
              !DO ig = 1,npwq
              DO ig = 1,ngm
                 pertg(ig) = pertg(ig) * pot3D%sqvc(ig) 
              ENDDO
              !sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - peso * DDOT( 2*npwq, pertg(1), 1, pertg(1), 1) / omega
              sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - peso*DDOT(2*ngm, pertg(1), 1, pertg(1), 1)/omega * q_grid%weight(iq)
              !IF(gstart==2) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) + REAL( pertg(1), KIND = DP )**2 / omega
              IF( ib == iv .AND. gstart == 2 .AND. l_gammaq ) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - pot3D%div
              !
              ! -- < SXX >
              !
              IF (gamma_only) THEN
                 CALL pdep_db_read(westpp_n_pdep_eigen_to_use)
              ELSE
                 CALL pdep_db_read(westpp_n_pdep_eigen_to_use,iq,l_print_pdep_read)
              ENDIF
              !
              IF( gamma_only ) THEN  
                 CALL glbrak_gamma( pertg, dvg, dproj, npwq, npwqx, 1, pert%nloc, 1, npol)
                 CALL mp_sum( dproj, intra_bgrp_comm )
                 DO ip = 1, pert%nloc
                    ip_glob = pert%l2g(ip)
                    sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - dproj(1,ip)**2 * (ev(ip_glob)/(1._DP-ev(ip_glob))) / omega
                 ENDDO
                 IF( ib == iv ) sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - (1._DP/westpp_epsinfty-1._DP) * pot3D%div
              ELSE
                 CALL glbrak_k( pertg, dvg, zproj, npwq, npwqx, 1, pert%nloc, 1, npol)
                 CALL mp_sum( zproj, intra_bgrp_comm )
                 DO ip = 1, pert%nloc
                    ip_glob = pert%l2g(ip)
                    sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - REAL(zproj(1,ip)*CONJG(zproj(1,ip)),KIND=DP) &
                                         & * (ev(ip_glob)/(1._DP-ev(ip_glob))) / omega * q_grid%weight(iq)
                 ENDDO
                 !
                 IF( ib == iv ) sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - (1._DP/westpp_epsinfty-1._DP) * pot3D%div
              ENDIF 
              ! -- </ SXX > 
              !
           ENDDO ! iv
           !
        ENDDO ! iq
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO ! ib
     !
  ENDDO ! iks
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  CALL mp_sum( sigma_exx, intra_bgrp_comm )
  CALL mp_sum( sigma_sxx, inter_image_comm )
  !
  sigma_sxx = sigma_exx + sigma_sxx 
  !
  DEALLOCATE( pertg ) 
  IF( noncolin ) THEN 
    DEALLOCATE( pertr_nc ) 
  ELSE
    DEALLOCATE( pertr ) 
  ENDIF
  IF( gamma_only ) THEN
     DEALLOCATE( dproj )
  ELSE
     DEALLOCATE( zproj )
     DEALLOCATE( evckmq )
     DEALLOCATE( phase )
  ENDIF  
  !
  ! Output it per k-point
  !
  IF(mpime==root) THEN
     CALL json%initialize()
     CALL json%load_file(filename=TRIM(logfile))
  ENDIF
  !
  ! STDOUT
  !
  WRITE(stdout,"(5X)")  
  CALL io_push_bar()
  WRITE(stdout,"(5X,a,1X,a,1X,a,1X,a,1X,a,1X,a)") &
  & 'K     ', 'B     ', '      Eks [eV]', '       Sx [eV]', '      Sxx [eV]', '        Sxx/Sx'
  CALL io_push_bar()
  !
  ALLOCATE(out_tab(westpp_range(2)-westpp_range(1)+1,5))
  !
  DO iks=1,k_grid%nps
     DO ib = westpp_range(1), westpp_range(2)
        out_tab( ib - westpp_range(1) + 1, 1) = REAL( ib, KIND=DP) 
        out_tab( ib - westpp_range(1) + 1, 2) = et(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 3) = sigma_exx(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 4) = sigma_sxx(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 5) = sigma_sxx(ib,iks) / sigma_exx(ib,iks)
        WRITE(stdout,"(5X,i5.5,1X,i6.6,1X,1f14.6,1X,1f14.6,1X,1f14.6,1X,1f14.6)") &
        & iks, ib, et(ib,iks)*rytoev, sigma_exx(ib,iks)*rytoev, sigma_sxx(ib,iks)*rytoev, sigma_sxx(ib,iks)/sigma_exx(ib,iks)
     ENDDO
     IF (k_grid%nps>1.AND.iks<k_grid%nps) CALL io_push_bar()
     WRITE(myglobk,'(i5.5)') iks_l2g(iks)
     !
     !CALL serial_table_output(mpime==root,4000,'sxx_K'//myglobk,out_tab,&
     !& westpp_range(2)-westpp_range(1)+1,5,&
     !& (/'      band','    E0[eV]','    Sx[eV]','   Sxx[eV]','    Sxx/Sx'/))
     !
     IF(mpime==root) THEN 
        !
        CALL json%add('output.S.K'//TRIM(myglobk)//'.bandmap',out_tab(:,1))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.eks',out_tab(:,2))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.sx',out_tab(:,3))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.sxx',out_tab(:,4))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.fraction',out_tab(:,5))
        !
     ENDIF
  ENDDO
  DEALLOCATE(out_tab)
  CALL io_push_bar()
  !
  IF( mpime == root ) THEN 
     !
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print_file( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
  DEALLOCATE( sigma_exx )
  DEALLOCATE( sigma_sxx ) 
  !
END SUBROUTINE
