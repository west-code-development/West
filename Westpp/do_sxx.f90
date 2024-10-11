!
! Copyright (C) 2015-2024 M. Govoni
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
  USE io_global,             ONLY : stdout
  USE pwcom,                 ONLY : igk_k,npw,npwx,current_k,ngk,et
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE westcom,               ONLY : iuwfc,lrwfc,westpp_range,nbnd_occ,westpp_epsinfty,dvg,ev,&
                                  & npwq,npwqx,fftdriver,logfile,ngq,westpp_n_pdep_eigen_to_use
  USE mp_global,             ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_base,              ONLY : dffts
  USE wvfct,                 ONLY : nbnd
  USE buffers,               ONLY : get_buffer
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k,single_fwfft_k
  USE distribution_center,   ONLY : pert
  USE class_idistribute,     ONLY : idistribute
  USE control_flags,         ONLY : gamma_only
  USE gvect,                 ONLY : gstart,ngm
  USE cell_base,             ONLY : omega
  USE noncollin_module,      ONLY : noncolin,npol
  USE mp_world,              ONLY : mpime,root
  USE constants,             ONLY : rytoev
  USE json_module,           ONLY : json_file
  USE pdep_db,               ONLY : pdep_db_read
  USE types_bz_grid,         ONLY : k_grid,q_grid,compute_phase
  USE types_coulomb,         ONLY : pot3D
  USE wavefunctions,         ONLY : evc,psic,psic_nc
#if defined(__CUDA)
  USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ir, ip, ig, iks, ib, iv, ip_glob, ik, is, ikqs, ikq, iq, nbndval, npwkq
  INTEGER :: dffts_nnr
  COMPLEX(DP),ALLOCATABLE :: pertg(:), pertr(:), pertr_nc(:,:)
  COMPLEX(DP),ALLOCATABLE :: evckmq(:,:), phase(:)
  LOGICAL :: l_gammaq
  REAL(DP) :: g0(3)
  REAL(DP) :: dot_tmp
  TYPE(bar_type) :: barra
  REAL(DP),ALLOCATABLE :: sigma_exx(:,:)
  REAL(DP),ALLOCATABLE :: sigma_sxx(:,:)
  REAL(DP) :: peso
  CHARACTER(LEN=6) :: label_k
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  COMPLEX(DP),ALLOCATABLE :: zproj(:,:)
  REAL(DP),ALLOCATABLE :: dproj(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  LOGICAL :: l_print_pdep_read
  !
  IF(westpp_n_pdep_eigen_to_use < 1) CALL errore('do_sxx','westpp_n_pdep_eigen_to_use < 1',1)
  IF(westpp_range(2) > nbnd) CALL errore('do_sxx','westpp_range(2) > nbnd',1)
  !
  pert = idistribute()
  CALL pert%init(westpp_n_pdep_eigen_to_use,'i','npdep',.TRUE.)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( sigma_exx( westpp_range(1):westpp_range(2), k_grid%nps) )
  ALLOCATE( sigma_sxx( westpp_range(1):westpp_range(2), k_grid%nps) )
  !
  sigma_exx = 0._DP
  sigma_sxx = 0._DP
  !
  ALLOCATE( pertg(ngm) )
  !$acc enter data create(pertg)
  IF(noncolin) THEN
     ALLOCATE( pertr_nc( dffts%nnr, npol ) )
     !$acc enter data create(pertr_nc)
  ELSE
     ALLOCATE( pertr( dffts%nnr ) )
     !$acc enter data create(pertr)
  ENDIF
  !
  IF( gamma_only ) THEN
     CALL pdep_db_read(westpp_n_pdep_eigen_to_use)
     !
     peso = 2._DP
     ALLOCATE( dproj( 1, pert%nloc ) )
     !$acc enter data create(dproj) copyin(dvg)
  ELSE
     peso = 1._DP
     ALLOCATE( zproj( 1, pert%nloc ) )
     ALLOCATE( phase(dffts%nnr) )
     ALLOCATE( evckmq(npwx*npol,nbnd) )
     !$acc enter data create(zproj,phase,evckmq)
  ENDIF
  !
  CALL io_push_title('(S)creened eXact eXchange')
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
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(k_grid%nps > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     DO ib = 1, nbnd
        !
        IF( ib < westpp_range(1) .OR. ib > westpp_range(2) ) CYCLE
        !
        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ib),psic,'Wave')
        ELSEIF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1:npwx,ib),psic_nc(:,1),'Wave',igk_k(:,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx:npwx*2,ib),psic_nc(:,2),'Wave',igk_k(:,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(:,ib),psic,'Wave',igk_k(:,current_k))
        ENDIF
        !
        DO iq = 1, q_grid%np
           !
           IF(.NOT. gamma_only) THEN
              CALL pdep_db_read(westpp_n_pdep_eigen_to_use,iq,l_print_pdep_read)
              !$acc enter data copyin(dvg)
           ENDIF
           !
           IF ( gamma_only ) THEN
              !
              l_gammaq = .TRUE.
              nbndval = nbnd_occ(iks)
              CALL pot3D%init('Rho',.FALSE.,'gb')
              !
           ELSE
              !
              l_gammaq = q_grid%l_pIsGamma(iq)
              !
              CALL k_grid%find( k_grid%p_cart(:,ik) - q_grid%p_cart(:,iq), 'cart', ikq, g0 )
              ikqs = k_grid%ipis2ips(ikq,is)
              CALL compute_phase( g0, 'cart', phase )
              !
              nbndval = nbnd_occ(ikqs)
              npwq = ngq(iq)
              npwkq = ngk(ikqs)
              !
              CALL pot3D%init('Rho',.FALSE.,'gb',iq)
              !
              IF ( my_image_id == 0 ) CALL get_buffer( evckmq, lrwfc, iuwfc, ikqs )
              CALL mp_bcast( evckmq, 0, inter_image_comm )
              !
              !$acc update device(evckmq,phase)
              !
              IF (iks==1 .AND. ib==1 .AND. iq==1) THEN
                 l_print_pdep_read = .TRUE.
              ELSE
                 l_print_pdep_read = .FALSE.
              ENDIF
              !
           ENDIF
           !
           !$acc enter data copyin(pot3D)
           !$acc enter data copyin(pot3D%sqvc)
           !
           DO iv = 1, nbndval
              !
              ! Bring it to R-space
              IF(gamma_only) THEN
                 !
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(:,iv),pertr,'Wave')
                 !
                 !$acc parallel loop present(pertr)
                 DO ir = 1, dffts_nnr
                    pertr(ir)=psic(ir)*pertr(ir)
                 ENDDO
                 !$acc end parallel
                 !
                 CALL single_fwfft_gamma(dffts,npwq,npwqx,pertr,pertg,TRIM(fftdriver))
                 !
              ELSEIF(noncolin) THEN
                 !
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1:npwx,iv),pertr_nc(:,1),'Wave',igk_k(:,ikqs))
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1+npwx:npwx*2,iv),pertr_nc(:,2),'Wave',igk_k(:,ikqs))
                 !
                 !$acc parallel loop present(pertr_nc,phase)
                 DO ir = 1, dffts_nnr
                    pertr_nc(ir,1)=CONJG(pertr_nc(ir,1)*phase(ir))*psic_nc(ir,1)+CONJG(pertr_nc(ir,2)*phase(ir))*psic_nc(ir,2)
                 ENDDO
                 !$acc end parallel
                 !
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr_nc(:,1),pertg,'Rho') ! no igk
                 !
              ELSE
                 !
                 CALL single_invfft_k(dffts,npwkq,npwx,evckmq(:,iv),pertr,'Wave',igk_k(:,ikqs))
                 !
                 !$acc parallel loop present(pertr,phase)
                 DO ir = 1, dffts_nnr
                    pertr(ir)=CONJG(pertr(ir)*phase(ir)) * psic(ir)
                 ENDDO
                 !$acc end parallel
                 !
                 CALL single_fwfft_k(dffts,ngm,ngm,pertr,pertg,'Rho') ! no igk
                 !
              ENDIF
              !
              !$acc parallel loop present(pertg,pot3D,pot3D%sqvc)
              DO ig = 1, ngm
                 pertg(ig) = pertg(ig) * pot3D%sqvc(ig)
              ENDDO
              !$acc end parallel
              !
              dot_tmp = 0._DP
              !$acc parallel loop reduction(+:dot_tmp) present(pertg) copy(dot_tmp)
              DO ig = 1, ngm
                 dot_tmp = dot_tmp+REAL(pertg(ig),KIND=DP)**2+AIMAG(pertg(ig))**2
              ENDDO
              !$acc end parallel
              !
              sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - peso*dot_tmp/omega * q_grid%weight(iq)
              IF( ib == iv .AND. gstart == 2 .AND. l_gammaq ) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - pot3D%div
              !
              ! -- < SXX >
              !
              IF( gamma_only ) THEN
                 CALL glbrak_gamma( pertg, dvg, dproj, npwq, npwqx, 1, pert%nloc, 1, npol)
                 !$acc update host(dproj)
                 CALL mp_sum( dproj, intra_bgrp_comm )
                 DO ip = 1, pert%nloc
                    ip_glob = pert%l2g(ip)
                    sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - dproj(1,ip)**2 * (ev(ip_glob)/(1._DP-ev(ip_glob))) / omega
                 ENDDO
                 IF( ib == iv ) sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - (1._DP/westpp_epsinfty-1._DP) * pot3D%div
              ELSE
                 CALL glbrak_k( pertg, dvg, zproj, npwq, npwqx, 1, pert%nloc, 1, npol)
                 !$acc update host(zproj)
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
           IF(.NOT. gamma_only) THEN
              !$acc exit data delete(dvg)
           ENDIF
           !
           !$acc exit data delete(pot3D%sqvc)
           !$acc exit data delete(pot3D)
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
  !$acc exit data delete(pertg)
  DEALLOCATE( pertg )
  IF( noncolin ) THEN
     !$acc exit data delete(pertr_nc)
     DEALLOCATE( pertr_nc )
  ELSE
     !$acc exit data delete(pertr)
     DEALLOCATE( pertr )
  ENDIF
  IF( gamma_only ) THEN
     !$acc exit data delete(dproj,dvg)
     DEALLOCATE( dproj )
  ELSE
     !$acc exit data delete(zproj,phase,evckmq)
     DEALLOCATE( zproj )
     DEALLOCATE( phase )
     DEALLOCATE( evckmq )
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
  ! Output it per k-point
  !
  IF(mpime==root) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
  !
  ! STDOUT
  !
  WRITE(stdout,*)
  CALL io_push_bar()
  WRITE(stdout,'(5X,a)') 'K       B            Eks [eV]        Sx [eV]       Sxx [eV]         Sxx/Sx'
  CALL io_push_bar()
  !
  ALLOCATE(out_tab(westpp_range(2)-westpp_range(1)+1,5))
  !
  DO iks = 1, k_grid%nps
     DO ib = westpp_range(1), westpp_range(2)
        out_tab( ib - westpp_range(1) + 1, 1) = REAL( ib, KIND=DP)
        out_tab( ib - westpp_range(1) + 1, 2) = et(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 3) = sigma_exx(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 4) = sigma_sxx(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 5) = sigma_sxx(ib,iks) / sigma_exx(ib,iks)
        WRITE(stdout,'(5X,i6.6,2X,i6.6,X,f14.6,X,f14.6,X,f14.6,X,f14.6)') &
        & iks, ib, et(ib,iks)*rytoev, sigma_exx(ib,iks)*rytoev, sigma_sxx(ib,iks)*rytoev, sigma_sxx(ib,iks)/sigma_exx(ib,iks)
     ENDDO
     IF (k_grid%nps>1.AND.iks<k_grid%nps) CALL io_push_bar()
     WRITE(label_k,'(i6.6)') iks
     !
     IF(mpime==root) THEN
        !
        CALL json%add('output.S.K'//label_k//'.bandmap',out_tab(:,1))
        CALL json%add('output.S.K'//label_k//'.eks',out_tab(:,2))
        CALL json%add('output.S.K'//label_k//'.sx',out_tab(:,3))
        CALL json%add('output.S.K'//label_k//'.sxx',out_tab(:,4))
        CALL json%add('output.S.K'//label_k//'.fraction',out_tab(:,5))
        !
     ENDIF
  ENDDO
  DEALLOCATE(out_tab)
  CALL io_push_bar()
  !
  IF( mpime == root ) THEN
     !
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print( iunit )
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
