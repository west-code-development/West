
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
SUBROUTINE calc_exx2(sigma_exx)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_exx(n,iks) = < n,iks | V_exx | n,iks >     n = qp_bands(1):qp_bands(SIZE(qp_bands))
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout
  USE gvect,                ONLY : gstart,ngm
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : npw,npwx,nbnd,igk_k,current_k,ngk
  USE fft_at_gamma,         ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,             ONLY : single_invfft_k,single_fwfft_k
  USE wavefunctions,        ONLY : evc,psic,psic_nc
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,occupation,qp_bands
  USE westcom,              ONLY : l_enable_off_diagonal,sigma_exx_full,ijpmap,npair
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE distribution_center,  ONLY : band_group,kpt_pool
  USE class_idistribute,    ONLY : idistribute
  USE types_bz_grid,        ONLY : k_grid,q_grid,compute_phase
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(OUT) :: sigma_exx(SIZE(qp_bands),k_grid%nps)
  !
  ! Workspace
  !
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  COMPLEX(DP),ALLOCATABLE :: psic1(:), pertr1(:), pertg1(:)
  COMPLEX(DP), ALLOCATABLE :: evckmq(:,:), phase(:)
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  INTEGER :: ib,iv,ir,iks,ik,is,ig,ivloc,ibloc,iq,ikqs,ikq,iks_g,jb,ib_index,jb_index,index
  INTEGER :: nbndval
  INTEGER :: npwkq
  TYPE(idistribute) :: vband
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  LOGICAL :: l_gammaq
  REAL(DP) :: g0(3), peso
  COMPLEX(DP) :: braket
  !
  WRITE(stdout,'(5x,a)') ''
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') '(X)-Sigma'
  CALL io_push_bar()
  !
  ALLOCATE(pertg(ngm))
  IF (gamma_only) THEN
     peso = 2._DP
     IF (l_enable_off_diagonal) THEN
        ALLOCATE( psic1(dffts%nnr) )
        ALLOCATE( pertr1(dffts%nnr) )
        ALLOCATE( pertg1(ngm) )
     ENDIF
  ELSE
     peso = 1._DP
     ALLOCATE(phase(dffts%nnr))
     ALLOCATE(evckmq(npwx*npol,nbnd))
  ENDIF
  IF(noncolin) THEN
     ALLOCATE(pertr_nc(dffts%nnr,npol))
  ELSE
     ALLOCATE(pertr(dffts%nnr))
  ENDIF
  !
  !
  ! Set to zero
  !
  sigma_exx = 0._DP
  !
  CALL band_group%init(SIZE(qp_bands),'b','band_group',.FALSE.)
  !
  IF (l_enable_off_diagonal) THEN
      barra_load = kpt_pool%nloc*band_group%nloc*band_group%nglob
  ELSE
      barra_load = kpt_pool%nloc*band_group%nloc
  ENDIF
  CALL start_bar_type(barra,'sigmax',barra_load)
  !
  ! LOOP
  !
  DO iks = 1,kpt_pool%nloc ! KPOINT-SPIN
     !
     iks_g = kpt_pool%l2g(iks)
     ik = k_grid%ip(iks_g)
     is = k_grid%is(iks_g)
     current_k = iks
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     DO ibloc = 1,band_group%nloc
        !
        ib_index = band_group%l2g(ibloc)
        ib = qp_bands(ib_index)
        !
        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ib),psic,'Wave')
        ELSEIF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic_nc(1,1),'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1,ib),psic_nc(1,2),'Wave',igk_k(1,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic,'Wave',igk_k(1,current_k))
        ENDIF
        !
        DO jb_index = 1, qp_bands
           !
           jb = qp_bands(jb_index)
           !
           IF ( l_enable_off_diagonal ) THEN
              IF (jb > ib) CYCLE
              index = ijpmap(jb_index,ib_index)
              !
              IF (jb < ib) THEN 
                  IF (gamma_only) THEN
                     CALL single_invfft_gamma(dffts,npw,npwx,evc(1,jb),psic1,'Wave')
                  ENDIF
               ENDIF
              !
           ELSE 
              IF (jb /= ib) CYCLE
           ENDIF
           !
           DO iq = 1, q_grid%np
              !
              IF (gamma_only) THEN
                 l_gammaq = .TRUE.
                 CALL pot3D%init('Rho',.FALSE.,'gb')
                 nbndval = nbnd_occ(iks)
              ELSE
                 l_gammaq = q_grid%l_pIsGamma(iq)
                 CALL pot3D%init('Rho',.FALSE.,'gb',iq)
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
              DO ivloc = 1,vband%nloc
                 !
                 iv = vband%l2g(ivloc)
                 !
                 ! Bring it to R-space
                 IF (gamma_only) THEN
                    CALL single_invfft_gamma(dffts,npw,npwx,evc(1,iv),pertr,'Wave')
                    DO ir=1,dffts%nnr
                       IF ( l_enable_off_diagonal .AND. jb < ib ) pertr1(ir)=psic1(ir)*pertr(ir)
                       pertr(ir)=psic(ir)*pertr(ir)
                    ENDDO
                    CALL single_fwfft_gamma(dffts,ngm,ngm,pertr,pertg,'Rho')
                    IF ( l_enable_off_diagonal .AND. jb < ib ) CALL single_fwfft_gamma(dffts,ngm,ngm,pertr1,pertg1,'Rho')
                 ELSEIF(noncolin) THEN
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1     ,iv),pertr_nc(1,1),'Wave',igk_k(1,ikqs))
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1+npwx,iv),pertr_nc(1,2),'Wave',igk_k(1,ikqs))
                    DO ir=1,dffts%nnr 
                       pertr_nc(ir,1)=DCONJG(pertr_nc(ir,1)*phase(ir))*psic_nc(ir,1)+DCONJG(pertr_nc(ir,2)*phase(ir))*psic_nc(ir,2)
                    ENDDO
                    CALL single_fwfft_k(dffts,ngm,ngm,pertr_nc(1,1),pertg,'Rho') ! no igk
                 ELSE
                    CALL single_invfft_k(dffts,npwkq,npwx,evckmq(1,iv),pertr,'Wave',igk_k(1,ikqs))
                    DO ir=1,dffts%nnr 
                       pertr(ir)=DCONJG(pertr(ir)*phase(ir)) * psic(ir)
                    ENDDO
                    CALL single_fwfft_k(dffts,ngm,ngm,pertr,pertg,'Rho') ! no igk
                 ENDIF 
                 !
                 DO ig = 1,ngm
                    pertg(ig) = pertg(ig) * pot3D%sqvc(ig) 
                    IF ( l_enable_off_diagonal .AND. jb < ib ) pertg1(ig) = pertg1(ig) * pot3D%sqvc(ig)
                 ENDDO
                 !
                 IF ( l_enable_off_diagonal .AND. jb < ib ) THEN
                    braket = peso*REAL( ZDOTC( ngm, pertg(1), 1, pertg1(1), 1 ) )/omega*q_grid%weight(iq)
                    sigma_exx_full( index, iks_g ) = sigma_exx_full( index, iks_g ) - occupation(iv,iks)*braket
                 ELSEIF ( jb == ib ) THEN
                    braket = peso*DDOT( 2*ngm, pertg(1), 1, pertg(1), 1)/omega*q_grid%weight(iq)
                    sigma_exx( ib_index, iks_g ) = sigma_exx( ib_index, iks_g ) - occupation(iv,iks)*braket
                    IF ( l_enable_off_diagonal ) sigma_exx_full( index, iks_g ) = sigma_exx_full( index, iks_g ) &
                    & - occupation(iv,iks)*braket
                    IF( ib == iv .AND. gstart == 2 .AND. l_gammaq ) THEN
                       sigma_exx( ib_index, iks_g ) = sigma_exx( ib_index, iks_g ) - occupation(iv,iks)*pot3D%div
                       IF ( l_enable_off_diagonal ) sigma_exx_full( index, iks_g ) = sigma_exx_full( index, iks_g )&
                       & - occupation(iv,iks)*pot3D%div
                    ENDIF 
                 ENDIF
                 !
              ENDDO ! iv
              !
            ENDDO ! iq
            !
        ENDDO ! jb
        !
        CALL update_bar_type( barra, 'sigmax', nb2-nb1+1  )
        !
     ENDDO ! ibloc
     !
  ENDDO ! iks
  !
  CALL stop_bar_type(barra,'sigmax')
  !
  CALL mp_sum(sigma_exx,intra_bgrp_comm)
  CALL mp_sum(sigma_exx,inter_bgrp_comm)
  CALL mp_sum(sigma_exx,inter_pool_comm)
  CALL mp_sum(sigma_exx,inter_image_comm) 
  IF (l_enable_off_diagonal) THEN
     CALL mp_sum(sigma_exx_full,intra_bgrp_comm)
     CALL mp_sum(sigma_exx_full,inter_bgrp_comm)
     CALL mp_sum(sigma_exx_full,inter_pool_comm)
     CALL mp_sum(sigma_exx_full,inter_image_comm)
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
  IF(.NOT. gamma_only) THEN
     DEALLOCATE(phase)
     DEALLOCATE(evckmq)
  ENDIF
  !
END SUBROUTINE