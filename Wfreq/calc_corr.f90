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
SUBROUTINE calc_corr_gamma( sigma_corr, energy, l_verbose, l_full)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_corr(n,iks) = < ib,iks | S_c(energy(ib,iks))  | ib,iks >
  ! ... ib = qp_bandrange(1):qp_bandrange(2)
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_image_comm,inter_pool_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE constants,            ONLY : pi
  USE pwcom,                ONLY : et
  USE westcom,              ONLY : qp_bandrange,l_enable_lanczos,n_lanczos,l_macropol,&
                                 & d_head_ifr,z_head_rfr,d_body1_ifr,d_body2_ifr,d_diago,&
                                 & z_body_rfr,l_enable_off_diagonal,ijpmap,npair, d_body1_ifr_full,&
                                 & d_body2_ifr_full,d_diago_full,z_body_rfr_full,sigma_corr_full,&
                                 & l_frac_occ,occupation,nbnd_occ,nbnd_occ_full
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar,io_push_title
  USE distribution_center,  ONLY : pert,ifr,rfr,aband,kpt_pool
  USE types_coulomb,        ONLY : pot3D
  USE types_bz_grid,        ONLY : k_grid
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP),INTENT(OUT) :: sigma_corr( qp_bandrange(1):qp_bandrange(2), k_grid%nps )  ! The correlation self-energy, imaginary part is lifetime.
  REAL(DP),INTENT(IN) :: energy( qp_bandrange(1):qp_bandrange(2), k_grid%nps )          ! The energy variable
  LOGICAL,INTENT(IN) :: l_verbose
  LOGICAL,INTENT(IN) :: l_full
  !
  ! Workspace
  !
  INTEGER :: iks,ib,ifreq,glob_ifreq,il,im,glob_im,ip,iks_g,jb,index
  INTEGER :: nbndval,nbndval_full
  REAL(DP) :: peso
  !
  REAL(DP),EXTERNAL :: integrate_imfreq
  !
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  !
  REAL(DP) :: partial_b,partial_h
  REAL(DP) :: segno, enrg, enrg1 
  COMPLEX(DP) :: residues_b,residues_h
  LOGICAL :: this_is_a_pole
  !
  ! PRINT TITLE of CALC
  !
  IF(l_verbose) CALL io_push_title('Sigma-C')
  !
  ! ZERO
  !
  sigma_corr = 0._DP
  !
  CALL pot3D%init('Wave',.FALSE.,'default')
  !
  ! -----------------------------------
  ! The part with imaginary integration
  ! -----------------------------------
  !
  IF(l_verbose) WRITE(stdout,'(5x,"Integrating along the IM axis...")')
  IF(l_verbose) CALL io_push_bar
  !
  IF(l_verbose) THEN
     IF (l_enable_off_diagonal .AND. l_full) THEN
        barra_load = kpt_pool%nloc * npair
     ELSE
        barra_load = kpt_pool%nloc * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
     ENDIF
  ENDIF
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_i', barra_load )
  !
  ! LOOP
  !
  DO iks = 1, kpt_pool%nloc ! KPOINT-SPIN
     !
     iks_g = kpt_pool%l2g(iks)
     !
     nbndval = nbnd_occ(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        DO jb = qp_bandrange(1), qp_bandrange(2)
           !
           IF (l_enable_off_diagonal .AND. l_full .AND. jb <= ib .OR. &
           & l_enable_off_diagonal .AND. .NOT. l_full .AND. jb == ib) THEN
              index = ijpmap(jb,ib)
           ELSEIF ( .NOT. l_enable_off_diagonal .AND. jb == ib ) THEN
              CONTINUE
           ELSE
              CYCLE
           ENDIF
           !
           ! HEAD PART
           !
           partial_h = 0._DP
           ! 
           IF(l_macropol .AND. jb == ib) THEN
              !
              DO ifreq = 1,ifr%nloc
                 enrg = et(ib,iks) - energy(ib,iks_g)
                 partial_h = partial_h + d_head_ifr(ifreq)*integrate_imfreq(ifreq,enrg)
              ENDDO
              !
           ENDIF
           !
           ! BODY 1st part : poles of H
           !
           partial_b = 0._DP 
           !
           DO ifreq = 1,ifr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 enrg = et(glob_im,iks) - energy(ib,iks_g)
                 IF (l_enable_off_diagonal .AND. l_full .AND. jb <= ib .OR. &
                 & l_enable_off_diagonal .AND. .NOT. l_full .AND. jb == ib) THEN
                    enrg1 = et(glob_im,iks) - energy(jb,iks_g)
                    partial_b = partial_b + d_body1_ifr_full(im,ifreq,index,iks_g)*0.5_DP*&
                    &(integrate_imfreq(ifreq,enrg) + integrate_imfreq(ifreq,enrg1))
                 ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                    partial_b = partial_b + d_body1_ifr(im,ifreq,ib,iks_g)*integrate_imfreq(ifreq,enrg)
                 ENDIF
              ENDDO
           ENDDO
           !
           ! BODY 2nd part : Lanczos
           !
           IF( l_enable_lanczos ) THEN
              !
              DO ifreq = 1,ifr%nloc
                 DO ip = 1, pert%nloc
                    DO il = 1, n_lanczos
                       IF (l_enable_off_diagonal .AND. l_full .AND. jb <= ib .OR. &
                       & l_enable_off_diagonal .AND. .NOT. l_full .AND. jb == ib) THEN
                          enrg = d_diago_full(il,ip,index,iks_g) - energy(ib,iks_g)
                          enrg1 = d_diago_full(il,ip,index,iks_g) - energy(jb,iks_g)
                          partial_b = partial_b + d_body2_ifr_full(il,ip,ifreq,index,iks_g)* 0.5_DP &
                          &* (integrate_imfreq(ifreq,enrg) + integrate_imfreq(ifreq,enrg1))
                       ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                          enrg = d_diago(il,ip,ib,iks_g) - energy(ib,iks_g)
                          partial_b = partial_b + d_body2_ifr(il,ip,ifreq,ib,iks_g)*integrate_imfreq(ifreq,enrg)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
              !
           ENDIF
           !
           CALL mp_sum( partial_h, intra_bgrp_comm) 
           CALL mp_sum( partial_b, intra_bgrp_comm) 
           CALL mp_sum( partial_b, inter_image_comm) 
           !
           IF (jb == ib) sigma_corr(ib,iks_g) = sigma_corr(ib,iks_g) &
           & + CMPLX( partial_b/omega/pi + partial_h*pot3D%div/pi, 0._DP, KIND=DP ) 
           IF (l_enable_off_diagonal .AND. l_full ) sigma_corr_full(index,iks_g) = sigma_corr_full(index,iks_g) &
           & + CMPLX( partial_b/omega/pi + partial_h*pot3D%div/pi, 0._DP, KIND=DP ) 
           !
           IF(l_verbose) CALL update_bar_type( barra, 'sigmac_i', 1 )
           !
        ENDDO !jb
        !
     ENDDO ! ib
     !
  ENDDO ! iks
  !
  IF(l_verbose) CALL stop_bar_type( barra, 'sigmac_i' )
  !
  ! -------------------
  ! The part with poles
  ! -------------------
  !
  IF(l_verbose) WRITE(stdout,'(5x,"Residues along the RE axis...")')
  IF(l_verbose) CALL io_push_bar
  !
  IF(l_verbose) THEN
     IF (l_enable_off_diagonal .AND. l_full) THEN
        barra_load = kpt_pool%nloc * npair
     ELSE
        barra_load = kpt_pool%nloc * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
     ENDIF
  ENDIF  
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_r', barra_load )
  !
  ! LOOP
  !
  DO iks = 1, kpt_pool%nloc
     !
     iks_g = kpt_pool%l2g(iks)
     nbndval = nbnd_occ(iks)
     IF (l_frac_occ) nbndval_full = nbnd_occ_full(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        enrg = energy(ib,iks_g)
        !
        DO jb = qp_bandrange(1), qp_bandrange(2)
           !
           IF (l_enable_off_diagonal .AND. l_full .AND. jb <= ib .OR. &
           & l_enable_off_diagonal .AND. .NOT. l_full .AND. jb == ib) THEN
              index = ijpmap(jb,ib)
           ELSEIF ( .NOT. l_enable_off_diagonal .AND. jb == ib ) THEN
              CONTINUE
           ELSE
              CYCLE
           ENDIF
           !
           enrg1 = energy(jb,iks_g)
           !
           residues_b = 0._DP
           residues_h = 0._DP
           !
           DO im = 1,aband%nloc
              !
              glob_im = aband%l2g(im)
              !
              this_is_a_pole=.FALSE.
              !
              IF (l_frac_occ) THEN
                 IF( glob_im > nbndval_full .AND. glob_im <= nbndval ) THEN
                    peso = occupation(glob_im,iks)
                 ELSE
                    peso = 1._DP
                 ENDIF
              ELSE
                 peso = 1._DP
              ENDIF
              !
              IF( glob_im <= nbndval ) THEN
                 segno = -1._DP
                 IF( et(glob_im,iks) - enrg >  0.00001_DP ) this_is_a_pole=.TRUE.
              ELSE
                 segno = 1._DP
                 IF( et(glob_im,iks) - enrg < -0.00001_DP ) this_is_a_pole=.TRUE.
              ENDIF  
              !
              IF( this_is_a_pole ) THEN
                 !
                 CALL retrieve_glob_freq( et(glob_im,iks) - enrg, glob_ifreq )
                 !
                 DO ifreq = 1, rfr%nloc 
                    !
                    IF( rfr%l2g(ifreq) .NE. glob_ifreq ) CYCLE
                    !
                    IF(glob_im==ib.AND.l_macropol.AND.jb==ib) residues_h = residues_h &
                    & + peso * segno * z_head_rfr(ifreq)
                    !
                    IF (l_enable_off_diagonal .AND. l_full .AND. jb <= ib .OR. &
                    & l_enable_off_diagonal .AND. .NOT. l_full .AND. jb == ib) THEN
                       residues_b = residues_b + 0.5_DP * peso * segno &
                       & * z_body_rfr_full( im, ifreq, index, iks_g )
                    ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                       residues_b = residues_b + peso * segno * z_body_rfr( im, ifreq, ib, iks_g )
                    ENDIF
                    ! 
                 ENDDO
                 !
              ENDIF 
              !
              IF (l_frac_occ) THEN
                 !
                 this_is_a_pole=.FALSE.
                 IF( glob_im > nbndval_full .AND. glob_im <= nbndval ) THEN
                    segno = 1._DP
                    IF( et(glob_im,iks) - enrg < -0.00001_DP ) this_is_a_pole=.TRUE.
                 ENDIF
                 !
                 IF( this_is_a_pole ) THEN
                    !
                    CALL retrieve_glob_freq( et(glob_im,iks) - enrg, glob_ifreq )
                    !
                    DO ifreq = 1, rfr%nloc
                       !
                       IF( rfr%l2g(ifreq) /= glob_ifreq ) CYCLE
                       !
                       IF(glob_im==ib.AND.l_macropol) residues_h = residues_h &
                       & + (1.0 - peso) * segno * z_head_rfr(ifreq)
                       !
                       IF (l_enable_off_diagonal .AND. l_full .AND. jb <= ib .OR. &
                       & l_enable_off_diagonal .AND. .NOT. l_full .AND. jb == ib) THEN    
                          residues_b = residues_b + (1._DP - peso) * segno &
                          & * z_body_rfr_full( im, ifreq, index, iks_g )
                       ELSEIF (.NOT. l_enable_off_diagonal .AND. jb == ib) THEN
                          residues_b = residues_b + (1._DP - peso) * segno &
                          & * z_body_rfr( im, ifreq, ib, iks_g )
                       ENDIF
                       !
                    ENDDO
                    !
                 ENDIF
                 !
              ENDIF
              !
           ENDDO ! im
           !
           IF (l_enable_off_diagonal .AND. l_full .AND. jb <= ib .OR. &
           & l_enable_off_diagonal .AND. .NOT. l_full .AND. jb == ib) THEN
              !
              DO im = 1,aband%nloc
                 !
                 glob_im = aband%l2g(im)
                 !
                 this_is_a_pole=.false.
                 !
                 IF (l_frac_occ) THEN
                    IF( glob_im > nbndval_full .AND. glob_im <= nbndval ) THEN
                       peso = occupation(glob_im,iks)
                    ELSE
                       peso = 1._DP
                    ENDIF
                 ELSE
                    peso = 1._DP
                 ENDIF
                 !
                 IF( glob_im <= nbndval ) THEN
                    segno = -1._DP
                    IF( et(glob_im,iks) - enrg1 >  0.00001_DP ) this_is_a_pole=.TRUE.
                 ELSE
                    segno = 1._DP
                    IF( et(glob_im,iks) - enrg1 < -0.00001_DP ) this_is_a_pole=.TRUE.
                 ENDIF  
                 !
                 IF( this_is_a_pole ) THEN
                    !
                    CALL retrieve_glob_freq( et(glob_im,iks) - enrg1, glob_ifreq )
                    !
                    DO ifreq = 1, rfr%nloc 
                       !
                       IF( rfr%l2g(ifreq) .NE. glob_ifreq ) CYCLE
                       !
                       residues_b = residues_b + 0.5_DP * peso * segno &
                       & * z_body_rfr_full( im, ifreq, index, iks_g )
                       ! 
                    ENDDO
                    !
                 ENDIF 
                 !
                 IF (l_frac_occ) THEN
                    !
                    this_is_a_pole=.FALSE.
                    IF( glob_im > nbndval_full .AND. glob_im <= nbndval ) THEN
                       segno = 1._DP
                       IF( et(glob_im,iks) - enrg1 < -0.00001_DP ) this_is_a_pole=.TRUE.
                    ENDIF
                    !
                    IF( this_is_a_pole ) THEN
                       !
                       CALL retrieve_glob_freq( et(glob_im,iks) - enrg1, glob_ifreq )
                       !
                       DO ifreq = 1, rfr%nloc 
                          !
                          IF( rfr%l2g(ifreq) .NE. glob_ifreq ) CYCLE
                          !
                          residues_b = residues_b + 0.5_DP * peso * segno &
                          & * z_body_rfr_full( im, ifreq, index, iks_g )
                          ! 
                       ENDDO
                       !
                    ENDIF 
                    !
                 ENDIF
                 !
              ENDDO ! im
              !
           ENDIF
           !
           CALL mp_sum( residues_h, intra_bgrp_comm )
           CALL mp_sum( residues_h, inter_image_comm )
           CALL mp_sum( residues_b, intra_bgrp_comm )
           CALL mp_sum( residues_b, inter_image_comm )
           !
           IF (jb == ib) sigma_corr(ib,iks_g) = sigma_corr(ib,iks_g) + residues_b/omega + residues_h*pot3D%div
           IF (l_enable_off_diagonal .AND. l_full) sigma_corr_full(index,iks_g) = sigma_corr_full(index,iks_g) &
           & + residues_b/omega + residues_h*pot3D%div
           !
           IF(l_verbose) CALL update_bar_type( barra, 'sigmac_r', 1 )
           !
        ENDDO ! jb
        !
     ENDDO ! ib
     !
  ENDDO ! ik
  !
  IF(l_verbose) CALL stop_bar_type( barra, 'sigmac_r' )
  !
  CALL mp_sum( sigma_corr, inter_pool_comm )
  IF (l_enable_off_diagonal .AND. l_full) CALL mp_sum( sigma_corr_full, inter_pool_comm )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE calc_corr_k( sigma_corr, energy, l_verbose)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_corr(n,iks) = < ib,iks | S_c(energy(ib,iks))  | ib,iks >
  ! ... ib = qp_bandrange(1):qp_bandrange(2)
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE constants,            ONLY : pi
  USE pwcom,                ONLY : et
  USE westcom,              ONLY : qp_bandrange,nbnd_occ,l_enable_lanczos,n_lanczos,l_macropol,&
                                 & z_head_ifr,z_head_rfr,z_body1_ifr_q,z_body2_ifr_q,d_diago_q,&
                                 & z_body_rfr_q
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar,io_push_title
  USE distribution_center,  ONLY : pert,ifr,rfr,aband
  USE types_bz_grid,        ONLY : k_grid,q_grid
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP),INTENT(OUT) :: sigma_corr( qp_bandrange(1):qp_bandrange(2), k_grid%nps )  ! The correlation self-energy,
                                                                                        ! imaginary part is lifetime.
  REAL(DP),INTENT(IN) :: energy( qp_bandrange(1):qp_bandrange(2), k_grid%nps )          ! The energy variable
  LOGICAL,INTENT(IN) :: l_verbose
  !
  ! Workspace
  !
  INTEGER :: ik,ikk,iks,ikks,iq,ib,ifreq,glob_ifreq,il,im,glob_im,ip,is,iss
  INTEGER :: nbndval
  !
  REAL(DP),EXTERNAL :: integrate_imfreq
  !
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  !
  COMPLEX(DP) :: partial_b,partial_h
  REAL(DP) :: segno, enrg
  REAL(DP) :: g0(3)
  COMPLEX(DP) :: residues_b,residues_h
  LOGICAL :: this_is_a_pole
  LOGICAL :: l_gammaq
  !
  ! PRINT TITLE of CALC
  !
  IF(l_verbose) CALL io_push_title('Sigma_C')
  !
  ! ZERO
  !
  sigma_corr = 0._DP
  !
  ! -----------------------------------
  ! The part with imaginary integration
  ! -----------------------------------
  !
  IF(l_verbose) WRITE(stdout,'(5x,"Integrating along the IM axis...")')
  IF(l_verbose) CALL io_push_bar
  !
  IF(l_verbose) barra_load = k_grid%nps * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_i', barra_load )
  !
  ! LOOP
  !
  DO iks = 1, k_grid%nps   ! KPOINT-SPIN (MATRIX ELEMENT)
     !
     ik = k_grid%ip(iks)
     is = k_grid%is(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        partial_h = 0._DP
        partial_b = 0._DP
        !
!        DO iq = 1, q_grid%np   ! Q-POINT
        DO ikks = 1, k_grid%nps    ! KPOINT-SPIN (INTEGRAL OVER K')
           !
           ikk = k_grid%ip(ikks)
           iss = k_grid%is(ikks)
           IF( is /= iss ) CYCLE
           !
           !CALL k_grid%find( k_grid%p_cart(:,ik) - k_grid%p_cart(:,ikk), 1, 'cart', iq, g0 ) !M
           CALL k_grid%find( k_grid%p_cart(:,ik) - k_grid%p_cart(:,ikk), 'cart', iq, g0 )
           l_gammaq = q_grid%l_pIsGamma(iq)
           nbndval = nbnd_occ(ikks)
           !
           CALL pot3D%init('Wave',.TRUE.,'default',iq)
           !
           ! HEAD PART
           !
           IF(l_macropol .AND. l_gammaq) THEN
              !
              DO ifreq = 1,ifr%nloc
                 enrg = et(ib,iks) - energy(ib,iks)
                 partial_h = partial_h + z_head_ifr(ifreq)*integrate_imfreq(ifreq,enrg)
              ENDDO
              !
           ENDIF
           !
           ! BODY 1st part : poles of H
           !
           DO ifreq = 1,ifr%nloc
              DO im = 1, aband%nloc
                 glob_im = aband%l2g(im)
                 enrg = et(glob_im,ikks) - energy(ib,iks)
                 partial_b = partial_b + z_body1_ifr_q(im,ifreq,ib,iks,iq)*integrate_imfreq(ifreq,enrg)*q_grid%weight(iq)
              ENDDO
           ENDDO
           !
           !
           ! BODY 2nd part : Lanczos
           !
           IF( l_enable_lanczos ) THEN
              !
              DO ifreq = 1,ifr%nloc
                 DO ip = 1, pert%nloc
                    DO il = 1, n_lanczos
                       enrg = d_diago_q(il,ip,ib,iks,iq) - energy(ib,iks)
                       partial_b = partial_b + z_body2_ifr_q(il,ip,ifreq,ib,iks,iq)*integrate_imfreq(ifreq,enrg)*q_grid%weight(iq)
                    ENDDO
                 ENDDO
              ENDDO
              !
           ENDIF
           !
        ENDDO ! ikks
        !
        CALL mp_sum( partial_h, intra_bgrp_comm)
        CALL mp_sum( partial_b, intra_bgrp_comm)
        CALL mp_sum( partial_b, inter_image_comm)
        !
        sigma_corr(ib,iks) = sigma_corr(ib,iks) + partial_b/omega/pi + partial_h*pot3D%div/pi
        !
        IF(l_verbose) CALL update_bar_type( barra, 'sigmac_i', 1 )
        !
     ENDDO ! ib
     !
  ENDDO ! iks
  !
  IF(l_verbose) CALL stop_bar_type( barra, 'sigmac_i' )
  !
  ! -------------------
  ! The part with poles
  ! -------------------
  !
  IF(l_verbose) WRITE(stdout,'(5x,"Residues along the RE axis...")')
  IF(l_verbose) CALL io_push_bar
  !
  IF(l_verbose) barra_load = k_grid%nps * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_r', barra_load )
  !
  ! LOOP
  !
  DO iks = 1, k_grid%nps   ! KPOINT-SPIN (MATRIX ELEMENT)
     !
     ik = k_grid%ip(iks)
     is = k_grid%is(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        enrg = energy(ib,iks)
        !
        residues_b = 0._DP
        residues_h = 0._DP
        !
!        DO iq = 1, q_grid%np   ! Q-POINT
        DO ikks = 1, k_grid%nps   ! KPOINT-SPIN (INTEGRAL OVER K')
           !
           ikk = k_grid%ip(ikks)
           iss = k_grid%is(ikks)
           !
           IF( is /= iss ) CYCLE
           !
           !CALL k_grid%find( k_grid%p_cart(:,ik) - k_grid%p_cart(:,ikk), 1, 'cart', iq, g0 )  !M
           CALL k_grid%find( k_grid%p_cart(:,ik) - k_grid%p_cart(:,ikk), 'cart', iq, g0 )
           l_gammaq = q_grid%l_pIsGamma(iq)
           nbndval = nbnd_occ(ikks)
           !
           DO im = 1,aband%nloc
              !
              glob_im = aband%l2g(im)
              !
              this_is_a_pole=.FALSE.
              IF( glob_im <= nbndval ) THEN ! poles inside G+
                 segno = -1._DP
                 IF( et(glob_im,ikks) - enrg >  0.00001_DP ) this_is_a_pole=.TRUE.
              ELSE ! poles inside G-
                 segno = 1._DP
                 IF( et(glob_im,ikks) - enrg < -0.00001_DP ) this_is_a_pole=.TRUE.
              ENDIF
              !
              IF( this_is_a_pole ) THEN
                 !
                 CALL retrieve_glob_freq( et(glob_im,ikks) - enrg, glob_ifreq )
                 !
                 DO ifreq = 1, rfr%nloc
                    !
                    IF( rfr%l2g(ifreq) .NE. glob_ifreq ) CYCLE
                    !
                    IF(glob_im==ib.AND.l_macropol.AND.l_gammaq) residues_h = residues_h + segno * z_head_rfr(ifreq)
                    !
                    residues_b = residues_b + segno * z_body_rfr_q( im, ifreq, ib, iks, iq )*q_grid%weight(iq)
                    !
                 ENDDO
                 !
              ENDIF
              !
           ENDDO ! im
           !
        ENDDO ! ikks
        !
        CALL mp_sum( residues_h, intra_bgrp_comm )
        CALL mp_sum( residues_h, inter_image_comm )
        CALL mp_sum( residues_b, intra_bgrp_comm )
        CALL mp_sum( residues_b, inter_image_comm )
        !
        sigma_corr(ib,iks) = sigma_corr(ib,iks) + residues_b/omega + residues_h*pot3D%div
        !
        IF(l_verbose) CALL update_bar_type( barra, 'sigmac_r', 1 )
        !
     ENDDO ! ib
     !
  ENDDO ! iks
  !
  IF(l_verbose) CALL stop_bar_type( barra, 'sigmac_r' )
  !
END SUBROUTINE
!
FUNCTION integrate_imfreq( ifreq, c )
  !
  ! = int_{a}^{b}  dz  \frac{ c }{ c^2 + z^2 }
  !
  USE kinds,     ONLY : DP
  USE westcom,   ONLY : imfreq_list_integrate
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ifreq
  REAL(DP),INTENT(IN) :: c
  REAL(DP) :: integrate_imfreq
  !
  REAL(DP) :: a, b
  !
  integrate_imfreq = 0._DP
  !
  IF( ABS( c ) < 0.000001_DP ) RETURN
  !
  a = imfreq_list_integrate( 1, ifreq )
  b = imfreq_list_integrate( 2, ifreq )
  !
  integrate_imfreq = ATAN( c * (b-a) / (c*c+a*b) )
  !
END FUNCTION
!
SUBROUTINE retrieve_glob_freq( freq, glob_ifreq )
  !
  USE kinds,               ONLY : DP
  USE westcom,             ONLY : n_refreq, ecut_refreq
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(IN) :: freq
  INTEGER, INTENT(OUT) :: glob_ifreq
  !
  ! Find glob_ifreq
  !
  glob_ifreq = 1 + NINT( REAL(n_refreq-1,KIND=DP) * ABS(freq) / ecut_refreq )
  glob_ifreq = MIN( n_refreq, glob_ifreq )
  glob_ifreq = MAX( 1, glob_ifreq )
  !
END SUBROUTINE
