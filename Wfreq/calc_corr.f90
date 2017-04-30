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
SUBROUTINE calc_corr_gamma( sigma_corr, energy, l_verbose)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_corr(n,iks) = < ib,iks | S_c(energy(ib,iks))  | ib,iks >     
  ! ... ib = qp_bandrange(1):qp_bandrange(2)
  !
  USE kinds,                ONLY : DP 
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm,world_comm
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE constants,            ONLY : tpi,fpi,rytoev,e2,pi
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,g2kin,nspin,current_k,wk
  USE westcom,              ONLY : qp_bandrange,isz,&
                                 & nbnd_occ,l_enable_lanczos,&
                                 & n_lanczos,iks_l2g,l_macropol,&
                                 & d_head_ifr,z_head_rfr,d_body1_ifr,d_body2_ifr,d_diago,z_body_rfr
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar,io_push_value,io_push_title
  USE distribution_center,  ONLY : pert,ifr,rfr,aband
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP),INTENT(OUT) :: sigma_corr( qp_bandrange(1):qp_bandrange(2), nks )  ! The correlation self-energy, imaginary part is lifetime.  
  REAL(DP),INTENT(IN) :: energy( qp_bandrange(1):qp_bandrange(2), nks )          ! The energy variable
  LOGICAL,INTENT(IN) :: l_verbose
  !
  ! Workspace
  !
  INTEGER :: iks,ib,ifreq,glob_ifreq,il,im,glob_im,ip
  INTEGER :: nbndval
  !
  REAL(DP),EXTERNAL :: integrate_imfreq
  !
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  !
  REAL(DP) :: partial_b,partial_h
  REAL(DP) :: segno, enrg 
  COMPLEX(DP) :: residues_b,residues_h
  LOGICAL :: this_is_a_pole
  !
  ! PRINT TITLE of CALC
  !
  IF(l_verbose) CALL io_push_title('Sigma-C')
  !
  ! BARRIER : ALL 
  !
  CALL mp_barrier( world_comm )
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
  IF(l_verbose) barra_load = nks * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_i', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     nbndval = nbnd_occ(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        ! HEAD PART
        !
        partial_h = 0._DP
        ! 
        IF(l_macropol) THEN
           !
           DO ifreq = 1,ifr%nloc
              enrg = et(ib,iks) - energy(ib,iks)
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
              enrg = et(glob_im,iks) - energy(ib,iks)
              partial_b = partial_b + d_body1_ifr(im,ifreq,ib,iks)*integrate_imfreq(ifreq,enrg)
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
                    enrg = d_diago(il,ip,ib,iks) - energy(ib,iks)
                    partial_b = partial_b + d_body2_ifr(il,ip,ifreq,ib,iks)*integrate_imfreq(ifreq,enrg)
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
        sigma_corr(ib,iks) = sigma_corr(ib,iks) + CMPLX( partial_b/omega/pi + partial_h*isz/pi, 0._DP, KIND=DP ) 
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
  IF(l_verbose) barra_load = nks * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_r', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, nks
     !
     nbndval = nbnd_occ(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        enrg = energy(ib,iks)
        !
        residues_b = 0._DP
        residues_h = 0._DP
        !
        DO im = 1,aband%nloc
           !
           glob_im = aband%l2g(im)
           !
           this_is_a_pole=.false.
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
                 IF(glob_im==ib.AND.l_macropol) residues_h = residues_h + segno * z_head_rfr(ifreq)
                 !
                 residues_b = residues_b + segno * z_body_rfr( im, ifreq, ib, iks )
                 ! 
              ENDDO
              !
           ENDIF 
           !
        ENDDO ! im
        !
        CALL mp_sum( residues_h, intra_bgrp_comm )
        CALL mp_sum( residues_h, inter_image_comm )
        CALL mp_sum( residues_b, intra_bgrp_comm )
        CALL mp_sum( residues_b, inter_image_comm )
        !
        sigma_corr(ib,iks) = sigma_corr(ib,iks) + residues_b/omega + residues_h*isz
        !
        IF(l_verbose) CALL update_bar_type( barra, 'sigmac_r', 1 )
        !
     ENDDO ! ib  
     !
  ENDDO ! ik
  !
  IF(l_verbose) CALL stop_bar_type( barra, 'sigmac_r' )
  !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE calc_corr_k( sigma_corr, energy, l_verbose)
  !-----------------------------------------------------------------------
  !
  ! store in sigma_corr(n,iks) = < ib,iks | S_c(energy(ib,iks))  | ib,iks >     
  ! ... ib = qp_bandrange(1):qp_bandrange(2)
  !
  USE kinds,                ONLY : DP 
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm,world_comm
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE constants,            ONLY : tpi,fpi,rytoev,e2,pi
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,g2kin,nspin,current_k,wk
  USE westcom,              ONLY : qp_bandrange,isz,&
                                 & nbnd_occ,l_enable_lanczos,&
                                 & n_lanczos,iks_l2g,l_macropol,&
                                 & z_head_ifr,z_head_rfr,z_body1_ifr,z_body2_ifr,d_diago,z_body_rfr
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar,io_push_value,io_push_title
  USE distribution_center,  ONLY : pert,ifr,rfr,aband
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP),INTENT(OUT) :: sigma_corr( qp_bandrange(1):qp_bandrange(2), nks )  ! The correlation self-energy, imaginary part is lifetime.  
  REAL(DP),INTENT(IN) :: energy( qp_bandrange(1):qp_bandrange(2), nks )          ! The energy variable
  LOGICAL,INTENT(IN) :: l_verbose
  !
  ! Workspace
  !
  INTEGER :: iks,ib,ifreq,glob_ifreq,il,im,glob_im,ip
  INTEGER :: nbndval
  !
  REAL(DP),EXTERNAL :: integrate_imfreq
  !
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  !
  COMPLEX(DP) :: partial_b,partial_h
  REAL(DP) :: segno, enrg 
  COMPLEX(DP) :: residues_b,residues_h
  LOGICAL :: this_is_a_pole
  !
  ! PRINT TITLE of CALC
  !
  IF(l_verbose) CALL io_push_title('Sigma_C')
  !
  ! BARRIER : ALL 
  !
  CALL mp_barrier( world_comm )
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
  IF(l_verbose) barra_load = nks * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_i', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     nbndval = nbnd_occ(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        ! HEAD PART
        !
        partial_h = 0._DP
        ! 
        IF(l_macropol) THEN
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
        partial_b = 0._DP 
        !
        DO ifreq = 1,ifr%nloc
           DO im = 1, aband%nloc
              glob_im = aband%l2g(im)
              enrg = et(glob_im,iks) - energy(ib,iks)
              partial_b = partial_b + z_body1_ifr(im,ifreq,ib,iks)*integrate_imfreq(ifreq,enrg)
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
                    enrg = d_diago(il,ip,ib,iks) - energy(ib,iks)
                    partial_b = partial_b + z_body2_ifr(il,ip,ifreq,ib,iks)*integrate_imfreq(ifreq,enrg)
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
        sigma_corr(ib,iks) = sigma_corr(ib,iks) + partial_b/omega/pi + partial_h*isz/pi 
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
  IF(l_verbose) barra_load = nks * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  IF(l_verbose) CALL start_bar_type( barra, 'sigmac_r', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, nks
     !
     nbndval = nbnd_occ(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        enrg = energy(ib,iks)
        !
        residues_b = 0._DP
        residues_h = 0._DP
        !
        DO im = 1,aband%nloc
           !
           glob_im = aband%l2g(im)
           !
           this_is_a_pole=.false.
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
                 IF(glob_im==ib.AND.l_macropol) residues_h = residues_h + segno * z_head_rfr(ifreq)
                 !
                 residues_b = residues_b + segno * z_body_rfr( im, ifreq, ib, iks )
                 ! 
              ENDDO
              !
           ENDIF 
           !
        ENDDO ! im
        !
        CALL mp_sum( residues_h, intra_bgrp_comm )
        CALL mp_sum( residues_h, inter_image_comm )
        CALL mp_sum( residues_b, intra_bgrp_comm )
        CALL mp_sum( residues_b, inter_image_comm )
        !
        sigma_corr(ib,iks) = sigma_corr(ib,iks) + residues_b/omega + residues_h*isz
        !
        IF(l_verbose) CALL update_bar_type( barra, 'sigmac_r', 1 )
        !
     ENDDO ! ib  
     !
  ENDDO ! ik
  !
  IF(l_verbose) CALL stop_bar_type( barra, 'sigmac_r' )
  !
END SUBROUTINE
!
!
!
!
!
!
FUNCTION integrate_imfreq( ifreq, c )
  !
  ! = int_{a}^{b}  dz  \frac{ c }{ c^2 + z^2 }  
  !
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : et
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
!
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
