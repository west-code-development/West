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
SUBROUTINE solve_eri(ifreq,l_isFreqReal)
  !-----------------------------------------------------------------------
  !
  USE westcom,              ONLY : n_pdep_eigen_to_use,d_epsm1_ifr_a,z_epsm1_rfr_a,imfreq_list,&
                                 & refreq_list,z_head_rfr_a,d_head_ifr_a,n_pairs,eri_w
  USE distribution_center,  ONLY : pert,ifr,rfr
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : intra_bgrp_comm,me_bgrp
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE wfreq_db,             ONLY : qdet_db_write_eri
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ifreq
  LOGICAL,INTENT(IN) :: l_isFreqReal
  !
  INTEGER :: who, iloc
  REAL(DP) :: freq
  COMPLEX(DP) :: chi_head
  !
  REAL(DP),ALLOCATABLE :: braket(:,:,:)
  REAL(DP),ALLOCATABLE :: eri_vc(:,:,:,:)
  COMPLEX(DP),ALLOCATABLE :: chi_body(:,:)
  !
  ! TODO: check if G=0 component is double counted in any dot products
  !
  ! Compute 4-center integrals of W (screened electron repulsion integrals, eri)
  ! (ij|W|kl) = int dx dx' f_i(x) f_j(x) W(r,r') f_k(x') f_l(x')
  ! f_i is the Fourier transform to direct space of westcom/proj_c_i
  !
  CALL io_push_title('ERI')
  !
  ALLOCATE( chi_body( pert%nglob, pert%nloc ) )
  !
  ! Load chi at given frequency (index and logical from input)
  !
  IF ( l_isFreqReal ) THEN
     !
     CALL rfr%g2l(ifreq,iloc,who)
     IF ( me_bgrp == who ) THEN
        !
        freq = refreq_list(iloc)
        chi_head = z_head_rfr_a(iloc)
        chi_body(:,:) = z_epsm1_rfr_a(:,:,iloc)
        !
     ENDIF
     !
  ELSE
     !
     CALL ifr%g2l(ifreq,iloc,who)
     IF ( me_bgrp == who ) THEN
        !
        freq = imfreq_list(iloc)
        chi_head = CMPLX(d_head_ifr_a(ifreq),KIND=DP)
        chi_body(:,:) = CMPLX(d_epsm1_ifr_a(:,:,ifreq),KIND=DP)
        !
     ENDIF
     !
  ENDIF
  !
  CALL mp_bcast( freq, who, intra_bgrp_comm )
  CALL mp_bcast( chi_head, who, intra_bgrp_comm )
  CALL mp_bcast( chi_body, who, intra_bgrp_comm )
  !
  ! Compute ERI (Electron Repulsion Integrals)
  !
  ALLOCATE( eri_vc(n_pairs,n_pairs,nspin,nspin) )
  ALLOCATE( eri_w(n_pairs,n_pairs,nspin,nspin) )
  !
  ! 4-center integrals of Vc
  !
  CALL compute_eri_vc(eri_vc)
  !
  ! 4-center integrals of Wp
  !
  ALLOCATE( braket(n_pairs,nspin,n_pdep_eigen_to_use) )
  !
  CALL compute_braket(braket)
  CALL compute_eri_wp(braket, chi_head, chi_body, eri_w)
  !
  ! 4-center integrals of W
  !
  eri_w(:,:,:,:) = eri_vc + eri_w
  !
  CALL qdet_db_write_eri(eri_vc,eri_w)
  !
  CALL io_push_bar()
  !
  DEALLOCATE( chi_body )
  DEALLOCATE( braket )
  DEALLOCATE( eri_vc )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_braket(braket)
  !-----------------------------------------------------------------------
  !
  ! braket_{ij(s)m} = < ij(s) | pdep_m * V_c^0.5 >
  ! ij is a pair of functions taken from westcom/proj_c
  ! s is the spin index
  ! m is the PDEP index
  ! V_c is the bare Coulomb potential
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin
  USE gvect,                ONLY : gstart
  USE westcom,              ONLY : wstat_save_dir,npwq,npwqx,n_pdep_eigen_to_use,fftdriver,&
                                 & proj_c,n_pairs,pijmap
  USE mp_global,            ONLY : intra_bgrp_comm,inter_image_comm
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,double_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE mp,                   ONLY : mp_sum
  USE types_coulomb,        ONLY : pot3D
  USE io_push,              ONLY : io_push_title
  USE distribution_center,  ONLY : macropert
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute
  USE wavefunctions,        ONLY : psic
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(OUT) :: braket(n_pairs,nspin,n_pdep_eigen_to_use)
  !
  COMPLEX(DP),ALLOCATABLE :: phi(:), rho_r(:), rho_g(:)
  !
  INTEGER :: s, m, p1, i, j, ig, ir, mloc
  CHARACTER(LEN=25) :: filepot
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  !
  REAL(DP),EXTERNAL :: DDOT
  !
  CALL io_push_title('braket')
  !
  ALLOCATE( phi(npwqx) )
  ALLOCATE( rho_r(dffts%nnr) )
  ALLOCATE( rho_g(npwqx) )
  !
  CALL pot3D%init(fftdriver,.FALSE.,'default')
  !
  barra_load = 0
  DO mloc = 1, macropert%nloc
     m = macropert%l2g(mloc)
     IF (m > n_pdep_eigen_to_use) CYCLE
     barra_load = barra_load+1
  ENDDO
  barra_load = barra_load*nspin
  CALL start_bar_type ( barra, 'eri_brak', barra_load )
  !
  braket = 0._DP
  !
  DO s = 1, nspin
     !
     DO mloc = 1, macropert%nloc
        !
        m = macropert%l2g(mloc)
        !
        IF (m > n_pdep_eigen_to_use) CYCLE ! skip the head, use only the body
        !
        ! Read the m-th PDEP
        !
        CALL generate_pdep_fname( filepot, m )
        fname = TRIM( wstat_save_dir ) // '/'// filepot
        CALL pdep_read_G_and_distribute(fname,phi(:))
        !
        ! Mulitply by V_c^0.5
        !
        DO ig = 1, npwq
           phi(ig) = pot3D%sqvc(ig) * phi(ig)
        ENDDO
        !
        ! Compute the braket for each pair of functions
        !
        DO p1 = 1, n_pairs
           !
           i = pijmap(1,p1)
           j = pijmap(2,p1)
           !
           CALL double_invfft_gamma(dffts,npwq,npwqx,proj_c(:,i,s),proj_c(:,j,s),psic,'Wave')
           !
           DO ir = 1, dffts%nnr
              rho_r(ir) = REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
           ENDDO
           !
           CALL single_fwfft_gamma(dffts,npwq,npwqx,rho_r,rho_g,TRIM(fftdriver))
           !
           ! assume Gamma-only case
           !
           braket(p1,s,m) = 2._DP * DDOT(2*npwq,rho_g,1,phi,1)
           IF(gstart == 2) braket(p1,s,m) = braket(p1,s,m) - REAL(rho_g(1),KIND=DP)*REAL(phi(1),KIND=DP)
           !
        ENDDO
        !
        CALL update_bar_type( barra,'eri_brak', 1 )
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(braket, inter_image_comm)
  CALL mp_sum(braket, intra_bgrp_comm)
  !
  CALL stop_bar_type( barra,'eri_brak' )
  !
  DEALLOCATE( phi )
  DEALLOCATE( rho_r )
  DEALLOCATE( rho_g )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_eri_vc(eri_vc)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin
  USE westcom,              ONLY : n_pairs,pijmap,proj_c,npwq,npwqx
  USE mp_global,            ONLY : intra_bgrp_comm,inter_image_comm
  USE distribution_center,  ONLY : bandpair
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,double_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE mp,                   ONLY : mp_sum
  USE types_coulomb,        ONLY : pot3D
  USE io_push,              ONLY : io_push_title
  USE gvect,                ONLY : gstart,ngm
  USE cell_base,            ONLY : omega
  USE wavefunctions,        ONLY : psic
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(OUT) :: eri_vc(n_pairs,n_pairs,nspin,nspin)
  !
  COMPLEX(DP),ALLOCATABLE :: rho_g1(:), rho_g2(:), rho_r(:)
  !
  INTEGER :: i, j, k, l, p1, p1loc, p2, s1, s2
  INTEGER :: ir, ig
  TYPE(bar_type) :: barra
  !
  COMPLEX(DP),EXTERNAL :: ZDOTC
  REAL(DP),EXTERNAL :: DDOT
  !
  CALL io_push_title('pair density')
  !
  ! Distribute pairs over images
  !
  CALL bandpair%init(n_pairs,'i','n_pairs',.FALSE.)
  !
  CALL pot3D%init('Rho',.FALSE.,'gb')
  !
  ALLOCATE( rho_r(dffts%nnr) )
  ALLOCATE( rho_g1(ngm) )
  ALLOCATE( rho_g2(ngm) )
  !
  eri_vc(:,:,:,:) = 0._DP
  !
  ! eri_vc with SAME SPIN --> compute only p2>=p1, obtain p2<p1 by symmetry
  !
  CALL start_bar_type ( barra, 'eri_vc_s1s1', nspin*bandpair%nloc )
  !
  DO s1 = 1, nspin
     DO p1loc = 1, bandpair%nloc
        !
        p1 = bandpair%l2g(p1loc)
        i = pijmap(1,p1)
        j = pijmap(2,p1)
        !
        CALL double_invfft_gamma(dffts,npwq,npwqx,proj_c(:,i,s1),proj_c(:,j,s1),psic,'Wave')
        !
        DO ir = 1, dffts%nnr
           rho_r(ir) = REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,ngm,ngm,rho_r,rho_g1,'Rho')
        !
        DO ig = 1, ngm
           rho_g1(ig) = pot3D%sqvc(ig) * rho_g1(ig)
        ENDDO
        !
        ! (s1 = s2); (p2 = p1)
        !
        eri_vc(p1,p1,s1,s1) = eri_vc(p1,p1,s1,s1) &
        & + 2._DP * DDOT(2*(ngm-gstart+1),rho_g1(gstart:ngm),1,rho_g1(gstart:ngm),1)/omega
        IF(i==j .AND. gstart==2) THEN
           eri_vc(p1,p1,s1,s1) = eri_vc(p1,p1,s1,s1) + pot3D%div
        ENDIF
        !
        ! (s1 = s2); (p2 > p1)
        !
        DO p2 = p1+1, n_pairs
           !
           k = pijmap(1,p2)
           l = pijmap(2,p2)
           !
           CALL double_invfft_gamma(dffts,npwq,npwqx,proj_c(:,k,s1),proj_c(:,l,s1),psic,'Wave')
           !
           DO ir = 1, dffts%nnr
              rho_r(ir) = REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
           ENDDO
           !
           CALL single_fwfft_gamma(dffts,ngm,ngm,rho_r,rho_g2,'Rho')
           !
           DO ig = 1, ngm
              rho_g2(ig) = pot3D%sqvc(ig) * rho_g2(ig)
           ENDDO
           !
           eri_vc(p1,p2,s1,s1) = eri_vc(p1,p2,s1,s1) &
           & + 2._DP * DDOT(2*(ngm-gstart+1),rho_g1(gstart:ngm),1,rho_g2(gstart:ngm),1)/omega
           IF(i==j .AND. k==l .AND. gstart==2) THEN
              eri_vc(p1,p2,s1,s1) = eri_vc(p1,p2,s1,s1) + pot3D%div
           ENDIF
           !
           ! Apply symmetry
           !
           eri_vc(p2,p1,s1,s1) = eri_vc(p1,p2,s1,s1)
           !
        ENDDO
        !
        CALL update_bar_type ( barra, 'eri_vc_s1s1', 1 )
        !
     ENDDO
  ENDDO
  !
  CALL stop_bar_type ( barra, 'eri_vc_s1s1' )
  !
  ! eri_vc with DIFFERENT SPIN --> compute only s2>s1, obtain s2<s1 by symmetry
  !
  IF (nspin==2) THEN
     !
     CALL start_bar_type ( barra, 'eri_vc_s1s2', bandpair%nloc )
     !
     s1 = 1
     s2 = 2
     !
     DO p1loc = 1, bandpair%nloc
        !
        p1 = bandpair%l2g(p1loc)
        i = pijmap(1,p1)
        j = pijmap(2,p1)
        !
        CALL double_invfft_gamma(dffts,npwq,npwqx,proj_c(:,i,s1),proj_c(:,j,s2),psic,'Wave')
        !
        DO ir = 1, dffts%nnr
           rho_r(ir) = REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,ngm,ngm,rho_r,rho_g1,'Rho')
        !
        DO ig = 1, ngm
           rho_g1(ig) = pot3D%sqvc(ig) * rho_g1(ig)
        ENDDO
        !
        ! (s1 < s2)
        !
        DO p2 = 1, n_pairs
           !
           k = pijmap(1,p2)
           l = pijmap(2,p2)
           !
           CALL double_invfft_gamma(dffts,npwq,npwqx,proj_c(:,k,s1),proj_c(:,l,s2),psic,'Wave')
           !
           DO ir = 1, dffts%nnr
              rho_r(ir) = REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
           ENDDO
           !
           CALL single_fwfft_gamma(dffts,ngm,ngm,rho_r,rho_g2,'Rho')
           !
           DO ig = 1, ngm
              rho_g2(ig) = pot3D%sqvc(ig) * rho_g2(ig)
           ENDDO
           !
           eri_vc(p1,p2,s1,s2) = eri_vc(p1,p2,s1,s2) &
           & + 2._DP * DDOT(2*(ngm-gstart+1),rho_g1(gstart:ngm),1,rho_g2(gstart:ngm),1)/omega
           IF( i==j .AND. k==l .AND. gstart==2 ) THEN
              eri_vc(p1,p2,s1,s2) = eri_vc(p1,p2,s1,s2) + pot3D%div
           ENDIF
           !
           ! Apply symmetry
           !
           eri_vc(p2,p1,s2,s1) = eri_vc(p1,p2,s1,s2)
           !
        ENDDO
        !
        CALL update_bar_type ( barra, 'eri_vc_s1s2', 1 )
        !
     ENDDO
     !
     CALL stop_bar_type ( barra, 'eri_vc_s1s2' )
     !
  ENDIF
  !
  CALL mp_sum( eri_vc, intra_bgrp_comm )
  CALL mp_sum( eri_vc, inter_image_comm )
  !
  DEALLOCATE( rho_r )
  DEALLOCATE( rho_g1 )
  DEALLOCATE( rho_g2 )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_eri_wp(braket, chi_head, chi_body, eri_wp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE distribution_center,  ONLY : pert,macropert
  USE pwcom,                ONLY : nspin
  USE westcom,              ONLY : n_pdep_eigen_to_use,fftdriver,n_pairs,pijmap
  USE types_coulomb,        ONLY : pot3D
  USE mp_global,            ONLY : inter_image_comm
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE mp,                   ONLY : mp_sum
  USE io_push,              ONLY : io_push_title
  USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(IN) :: braket(n_pairs,nspin,n_pdep_eigen_to_use)
  COMPLEX(DP),INTENT(IN) :: chi_head, chi_body(pert%nglob,pert%nloc)
  COMPLEX(DP),INTENT(OUT) :: eri_wp(n_pairs,n_pairs,nspin,nspin)
  !
  INTEGER :: s1, s2, p1, p2, i, j, k, l, m, n, nloc, nloc_max
  TYPE(bar_type) :: barra
  !
  CALL io_push_title('Wp Integrals (PDEP)')
  !
  DO nloc = 1, macropert%nloc
     n = macropert%l2g(nloc)
     IF (n > n_pdep_eigen_to_use) EXIT
     nloc_max = nloc
  ENDDO
  !
  ! Workload too small for progress bar
  !
  CALL start_bar_type(barra, 'Wp', 1)
  !
  CALL pot3D%init(fftdriver,.FALSE.,'default')
  !
  eri_wp(:,:,:,:) = 0._DP
  !
  DO nloc = 1, nloc_max ! iterate over m, n
     !
     n = macropert%l2g(nloc)
     !
     DO m = 1, n_pdep_eigen_to_use
        DO s2 = 1, nspin
           DO s1 = 1, nspin
              DO p2 = 1, n_pairs
                 !
                 k = pijmap(1,p2)
                 l = pijmap(2,p2)
                 !
                 DO p1 = 1, n_pairs
                    !
                    i = pijmap(1,p1)
                    j = pijmap(2,p1)
                    !
                    eri_wp(p1,p2,s1,s2) = eri_wp(p1,p2,s1,s2) &
                    & + braket(p1,s1,m)*chi_body(m,nloc)*braket(p2,s2,n)/omega
                    IF( i==j .AND. k==l ) THEN
                       eri_wp(p1,p2,s1,s2) = eri_wp(p1,p2,s1,s2) + chi_head*pot3D%div
                    ENDIF
                    !
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO ! iterate over m, n
  !
  CALL mp_sum(eri_wp, inter_image_comm)
  !
  CALL update_bar_type(barra, 'Wp', 1)
  CALL stop_bar_type(barra, 'Wp')
  !
END SUBROUTINE
