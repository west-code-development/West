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
SUBROUTINE solve_eri(ifreq,real_freq)
  !-----------------------------------------------------------------------
  !
  USE westcom,              ONLY : n_pdep_eigen_to_use,d_epsm1_ifr_a,z_epsm1_rfr_a,iuwfc,lrwfc,&
                                 & imfreq_list,refreq_list,ijpmap,pijmap,braket,eri,&
                                 & wfreq_save_dir,z_head_rfr_a,d_head_ifr_a,n_bands,n_pairs
  USE distribution_center,  ONLY : pert,macropert,ifr,rfr
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE mp_global,            ONLY : intra_bgrp_comm,me_bgrp,inter_image_comm,my_image_id
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title, io_push_bar
  USE wfreq_db,             ONLY : qdet_db_write
  !
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ifreq
  LOGICAL,INTENT(IN) :: real_freq
  !
  INTEGER  :: i, j, p1
  INTEGER  :: who, iloc, iunit
  COMPLEX(DP) :: freq, chi_head
  REAL(DP)    :: ry_to_ha = 0.5_DP
  REAL(DP), ALLOCATABLE :: bare_eri(:,:,:,:), screened_eri(:,:,:,:)
  !
  COMPLEX(DP),ALLOCATABLE  :: chi_body(:,:)
  !
  ! TODO: check if G=0 component is double counted in any dot products
  !
  ! Compute 4-center integrals of W (screened electron repulsion integrals, eri)
  ! (ij|W|kl) = int dx dx' f*_i(x) f_j(x) W(r,r') f*_k(x') f_l(x'), f are projectors in proj_r
  !
  CALL io_push_title("ERI")
  !
  ALLOCATE( chi_body( pert%nglob, pert%nloc ) )
  !
  ! Load chi at given frequency
  !
  freq = 0._DP
  chi_head = 0._DP
  chi_body = 0._DP
  !
  IF ( real_freq ) THEN
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
        chi_head = d_head_ifr_a(ifreq)
        chi_body(:,:) = d_epsm1_ifr_a(:,:,ifreq)
        !
     ENDIF
     !
  ENDIF
  !
  CALL mp_sum( freq, intra_bgrp_comm )
  CALL mp_sum( chi_head, intra_bgrp_comm )
  CALL mp_sum( chi_body, intra_bgrp_comm )
  !
  ! Compute ERI
  !
  ALLOCATE( braket(n_pairs,nspin,n_pdep_eigen_to_use) )
  !  
  CALL compute_braket()
  !
  ALLOCATE( screened_eri(n_pairs,n_pairs,nspin,nspin) )
  ALLOCATE( bare_eri(n_pairs,n_pairs,nspin,nspin) )
  ALLOCATE( eri(n_pairs,n_pairs,nspin,nspin) )
  !
  bare_eri = 0._DP
  eri = 0._DP
  !
  ! 4-center integrals of Vc
  !
  CALL compute_bv_direct(bare_eri)
  !
  CALL compute_hv(bare_eri)
  !
  ! 4-center integrals of Wp
  screened_eri = 0._DP
  !
  CALL compute_wp_pdep(chi_head, chi_body, screened_eri)
  !
  ! calculate total 4-center integrals
  eri = bare_eri + REAL(screened_eri)
  !
  CALL qdet_db_write()
  !
  CALL io_push_bar()
  !
  DEALLOCATE( eri, braket, chi_body )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_braket()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx
  USE westcom,              ONLY : wstat_save_dir,npwq,n_pdep_eigen_to_use,npwqx,fftdriver,iuwfc,lrwfc,&
                                 & n_bands,n_pairs,proj_r,braket,pijmap
  USE mp_global,            ONLY : intra_bgrp_comm,me_bgrp,inter_image_comm,my_image_id
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE types_coulomb,        ONLY : pot3D
  USE io_push,              ONLY : io_push_title
  USE distribution_center,  ONLY : macropert,kpt_pool
  USE pdep_db,              ONLY : generate_pdep_fname
  USE pdep_io,              ONLY : pdep_read_G_and_distribute 
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),ALLOCATABLE :: phi(:), rho_r(:), rho_g(:), psi1(:), psi2(:)
  !
  INTEGER :: s, m, p1, i, j, ig, ir, mloc
  INTEGER :: iunit
  CHARACTER(LEN=25) :: filepot
  CHARACTER(LEN=:),ALLOCATABLE :: fname
  TYPE(bar_type) :: barra
  !
  REAL(DP), EXTERNAL    :: DDOT
  !
  CALL io_push_title('braket')
  !
  ALLOCATE( phi(npwqx), rho_r(dffts%nnr), rho_g(npwqx) )
  !
  CALL pot3D%init(fftdriver,.FALSE.,"default")
  !
  CALL start_bar_type ( barra, 'braket', nspin * macropert%nloc )
  !
  braket = 0._DP
  !
  DO s = 1, nspin
     !
     DO mloc=1,macropert%nloc
        !
        m = macropert%l2g(mloc)
        !
        IF (m > n_pdep_eigen_to_use) CYCLE
        !
        CALL generate_pdep_fname( filepot, m ) 
        fname = TRIM( wstat_save_dir ) // "/"// filepot
        CALL pdep_read_G_and_distribute(fname,phi(:))
        !
        DO ig = 1, npwq
           phi(ig) = pot3D%sqvc(ig) * phi(ig)
        ENDDO
        !
        DO p1 = 1, n_pairs
           !
           i = pijmap(1,p1)
           j = pijmap(2,p1)
           !
           DO ir=1,dffts%nnr
              rho_r(ir)=proj_r(ir,i,s)*proj_r(ir,j,s)
           ENDDO
           !
           CALL single_fwfft_gamma(dffts,npwq,npwqx,rho_r,rho_g,TRIM(fftdriver))
           !
           braket(p1,s,m) = CMPLX(2.0_DP * DDOT(2*npwq,rho_g,1,phi,1), 0._DP)  ! assume Gamma-only case
           !
        ENDDO
        !
        CALL update_bar_type( barra,'braket', 1 )
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(braket, inter_image_comm)
  CALL mp_sum(braket, intra_bgrp_comm)
  !
  CALL stop_bar_type( barra,'braket' )
  !
  DEALLOCATE( phi, rho_r, rho_g )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_bv_direct(bare)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin,npw,npwx
  USE westcom,              ONLY : iuwfc,lrwfc,eri,n_pdep_eigen_to_use,n_bands,n_pairs,proj_r,&
                                 & pijmap
  USE mp_global,            ONLY : intra_bgrp_comm,me_bgrp,inter_image_comm,my_image_id
  USE distribution_center,  ONLY : bandpair,kpt_pool
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_fwfft_gamma
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE mp,                   ONLY : mp_sum
  USE types_coulomb,        ONLY : pot3D
  USE io_push,              ONLY : io_push_title
  USE gvect,                ONLY : ngm
  USE cell_base,            ONLY : omega
  USE class_idistribute,    ONLY : idistribute
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: bare(n_pairs,n_pairs,nspin,nspin)
  COMPLEX(DP),ALLOCATABLE  :: rho_gs(:,:,:), rho_g(:), rho_r(:), psi1(:), psi2(:)
  !
  COMPLEX(DP) :: Hv,Bv
  INTEGER     :: i, j, k, l, p1, p1loc, p2, s1, s2
  INTEGER     :: ir, ig
  TYPE(bar_type) :: barra
  !
  COMPLEX(DP), EXTERNAL :: ZDOTC
  REAL(DP), EXTERNAL    :: DDOT
  !
  CALL io_push_title('pair density')
  !
  ! Distribute pairs over images
  !
  CALL bandpair%init(n_pairs,'i','n_pairs',.FALSE.)
  !
  CALL pot3D%init('Rho',.FALSE.,"gb")
  !
  ALLOCATE( rho_gs(ngm,bandpair%nloc,nspin) )
  ALLOCATE( psi1(dffts%nnr), psi2(dffts%nnr), rho_r(dffts%nnr), rho_g(ngm) )
  !
  CALL start_bar_type ( barra, 'rho', nspin * bandpair%nloc )
  !
  rho_gs = 0._DP
  !
  ! Iterate over i, j, store results in rho_gs
  !
  DO s1 = 1, nspin
     !
     DO p1loc = 1, bandpair%nloc
        !
        p1 = bandpair%l2g(p1loc)
        i = pijmap(1,p1)
        j = pijmap(2,p1)
        !
        DO ir=1,dffts%nnr
           rho_r(ir)=proj_r(ir,i,s1)*proj_r(ir,j,s1)
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,ngm,ngm,rho_r,rho_gs(:,p1loc,s1),'Rho')
        !
        DO ig = 1, ngm
           rho_gs(ig,p1loc,s1) = pot3D%sqvc(ig)**2 * rho_gs(ig,p1loc,s1)
        ENDDO
        !
        CALL update_bar_type ( barra, 'rho', 1 )
        !
     ENDDO
  ENDDO  ! Iterate over i, j
  !
  CALL stop_bar_type ( barra, 'rho' )
  !
  CALL io_push_title('Vc Integrals (Direct)')
  CALL start_bar_type ( barra, 'Bv', nspin * n_pairs )
  !
  ! Iterate over k, l
  !
  DO s2 = 1, nspin
     !
     DO p2 = 1, n_pairs
        !
        k = pijmap(1,p2)
        l = pijmap(2,p2)
        !
        DO ir=1,dffts%nnr
           rho_r(ir)=proj_r(ir,k,s2)*proj_r(ir,l,s2)
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,ngm,ngm,rho_r,rho_g,'Rho')
        !
        DO s1 = 1, nspin
           DO p1loc = 1, bandpair%nloc
              !
              p1 = bandpair%l2g(p1loc)
              !
              Bv = (1/omega) * 2.0_DP * DDOT(2*ngm,rho_gs(:,p1loc,s1),1,rho_g,1)
              bare(p2,p1,s2,s1) = bare(p2,p1,s2,s1) + REAL(Bv)
           ENDDO
        ENDDO
        !
        CALL update_bar_type ( barra, 'Bv', 1 )
        !
     ENDDO
  ENDDO  ! Iterate over k, l
  !
  CALL mp_sum( bare, intra_bgrp_comm )
  CALL mp_sum( bare, inter_image_comm )
  !
  CALL stop_bar_type ( barra, 'Bv' )
  !
  DEALLOCATE( rho_gs, psi1, psi2, rho_r, rho_g )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_bv_pdep(bare)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : nspin
  USE westcom,              ONLY : eri,n_pdep_eigen_to_use,n_bands,n_pairs,fftdriver,braket
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE mp,                   ONLY : mp_sum
  USE types_coulomb,        ONLY : pot3D
  USE io_push,              ONLY : io_push_title
  USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: bare(n_pairs,n_pairs,nspin,nspin)
  COMPLEX(DP) :: Hv,Bv
  INTEGER     :: s1, s2, p1, p2
  TYPE(bar_type) :: barra
  !
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL io_push_title('Vc Integrals (PDEP)')
  !
  CALL pot3D%init(fftdriver,.FALSE.,"default")
  !
  CALL start_bar_type ( barra, 'Bv', nspin * nspin * n_pairs )
  !
  DO s1 = 1, nspin
     DO s2 = 1, nspin
        DO p1 = 1, n_pairs
           DO p2 = 1, n_pairs
              !
              Bv = (1/omega) * ZDOTC(n_pdep_eigen_to_use,braket(p2,s2,:),1,braket(p1,s1,:),1)
              !
              bare(p2,p1,s2,s1) = bare(p2,p1,s2,s1) + REAL(Bv)
              !
           ENDDO
           !
           CALL update_bar_type ( barra, 'Bv', 1 )
           !
        ENDDO
     ENDDO
  ENDDO
  !
  CALL stop_bar_type ( barra, 'Bv' )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_hv(bare)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : n_bands,eri,ijpmap,n_pairs
  USE pwcom,                ONLY : nspin
  USE types_coulomb,        ONLY : pot3D
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: bare(n_pairs,n_pairs,nspin,nspin)
  COMPLEX(DP) :: Hv
  INTEGER     :: s1, s2, p1, p2, i, k
  !
  Hv = pot3D%div
  !
  DO s1 = 1, nspin
     DO s2 = 1, nspin
        DO i = 1, n_bands
           DO k = 1, n_bands
              !
              p1 = ijpmap(i,i)
              p2 = ijpmap(k,k)
              !
              bare(p2,p1,s2,s1) = bare(p2,p1,s2,s1) + REAL(Hv)
              !
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_wp_pdep(chi_head, chi_body, screened)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE distribution_center,  ONLY : pert,macropert
  USE pwcom,                ONLY : nspin
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_bands,fftdriver,n_pairs,braket,ijpmap
  USE types_coulomb,        ONLY : pot3D
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE mp,                   ONLY : mp_sum
  USE io_push,              ONLY : io_push_title
  USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: chi_head, chi_body(pert%nglob,pert%nloc)
  COMPLEX(DP), INTENT(INOUT)  :: screened(n_pairs,n_pairs,nspin,nspin)
  !
  COMPLEX(DP) :: Hp,Bp
  INTEGER     :: s1, s2, p1, p2, i, k, m, nloc, n
  TYPE(bar_type) :: barra
  !
  CALL io_push_title('Wp Integrals (PDEP)')
  !
  CALL start_bar_type(barra, 'Wp', macropert%nloc)
  !
  ! Body of Wp (Bp)
  !
  DO nloc = 1, macropert%nloc  ! Iterate over m, n
     n = macropert%l2g(nloc)
     IF (n > n_pdep_eigen_to_use) CYCLE
     DO m = 1, macropert%nglob
        IF (m > n_pdep_eigen_to_use) CYCLE
        !
        DO s1 = 1, nspin
           DO s2 = 1, nspin
              DO p1 = 1, n_pairs
                 DO p2 = 1, n_pairs
                    !
                    Bp = (1/omega) * chi_body(m,nloc) * braket(p1,s1,m) * CONJG(braket(p2,s2,n))
                    screened(p2,p1,s2,s1) = screened(p2,p1,s2,s1) + Bp
                    !
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        !
     ENDDO
     !
     CALL update_bar_type ( barra, 'Wp', 1)
     !
  ENDDO ! Iterate over m, n
  !
  CALL mp_sum(screened, inter_image_comm)
  !
  CALL stop_bar_type ( barra, 'Wp' )
  !
  ! head of Wp (Hp)
  !
  CALL pot3D%init(fftdriver,.FALSE.,"default")
  !
  Hp = chi_head * pot3D%div
  !
  DO s1 = 1, nspin
     DO s2 = 1, nspin
        DO i = 1, n_bands
           DO k = 1, n_bands
              !
              p1 = ijpmap(i,i)
              p2 = ijpmap(k,k)
              !
              screened(p2,p1,s2,s1) = screened(p2,p1,s2,s1) + Hp
              !
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
