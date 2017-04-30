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
SUBROUTINE gradcorr_x( rho, rhog, rho_core, rhog_core, etx, vtx, v )
  !-----------------------------------------------------------------------
  !
  ! [NO CORRELATION!!!!!!!]
  ! Modified from PW/src/gradcorr.f90 to exclude correlation 
  ! Ideally this routine should be committed directly in QE.
  !
  USE constants,            ONLY : e2
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nl, ngm, g
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : omega, alat
  USE funct,                ONLY : gcxc, gcx_spin, gcc_spin, igcc_is_lyp, &
                                   gcc_spin_more, dft_is_gradient, get_igcc
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : ux
  USE wavefunctions_module, ONLY : psic
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: rho(dfftp%nnr,nspin), rho_core(dfftp%nnr)
  COMPLEX(DP), INTENT(IN)    :: rhog(ngm,nspin), rhog_core(ngm)
  REAL(DP),    INTENT(INOUT) :: v(dfftp%nnr,nspin)
  REAL(DP),    INTENT(INOUT) :: vtx, etx
  !
  INTEGER :: k, ipol, is, nspin0, ir, jpol
  !
  REAL(DP),    ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP),    ALLOCATABLE :: rhoout(:,:), segni(:), vgg(:,:), vsave(:,:)
  REAL(DP),    ALLOCATABLE :: gmag(:,:,:)

  COMPLEX(DP), ALLOCATABLE :: rhogsum(:,:)
  !
  REAL(DP) :: grho2(2), sx, sc, v1x, v2x, v1c, v2c, &
              v1xup, v1xdw, v2xup, v2xdw, &
              etxgc, vtxgc, segno, arho, fac, rh, amag 
  REAL(DP) :: rup, rdw, grup, grdw, seg
  !
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  !
  !
  IF ( .NOT. dft_is_gradient() ) RETURN
  !
  etxgc = 0.D0
  vtxgc = 0.D0
  !
  nspin0=nspin
  if (nspin==4) nspin0=1
  if (nspin==4.and.domag) nspin0=2
  fac = 1.D0 / DBLE( nspin0 )
  !
  ALLOCATE(    h( 3, dfftp%nnr, nspin0) )
  ALLOCATE( grho( 3, dfftp%nnr, nspin0) )
  ALLOCATE( rhoout( dfftp%nnr, nspin0) )
  IF (nspin==4.AND.domag) THEN
     ALLOCATE( vgg( dfftp%nnr, nspin0 ) )
     ALLOCATE( vsave( dfftp%nnr, nspin ) )
     ALLOCATE( segni( dfftp%nnr ) )
     vsave=v
     v=0.d0
  ENDIF
  !
  ALLOCATE( rhogsum( ngm, nspin0 ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF ( nspin == 4 .AND. domag ) THEN
     !
     CALL compute_rho(rho,rhoout,segni,dfftp%nnr)
     !
     ! ... bring starting rhoout to G-space
     !
     DO is = 1, nspin0
        !
        psic(:) = rhoout(:,is)
        !
        CALL fwfft ('Dense', psic, dfftp)
        !
        rhogsum(:,is) = psic(nl(:))
        !
     END DO
  ELSE
     !
     rhoout(:,1:nspin0)  = rho(:,1:nspin0)
     rhogsum(:,1:nspin0) = rhog(:,1:nspin0)
     !
  ENDIF
  DO is = 1, nspin0
     !
     rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)
     rhogsum(:,is) = fac * rhog_core(:) + rhogsum(:,is)
     !
     CALL gradrho( dfftp%nnr, rhogsum(1,is), ngm, g, nl, grho(1,1,is) )
     !
  END DO
  !
  DEALLOCATE( rhogsum )
  !
  IF ( nspin0 == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     !
     DO k = 1, dfftp%nnr
        !
        arho = ABS( rhoout(k,1) )
        !
        IF ( arho > epsr ) THEN
           !
           grho2(1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
           !
           IF ( grho2(1) > epsg ) THEN
              !
              segno = SIGN( 1.D0, rhoout(k,1) )
              !
              CALL gcxc( arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c )
              !
              ! ... first term of the gradient correction : D(rho*Ex)/D(rho)
              !
              v(k,1) = v(k,1) + e2 * v1x 
              !
              ! ... h contains :
              !
              ! ...    D(rho*Ex) / D(|grad rho|) * (grad rho) / |grad rho|
              !
              h(:,k,1) = e2 * v2x * grho(:,k,1)
              !
              vtxgc = vtxgc+e2* v1x * ( rhoout(k,1) - rho_core(k) )
              etxgc = etxgc+e2* sx * segno
              !
           ELSE
              h(:,k,1)=0.D0
           END IF
           !
        ELSE
           !
           h(:,k,1) = 0.D0
           !
        END IF
        !
     END DO
     !
  ELSE
     !
     ! ... spin-polarised case
     !
!$omp parallel do private( rh, grho2, sx, v1xup, v1xdw, v2xup, v2xdw, rup, rdw, grup, grdw  ), &
!$omp             reduction(+:etxgc,vtxgc)
     DO k = 1, dfftp%nnr
        !
        rh = rhoout(k,1) + rhoout(k,2)
        !
        grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
        !
        CALL gcx_spin( rhoout(k,1), rhoout(k,2), grho2(1), &
                       grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
        !
        ! ... first term of the gradient correction : D(rho*Ex)/D(rho)
        !
        v(k,1) = v(k,1) + e2 * v1xup
        v(k,2) = v(k,2) + e2 * v1xdw
        !
        ! ... h contains D(rho*Ex)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        DO ipol = 1, 3
           !
           grup = grho(ipol,k,1)
           grdw = grho(ipol,k,2)
           h(ipol,k,1) = e2 * ( v2xup  * grup )
           h(ipol,k,2) = e2 * ( v2xdw  * grdw )
           !
        END DO
        !
        vtxgc = vtxgc + &
                 e2 * v1xup * ( rhoout(k,1) - rho_core(k) * fac )
        vtxgc = vtxgc + &
                 e2 * v1xdw * ( rhoout(k,2) - rho_core(k) * fac )
        etxgc = etxgc + e2 * sx 
        !
     END DO
!$omp end parallel do
     !
  END IF
  !
  DO is = 1, nspin0
     !
     rhoout(:,is) = rhoout(:,is) - fac * rho_core(:)
     !
  END DO
  !
  DEALLOCATE( grho )
  !
  ALLOCATE( dh( dfftp%nnr ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin0
     !
     CALL grad_dot( dfftp%nnr, h(1,1,is), ngm, g, nl, alat, dh )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     vtxgc = vtxgc - SUM( dh(:) * rhoout(:,is) )
     !
  END DO
  !
  vtx = vtx + omega * vtxgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  etx = etx + omega * etxgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )

  IF (nspin==4.AND.domag) THEN
     DO is=1,nspin0
        vgg(:,is)=v(:,is)
     ENDDO
     v=vsave
     DO k=1,dfftp%nnr
        v(k,1)=v(k,1)+0.5d0*(vgg(k,1)+vgg(k,2))
        amag=sqrt(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
        IF (amag.GT.1.d-12) THEN
           v(k,2)=v(k,2)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,2)/amag
           v(k,3)=v(k,3)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,3)/amag
           v(k,4)=v(k,4)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,4)/amag
        ENDIF
     ENDDO
  ENDIF
  !
  DEALLOCATE( dh )
  DEALLOCATE( h )
  DEALLOCATE( rhoout )
  IF (nspin==4.and.domag) THEN
     DEALLOCATE( vgg )
     DEALLOCATE( vsave )
     DEALLOCATE( segni )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
