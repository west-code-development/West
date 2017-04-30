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
SUBROUTINE v_x( rho, rho_core, rhog_core, etx, vtx, v )
  !-----------------------------------------------------------------------
  !
  ! [NO CORRELATION!!!!!!!]
  ! Modified from PW/src/v_of_rho.f90 to exclude correlation 
  ! Ideally this routine should be committed directly in QE.
  !
  ! ... Exchange potential Vx(r) from n(r) [NO CORRELATION]
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : xc, xc_spin, nlc, dft_is_nonlocc
  USE scf,              ONLY : scf_type
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin), vtx, etx
    ! V_x potential
    ! integral V_x * rho
    ! E_x energy
  !
  ! ... local variables
  !
  REAL(DP) :: rhox, arhox, zeta, amag, vs, ex, ec, vx(2), vc(2), rhoneg(2)
    ! the total charge in each point
    ! the absolute value of the charge
    ! the absolute value of the charge
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir, ipol
    ! counter on mesh points
    ! counter on nspin
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag    = 1.D-20
  !
  !
  CALL start_clock( 'v_x' )
  !
  etx    = 0.D0
  vtx    = 0.D0
  v(:,:) = 0.D0
  rhoneg = 0.D0
  !
  IF ( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
     !
     ! ... spin-unpolarized case
     !
!$omp parallel do private( rhox, arhox, ex, ec, vx, vc ), &
!$omp             reduction(+:etx,vtx), reduction(-:rhoneg)
     DO ir = 1, dfftp%nnr
        !
        rhox = rho%of_r(ir,1) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           CALL xc( arhox, ex, ec, vx(1), vc(1) )
           !
           v(ir,1) = e2 * vx(1) 
           !
           etx = etx + e2 * ex * rhox
           !
           vtx = vtx + v(ir,1) * rho%of_r(ir,1)
           !
        ENDIF
        !
        IF ( rho%of_r(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        !
     END DO
!$omp end parallel do
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ! ... spin-polarized case
     !
!$omp parallel do private( rhox, arhox, zeta, ex, ec, vx, vc ), &
!$omp             reduction(+:etx,vtx), reduction(-:rhoneg)
     DO ir = 1, dfftp%nnr
        !
        rhox = rho%of_r(ir,1) + rho%of_r(ir,2) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = ( rho%of_r(ir,1) - rho%of_r(ir,2) ) / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) zeta = SIGN( 1.D0, zeta )
           !
           IF ( rho%of_r(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
           IF ( rho%of_r(ir,2) < 0.D0 ) rhoneg(2) = rhoneg(2) - rho%of_r(ir,2)
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           v(ir,:) = e2 * vx(:)
           !
           etx = etx + e2* ex * rhox
           !
           vtx = vtx + ( v(ir,1)*rho%of_r(ir,1) + v(ir,2)*rho%of_r(ir,2) )
           !
        END IF
        !
     END DO
!$omp end parallel do
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     ! ... noncolinear case
     !
     DO ir = 1,dfftp%nnr
        !
        amag = SQRT( rho%of_r(ir,2)**2 + rho%of_r(ir,3)**2 + rho%of_r(ir,4)**2 )
        !
        rhox = rho%of_r(ir,1) + rho_core(ir)
        !
        IF ( rho%of_r(ir,1) < 0.D0 )  rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = amag / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) THEN
              !
              rhoneg(2) = rhoneg(2) + 1.D0 / omega
              !
              zeta = SIGN( 1.D0, zeta )
              !
           END IF
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           vs = 0.5D0*( vx(1) - vx(2) )
           !
           v(ir,1) = e2*( 0.5D0*( vx(1) + vx(2) ) )
           !
           IF ( amag > vanishing_mag ) THEN
              !
              DO ipol = 2, 4
                 !
                 v(ir,ipol) = e2 * vs * rho%of_r(ir,ipol) / amag
                 !
                 vtx = vtx + v(ir,ipol) * rho%of_r(ir,ipol)
                 !
              END DO
              !
           END IF
           !
           etx = etx + e2*( ex ) * rhox
           vtx = vtx + v(ir,1) * rho%of_r(ir,1)
           !
        END IF
        !
     END DO
     !
  END IF
  !
  CALL mp_sum(  rhoneg , intra_bgrp_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2ES10.3)') rhoneg
  !
  ! ... energy terms, local-density contribution
  !
  vtx = omega * vtx / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  etx = omega * etx / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  ! ... add gradient corrections (if any)
  !
  CALL gradcorr_x( rho%of_r, rho%of_g, rho_core, rhog_core, etx, vtx, v )
  !
  CALL mp_sum(  vtx , intra_bgrp_comm )
  CALL mp_sum(  etx , intra_bgrp_comm )
  !
  CALL stop_clock( 'v_x' )
  !
  RETURN
  !
END SUBROUTINE
