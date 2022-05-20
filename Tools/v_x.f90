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
SUBROUTINE v_x( rho, rho_core, rhog_core, etx, vtx, v )
  !----------------------------------------------------------------------------
  !! Exchange potential Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE noncollin_module, ONLY : domag
  USE funct,            ONLY : nlc, dft_is_nonlocc
  USE scf,              ONLY : scf_type
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin)
  !! V_xc potential
  REAL(DP), INTENT(OUT) :: vtx
  !! integral V_xc * rho
  REAL(DP), INTENT(OUT) :: etx
  !! E_xc energy
  !
  ! ... local variables
  !
  REAL(DP) :: rhoneg(2), vs
  !
  REAL(DP) :: arho, amag
  REAL(DP) :: rhoup2, rhodw2
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: vx(:,:), vc(:,:)
  ! In order:
    ! the absolute value of the total charge
    ! the absolute value of the magnetization
    ! zeta = amag / arhox
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir
    ! counter on mesh points
    ! counter on nspin
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.E-10_DP, &
                         vanishing_mag    = 1.E-20_DP
  !
  !
  CALL start_clock( 'v_x' )
  !
  ALLOCATE( ex(dfftp%nnr) )
  ALLOCATE( ec(dfftp%nnr) )
  ALLOCATE( vx(dfftp%nnr,nspin) )
  ALLOCATE( vc(dfftp%nnr,nspin) )
  !
  etx   = 0._DP
  vtx   = 0._DP
  v(:,:) = 0._DP
  rhoneg = 0._DP
  !
  !
  rho%of_r(:,1) = rho%of_r(:,1) + rho_core(:)
  !
  IF ( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
     ! ... spin-unpolarized case
     !
     CALL xc( dfftp%nnr, 1, 1, rho%of_r(:,1), ex, ec, vx(:,1), vc(:,1) )
     !
     DO ir = 1, dfftp%nnr
        v(ir,1) = e2*( vx(ir,1) )
        etx = etx + e2*( ex(ir) )*rho%of_r(ir,1)
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        vtx = vtx + v(ir,1)*rho%of_r(ir,1)
        IF (rho%of_r(ir,1) < 0._DP) rhoneg(1) = rhoneg(1)-rho%of_r(ir,1)
     ENDDO
     !
     !
  ELSEIF ( nspin == 2 ) THEN
     ! ... spin-polarized case
     !
     CALL xc( dfftp%nnr, 2, 2, rho%of_r, ex, ec, vx, vc )
     !
     DO ir = 1, dfftp%nnr   !OMP ?
        v(ir,:) = e2*( vx(ir,:) )
        etx = etx + e2*( (ex(ir) )*rho%of_r(ir,1) )
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        vtx = vtx + ( ( v(ir,1) + v(ir,2) )*rho%of_r(ir,1) + &
                        ( v(ir,1) - v(ir,2) )*rho%of_r(ir,2) )
        !
        rhoup2 = rho%of_r(ir,1)+rho%of_r(ir,2)
        rhodw2 = rho%of_r(ir,1)-rho%of_r(ir,2)
        IF (rhoup2 < 0._DP) rhoneg(1) = rhoneg(1) + rhoup2
        IF (rhodw2 < 0._DP) rhoneg(2) = rhoneg(2) + rhodw2
     ENDDO
     !
     vtx   = 0.5_DP * vtx
     rhoneg = 0.5_DP * rhoneg
     !
     !
  ELSE IF ( nspin == 4 ) THEN
     ! ... noncolinear case
     !
     CALL xc( dfftp%nnr, 4, 2, rho%of_r, ex, ec, vx, vc )
     !
     DO ir = 1, dfftp%nnr  !OMP ?
        arho = ABS( rho%of_r(ir,1) )
        IF ( arho < vanishing_charge ) CYCLE
        vs = 0.5_DP*( vx(ir,1) - vx(ir,2)  )
        v(ir,1) = e2*( 0.5_DP*( vx(ir,1) + vx(ir,2) ) )
        !
        amag = SQRT( SUM( rho%of_r(ir,2:4)**2 ) )
        IF ( amag > vanishing_mag ) THEN
           v(ir,2:4) = e2 * vs * rho%of_r(ir,2:4) / amag
           vtx = vtx + SUM( v(ir,2:4) * rho%of_r(ir,2:4) )
        ENDIF
        etx = etx + e2*( ex(ir) ) * arho
        !
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        IF ( rho%of_r(ir,1) < 0._DP )  rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        IF ( amag / arho > 1._DP )  rhoneg(2) = rhoneg(2) + 1._DP/omega
        vtx = vtx + v(ir,1) * rho%of_r(ir,1)
     ENDDO
     !
     !
  ENDIF
  !
  DEALLOCATE( ex, ec )
  DEALLOCATE( vx, vc )
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
  CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etx, vtx, v )
  !
  ! ... add non local corrections (if any)
  !
  IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etx, vtx, v )
  !
  CALL mp_sum(  vtx , intra_bgrp_comm )
  CALL mp_sum(  etx , intra_bgrp_comm )
  !
  CALL stop_clock( 'v_x' )
  !
END SUBROUTINE
