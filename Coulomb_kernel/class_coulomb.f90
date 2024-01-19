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
! Marco Govoni
!
!-----------------------------------------------------------------------
MODULE class_coulomb
!-----------------------------------------------------------------------
   !
   USE kinds,            ONLY : DP
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   TYPE, PUBLIC :: coulomb
      !
      REAL(DP) :: div                              ! divergence
      CHARACTER(LEN=7) :: singularity_removal_mode ! singularity_removal_mode
      CHARACTER(LEN=5) :: cdriver                  ! FFT driver = 'Wave', 'Rho'
      LOGICAL :: l_use_igq                         ! use igq map
      INTEGER :: iq                                ! q-point
      REAL(DP), ALLOCATABLE :: sqvc(:)             ! square root of Coulomb potential in PW
      REAL(DP) :: mya, myb, mymu                   ! mya and mya + myb denote the factor of the short-range and of the
                                                   ! long-range parts of Fock exchange, respectively
                                                   ! mymu is inverse screening length in bohr^{-1}
                                                   ! See Eq. (3) of Phys. Rev. Mater. 2, 073803 (2018) for details
                                                   ! For HSE functional, mya = 1, myb = -1, mymu = 0.106
                                                   ! For PBE functional, mya = 1, myb = 0, mymu = c
      !
      CONTAINS
      !
      PROCEDURE :: init => sqvc_init
      PROCEDURE :: compute_divergence => compute_divergence
      PROCEDURE :: print_divergence => print_divergence
      !
   END TYPE coulomb
   !
   CONTAINS
   !
   !-----------------------------------------------------------------------
   SUBROUTINE sqvc_init(this,cdriver,l_use_igq,singularity_removal_mode,iq,mya,myb,mymu)
      !-----------------------------------------------------------------------
      !
      ! This routine computes results of applying sqVc to a potential
      ! associated with vector q. Coulomb cutoff technique can be used
      !
      USE kinds,                ONLY : DP
      USE constants,            ONLY : fpi,e2,eps8
      USE cell_base,            ONLY : at,tpiba2
      USE gvect,                ONLY : g,ngm
      USE westcom,              ONLY : igq_q,npwqx,npwq
      USE types_bz_grid,        ONLY : q_grid
      USE control_flags,        ONLY : gamma_only
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(coulomb) :: this
      CHARACTER(LEN=*), INTENT(IN) :: cdriver
      LOGICAL, INTENT(IN) :: l_use_igq
      CHARACTER(LEN=*), INTENT(IN) :: singularity_removal_mode
      INTEGER, INTENT(IN), OPTIONAL :: iq
      REAL(DP), INTENT(IN), OPTIONAL :: mya,myb,mymu
      !
      ! Workspace
      !
      LOGICAL :: l_abmu
      REAL(DP) :: qgnorm2, qg(3), x, einv
      INTEGER :: numg, numgx
      INTEGER :: ig, ipol
      LOGICAL :: on_double_grid
      REAL(DP) :: grid_factor
      !
      CALL start_clock('sqvc_init')
      !
      this%cdriver = cdriver
      this%l_use_igq = l_use_igq
      IF ( PRESENT(iq) ) THEN
         this%iq = iq
      ELSE
         this%iq = 1   ! gamma-only
      ENDIF
      !
      IF ( PRESENT(mya) .AND. PRESENT(myb) .AND. PRESENT(mymu) ) THEN
         this%mya = mya
         this%myb = myb
         this%mymu = mymu
         l_abmu = .TRUE.
      ELSEIF ( PRESENT(mya) .OR. PRESENT(myb) .OR. PRESENT(mymu) ) THEN
         CALL errore('sqvc_init','mya, myb, mymu must be specified together',1)
      ELSE
         l_abmu = .FALSE.
      ENDIF
      !
      ! ... Check compatibility between singularity removal mode and FFT driver
      !
      SELECT CASE ( singularity_removal_mode )
      CASE ( 'vcut_spherical', 'default' )
         this%singularity_removal_mode = 'default'
      CASE ( 'gygi-baldereschi', 'gygi-bald', 'g-b', 'gb' )
         this%singularity_removal_mode = 'gb'
         IF ( this%cdriver == 'Wave' ) CALL errore('sqvc_init','gb singularity removal mode requires Rho grid',1)
      CASE DEFAULT
         CALL errore('sqvc_init','singularity removal mode not supported, supported only default and gb',1)
      END SELECT
      !
      ! ... Check compatibility between FFT driver and use of igq map
      !
      SELECT CASE ( this%cdriver )
      CASE ( 'Wave' )
         IF (.NOT.gamma_only .AND. .NOT.this%l_use_igq) THEN
            CALL errore('sqvc_init','q-points case requires igq map when using Wave grid',1)
         ELSEIF (gamma_only .AND. this%l_use_igq) THEN
            CALL errore('sqvc_init','igq map not needed in gamma-only case',1)
         ENDIF
         numg = npwq
         numgx = npwqx
      CASE ( 'Rho' )
         IF (this%l_use_igq) CALL errore('sqvc_init','igq map not used with Rho grid',1)
         numg = ngm
         numgx = ngm
      CASE DEFAULT
         CALL errore('sqvc_init','cdriver value not supported, supported only Wave and Rho',1)
      END SELECT
      !
      IF( ALLOCATED(this%sqvc) )  DEALLOCATE( this%sqvc )
      ALLOCATE( this%sqvc( numgx ) )
      !
      this%sqvc = 0._DP
      DO ig = 1,numg
         !
         IF ( this%l_use_igq ) THEN
            qg(:) = g(:,igq_q(ig,this%iq)) + q_grid%p_cart(:,this%iq)
         ELSE
            qg(:) = g(:,ig) + q_grid%p_cart(:,this%iq)
         ENDIF
         !
         qgnorm2 = SUM( qg(:)**2 ) * tpiba2
         !
         IF( qgnorm2 < eps8 ) CYCLE ! don't touch sqvc_tmp of |q+G|=0
         !
         grid_factor = 1._DP
         on_double_grid = .FALSE.
         !
         IF( this%singularity_removal_mode == 'gb' ) THEN
            !
            ! In this case we use Gygi-Baldereschi method
            !
            grid_factor = 8._DP/7._DP
            on_double_grid = .TRUE.
            DO ipol = 1,3
               x = 0.5_DP*( qg(1)*at(1,ipol)+qg(2)*at(2,ipol)+qg(3)*at(3,ipol) )*REAL(q_grid%ngrid(ipol),KIND=DP)
               on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps8)
            ENDDO
            !
         ENDIF
         !
         IF( on_double_grid ) CYCLE
         !
         IF( l_abmu ) THEN
            einv = mya + myb * EXP(- qgnorm2 / 4._DP / mymu**2)
            this%sqvc(ig) = SQRT(e2*fpi*grid_factor/qgnorm2 * einv)
         ELSE
            this%sqvc(ig) = SQRT(e2*fpi*grid_factor/qgnorm2)
         ENDIF
         !
      ENDDO
      !
      IF ( q_grid%l_pIsGamma(this%iq) ) THEN
         IF ( l_abmu ) THEN
            CALL this%compute_divergence(singularity_removal_mode,mya,myb,mymu)
         ELSE
            CALL this%compute_divergence()
         ENDIF
      ENDIF
      !
      CALL stop_clock('sqvc_init')
      !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE compute_divergence(this,singularity_removal_mode,mya,myb,mymu)
      !-----------------------------------------------------------------------
      !
      USE constants,            ONLY : pi,tpi,fpi,e2,eps8
      USE cell_base,            ONLY : omega,at,bg,tpiba2
      USE mp,                   ONLY : mp_sum
      USE mp_global,            ONLY : nimage,my_image_id,nproc_image,inter_image_comm,&
                                     & intra_bgrp_comm,me_bgrp,nproc_bgrp
      USE control_flags,        ONLY : gamma_only
      USE gvecw,                ONLY : ecutwfc
      USE random_numbers,       ONLY : randy
      USE gvect,                ONLY : g,ngm
      USE types_bz_grid,        ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(coulomb) :: this
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: singularity_removal_mode
      REAL(DP), INTENT(IN), OPTIONAL :: mya, myb, mymu
      !
      ! Workspace
      !
      LOGICAL :: l_abmu
      CHARACTER(LEN=7) :: singularity_removal_mode_use
      REAL(DP) :: div
      LOGICAL :: i_am_ort, on_double_grid
      REAL(DP) :: qg(3), qgnorm2, alpha, peso
      INTEGER :: i1, i2, i3, iq, ig, ipol
      REAL(DP) :: prod(3,3), qhelp, edge(3), qbz(3), rand, qmo, vbz, vhelp, intcounter, x
      REAL(DP) :: q_, qq, aa, dq1
      LOGICAL, PARAMETER :: try_ort_div = .TRUE.
      INTEGER, PARAMETER :: nqq = 100000
      REAL(DP), PARAMETER :: grid_factor = 8._DP/7._DP
      !
      IF ( PRESENT(singularity_removal_mode) ) THEN
         singularity_removal_mode_use = singularity_removal_mode
      ELSE
         singularity_removal_mode_use = this%singularity_removal_mode
      ENDIF
      !
      IF ( PRESENT(mya) .AND. PRESENT(myb) .AND. PRESENT(mymu) ) THEN
         l_abmu = .TRUE.
      ELSEIF ( PRESENT(mya) .OR. PRESENT(myb) .OR. PRESENT(mymu) ) THEN
         CALL errore('sqvc_init','mya, myb, mymu must be specified together',1)
      ELSE
         l_abmu = .FALSE.
      ENDIF
      !
      div = 0._DP
      !
      SELECT CASE( singularity_removal_mode_use )
      !
      CASE( 'vcut_spherical', 'default' )
         !
         ! In this case we use the spherical region
         !
         div = ( (6._DP * pi * pi / ( omega*REAL(q_grid%np,KIND=DP) ) )**(1._DP/3._DP) ) / ( 2._DP * pi * pi ) * fpi * e2
         !
         ! If the angles are all 90 deg then overwrite div as follows:
         !
         IF( try_ort_div ) THEN
            !
            ! prod( i, j ) = (b_i)^t * b_j    (if off-diagonal prods are all zero --> the angles are all 90 deg --> cell is orthorombic)
            !
            prod = 0._DP
            DO i1 = 1, 3
               DO i2 = 1, 3
                  DO i3 = 1, 3
                     prod(i1,i2) = prod(i1,i2) + bg(i3,i1) * bg(i3,i2) * tpiba2
                  ENDDO
               ENDDO
            ENDDO
            !
            ! check if the off-diagonal prods are zero
            !
            i_am_ort = .TRUE.
            DO i1 = 1, 3
               DO i2 = 1, 3
                  IF( i1 == i2 ) CYCLE
                  IF( ABS( prod(i1,i2)) > eps8 ) THEN
                     i_am_ort = .FALSE.
                  ENDIF
               ENDDO
            ENDDO
            !
            ! if the system is not ort, spherical div will remain
            !
            IF ( i_am_ort ) THEN
               !
               edge(1) = SQRT(prod(1,1)) / 2._DP
               edge(2) = SQRT(prod(2,2)) / 2._DP
               edge(3) = SQRT(prod(3,3)) / 2._DP
               edge(:) = edge(:) / REAL(q_grid%ngrid(:),KIND=DP)
               !
               qhelp = MIN( edge(1),edge(2),edge(3)  )
               vbz = tpi**3 / ( omega * REAL(q_grid%np,KIND=DP) )
               vhelp = fpi / 3._DP * qhelp**3
               !
               rand = randy(my_image_id*nproc_image+me_bgrp)
               div = 0._DP
               intcounter = 0
               !
               DO i1 = 1, nqq
                  qmo=0._DP
                  DO i2 = 1, 3
                     qbz(i2) = randy() * edge(i2)
                     qmo = qmo + qbz(i2)**2
                  ENDDO
                  qmo = SQRT( qmo )
                  IF( qmo < qhelp ) CYCLE
                  div = div + 1._DP/qmo/qmo
                  intcounter = intcounter + 1._DP
               ENDDO
               !
               div = div * ( vbz - vhelp  ) / intcounter
               div = div + fpi * qhelp
               div = div * fpi * e2 / ( tpi * tpi * tpi )
               !
               div = div / REAL(nimage*nproc_bgrp,KIND=DP)
               !
               ! Cannot sum over world_comm, inter_pool_comm, or inter_bgrp_comm,
               ! because not all pools or band groups may enter this routine
               !
               CALL mp_sum(div,intra_bgrp_comm)
               CALL mp_sum(div,inter_image_comm)
               !
            ENDIF
            !
         ENDIF
         !
         IF ( l_abmu ) div = div * (mya + myb) - e2 * pi * myb/(mymu**2)
         !
      CASE( 'gygi-baldereschi', 'gygi-bald', 'g-b', 'gb' )
         !
         ! In this case we use Gygi-Baldereschi method
         !
         alpha = 10._DP / ecutwfc  ! DEFINITION OF ALPHA
         !
         div = 0._DP
         !
         DO iq = 1, q_grid%np
            !
            DO ig = 1,ngm
               qg(:) = q_grid%p_cart(:,iq) + g(:,ig)
               qgnorm2 = SUM( qg(:)**2 ) * tpiba2
               on_double_grid = .TRUE.
               DO ipol = 1,3
                  x = 0.5_DP*( qg(1)*at(1,ipol)+qg(2)*at(2,ipol)+qg(3)*at(3,ipol) )*REAL(q_grid%ngrid(ipol),KIND=DP)
                  on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps8)
               ENDDO
               IF( .NOT.on_double_grid .AND. qgnorm2 > eps8 ) THEN
                  !
                  ! with range separation
                  ! need to be fixed for mya != myb
                  !
                  IF ( l_abmu ) THEN
                     IF ( ABS(myb) > 0._DP ) THEN
                        div = div - EXP( -alpha * qgnorm2 ) / qgnorm2 &
                        & * (1._DP - EXP(-qgnorm2/4._DP/mymu**2))
                     ELSE
                        div = div - EXP( -alpha * qgnorm2 ) / qgnorm2
                     ENDIF
                  ELSE
                     div = div - EXP( -alpha * qgnorm2 ) / qgnorm2
                  ENDIF
               ENDIF
            ENDDO
            !
         ENDDO
         !
         CALL mp_sum( div, intra_bgrp_comm )
         !
         IF( gamma_only ) THEN
            peso = 2._DP
         ELSE
            peso = 1._DP
         ENDIF
         !
         IF( l_abmu ) THEN
            div = div * grid_factor * e2 * fpi / (omega * REAL(q_grid%np,KIND=DP)) * peso
            !
            dq1 = 5._DP / SQRT(alpha) / nqq
            aa = 0._DP
            IF ( ABS(myb) > 0._DP ) THEN
               DO iq = 0, nqq
                  q_ = dq1 * (iq+0.5_DP)
                  qq = q_**2
                  aa = aa - EXP(-alpha * qq) * EXP(-qq/4._DP/mymu**2) * dq1
               ENDDO
            ENDIF
            aa = aa * 8._DP / fpi
            aa = aa + 1._DP / SQRT(alpha * pi)
            !
            div = div + e2 * aa
         ELSE
            div = div * grid_factor * e2 * fpi / (omega * REAL(q_grid%np,KIND=DP)) * peso + e2 / SQRT( alpha * pi )
         ENDIF
         !
      END SELECT
      !
      this%div = div
      !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE print_divergence( this )
      !-----------------------------------------------------------------------
      !
      USE io_global,            ONLY : stdout
      USE types_bz_grid,        ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(coulomb) :: this
      !
      IF ( .NOT. q_grid%l_pIsGamma(this%iq) ) RETURN
      WRITE(stdout,"(5X,'Divergence = ',es14.6)") this%div
      !
   END SUBROUTINE
   !
END MODULE
