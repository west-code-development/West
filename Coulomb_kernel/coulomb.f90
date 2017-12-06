!
! Copyright (C) 2015-2017 M. Govoni 
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
MODULE coulomb
   !--------------------------------------------------------------------
   !
   USE kinds,                ONLY : DP
   USE constants,            ONLY : pi, tpi, fpi, e2, eps8
   USE cell_base,            ONLY : alat, omega, at, bg, tpiba, tpiba2
   USE gvect,                ONLY : g
   USE westcom,              ONLY : igq_q
   USE class_bz_grid,        ONLY : bz_grid
   USE types_bz_grid,        ONLY : q_grid
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   PUBLIC :: store_sqvc, compute_divergence, store_sqvc_sphcut
   !
   CONTAINS
   !
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE store_sqvc(sqvc_tmp,numg,singularity_removal_mode,iq,l_use_igq,div,l_printout_div)
      !-----------------------------------------------------------------------
      !
      ! This routine computes results of applying sqVc to a potential
      ! associated with vector q. Coulomb cutoff technique can be used
      !
      !USE kinds,                ONLY : DP
      !USE constants,            ONLY : fpi, e2, pi, eps8
      !USE cell_base,            ONLY : tpiba2,tpiba,at,alat
      !USE gvect,                ONLY : g
      USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                       vcut_get,  vcut_spheric_get, vcut_destroy
      !USE westcom,              ONLY : igq_q
      !USE class_bz_grid,        ONLY : bz_grid
      !USE types_bz_grid,        ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: numg
      REAL(DP),INTENT(OUT) :: sqvc_tmp(numg)
      CHARACTER(LEN=*), INTENT(IN) :: singularity_removal_mode
      INTEGER, INTENT(IN) :: iq
      LOGICAL, INTENT(IN) :: l_use_igq
      REAL(DP), INTENT(OUT) :: div
      LOGICAL,INTENT(IN), OPTIONAL :: l_printout_div
      !
      ! Workspace
      !
      REAL(DP) :: qgnorm2,qg(3),x,ecutvcut,atws(3,3)
      INTEGER :: ig, ipol
      LOGICAL :: on_double_grid, l_print
      REAL(DP) :: grid_factor = 8.d0/7.d0
      TYPE(vcut_type) :: vcut
      !
      CALL start_clock('storesqvc')
      !
      ! ======
      !  BODY
      ! ======
      !
      SELECT CASE(singularity_removal_mode)
         !
      CASE("spherical")
         !
         ! In this case we use the spherical region 
         !
         sqvc_tmp = 0._DP
         DO ig = 1,numg
            IF ( l_use_igq ) THEN
               qg(:) = q_grid%p_cart(:,iq) + g(:,igq_q(ig,iq)) 
            ELSE
               qg(:) = q_grid%p_cart(:,iq) + g(:,ig)
            ENDIF
            qgnorm2 = SUM( qg(:)**2 ) * tpiba2
            IF ( qgnorm2 > eps8 ) sqvc_tmp(ig) = SQRT(e2*fpi/qgnorm2)
         ENDDO
         !
      CASE("gb")
         !
         ! In this case we use Gygi-Baldereschi method
         !
         sqvc_tmp = 0._DP
         DO ig = 1,numg
            IF ( l_use_igq ) THEN
               qg(:) = q_grid%p_cart(:,iq) + g(:,igq_q(ig,iq)) 
            ELSE
               qg(:) = q_grid%p_cart(:,iq) + g(:,ig)
            ENDIF
            qgnorm2 = SUM( qg(:)**2 ) * tpiba2
            on_double_grid = .TRUE.
            DO ipol = 1,3
               x = 0.5_DP*( qg(1)*at(1,ipol)+qg(2)*at(2,ipol)+qg(3)*at(3,ipol) )*DBLE(q_grid%ngrid(ipol))
               on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps8)
            ENDDO
            IF( on_double_grid .OR. qgnorm2 < eps8 ) CYCLE
            sqvc_tmp(ig)=SQRT(e2*fpi*grid_factor/qgnorm2)
         ENDDO
         !
      CASE("ws")
         !
         ! In this case we use CUT_WS
         ! 
         ecutvcut = 0.7_DP
         !
         ! build the superperiodicity direct lattice
         !
         atws = alat * at
         !
         DO ipol=1,3
            atws(:,ipol) = atws(:,ipol) * DBLE(q_grid%ngrid(ipol))
         ENDDO
         !
         CALL vcut_init( vcut, atws, ecutvcut )
         !
         DO ig = 1,numg
            !
            IF ( l_use_igq ) THEN
               qg(:) = ( q_grid%p_cart(:,iq) + g(:,igq_q(ig,iq)) ) * tpiba
            ELSE
               qg(:) = ( q_grid%p_cart(:,iq) + g(:,ig) ) * tpiba
            ENDIF
            !
            sqvc_tmp( ig ) = DSQRT( vcut_get(vcut,qg) )
            !
         ENDDO
         !
         CALL vcut_destroy(vcut)
         !
      CASE DEFAULT 
         !
         CALL errore("store_sqvc", "singularity_removal_mode not supported, supported only spherical, gb, ws", 1 )
         !
      END SELECT
      !
      IF ( PRESENT(l_printout_div) ) THEN
         l_print = l_printout_div
      ELSE
         l_print = .FALSE.
      ENDIF
      !
      IF ( q_grid%l_pIsGamma(iq) ) CALL compute_divergence( numg, singularity_removal_mode, div, l_print )
      !
   END SUBROUTINE
   !
   !
   !
   SUBROUTINE compute_divergence( numg, singularity_removal_mode, div, l_printout_div )
      !    
      !USE kinds,                ONLY : DP
      !USE constants,            ONLY : pi, tpi, fpi, e2, eps8
      !USE cell_base,            ONLY : omega, at, bg, tpiba2
      !USE gvect,                ONLY : g
      USE io_global,            ONLY : stdout
      USE mp,                   ONLY : mp_sum
      USE mp_global,            ONLY : intra_bgrp_comm
      USE mp_world,             ONLY : mpime, world_comm, nproc
      USE control_flags,        ONLY : gamma_only
      USE gvecw,                ONLY : ecutwfc
      USE random_numbers,       ONLY : randy
      !USE class_bz_grid,        ONLY : bz_grid
      !USE types_bz_grid,        ONLY : q_grid
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: numg
      CHARACTER(LEN=*), INTENT(IN) :: singularity_removal_mode
      REAL(DP), INTENT(INOUT) :: div
      LOGICAL, INTENT(IN) :: l_printout_div
      !
      ! Workspace
      !
      LOGICAL :: try_ort_div=.TRUE., i_am_ort, on_double_grid
      REAL(DP) :: qg(3), qgnorm2, alpha, peso
      REAL(DP) :: grid_factor = 8.d0/7.d0
      INTEGER :: i1, i2, i3, iq, ig, ipol
      REAL(DP) :: prod(3,3), qhelp, edge(3), qbz(3), rand, qmo, vbz, vhelp, intcounter, x
      !
      SELECT CASE( singularity_removal_mode )
      !
      CASE("spherical")
         !
         ! In this case we use the spherical region 
         !
         div = ( (6._DP * pi * pi / ( omega*DBLE(q_grid%np) ) )**(1._DP/3._DP) ) / ( 2._DP * pi * pi ) * fpi * e2
         !
         IF ( l_printout_div ) WRITE(stdout,"(5X,'Spherical div        = ',es14.6)") div 
         !
         IF( try_ort_div ) THEN
            !  
            ! Redefine div according to orthorombic cell
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
            i_am_ort = .TRUE.
            DO i1 = 1, 3
               DO i2 = 1, 3
                  IF( i1 == i2 ) CYCLE
                  IF( ABS( prod(i1,i2)) > eps8 ) THEN 
                     IF ( l_printout_div ) WRITE(stdout,"(5X,'Warning: non orthorombic RL')") 
                     i_am_ort = .FALSE.
                  ENDIF
               ENDDO
            ENDDO 
            !
            IF ( i_am_ort ) THEN
               !
               edge(1) = SQRT(prod(1,1)) / 2._DP 
               edge(2) = SQRT(prod(2,2)) / 2._DP 
               edge(3) = SQRT(prod(3,3)) / 2._DP 
               edge(:) = edge(:) / DBLE(q_grid%ngrid(:))
               !
               qhelp = MIN( edge(1),edge(2),edge(3)  )  
               vbz = tpi**3 / ( omega * DBLE(q_grid%np) )
               vhelp = fpi / 3._DP * qhelp**3
               !
               rand=randy(mpime)
               div = 0._DP
               intcounter = 0 
               !
               DO i1 = 1, 100000
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
               div = div / REAL(nproc,KIND=DP)
               CALL mp_sum(div,world_comm)
               !
               IF ( l_printout_div ) WRITE(stdout,"(5X,'Orthorombic div      = ',es14.6)") div 
               !
            ENDIF
            !
         ENDIF
         !
      CASE("gb")
         !
         ! In this case we use Gygi-Baldereschi method
         !
         alpha = 10._DP / ecutwfc
         !
         div = 0._DP
         !
         DO iq = 1, q_grid%np
            ! 
            DO ig = 1,numg
               qg(:) = q_grid%p_cart(:,iq) + g(:,ig)
               qgnorm2 = SUM( qg(:)**2 ) * tpiba2
               on_double_grid = .TRUE.
               DO ipol = 1,3
                  x = 0.5_DP*( qg(1)*at(1,ipol)+qg(2)*at(2,ipol)+qg(3)*at(3,ipol) )*DBLE(q_grid%ngrid(ipol))
                  on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps8)
               ENDDO
               IF( .NOT.on_double_grid .AND. qgnorm2 > eps8 ) THEN 
                  div = div - EXP( -alpha * qgnorm2 ) / qgnorm2
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
         div = div * grid_factor * e2 * fpi / (omega * DBLE(q_grid%np)) * peso
         ! 
         div = div + e2 / SQRT( alpha * pi )
         !
      CASE("ws")
         !
         ! In this case we use CUT_WS
         ! 
         div = 0._DP
         !
      END SELECT
      !
      CALL stop_clock('storesqvc')
      !
   END SUBROUTINE
   !
   !
   !
   SUBROUTINE store_sqvc_sphcut(sqvc_tmp,numg,iq,l_use_igq,rcut)
      !
      ! This routine computes results of applying sqVc to a potential
      ! associated with vector q. Coulomb cutoff technique can be used
      !
      !USE kinds,            ONLY : DP
      !USE constants,        ONLY : fpi, e2, tpi, eps8
      !USE cell_base,        ONLY : tpiba2
      !USE gvect,            ONLY : g !, gstart
      !USE westcom,          ONLY : igq_q
      !USE class_bz_grid,    ONLY : bz_grid
      !USE types_bz_grid,    ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: numg
      REAL(DP),INTENT(OUT) :: sqvc_tmp(numg)
      INTEGER, INTENT(IN) :: iq
      LOGICAL, INTENT(IN) :: l_use_igq
      REAL(DP),INTENT(IN) :: rcut
      !
      ! Workspace
      !
      REAL(DP) :: qg(3), qgnorm2, qgnorm
      INTEGER :: ig
      !
      !
      CALL start_clock('storesqvc')
      !
      !
      !IF(gstart==2) sqvc_tmp(1)=SQRT(tpi*e2*rcut*rcut)
      !
!!$OMP PARALLEL private(ig,qgnorm2,qgnorm,qg)
!!$OMP DO
      sqvc_tmp = 0._DP
      DO ig = 1,numg
         IF ( l_use_igq ) THEN
            qg(:) = q_grid%p_cart(:,iq) + g(:,igq_q(ig,iq))
         ELSE
            qg(:) = q_grid%p_cart(:,iq) + g(:,ig)
         ENDIF
         qgnorm2 = ( SUM(qg(:))**2 ) * tpiba2
         IF ( qgnorm2 < eps8 ) THEN
            sqvc_tmp(ig)=SQRT(tpi*e2*rcut*rcut)
         ELSE
            qgnorm  = SQRT(qgnorm2)
            sqvc_tmp(ig) = SQRT(e2*fpi/qgnorm2 * (1._DP-COS(qgnorm*rcut)))
         ENDIF
      ENDDO
!!$OMP ENDDO
!!$OMP END PARALLEL
      !
      CALL stop_clock('storesqvc')
      !
   END SUBROUTINE
   !
   !
   !
   !!---------------------------------------------------------------------------------
   !SUBROUTINE store_sqvc_q(sqvc_tmp,numg,singularity_removal_mode,iq,l_use_igq)
   !  !-------------------------------------------------------------------------------
   !  !
   !  ! This routine computes results of applying sqVc to a potential
   !  ! associated with vector q.
   !  !
   !  USE kinds,                ONLY : DP
   !  USE io_global,            ONLY : stdout
   !  USE constants,            ONLY : fpi, e2, tpi, pi
   !  USE cell_base,            ONLY : tpiba2,tpiba,omega,at,alat,bg
   !  USE gvect,                ONLY : g, gstart
   !  USE mp,                   ONLY : mp_sum
   !  USE mp_global,            ONLY : intra_bgrp_comm
   !  USE mp_world,             ONLY : mpime,world_comm,nproc
   !  USE gvecw,                ONLY : ecutwfc
   !  USE gvect,                ONLY : g, gstart
   !  USE westcom,              ONLY : igq_q
   !  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
   !                                  vcut_get,  vcut_spheric_get, vcut_destroy
   !  USE random_numbers,       ONLY : randy
   !  USE class_bz_grid,        ONLY : bz_grid
   !  USE types_bz_grid,        ONLY : q_grid
   !  USE constants,            ONLY : eps8
   !  !
   !  IMPLICIT NONE
   !  !
   !  ! I/O
   !  !
   !  INTEGER, INTENT(IN) :: numg
   !  REAL(DP), INTENT(OUT) :: sqvc_tmp(numg)
   !  INTEGER, INTENT(IN) :: singularity_removal_mode
   !  INTEGER, INTENT(IN) :: iq
   !  LOGICAL, INTENT(IN) :: l_use_igq
   !  !
   !  ! Workspace
   !  !
   !  REAL(DP) :: gnorm2,nq(3),q(3),qq(3),ecutvcut,atws(3,3),alpha,x
   !  INTEGER :: ig, ipol, i1, i2, i3
   !  LOGICAL :: on_double_grid
   !  REAL(DP) :: grid_factor = 8.d0/7.d0
   !! REAL(DP) :: grid_factor = 1.0_DP
   !  TYPE(vcut_type)   :: vcut
   !  LOGICAL :: try_ort_div=.TRUE.,i_am_ort 
   !  REAL(DP) :: prod(3,3), qhelp, edge(3), qbz(3), rand, qmo, vbz, vhelp 
   !  REAL(DP) :: intcounter
   !  !
   !  CALL start_clock( 'storesqvcq' )
   !  !
   !  nq(1) = REAL( q_grid%ngrid(1), KIND=DP )
   !  nq(2) = REAL( q_grid%ngrid(2), KIND=DP )
   !  nq(3) = REAL( q_grid%ngrid(3), KIND=DP )
   !  !
   !  q(:) = q_grid%p_cart(:,iq)
   !  !
   !  ! =======
   !  !  BODY
   !  ! =======
   !  !
   !  SELECT CASE(singularity_removal_mode)
   !     !
   !  CASE(1)
   !     !
   !     ! In this case we use the spherical region 
   !     !
   !     DO ig = 1,numg
   !        IF (l_use_igq) THEN
   !           qq(:) = q(:) + g(:,igq_q(ig,iq))
   !        ELSE
   !           qq(:) = q(:) + g(:,ig)
   !        ENDIF
   !        gnorm2 = SUM(qq(:)**2) * tpiba2
   !        IF ( gnorm2 < eps8 ) THEN
   !           sqvc_tmp(ig) = 0._DP
   !        ELSE
   !           sqvc_tmp(ig) = SQRT(e2*fpi/gnorm2)
   !        ENDIF
   !     ENDDO
   !     !
   !  CASE(2)
   !     !
   !     ! In this case we use Gygi-Baldereschi method
   !     !
   !     DO ig = 1,numg
   !        IF (l_use_igq) THEN
   !           qq(:) = q(:) + g(:,igq_q(ig,iq))
   !        ELSE
   !           qq(:) = q(:) + g(:,ig)
   !        ENDIF
   !        gnorm2 = SUM(qq(:)**2) * tpiba2
   !        on_double_grid = .TRUE.
   !        DO ipol = 1,3
   !           x = 0.5_DP*( qq(1)*at(1,ipol)+qq(2)*at(2,ipol)+qq(3)*at(3,ipol) )*nq(ipol)
   !           on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps8)
   !        ENDDO
   !        IF( on_double_grid .OR. gnorm2 < eps8 ) THEN
   !           sqvc_tmp(ig)=0._DP
   !        ELSE
   !           sqvc_tmp(ig)=SQRT(e2*fpi*grid_factor/gnorm2)
   !        ENDIF
   !     ENDDO
   !     !
   !  CASE(3)
   !     !
   !     ! In this case we use CUT_WS
   !     ! 
   !     !
   !     ecutvcut = 0.7_DP
   !     !
   !     ! build the superperiodicity direct lattice
   !     !
   !     atws = alat * at
   !     !
   !     DO ipol=1,3
   !        atws(:,ipol) = atws(:,ipol) * nq(ipol)
   !     ENDDO
   !     !
   !     CALL vcut_init( vcut, atws, ecutvcut )
   !     !
   !     DO ig = 1,numg
   !        !
   !        IF (l_use_igq) THEN
   !           qq(:) = (q(:) + g(:,igq_q(ig,iq))) * tpiba
   !        ELSE
   !           qq(:) = (q(:) + g(:,ig)) * tpiba
   !        ENDIF
   !        !
   !        sqvc_tmp( ig ) = DSQRT( vcut_get(vcut,qq) )
   !        !
   !     ENDDO
   !     !
   !     CALL vcut_destroy(vcut)
   !     !
   !  END SELECT
   !  !
   !  CALL stop_clock( 'storesqvcq' )
   !  !
   !END SUBROUTINE
   !
   !
   !
END MODULE
