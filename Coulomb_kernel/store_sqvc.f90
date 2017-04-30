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
SUBROUTINE store_sqvc(sqvc_tmp,numg,singularity_removal_mode,div)
 !-----------------------------------------------------------------------
 !
 ! This routine computes results of applying sqVc to a potential
 ! associated with vector q. Coulomb cutoff technique can be used
 !
 USE kinds,                ONLY : DP
 USE io_global,            ONLY : stdout
 USE constants,            ONLY : fpi, e2, tpi, pi
 USE cell_base,            ONLY : tpiba2, tpiba,omega,at,alat,bg
 USE gvect,                ONLY : g, gstart
 USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                  vcut_get,  vcut_spheric_get, vcut_destroy
 USE mp,                   ONLY : mp_sum
 USE mp_global,            ONLY : intra_bgrp_comm
 USE mp_world,             ONLY : mpime,world_comm,nproc
 USE gvecw,                ONLY : ecutwfc
 USE control_flags,        ONLY : gamma_only
 USE random_numbers,       ONLY : randy
 !
 IMPLICIT NONE
 !
 ! I/O
 !
 INTEGER,INTENT(IN) :: numg
 REAL(DP),INTENT(OUT) :: sqvc_tmp(numg)
 INTEGER,INTENT(IN) :: singularity_removal_mode
 REAL(DP),INTENT(OUT) :: div
 !
 ! Workspace
 !
 REAL(DP) :: gnorm2,q(3),x,nq(3),ecutvcut,atws(3,3),dq(3),xq(3),alpha
 INTEGER :: ig, ipol, iq1, iq2, iq3, i1, i2, i3
 LOGICAL :: on_double_grid
 REAL(DP),PARAMETER :: eps=1.d-6
 REAL(DP) :: grid_factor = 8.d0/7.d0
 TYPE(vcut_type)   :: vcut
 LOGICAL :: try_ort_div=.TRUE.,i_am_ort 
 REAL(DP) :: prod(3,3), qhelp, edge(3), qbz(3), rand, qmo, vbz, vhelp 
 REAL(DP) :: intcounter
 !
 CALL start_clock('storesqvc')
 !
 ! ======
 !  BODY
 ! ======
 !
 SELECT CASE(singularity_removal_mode)
    !
 CASE(1)
    !
    ! In this case we use the spherical region 
    !
    sqvc_tmp(1)=0._DP
    DO ig = gstart,numg
       gnorm2 = ( g(1,ig)*g(1,ig) + g(2,ig)*g(2,ig) + g(3,ig)*g(3,ig) ) * tpiba2
       sqvc_tmp(ig) = SQRT(e2*fpi/gnorm2)
    ENDDO
    !
 CASE(2)
    !
    ! In this case we use Gygi-Baldereschi method
    !
    nq(1)=1._DP
    nq(2)=1._DP
    nq(3)=1._DP
    !
    sqvc_tmp(1)=0._DP
    DO ig = gstart,numg
       q(:) = g(:,ig)
       gnorm2 = SUM(q(:)**2) * tpiba2
       on_double_grid = .TRUE.
       DO ipol = 1,3
          x = 0.5_DP*( q(1)*at(1,ipol)+q(2)*at(2,ipol)+q(3)*at(3,ipol) )*nq(ipol)
          on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
       ENDDO
       IF(on_double_grid) THEN
          sqvc_tmp(ig)=0._DP
       ELSE
          sqvc_tmp(ig)=SQRT(e2*fpi*grid_factor/gnorm2)
       ENDIF
    ENDDO
    !
 CASE(3)
    !
    ! In this case we use CUT_WS
    ! 
    !
    nq(1)=1._DP
    nq(2)=1._DP
    nq(3)=1._DP
    ecutvcut = 0.7_DP
    !
    ! build the superperiodicity direct lattice
    !
    atws = alat * at
    !
    DO ipol=1,3
       atws(:,ipol) = atws(:,ipol) * nq(ipol)
    ENDDO
    !
    CALL vcut_init( vcut, atws, ecutvcut )
    !
    DO ig = 1,numg
       !
       q(:) = g(:,ig) * tpiba
       !
       sqvc_tmp( ig ) = DSQRT( vcut_get(vcut,q) )
       !
    ENDDO
    !
    CALL vcut_destroy(vcut)
    !
 END SELECT
 !
 ! ======
 !  HEAD
 ! ======
 !
 SELECT CASE(singularity_removal_mode)
    !
 CASE(1)
    !
    ! In this case we use the spherical region 
    !
    div = ( (6._DP * pi * pi / omega )**(1._DP/3._DP) ) / ( 2._DP * pi * pi ) * fpi * e2
    WRITE(stdout,"(5X,'Spherical div        = ',es14.6)") div 
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
             IF( ABS( prod(i1,i2)) > 1.d-8 ) THEN 
                WRITE(stdout,"(5X,'Warning: non orthorombic RL')") 
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
          !
          qhelp = MIN( edge(1),edge(2),edge(3)  )  
          vbz = tpi**3 / omega
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
          WRITE(stdout,"(5X,'Orthorombic div      = ',es14.6)") div 
          !
       ENDIF
       !
    ENDIF
    !
 CASE(2)
    !
    ! In this case we use Gygi-Baldereschi method
    !
    alpha = 10._DP / ecutwfc
    !
    nq(1)=1._DP
    nq(2)=1._DP
    nq(3)=1._DP
    !
    dq(1)=1._DP/nq(1)
    dq(2)=1._DP/nq(2)
    dq(3)=1._DP/nq(3)
    !
    div = 0._DP
    !
    DO iq1=1,1,INT(nq(1))
       DO iq2=1,1,INT(nq(2))
          DO iq3=1,1,INT(nq(3))
             xq(:) = bg(:,1)*(iq1-1)*dq(1) + bg(:,2)*(iq2-1)*dq(2) + bg(:,3)*(iq3-1)*dq(3)
             !
             DO ig = gstart,numg
                q(:) = g(:,ig) + xq(:)
                gnorm2 = SUM(q(:)**2) * tpiba2
                on_double_grid = .TRUE.
                DO ipol = 1,3
                   x = 0.5_DP*( q(1)*at(1,ipol)+q(2)*at(2,ipol)+q(3)*at(3,ipol) )*nq(ipol)
                   on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
                ENDDO
                IF(.NOT.on_double_grid) THEN 
                   div = div - EXP( -alpha * gnorm2 ) / gnorm2
                ENDIF
             ENDDO
             !
          ENDDO
       ENDDO
    ENDDO 
    !
    CALL mp_sum( div, intra_bgrp_comm )
    !
    IF( gamma_only ) THEN
       div = div * grid_factor * e2 * fpi / omega * 2._DP
    ELSE
       div = div * grid_factor * e2 * fpi / omega
    ENDIF
    ! 
    div = div + nq(1)*nq(2)*nq(3) * e2 / SQRT( alpha * pi )
    !
 CASE(3)
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
!------------------------------------------------------------------
SUBROUTINE store_vcspecial_H(vc_tmp,numg)
 !
 ! vcspecail_H
 !
 USE kinds,     ONLY : DP
 USE cell_base, ONLY : tpiba, at,omega
 USE pwcom,     ONLY : npw,npwx
 USE gvect,     ONLY : g,gstart
 USE constants, ONLY : fpi, e2, pi
 USE exx,       ONLY : exxdiv 
 !
 IMPLICIT NONE
 !
 ! I/O
 !
 INTEGER,INTENT(IN) :: numg
 REAL(DP),INTENT(OUT) :: vc_tmp(numg)
 ! 
 !Local variables
 INTEGER :: ig !Counters 
 REAL(DP) :: q(3), qq, x
 REAL(DP) :: grid_factor = 8.d0/7.d0
 REAL(DP) :: eps = 1.d-6 
 LOGICAL :: on_double_grid
 REAL(DP) :: kost
 !
 kost = e2*fpi
 !
 IF(gstart==2) vc_tmp(1) = -exxdiv 
 !
!$OMP PARALLEL SHARED(numg,g,tpiba,at,vc_tmp,kost,gstart) PRIVATE(ig,q,qq,on_double_grid,x)
!$OMP DO
 DO ig=gstart,numg
     !
     q(:)= g(:,ig)
     !
     q = q * tpiba
     !
     qq = SUM(q(:)**2) 
     !
     on_double_grid = .TRUE.
     x= 0.5_DP/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))!*nq1
     on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
     x= 0.5_DP/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))!*nq2
     on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
     x= 0.5_DP/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))!*nq3
     on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
     ! 
     IF( on_double_grid ) THEN
        !
        vc_tmp(ig) = -kost / qq
        !
     ELSE
        !
        vc_tmp(ig)= kost/qq / 7._DP
        !
     ENDIF
     !
  ENDDO
!$OMP ENDDO
!$OMP END PARALLEL 
  !
END SUBROUTINE 
!
!
!
SUBROUTINE store_sqvc_sphcut(sqvc_tmp,numg,rcut)
 !
 ! This routine computes results of applying sqVc to a potential
 ! associated with vector q. Coulomb cutoff technique can be used
 !
 USE kinds,            ONLY : DP
 USE io_global,        ONLY : stdout
 USE constants,        ONLY : fpi, e2, tpi
 USE cell_base,        ONLY : tpiba2, tpiba
 USE gvect,            ONLY : g, gstart
 !
 IMPLICIT NONE
 !
 ! I/O
 !
 INTEGER,INTENT(IN) :: numg
 REAL(DP),INTENT(OUT) :: sqvc_tmp(numg)
 REAL(DP),INTENT(IN) :: rcut
 !
 ! Workspace
 !
 REAL(DP) :: gnorm2,gnorm
 INTEGER :: ig
 !
 !
 CALL start_clock('storesqvc')
 !
 !
 IF(gstart==2) sqvc_tmp(1)=SQRT(tpi*e2*rcut*rcut)
 !
!$OMP PARALLEL private(ig,gnorm2,gnorm)
!$OMP DO
 DO ig = gstart,numg
    gnorm2 = ( g(1,ig)*g(1,ig) + g(2,ig)*g(2,ig) + g(3,ig)*g(3,ig) ) * tpiba2
    gnorm  = SQRT(gnorm2)
    sqvc_tmp(ig) = SQRT(e2*fpi/gnorm2 * (1._DP-COS(gnorm*rcut)))
 ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
 !
 CALL stop_clock('storesqvc')
 !
END SUBROUTINE
