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
SUBROUTINE commutator_Hx_psi (ik, m, ipol, psi, dpsi, l_skip_nlpp)
  !----------------------------------------------------------------------
  !
  ! On input : psi(m-bands)  = | psi_ik >
  ! On output: dpsi(m-bands) = | dpsi_ik > = [H,x_ipol] | psi_ik > in crystal axis 
  !            (projected on at(*,ipol) )
  !
  ! vkb, must be properly set for the appropriate k-point
  ! NB: here the last index of becp1 is missing, hence it refers 
  !     to a single k-point
  !
  !    CALL calbec (npw,  vkb, psi, becp1(:,:) ) 
  !    CALL calbec (npw, work, psi, becp2(:,:) )
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba, at
  USE ions_base,       ONLY : nat, ityp, ntyp => nsp
  USE io_global,       ONLY : stdout
  USE klist,           ONLY : xk
  USE gvect,           ONLY : g
  USE wvfct,           ONLY : npw, npwx, nbnd, g2kin, et
  USE lsda_mod,        ONLY : nspin
  USE noncollin_module,ONLY : noncolin, npol
  USE becmod,          ONLY : becp, bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE uspp,            ONLY : nkb, vkb
  USE uspp_param,      ONLY : nh, nhm
  USE control_flags,   ONLY : gamma_only
  USE pwcom,           ONLY : igk_k,current_k
  ! 
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: ik ! k-point index
  INTEGER,INTENT(IN) :: m ! number of bands to process
  INTEGER,INTENT(IN) :: ipol ! polarization index
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,m) ! input wavefunctions
  COMPLEX(DP), INTENT(OUT)    :: dpsi(npwx*npol,m) ! output wavefunctions
  LOGICAL, INTENT(IN) :: l_skip_nlpp ! skip NLPP
  !
  ! Workspace
  !
  TYPE(bec_type) :: becp1 ! dimensions ( nkb, m )
  TYPE(bec_type) :: becp2 ! dimensions ( nkb, m )
  REAL(DP),ALLOCATABLE :: gk(:,:)  
  INTEGER :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0, is, js, ijs
  COMPLEX(DP),ALLOCATABLE :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:), work (:,:), psc(:,:,:,:), deff_nc(:,:,:,:)
  REAL(DP),ALLOCATABLE :: deff(:,:,:)
  !
  CALL start_clock ('commutator_Hx_psi')
  !
  ALLOCATE( gk(3,npwx) )    
  !
  dpsi = 0._DP
!$OMP PARALLEL DEFAULT(none) SHARED(npw,gk,xk,ik,g,igk_k,tpiba,g2kin,current_k) PRIVATE(ig)
!$OMP DO
  DO ig = 1, npw
     gk (1:3, ig) = (xk (1:3, ik) + g (1:3, igk_k (ig,current_k) ) ) * tpiba
     g2kin (ig) = SUM(gk (1:3, ig) **2 )
  ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
  !
  ! this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  IF( noncolin ) THEN 
!$OMP PARALLEL DEFAULT(none) SHARED(m,npw,dpsi,at,ipol,gk,psi,npwx) PRIVATE(ibnd,ig)
!$OMP DO 
     DO ibnd = 1, m
        DO ig = 1, npw ! first spin-component
           dpsi(ig,ibnd) = SUM( at(1:3,ipol ) * gk(1:3,ig) ) * (0._DP,-2._DP) * psi(ig,ibnd)
        ENDDO
        DO ig = 1, npw ! second spin-component
           dpsi(ig+npwx,ibnd) = SUM( at(1:3,ipol) * gk(1:3,ig) ) * (0._DP,-2._DP) * psi(ig+npwx,ibnd)
        ENDDO
     ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ELSE
!$OMP PARALLEL DEFAULT(none) SHARED(m,npw,dpsi,at,ipol,gk,psi) PRIVATE(ibnd,ig)
!$OMP DO COLLAPSE(2)
     DO ibnd = 1, m
        DO ig = 1, npw 
           dpsi(ig,ibnd) = SUM(at(1:3,ipol)*gk(1:3,ig))*(0._DP,-2._DP)*psi(ig,ibnd)
        ENDDO
     ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF
  !
  ! Uncomment this goto and the continue below to calculate 
  ! the matrix elements of p without the commutator with the
  ! nonlocal potential.
  !
  IF( l_skip_nlpp ) GOTO 111
  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
  IF( nkb == 0 ) GOTO 111
  !
  ALLOCATE( work(npwx,nkb) )
  IF( noncolin ) THEN
     ALLOCATE( deff_nc(nhm, nhm, nat, nspin) )
  ELSE
     ALLOCATE( deff(nhm, nhm, nat) )
  ENDIF
  ALLOCATE( dvkb(npwx,nkb), dvkb1(npwx, nkb) )
  !
  dvkb   = 0._DP
  dvkb1  = 0._DP
  work   = 0._DP
  ! 
  CALL gen_us_dj (ik, dvkb)
  CALL gen_us_dy (ik, at (1, ipol), dvkb1)
  !
!$OMP PARALLEL DEFAULT(none) SHARED(npw,g2kin,gk) PRIVATE(ig)
!$OMP DO 
  DO ig = 1, npw
     IF (g2kin (ig) < 1.0d-10) THEN
        gk (1, ig) = 0._DP
        gk (2, ig) = 0._DP
        gk (3, ig) = 0._DP
     ELSE
        gk (1, ig) = gk (1, ig) / DSQRT (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / DSQRT (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / DSQRT (g2kin (ig) )
     ENDIF
  ENDDO
!$OMP END DO
!$OMP END PARALLEL 
  !
  jkb = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (nt == ityp (na)) THEN
           DO ikb = 1, nh (nt)
              jkb = jkb + 1
              DO ig = 1, npw
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                      (at (1, ipol) * gk (1, ig) + &
                       at (2, ipol) * gk (2, ig) + &
                       at (3, ipol) * gk (3, ig) )
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE( gk )
  DEALLOCATE( dvkb )
  DEALLOCATE( dvkb1 )

  !
  ! In the case of gamma point systems becp2 is real
  ! so we have to include a factor of i before calling
  ! calbec otherwise we would be stuck with the wrong component
  ! of becp2 later on.
  IF (gamma_only) work=(0.0_DP,1.0_DP)*work
  !
  CALL ALLOCATE_BEC_TYPE ( nkb, m, becp1 )
  CALL ALLOCATE_BEC_TYPE ( nkb, m, becp2 )
  !
  CALL calbec (npw,  vkb, psi, becp1, m)
  CALL calbec (npw, work, psi, becp2, m)
  !
  IF( noncolin ) THEN
     ALLOCATE( psc (nkb,npol,nbnd,2) )
     psc=0._DP
  ELSE
     ALLOCATE( ps2 (nkb,nbnd,2) )
     ps2=0._DP
  END IF
  !
  DO ibnd = 1, m
     IF( noncolin ) THEN
        CALL compute_deff_nc(deff_nc,et(ibnd,ik))
     ELSE
        CALL compute_deff(deff,et(ibnd,ik))
     ENDIF
     ijkb0 = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (nt == ityp (na)) THEN
              DO ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    IF( noncolin ) THEN
                       ijs=0
                       DO is=1, npol
                          DO js = 1, npol
                             ijs=ijs+1
                             psc(ikb,is,ibnd,1)=psc(ikb,is,ibnd,1)+  &
                                       (0._DP,-1._DP)*    &
                                  becp2%nc(jkb,js,ibnd)*deff_nc(ih,jh,na,ijs) 
                             psc(ikb,is,ibnd,2)=psc(ikb,is,ibnd,2)+ &
                                     (0._DP,-1._DP)* &
                                 becp1%nc(jkb,js,ibnd)*deff_nc(ih,jh,na,ijs) 
                          END DO
                       END DO
                    ELSEIF( gamma_only ) THEN
                       ! Note the different prefactors due to the factor 
                       ! of i introduced to work(:,:), as becp[1,2] are
                       ! real.
                       ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1) + becp2%r(jkb,ibnd) * &
                            (1._DP, 0._DP)*deff(ih,jh,na) 
                       ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) + becp1%r(jkb,ibnd) * &
                            (-1._DP, 0._DP)*deff(ih,jh,na)
                    ELSE
                       ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1) + becp2%k(jkb,ibnd) * &
                            (0._DP,-1._DP)*deff(ih,jh,na) 
                       ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) + becp1%k(jkb,ibnd) * &
                            (0._DP,-1._DP)*deff(ih,jh,na)
                    END IF
                 ENDDO
              ENDDO
              ijkb0=ijkb0+nh(nt)
           ENDIF
        ENDDO  ! na
     ENDDO  ! nt
  ENDDO ! nbnd
  IF (ikb /= nkb .OR. jkb /= nkb) call errore ('commutator_Hx_psi', 'unexpected error',1)
  !
  CALL deallocate_bec_type( becp1 )
  CALL deallocate_bec_type( becp2 )
  !
  IF( noncolin ) THEN
     CALL ZGEMM( 'N', 'N', npw, m*npol, nkb, &
          (1._DP,0._DP), vkb(1,1), npwx, psc(1,1,1,1), nkb, (1._DP,0._DP), &
          dpsi, npwx )
     CALL ZGEMM( 'N', 'N', npw, m*npol, nkb, &
          (1._DP,0._DP),work(1,1), npwx, psc(1,1,1,2), nkb, (1._DP,0._DP), &
          dpsi, npwx )
  ELSE
     CALL ZGEMM( 'N', 'N', npw, m, nkb, &
          (1._DP,0._DP), vkb(1,1), npwx, ps2(1,1,1), nkb, (1._DP,0._DP), &
          dpsi(1,1), npwx )
     CALL ZGEMM( 'N', 'N', npw, m, nkb, &
          (1._DP,0._DP),work(1,1), npwx, ps2(1,1,2), nkb, (1._DP,0._DP), &
          dpsi(1,1), npwx )
  ENDIF
  !
  IF( noncolin ) THEN
     DEALLOCATE( psc )
     DEALLOCATE( deff_nc )
  ELSE
     DEALLOCATE( ps2 )
     DEALLOCATE( deff )
  END IF
  DEALLOCATE( work )
  !
  111 CONTINUE
  !
  CALL stop_clock ('commutator_Hx_psi')
  RETURN
END SUBROUTINE commutator_Hx_psi
