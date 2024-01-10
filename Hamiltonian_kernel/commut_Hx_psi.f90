!
! Copyright (C) 2015-2023 M. Govoni
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
#if !defined(__CUDA)
!-----------------------------------------------------------------------
SUBROUTINE commut_Hx_psi(ik, m, ipol, psi, dpsi, l_skip_nlpp)
  !-----------------------------------------------------------------------
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
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : tpiba, at
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE klist,            ONLY : xk
  USE gvect,            ONLY : g
  USE wvfct,            ONLY : npw, npwx, g2kin, et
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : noncolin, npol
  USE becmod,           ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE uspp,             ONLY : nkb, vkb
  USE uspp_param,       ONLY : nh, nhm
  USE uspp_init,        ONLY : gen_us_dj, gen_us_dy
  USE control_flags,    ONLY : gamma_only
  USE pwcom,            ONLY : igk_k, current_k
  USE control_flags,    ONLY : offload_type
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: ik ! k-point index
  INTEGER,INTENT(IN) :: m ! number of bands to process
  INTEGER,INTENT(IN) :: ipol ! polarization index
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,m) ! input wavefunctions
  COMPLEX(DP), INTENT(OUT) :: dpsi(npwx*npol,m) ! output wavefunctions
  LOGICAL, INTENT(IN) :: l_skip_nlpp ! skip NLPP
  !
  ! Workspace
  !
  TYPE(bec_type) :: becp1 ! dimensions ( nkb, m )
  TYPE(bec_type) :: becp2 ! dimensions ( nkb, m )
  REAL(DP),ALLOCATABLE :: gk(:,:)
  INTEGER :: ig, na, ibnd, ikb, jkb, nt, ih, jh, ijkb0, is, js, ijs
  COMPLEX(DP),ALLOCATABLE :: ps2(:,:,:), dvkb (:,:), work (:,:), psc(:,:,:,:), deff_nc(:,:,:,:)
  REAL(DP),ALLOCATABLE :: deff(:,:,:)
  REAL(DP),PARAMETER :: tol = 1.E-10_DP
  !
  CALL start_clock ('commut_Hx_psi')
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
  ! this is the kinetic contribution to [H,x]: -2i (k+G)_ipol * psi
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
  IF( (.NOT. l_skip_nlpp) .AND. nkb > 0 ) THEN
     !
     ! this is the contribution from nonlocal pseudopotentials
     !
     ALLOCATE( work(npwx,nkb) )
     IF( noncolin ) THEN
        ALLOCATE( deff_nc(nhm, nhm, nat, nspin) )
     ELSE
        ALLOCATE( deff(nhm, nhm, nat) )
     ENDIF
     ALLOCATE( dvkb(npwx,nkb) )
     !
     work = 0._DP
     !
!$OMP PARALLEL DEFAULT(none) SHARED(npw,g2kin,gk) PRIVATE(ig)
!$OMP DO
     DO ig = 1, npw
        IF (g2kin (ig) < tol) THEN
           gk (1, ig) = 0._DP
           gk (2, ig) = 0._DP
           gk (3, ig) = 0._DP
        ELSE
           gk (1, ig) = gk (1, ig) / SQRT (g2kin (ig) )
           gk (2, ig) = gk (2, ig) / SQRT (g2kin (ig) )
           gk (3, ig) = gk (3, ig) / SQRT (g2kin (ig) )
        ENDIF
     ENDDO
!$OMP END DO
!$OMP END PARALLEL
     !
     CALL gen_us_dj (ik, dvkb)
     !
     jkb = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (nt == ityp (na)) THEN
              DO ikb = 1, nh (nt)
                 jkb = jkb + 1
                 DO ig = 1, npw
                    work (ig, jkb) = dvkb (ig, jkb) * &
                                    (at (1, ipol) * gk (1, ig) + &
                                     at (2, ipol) * gk (2, ig) + &
                                     at (3, ipol) * gk (3, ig) )
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     !
     CALL gen_us_dy (ik, at (1, ipol), dvkb)
     !
     jkb = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (nt == ityp (na)) THEN
              DO ikb = 1, nh (nt)
                 jkb = jkb + 1
                 DO ig = 1, npw
                    work (ig, jkb) = work (ig, jkb) + dvkb (ig, jkb)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     !
     DEALLOCATE( gk )
     DEALLOCATE( dvkb )
     !
     ! In the case of gamma point systems becp2 is real
     ! so we have to include a factor of i before calling
     ! calbec otherwise we would be stuck with the wrong component
     ! of becp2 later on.
     !
     IF (gamma_only) work=(0._DP,1._DP)*work
     !
     CALL allocate_bec_type ( nkb, m, becp1 )
     CALL allocate_bec_type ( nkb, m, becp2 )
     !
     CALL calbec (offload_type, npw,  vkb, psi, becp1, m)
     CALL calbec (offload_type, npw, work, psi, becp2, m)
     !
     IF( noncolin ) THEN
        ALLOCATE( psc (nkb,npol,m,2) )
        psc=0._DP
     ELSE
        ALLOCATE( ps2 (nkb,m,2) )
        ps2=0._DP
     ENDIF
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
                                psc(ikb,is,ibnd,1)=psc(ikb,is,ibnd,1)+ &
                                          (0._DP,-1._DP)* &
                                     becp2%nc(jkb,js,ibnd)*deff_nc(ih,jh,na,ijs)
                                psc(ikb,is,ibnd,2)=psc(ikb,is,ibnd,2)+ &
                                        (0._DP,-1._DP)* &
                                    becp1%nc(jkb,js,ibnd)*deff_nc(ih,jh,na,ijs)
                             ENDDO
                          ENDDO
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
                       ENDIF
                    ENDDO
                 ENDDO
                 ijkb0=ijkb0+nh(nt)
              ENDIF
           ENDDO ! na
        ENDDO ! nt
     ENDDO ! m
     IF (ikb /= nkb .OR. jkb /= nkb) CALL errore ('commut_Hx_psi', 'unexpected error',1)
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
     ENDIF
     DEALLOCATE( work )
     !
  ENDIF
  !
  CALL stop_clock ('commut_Hx_psi')
  !
END SUBROUTINE
!
#else
!-----------------------------------------------------------------------
SUBROUTINE commut_Hx_psi(ik, m, ipol, psi_d, dpsi_d, l_skip_nlpp)
  !-----------------------------------------------------------------------
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
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : tpiba, at
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE klist,            ONLY : xk
  USE gvect,            ONLY : g
  USE wvfct,            ONLY : npw, npwx, et, g2kin
  USE noncollin_module, ONLY : noncolin, npol
  USE becmod,           ONLY : calbec_cuf
  USE uspp,             ONLY : nkb, vkb
  USE uspp_param,       ONLY : nh
  USE uspp_init,        ONLY : gen_us_dj, gen_us_dy
  USE control_flags,    ONLY : gamma_only
  USE pwcom,            ONLY : igk_k, current_k
  USE control_flags,    ONLY : offload_type
  USE west_gpu,         ONLY : gk, dvkb, work, ps2, psc, deff, deff_nc, becp1, becp2
  USE cublas
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ik ! k-point index
  INTEGER, INTENT(IN) :: m ! number of bands to process
  INTEGER, INTENT(IN) :: ipol ! polarization index
  COMPLEX(DP), DEVICE, INTENT(IN) :: psi_d(npwx*npol,m) ! input wavefunctions
  COMPLEX(DP), DEVICE, INTENT(OUT) :: dpsi_d(npwx*npol,m) ! output wavefunctions
  LOGICAL, INTENT(IN) :: l_skip_nlpp ! skip NLPP
  !
  ! Workspace
  !
  INTEGER :: ig, na, ibnd, ikb, jkb, nt, ih, jh, ijkb0, is, js, ijs
  INTEGER :: nh_nt
  REAL(DP) :: xk1, xk2, xk3, at1, at2, at3
  COMPLEX(DP) :: reduce, reduce2, reduce3, reduce4
  REAL(DP), PARAMETER :: tol = 1.E-10_DP
  !
  CALL start_clock_gpu('commut_Hx_psi')
  !
  ! Initialize
  !
  xk1 = xk(1,ik)
  xk2 = xk(2,ik)
  xk3 = xk(3,ik)
  at1 = at(1,ipol)
  at2 = at(2,ipol)
  at3 = at(3,ipol)
  !
  !$acc parallel loop present(gk,g,igk_k,g2kin)
  DO ig = 1,npw
     gk(1,ig) = (xk1 + g(1,igk_k(ig,current_k))) * tpiba
     gk(2,ig) = (xk2 + g(2,igk_k(ig,current_k))) * tpiba
     gk(3,ig) = (xk3 + g(3,igk_k(ig,current_k))) * tpiba
     g2kin(ig) = gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !$acc end parallel
  !
  ! this is the kinetic contribution to [H,x]: -2i (k+G)_ipol * psi
  !
  !$acc parallel loop collapse(2) present(gk)
  DO ibnd = 1,m
     DO ig = 1,npwx
        IF(ig <= npw) THEN
           dpsi_d(ig,ibnd) = (at1*gk(1,ig) + at2*gk(2,ig) + at3*gk(3,ig)) * (0._DP,-2._DP) * psi_d(ig,ibnd)
           IF(noncolin) THEN
              dpsi_d(ig+npwx,ibnd) = (at1*gk(1,ig) + at2*gk(2,ig) + at3*gk(3,ig)) * (0._DP,-2._DP) * psi_d(ig+npwx,ibnd)
           ENDIF
        ELSE
           dpsi_d(ig,ibnd) = 0._DP
           IF(noncolin) THEN
              dpsi_d(ig+npwx,ibnd) = 0._DP
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  !$acc end parallel
  !
  IF((.NOT. l_skip_nlpp) .AND. nkb > 0) THEN
     !
     ! this is the contribution from nonlocal pseudopotentials
     !
     !$acc parallel loop present(g2kin,gk)
     DO ig = 1,npw
        IF(g2kin(ig) < tol) THEN
           gk(1,ig) = 0._DP
           gk(2,ig) = 0._DP
           gk(3,ig) = 0._DP
        ELSE
           gk(1,ig) = gk(1,ig)/SQRT(g2kin(ig))
           gk(2,ig) = gk(2,ig)/SQRT(g2kin(ig))
           gk(3,ig) = gk(3,ig)/SQRT(g2kin(ig))
        ENDIF
     ENDDO
     !$acc end parallel
     !
     CALL gen_us_dj(ik,dvkb)
     !
     !$acc kernels present(work)
     work(:,:) = 0._DP
     !$acc end kernels
     !
     ijkb0 = 0
     DO nt = 1,ntyp
        nh_nt = nh(nt)
        DO na = 1,nat
           IF(nt == ityp(na)) THEN
              !$acc parallel loop collapse(2) present(work,dvkb,gk)
              DO ikb = 1,nh_nt
                 DO ig = 1,npw
                    jkb = ijkb0+ikb
                    work(ig,jkb) = dvkb(ig,jkb) * (at1*gk(1,ig) + at2*gk(2,ig) + at3*gk(3,ig))
                 ENDDO
              ENDDO
              !$acc end parallel
              ijkb0 = ijkb0+nh(nt)
           ENDIF
        ENDDO
     ENDDO
     !
     CALL gen_us_dy(ik,at(1,ipol),dvkb)
     !
     ijkb0 = 0
     DO nt = 1,ntyp
        nh_nt = nh(nt)
        DO na = 1,nat
           IF(nt == ityp(na)) THEN
              !$acc parallel loop collapse(2) present(work,dvkb)
              DO ikb = 1,nh_nt
                 DO ig = 1,npw
                    jkb = ijkb0+ikb
                    work(ig,jkb) = work(ig,jkb) + dvkb(ig,jkb)
                 ENDDO
              ENDDO
              !$acc end parallel
              ijkb0 = ijkb0+nh(nt)
           ENDIF
        ENDDO
     ENDDO
     !
     ! In the case of gamma point systems becp2 is real
     ! so we have to include a factor of i before calling
     ! calbec otherwise we would be stuck with the wrong component
     ! of becp2 later on.
     !
     IF(gamma_only) THEN
        !$acc parallel loop collapse(2) present(work)
        DO ikb = 1,nkb
           DO ig = 1,npw
              work(ig,ikb) = (0.0_DP,1.0_DP)*work(ig,ikb)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     CALL calbec_cuf(offload_type,npw,vkb,psi_d,becp1,m)
     CALL calbec_cuf(offload_type,npw,work,psi_d,becp2,m)
     !
     IF(noncolin) THEN
        !$acc kernels present(psc)
        psc(:,:,:,:) = 0._DP
        !$acc end kernels
     ELSE
        !$acc kernels present(ps2)
        ps2(:,:,:) = 0._DP
        !$acc end kernels
     ENDIF
     !
     DO ibnd = 1,m
        IF(noncolin) THEN
           CALL compute_deff_nc(deff_nc,et(ibnd,ik))
        ELSE
           CALL compute_deff_real(deff,et(ibnd,ik))
        ENDIF
        ijkb0 = 0
        DO nt = 1,ntyp
           DO na = 1,nat
              IF(nt == ityp(na)) THEN
                 nh_nt = nh(nt)
                 IF(noncolin) THEN
                    IF(npol == 1) THEN
                       !$acc parallel present(deff_nc,psc,becp2%nc,becp1%nc)
                       !$acc loop
                       DO ih = 1,nh_nt
                          reduce = 0._DP
                          reduce2 = 0._DP
                          !$acc loop reduction(+:reduce,reduce2)
                          DO jh = 1,nh_nt
                             reduce = reduce + (0._DP,-1._DP) * becp2%nc(ijkb0+jh,1,ibnd) * deff_nc(ih,jh,na,1)
                             reduce2 = reduce2 + (0._DP,-1._DP) * becp1%nc(ijkb0+jh,1,ibnd) * deff_nc(ih,jh,na,1)
                          ENDDO
                          psc(ijkb0+ih,1,ibnd,1) = reduce
                          psc(ijkb0+ih,1,ibnd,2) = reduce2
                       ENDDO
                       !$acc end parallel
                    ELSEIF(npol == 2) THEN
                       !$acc parallel present(deff_nc,psc,becp2%nc,becp1%nc)
                       !$acc loop
                       DO ih = 1,nh_nt
                          reduce = 0._DP
                          reduce2 = 0._DP
                          reduce3 = 0._DP
                          reduce4 = 0._DP
                          !$acc loop reduction(+:reduce,reduce2,reduce3,reduce4)
                          DO jh = 1,nh_nt
                             ! is = 1, js = 1,2
                             reduce = reduce + (0._DP,-1._DP) * becp2%nc(ijkb0+jh,1,ibnd) * deff_nc(ih,jh,na,1)
                             reduce = reduce + (0._DP,-1._DP) * becp2%nc(ijkb0+jh,2,ibnd) * deff_nc(ih,jh,na,2)
                             reduce2 = reduce2 + (0._DP,-1._DP) * becp1%nc(ijkb0+jh,1,ibnd) * deff_nc(ih,jh,na,1)
                             reduce2 = reduce2 + (0._DP,-1._DP) * becp1%nc(ijkb0+jh,2,ibnd) * deff_nc(ih,jh,na,2)
                             ! is = 2, js = 1,2
                             reduce3 = reduce3 + (0._DP,-1._DP) * becp2%nc(ijkb0+jh,1,ibnd) * deff_nc(ih,jh,na,3)
                             reduce3 = reduce3 + (0._DP,-1._DP) * becp2%nc(ijkb0+jh,2,ibnd) * deff_nc(ih,jh,na,4)
                             reduce4 = reduce4 + (0._DP,-1._DP) * becp1%nc(ijkb0+jh,1,ibnd) * deff_nc(ih,jh,na,3)
                             reduce4 = reduce4 + (0._DP,-1._DP) * becp1%nc(ijkb0+jh,2,ibnd) * deff_nc(ih,jh,na,4)
                          ENDDO
                          psc(ijkb0+ih,1,ibnd,1) = reduce
                          psc(ijkb0+ih,1,ibnd,2) = reduce2
                          psc(ijkb0+ih,2,ibnd,1) = reduce3
                          psc(ijkb0+ih,2,ibnd,2) = reduce4
                       ENDDO
                       !$acc end parallel
                    ENDIF
                 ELSEIF(gamma_only) THEN
                    ! Note the different prefactors due to the factor of i introduced to work(:,:),
                    ! as becp[1,2] are real.
                    !
                    !$acc parallel present(deff,ps2,becp2%r,becp1%r)
                    !$acc loop
                    DO ih = 1,nh_nt
                       reduce = 0._DP
                       reduce2 = 0._DP
                       !$acc loop reduction(+:reduce,reduce2)
                       DO jh = 1,nh_nt
                          reduce = reduce + becp2%r(ijkb0+jh,ibnd) * (1._DP,0._DP) * deff(ih,jh,na)
                          reduce2 = reduce2 + becp1%r(ijkb0+jh,ibnd) * (-1._DP,0._DP) * deff(ih,jh,na)
                       ENDDO
                       ps2(ijkb0+ih,ibnd,1) = reduce
                       ps2(ijkb0+ih,ibnd,2) = reduce2
                    ENDDO
                    !$acc end parallel
                 ELSE
                    !$acc parallel present(deff,ps2,becp2%k,becp1%k)
                    !$acc loop
                    DO ih = 1,nh_nt
                       reduce = 0._DP
                       reduce2 = 0._DP
                       !$acc loop reduction(+:reduce,reduce2)
                       DO jh = 1,nh_nt
                          reduce = reduce + becp2%k(ijkb0+jh,ibnd) * (0._DP,-1._DP) * deff(ih,jh,na)
                          reduce2 = reduce2 + becp1%k(ijkb0+jh,ibnd) * (0._DP,-1._DP) * deff(ih,jh,na)
                       ENDDO
                       ps2(ijkb0+ih,ibnd,1) = reduce
                       ps2(ijkb0+ih,ibnd,2) = reduce2
                    ENDDO
                    !$acc end parallel
                 ENDIF
                 ijkb0 = ijkb0+nh(nt)
              ENDIF
           ENDDO ! na
        ENDDO ! nt
     ENDDO ! m
     !
     IF(noncolin) THEN
        !$acc host_data use_device(vkb,psc,work)
        CALL ZGEMM('N','N',npw,m*npol,nkb,(1._DP,0._DP),vkb,npwx,psc,nkb,(1._DP,0._DP),dpsi_d,npwx)
        CALL ZGEMM('N','N',npw,m*npol,nkb,(1._DP,0._DP),work,npwx,psc(1,1,1,2),nkb,(1._DP,0._DP),dpsi_d,npwx)
        !$acc end host_data
     ELSE
        !$acc host_data use_device(vkb,ps2,work)
        CALL ZGEMM('N','N',npw,m,nkb,(1._DP,0._DP),vkb,npwx,ps2,nkb,(1._DP,0._DP),dpsi_d,npwx)
        CALL ZGEMM('N','N',npw,m,nkb,(1._DP,0._DP),work,npwx,ps2(1,1,2),nkb,(1._DP,0._DP),dpsi_d,npwx)
        !$acc end host_data
     ENDIF
     !
  ENDIF
  !
  CALL stop_clock_gpu('commut_Hx_psi')
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE compute_deff_real(deff, et)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the effective value of the D-eS coefficients
  ! which appear often in many expressions in the US or PAW case.
  ! This routine is for the collinear case.
  !
  USE kinds,       ONLY : DP
  USE ions_base,   ONLY : nat
  USE uspp,        ONLY : deeq
  USE uspp_param,  ONLY : nhm
  USE lsda_mod,    ONLY : current_spin
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et
  ! The eigenvalues of the hamiltonian
  REAL(DP), INTENT(OUT) :: deff(nhm,nhm,nat)
  ! Effective values of the D-eS coefficients
  !
  ! ... local variables
  !
  INTEGER :: na, i, j
  !
  !$acc parallel loop collapse(3) present(deff,deeq)
  DO na = 1,nat
     DO i = 1,nhm
        DO j = 1,nhm
           deff(i,j,na) = deeq(i,j,na,current_spin)
        ENDDO
     ENDDO
  ENDDO
  !$acc end parallel
  !
END SUBROUTINE
#endif
