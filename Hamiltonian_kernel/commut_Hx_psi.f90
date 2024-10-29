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
  USE wvfct,            ONLY : npw, npwx, g2kin
  USE lsda_mod,         ONLY : current_spin
  USE pwcom,            ONLY : igk_k, current_k
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp,             ONLY : nkb, vkb, deeq, deeq_nc
  USE uspp_param,       ONLY : nh, nhm
  USE uspp_init,        ONLY : gen_us_dj, gen_us_dy
  USE control_flags,    ONLY : gamma_only, offload_type
  USE westcom,          ONLY : na_ikb, ijkb0_ikb
#if defined(__CUDA)
  USE becmod,           ONLY : calbec_cuf
  USE west_gpu,         ONLY : gk, dvkb, work, ps2, psc, becp1, becp2
  USE cublas
#else
  USE becmod,           ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: ik ! k-point index
  INTEGER, INTENT(IN) :: m ! number of bands to process
  INTEGER, INTENT(IN) :: ipol ! polarization index
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,m) ! input wavefunctions
  COMPLEX(DP), INTENT(OUT) :: dpsi(npwx*npol,m) ! output wavefunctions
  LOGICAL, INTENT(IN) :: l_skip_nlpp ! skip NLPP
  !
  ! Workspace
  !
  INTEGER :: ig, na, ibnd, ikb, jkb, nt, ih, jh, ijkb0
  INTEGER :: nh_nt
  REAL(DP) :: xk1, xk2, xk3, at1, at2, at3
  COMPLEX(DP) :: reduce, reduce2, reduce3, reduce4
#if !defined(__CUDA)
  TYPE(bec_type) :: becp1 ! dimensions ( nkb, m )
  TYPE(bec_type) :: becp2 ! dimensions ( nkb, m )
  REAL(DP), ALLOCATABLE :: gk(:,:)
  COMPLEX(DP), ALLOCATABLE :: ps2(:,:,:), dvkb(:,:), work(:,:), psc(:,:,:,:)
#endif
  REAL(DP), PARAMETER :: tol = 1.E-10_DP
  COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
  COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
  COMPLEX(DP), PARAMETER :: iota = (0._DP,1._DP)
  COMPLEX(DP), PARAMETER :: two_iota = (0._DP,2._DP)
  !
#if defined(__CUDA)
  CALL start_clock_gpu('commut_Hx_psi')
#else
  CALL start_clock('commut_Hx_psi')
#endif
  !
#if !defined(__CUDA)
  ALLOCATE(gk(3,npwx))
#endif
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
  ! This is the kinetic contribution to [H,x]: -2i (k+G)_ipol * psi .
  !
  !$acc parallel loop collapse(2) present(dpsi,gk,psi)
  DO ibnd = 1,m
     DO ig = 1,npwx
        IF(ig <= npw) THEN
           dpsi(ig,ibnd) = -two_iota * (at1*gk(1,ig) + at2*gk(2,ig) + at3*gk(3,ig)) * psi(ig,ibnd)
           IF(noncolin) THEN
              dpsi(ig+npwx,ibnd) = -two_iota * (at1*gk(1,ig) + at2*gk(2,ig) + at3*gk(3,ig)) * psi(ig+npwx,ibnd)
           ENDIF
        ELSE
           dpsi(ig,ibnd) = zero
           IF(noncolin) THEN
              dpsi(ig+npwx,ibnd) = zero
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  !$acc end parallel
  !
  IF((.NOT. l_skip_nlpp) .AND. nkb > 0) THEN
     !
     ! This is the contribution from nonlocal pseudopotentials.
     !
#if !defined(__CUDA)
     ALLOCATE(work(npwx,nkb))
     ALLOCATE(dvkb(npwx,nkb))
#endif
     !
     !$acc parallel loop present(g2kin,gk)
     DO ig = 1,npw
        IF(g2kin(ig) < tol) THEN
           gk(1,ig) = zero
           gk(2,ig) = zero
           gk(3,ig) = zero
        ELSE
           gk(1,ig) = gk(1,ig) / SQRT(g2kin(ig))
           gk(2,ig) = gk(2,ig) / SQRT(g2kin(ig))
           gk(3,ig) = gk(3,ig) / SQRT(g2kin(ig))
        ENDIF
     ENDDO
     !$acc end parallel
     !
     CALL gen_us_dj(ik,dvkb)
     !
     !$acc kernels present(work)
     work(:,:) = zero
     !$acc end kernels
     !
     !$acc parallel loop collapse(2) present(work,dvkb,gk)
     DO ikb = 1,nkb
        DO ig = 1,npw
           work(ig,ikb) = dvkb(ig,ikb) * (at1*gk(1,ig) + at2*gk(2,ig) + at3*gk(3,ig))
        ENDDO
     ENDDO
     !$acc end parallel
     !
     CALL gen_us_dy(ik,at(1,ipol),dvkb)
     !
     !$acc parallel loop collapse(2) present(work,dvkb)
     DO ikb = 1,nkb
        DO ig = 1,npw
           work(ig,ikb) = work(ig,ikb) + dvkb(ig,ikb)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     ! In the case of gamma point systems becp2 is real so we have to include a factor of i before
     ! calling calbec otherwise we would be stuck with the wrong component of becp2 later on.
     !
     IF(gamma_only) THEN
        !$acc parallel loop collapse(2) present(work)
        DO ikb = 1,nkb
           DO ig = 1,npw
              work(ig,ikb) = iota * work(ig,ikb)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
#if defined(__CUDA)
     CALL calbec_cuf(offload_type,npw,vkb,psi,becp1,m)
     CALL calbec_cuf(offload_type,npw,work,psi,becp2,m)
#else
     DEALLOCATE(gk)
     DEALLOCATE(dvkb)
     !
     CALL allocate_bec_type(nkb,m,becp1)
     CALL allocate_bec_type(nkb,m,becp2)
     !
     CALL calbec(offload_type,npw,vkb,psi,becp1,m)
     CALL calbec(offload_type,npw,work,psi,becp2,m)
#endif
     !
     IF(noncolin) THEN
#if !defined(__CUDA)
        ALLOCATE(psc(nkb,npol,m,2))
#endif
        !$acc kernels present(psc)
        psc(:,:,:,:) = zero
        !$acc end kernels
     ELSE
#if !defined(__CUDA)
        ALLOCATE(ps2(nkb,m,2))
#endif
        !$acc kernels present(ps2)
        ps2(:,:,:) = zero
        !$acc end kernels
     ENDIF
     !
     IF(.NOT. ALLOCATED(na_ikb)) THEN
        ALLOCATE(na_ikb(nkb))
        ALLOCATE(ijkb0_ikb(nkb))
        !
        ! Map ikb to na and ikb to ijkb0. This is to better parallelize the compute loop below.
        !
        ijkb0 = 0
        DO nt = 1,ntyp
           DO na = 1,nat
              IF(nt == ityp(na)) THEN
                 nh_nt = nh(nt)
                 DO ih = 1,nh_nt
                    ikb = ijkb0 + ih
                    na_ikb(ikb) = na
                    ijkb0_ikb(ikb) = ijkb0
                 ENDDO
                 ijkb0 = ijkb0 + nh_nt
              ENDIF
           ENDDO
        ENDDO
        !
        !$acc enter data copyin(na_ikb,ijkb0_ikb,ityp,nh)
     ENDIF
     !
     DO ibnd = 1,m
        IF(noncolin) THEN
           !$acc parallel loop present(na_ikb,ityp,nh,ijkb0_ikb,becp2%nc,deeq_nc,becp1%nc,psc)
           DO ikb = 1,nkb
              na = na_ikb(ikb)
              nt = ityp(na)
              nh_nt = nh(nt)
              ijkb0 = ijkb0_ikb(ikb)
              ih = ikb - ijkb0
              reduce = zero
              reduce2 = zero
              reduce3 = zero
              reduce4 = zero
              !$acc loop reduction(+:reduce,reduce2,reduce3,reduce4)
              DO jh = 1,nhm
                 IF(jh <= nh_nt) THEN
                    jkb = ijkb0 + jh
                    ! is = 1, js = 1,2
                    reduce = reduce - iota * becp2%nc(jkb,1,ibnd) * deeq_nc(ih,jh,na,1)
                    reduce = reduce - iota * becp2%nc(jkb,2,ibnd) * deeq_nc(ih,jh,na,2)
                    reduce2 = reduce2 - iota * becp1%nc(jkb,1,ibnd) * deeq_nc(ih,jh,na,1)
                    reduce2 = reduce2 - iota * becp1%nc(jkb,2,ibnd) * deeq_nc(ih,jh,na,2)
                    ! is = 2, js = 1,2
                    reduce3 = reduce3 - iota * becp2%nc(jkb,1,ibnd) * deeq_nc(ih,jh,na,3)
                    reduce3 = reduce3 - iota * becp2%nc(jkb,2,ibnd) * deeq_nc(ih,jh,na,4)
                    reduce4 = reduce4 - iota * becp1%nc(jkb,1,ibnd) * deeq_nc(ih,jh,na,3)
                    reduce4 = reduce4 - iota * becp1%nc(jkb,2,ibnd) * deeq_nc(ih,jh,na,4)
                 ENDIF
              ENDDO
              psc(ikb,1,ibnd,1) = reduce
              psc(ikb,1,ibnd,2) = reduce2
              psc(ikb,2,ibnd,1) = reduce3
              psc(ikb,2,ibnd,2) = reduce4
           ENDDO
           !$acc end parallel
        ELSEIF(gamma_only) THEN
           !$acc parallel loop present(na_ikb,ityp,nh,ijkb0_ikb,becp2%r,deeq,becp1%r,ps2)
           DO ikb = 1,nkb
              na = na_ikb(ikb)
              nt = ityp(na)
              nh_nt = nh(nt)
              ijkb0 = ijkb0_ikb(ikb)
              ih = ikb - ijkb0
              reduce = zero
              reduce2 = zero
              !$acc loop reduction(+:reduce,reduce2)
              DO jh = 1,nhm
                 IF(jh <= nh_nt) THEN
                    jkb = ijkb0 + jh
                    ! Note the different prefactors due to the factor of i introduced to work,
                    ! as becp[1,2] are real.
                    reduce = reduce + becp2%r(jkb,ibnd) * deeq(ih,jh,na,current_spin)
                    reduce2 = reduce2 - becp1%r(jkb,ibnd) * deeq(ih,jh,na,current_spin)
                 ENDIF
              ENDDO
              ps2(ikb,ibnd,1) = reduce
              ps2(ikb,ibnd,2) = reduce2
           ENDDO
           !$acc end parallel
        ELSE
           !$acc parallel loop present(na_ikb,ityp,nh,ijkb0_ikb,becp2%k,deeq,becp1%k,ps2)
           DO ikb = 1,nkb
              na = na_ikb(ikb)
              nt = ityp(na)
              nh_nt = nh(nt)
              ijkb0 = ijkb0_ikb(ikb)
              ih = ikb - ijkb0
              reduce = zero
              reduce2 = zero
              !$acc loop reduction(+:reduce,reduce2)
              DO jh = 1,nhm
                 IF(jh <= nh_nt) THEN
                    jkb = ijkb0 + jh
                    reduce = reduce - iota * becp2%k(jkb,ibnd) * deeq(ih,jh,na,current_spin)
                    reduce2 = reduce2 - iota * becp1%k(jkb,ibnd) * deeq(ih,jh,na,current_spin)
                 ENDIF
              ENDDO
              ps2(ikb,ibnd,1) = reduce
              ps2(ikb,ibnd,2) = reduce2
           ENDDO
           !$acc end parallel
        ENDIF
     ENDDO ! m
     !
#if !defined(__CUDA)
     CALL deallocate_bec_type(becp1)
     CALL deallocate_bec_type(becp2)
#endif
     !
     IF(noncolin) THEN
        !$acc host_data use_device(vkb,psc,dpsi,work)
        CALL ZGEMM('N','N',npw,m*npol,nkb,one,vkb,npwx,psc,nkb,one,dpsi,npwx)
        CALL ZGEMM('N','N',npw,m*npol,nkb,one,work,npwx,psc(1,1,1,2),nkb,one,dpsi,npwx)
        !$acc end host_data
     ELSE
        !$acc host_data use_device(vkb,ps2,dpsi,work)
        CALL ZGEMM('N','N',npw,m,nkb,one,vkb,npwx,ps2,nkb,one,dpsi,npwx)
        CALL ZGEMM('N','N',npw,m,nkb,one,work,npwx,ps2(1,1,2),nkb,one,dpsi,npwx)
        !$acc end host_data
     ENDIF
     !
#if !defined(__CUDA)
     IF(noncolin) THEN
        DEALLOCATE(psc)
     ELSE
        DEALLOCATE(ps2)
     ENDIF
     DEALLOCATE(work)
#endif
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('commut_Hx_psi')
#else
  CALL stop_clock('commut_Hx_psi')
#endif
  !
END SUBROUTINE
