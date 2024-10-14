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
SUBROUTINE linsolve_commut_Hx(iks,m,e,fin,fout)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : tr2_dfpt,l_skip_nl_part_of_hcomr
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : bg
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : noncolin,npol
  USE wvfct,                ONLY : g2kin
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : intra_bgrp_comm
#if defined(__CUDA)
  USE west_gpu,             ONLY : ep_pol,e_pol,phi
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: iks
  INTEGER,INTENT(IN) :: m
  REAL(DP),INTENT(IN) :: e
  COMPLEX(DP),INTENT(IN) :: fin(npwx*npol)
  COMPLEX(DP),INTENT(OUT) :: fout(npwx*npol,3)
  !
  ! Workspace
  !
  INTEGER :: ig,ipol
  INTEGER :: ierr
  REAL(DP) :: ep
#if !defined(__CUDA)
  REAL(DP),ALLOCATABLE :: ep_pol(:)
  REAL(DP),ALLOCATABLE :: e_pol(:)
  COMPLEX(DP),ALLOCATABLE :: phi(:,:)
#endif
  REAL(DP),PARAMETER :: factor = 1.35_DP
  !
#if defined(__CUDA)
  CALL start_clock_gpu('linHx')
#else
  CALL start_clock('linHx')
#endif
  !
#if !defined(__CUDA)
  ALLOCATE(ep_pol(3))
  ALLOCATE(e_pol(3))
  ALLOCATE(phi(npwx*npol,3))
#endif
  !
  !$acc enter data create(fout)
  !
  DO ipol = 1,3
     CALL commut_Hx_psi(iks,1,ipol,fin,fout(1,ipol),l_skip_nl_part_of_hcomr)
  ENDDO
  !
  !$acc parallel loop collapse(2) present(phi,fout,bg)
  DO ipol = 1,3
     DO ig = 1,npwx*npol
        phi(ig,ipol) = fout(ig,1)*bg(ipol,1)+fout(ig,2)*bg(ipol,2)+fout(ig,3)*bg(ipol,3)
     ENDDO
  ENDDO
  !$acc end parallel
  !
  CALL apply_alpha_pc_to_m_wfcs(m,3,phi,(1._DP,0._DP))
  !
  ep = 0._DP
  !$acc parallel loop reduction(+:ep) present(fin,g2kin) copy(ep)
  DO ig = 1,npw
     ep = ep+CONJG(fin(ig))*fin(ig)*g2kin(ig)
     IF(noncolin) THEN
        ep = ep+CONJG(fin(ig+npwx))*fin(ig+npwx)*g2kin(ig)
     ENDIF
  ENDDO
  !$acc end parallel
  !
  ep = ep*factor
  !
  CALL mp_sum(ep,intra_bgrp_comm)
  !
  !$acc kernels present(ep_pol)
  ep_pol(:) = ep
  !$acc end kernels
  !
  !$acc kernels present(e_pol)
  e_pol(:) = e
  !$acc end kernels
  !
  CALL precondition_m_wfcts(3,phi,fout,ep_pol)
  CALL linsolve_sternheimer_m_wfcts(m,3,phi,fout,e_pol,ep_pol,tr2_dfpt,ierr)
  !
  IF(ierr /= 0) WRITE(stdout,'(7X,"** WARNING : MACROPOL not converged, ierr = ",I8)') ierr
  !
#if !defined(__CUDA)
  DEALLOCATE(ep_pol)
  DEALLOCATE(e_pol)
  DEALLOCATE(phi)
#endif
  !
  !$acc exit data copyout(fout)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('linHx')
#else
  CALL stop_clock('linHx')
#endif
  !
END SUBROUTINE
