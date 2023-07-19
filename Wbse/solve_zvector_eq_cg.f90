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
! Yu Jin, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE solve_zvector_eq_cg(z_rhs, z_out)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE mp,                   ONLY : mp_bcast
  USE noncollin_module,     ONLY : npol
  USE pwcom,                ONLY : npwx,nspin
  USE westcom,              ONLY : forces_zeq_cg_tr,l_pre_shift,forces_inexact_krylov,&
                                 & forces_inexact_krylov_tr,do_inexact_krylov,dvg_exc_forces
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm
  USE wbse_bgrp,            ONLY : gather_bands
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: z_rhs(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: z_out(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! ... Local variables
  !
  INTEGER :: iks,ib,ig
  INTEGER :: kpt_pool_nloc,band_group_nlocx
  INTEGER :: iter
  INTEGER, PARAMETER :: max_cg_iters = 1000
  REAL(DP) :: residual_sq
  COMPLEX(DP) :: alpha,beta
  COMPLEX(DP), ALLOCATABLE :: residual_old(:),residual_new(:),dotp(:),rz_new(:),rz_old(:)
  COMPLEX(DP), ALLOCATABLE :: r_old(:,:,:),r_new(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: p(:,:,:),Ap(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z(:,:,:)
  !$acc declare device_resident(r_old,r_new,p,Ap,z)
  !
  REAL(DP) :: time_spent(2)
  REAL(DP), EXTERNAL :: get_clock
  CHARACTER(20), EXTERNAL :: human_readable_time
  !
#if defined(__CUDA)
  CALL start_clock_gpu('zvec_cg')
#else
  CALL start_clock('zvec_cg')
#endif
  !
  CALL io_push_title('Solve the Z vector equation using the CG algorithm')
  !
  kpt_pool_nloc = kpt_pool%nloc
  band_group_nlocx = band_group%nlocx
  !
  do_inexact_krylov = .FALSE.
  !
  ALLOCATE(r_new(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(r_old(npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(p    (npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(Ap   (npwx*npol, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(z    (npwx*npol, band_group%nlocx, kpt_pool%nloc))
  !
  ALLOCATE(residual_old(nspin))
  ALLOCATE(residual_new(nspin))
  ALLOCATE(dotp(nspin))
  ALLOCATE(rz_new(nspin))
  ALLOCATE(rz_old(nspin))
  !
  CALL wbse_dot(z_rhs,z_rhs,band_group%nlocx,dotp)
  !
  WRITE(stdout,"(5X,'Norm of z_rhs_vec    = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
  !
  ! Initial guess
  !
  !$acc kernels present(z_out,z_rhs)
  z_out(:,:,:) = z_rhs
  !$acc end kernels
  !
  ! no spin flip case for z-vector equation
  !
  CALL west_apply_liouvillian(z_out, r_new, .FALSE.)
  !
  DO iks = 1,kpt_pool%nloc
     CALL gather_bands( z_out(:,:,iks), dvg_exc_forces(:,:,iks) )
  ENDDO
  !
  CALL west_apply_liouvillian_btda(z_out, r_new, .FALSE.)
  !
  !$acc parallel loop collapse(3) present(r_new,z_rhs)
  DO iks = 1,kpt_pool_nloc
     DO ib = 1,band_group_nlocx
        DO ig = 1,npwx*npol
           r_new(ig,ib,iks) = z_rhs(ig,ib,iks) - r_new(ig,ib,iks)
        ENDDO
     ENDDO
  ENDDO
  !$acc end parallel
  !
  CALL wbse_dot(r_new,r_new,band_group%nlocx,residual_new)
  !
  WRITE(stdout,"(5X,'Initial residual     = ',ES15.8)") REAL(SUM(residual_new),KIND=DP)
  !
  CALL cg_precondition(r_new,z,l_pre_shift)
  !
  CALL wbse_dot(r_new,z,band_group%nlocx,rz_new)
  !
  !$acc kernels present(p,z)
  p(:,:,:) = z
  !$acc end kernels
  !
  !$acc kernels present(r_old,r_new)
  r_old(:,:,:) = r_new
  !$acc end kernels
  !
  residual_old(:) = residual_new
  rz_old(:) = rz_new
  !
  DO iter = 1, max_cg_iters
     !
     time_spent(1) = get_clock('zvec_cg')
     !
     CALL west_apply_liouvillian(p, Ap, .FALSE.)
     !
     DO iks = 1,kpt_pool%nloc
        CALL gather_bands( p(:,:,iks), dvg_exc_forces(:,:,iks) )
     ENDDO
     !
     CALL west_apply_liouvillian_btda(p, Ap, .FALSE.)
     !
     CALL wbse_dot(p,Ap,band_group%nlocx,dotp)
     !
     alpha = SUM(rz_old) / SUM(dotp)
     !
     !$acc parallel loop collapse(3) present(z_out,p)
     DO iks = 1,kpt_pool_nloc
        DO ib = 1,band_group_nlocx
           DO ig = 1,npwx*npol
              z_out(ig,ib,iks) = z_out(ig,ib,iks) + alpha * p(ig,ib,iks)
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc parallel loop collapse(3) present(r_new,r_old,Ap)
     DO iks = 1,kpt_pool_nloc
        DO ib = 1,band_group_nlocx
           DO ig = 1,npwx*npol
              r_new(ig,ib,iks) = r_old(ig,ib,iks) - alpha * Ap(ig,ib,iks)
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
     !
     CALL wbse_dot(r_new,r_new,band_group%nlocx,residual_new)
     !
     WRITE(stdout,"(/,5X,'                  *----------*              *-----------------*')")
     WRITE(stdout,"(  5X,'#     Iteration = | ', I8,' |','   ','Residual = | ', ES15.8,' |')") &
     & iter, REAL(SUM(residual_new),KIND=DP)
     WRITE(stdout,"(  5X,'                  *----------*              *-----------------*')")
     !
     residual_sq = ABS(SUM(residual_new))
     !
     CALL mp_bcast(residual_sq,0,inter_image_comm)
     !
     do_inexact_krylov = .FALSE.
     !
     IF(forces_inexact_krylov > 0) THEN
        IF(residual_sq < forces_inexact_krylov_tr) do_inexact_krylov = .TRUE.
     ENDIF
     !
     IF(residual_sq < forces_zeq_cg_tr) THEN
        !
        time_spent(2) = get_clock('zvec_cg')
        !
        WRITE(stdout,"(5X,'Time spent in last iteration ',a)") &
        & TRIM(human_readable_time(time_spent(2)-time_spent(1)))
        !
        EXIT
        !
     ENDIF
     !
     CALL cg_precondition(r_new,z,l_pre_shift)
     !
     CALL wbse_dot(r_new,z,band_group%nlocx,rz_new)
     !
     beta = SUM(rz_new) / SUM(rz_old)
     !
     !$acc parallel loop collapse(3) present(p,z)
     DO iks = 1,kpt_pool_nloc
        DO ib = 1,band_group_nlocx
           DO ig = 1,npwx*npol
              p(ig,ib,iks) = z(ig,ib,iks) + beta * p(ig,ib,iks)
           ENDDO
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc kernels present(r_old,r_new)
     r_old(:,:,:) = r_new
     !$acc end kernels
     !
     residual_old(:) = residual_new
     rz_old(:) = rz_new
     !
     time_spent(2) = get_clock('zvec_cg')
     !
     WRITE(stdout,"(5X,'Time spent in last iteration ',a)") &
     & TRIM(human_readable_time(time_spent(2)-time_spent(1)))
     !
  ENDDO ! iter
  !
  do_inexact_krylov = .FALSE.
  !
  DEALLOCATE(r_new)
  DEALLOCATE(r_old)
  DEALLOCATE(p)
  DEALLOCATE(Ap)
  DEALLOCATE(z)
  DEALLOCATE(residual_old)
  DEALLOCATE(residual_new)
  DEALLOCATE(dotp)
  DEALLOCATE(rz_old)
  DEALLOCATE(rz_new)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('zvec_cg')
#else
  CALL stop_clock('zvec_cg')
#endif
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE cg_precondition(x, px, turn_shift)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_bcast
  USE buffers,              ONLY : get_buffer
  USE pwcom,                ONLY : npwx
  USE westcom,              ONLY : nbnd_occ,lrwfc,iuwfc,n_trunc_bands
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d,evc_work=>evc_d
  USE wavefunctions,        ONLY : evc_host=>evc
  USE wvfct,                ONLY : g2kin
  USE wvfct_gpum,           ONLY : et=>et_d
  USE west_gpu,             ONLY : reallocate_ps_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc
  USE wvfct,                ONLY : g2kin,et
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: x(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: turn_shift
  COMPLEX(DP), INTENT(OUT) :: px(npwx,band_group%nlocx,kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: ig, ibnd, nbndval, nbnd_do, lbnd, iks, band_group_myoffset
  REAL(DP):: tmp,tmp_abs,tmp_sgn
  REAL(DP), ALLOCATABLE :: g2kin_save(:,:)
  REAL(DP), PARAMETER :: minimum = 1._DP
  !
#if defined(__CUDA)
  CALL start_clock_gpu('precd_cg')
#else
  CALL start_clock('precd_cg')
#endif
  !
  band_group_myoffset = band_group%myoffset
  !
  !$acc kernels present(px)
  px(:,:,:) = (0._DP, 0._DP)
  !$acc end kernels
  !
  ALLOCATE(g2kin_save(npwx,kpt_pool%nloc))
  !$acc enter data create(g2kin_save)
  !
  DO iks = 1,kpt_pool%nloc
     !
     CALL g2_kin(iks)
     !
     !$acc kernels present(g2kin_save,g2kin)
     g2kin_save(:,iks) = g2kin
     !$acc end kernels
     !
     nbndval = nbnd_occ(iks)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     !$acc parallel vector_length(1024) present(g2kin_save,px,x)
     !$acc loop
     DO lbnd = 1,nbnd_do
        DO ig = 1,npwx
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = band_group_myoffset+lbnd
           !
           IF(turn_shift) THEN
              tmp = g2kin_save(ig,iks)-et(ibnd+n_trunc_bands,iks)
           ELSE
              tmp = g2kin_save(ig,iks)
           ENDIF
           !
           ! Same as the following line but without thread divergence
           ! IF(ABS(tmp) < minimum) tmp = SIGN(minimum,tmp)
           !
           tmp_abs = MAX(ABS(tmp),minimum)
           tmp_sgn = SIGN(1._DP,tmp)
           tmp = tmp_sgn*tmp_abs
           !
           px(ig,lbnd,iks) = x(ig,lbnd,iks)/tmp
           !
        ENDDO
        !
     ENDDO
     !$acc end parallel
     !
     ! Pc[k]*px(k)
     !
     ! load evc from iks to apply Pc of the current spin channel
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
#if defined(__CUDA)
     CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
     !
     CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,px(:,:,iks),(1._DP,0._DP))
     !
  ENDDO
  !
  !$acc exit data delete(g2kin_save)
  DEALLOCATE(g2kin_save)
  !
#if defined(__CUDA)
  CALL stop_clock_gpu('precd_cg')
#else
  CALL stop_clock('precd_cg')
#endif
  !
END SUBROUTINE
