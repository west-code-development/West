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
SUBROUTINE set_npwq()
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE westcom,         ONLY : npwq,npwq_g,npwqx,ngq,ngq_g,igq_q,l_use_ecutrho,fftdriver
  USE mp,              ONLY : mp_max, mp_sum
  USE mp_global,       ONLY : intra_bgrp_comm, inter_bgrp_comm, nbgrp, inter_pool_comm, intra_pool_comm
  USE gvect,           ONLY : ig_l2g,ngm,ngmx,g
  USE gvecw,           ONLY : gcutw
  USE pwcom,           ONLY : npw,npwx
  USE control_flags,   ONLY : gamma_only
  USE class_bz_grid,   ONLY : bz_grid
  USE types_bz_grid,   ONLY : q_grid
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER, EXTERNAL :: n_plane_waves
  REAL(DP), ALLOCATABLE :: gq2kin(:)
  INTEGER :: iq, ig
!  INTEGER :: npwqx_g
!  INTEGER, ALLOCATABLE :: igq_l2g(:)
  !
  IF ( gamma_only ) THEN
     !
     IF( l_use_ecutrho ) THEN
        npwq      = ngm
        npwqx     = ngmx
        fftdriver = 'Dense'
     ELSE
        npwq      = npw
        npwqx     = npwx
        fftdriver = 'Wave'
     ENDIF
     !ALLOCATE(q0ig_l2g(npwq0))
     !q0ig_l2g(1:npwq0) = ig_l2g(1:npwq0)
     !npwq0_g=MAXVAL(q0ig_l2g(1:npwq0))
     npwq_g=MAXVAL(ig_l2g(1:npwq))
     CALL mp_max(npwq_g,intra_bgrp_comm)
     !
  ELSE
     !
     fftdriver = 'Wave'
     ! 'Dense' grid not yet implemented
     !
     npwqx = n_plane_waves( gcutw, q_grid%nps, q_grid%xp_cart, g, ngm )
     !
     ALLOCATE( gq2kin(npwqx) )
     ALLOCATE( ngq(q_grid%nps) )
     ALLOCATE( igq_q(npwqx,q_grid%nps) )
     !ALLOCATE( igq_l2g(npwqx,q_grid%nps) )
     igq_q(:,:) = 0
     !igq_l2g(:,:) = 0
     DO iq = 1, q_grid%nps
        CALL gq_sort( q_grid%xp_cart(:,iq), ngm, g, gcutw, ngq(iq), igq_q(:,iq), gq2kin )
        !CALL gq_l2gmap( ngm, ig_l2g(1), ngq(iq), igq_q(1,iq), igq_l2g(1,iq) )
     ENDDO
     !
     DEALLOCATE(gq2kin)
     !
     ! ... compute the global number of q+G vectors for each q-point
     !
     ALLOCATE( ngq_g(q_grid%nps) )
     !
     ngq_g = 0
     ngq_g(:) = ngq(:)
     !CALL mp_sum( ngq_g, inter_pool_comm )
     !CALL mp_sum( ngq_g, intra_pool_comm )
     CALL mp_sum( ngq_g, intra_bgrp_comm )
     !ngq_g = ngq_g / nbgrp
     !
     ! ... compute the maximum G vector index among all q+G in processors
     !
     npwq_g = 0 
     DO iq = 1, q_grid%nps
        DO ig = 1, ngq(iq)
           npwq_g = MAX( npwq_g, ig_l2g(igq_q(ig,iq)) )
        ENDDO
     ENDDO
     !npwq_g = MAXVAL( igq_l2g(:,:) )
     !
     !CALL mp_max( npwq_g, intra_pool_comm )
     CALL mp_max( npwq_g, intra_bgrp_comm )
     !
     ! ... compute the maximum number of G vectors among all q-points
     !
!     npwqx_g = MAXVAL( ngq_g(1:q_grid%nps) )
     !
     ! ... define a further l2g map
     !
!     ALLOCATE( igq_l2g_kdip(npwqx_g, q_grid%nps) )
!     igq_l2g_kdip(:,:) = 0
!     !
!     DO iq = 1, q_grid%nps
!        ALLOCATE( igq_l2g(ngq(iq)) )
!        DO ig = 1, ngq(iq)
!           igq_l2g(ig) = ig_l2g( igq_q(ig,iq) ) 
!        ENDDO 
!        CALL gq_l2gmap_kdip( npwq_g, ngq_g(iq), ngq(iq), igq_l2g, igq_l2g_kdip(1,iq) )
!        DEALLOCATE( igq_l2g ) 
!     ENDDO
     !
  ENDIF
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------
SUBROUTINE gq_sort( q, ngm, g, ecut, ngq, igq, gq )
   !----------------------------------------------------------------------------
   !
   ! ... sorts q+g in order of increasing magnitude, up to ecut
   ! ... NB: this version should yield the same ordering for different ecut
   ! ...     and the same ordering in all machines
   !
   USE kinds,     ONLY : DP
   USE constants, ONLY : eps8
   USE westcom,   ONLY : npwqx
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(in) :: q(3)       ! the q point
   INTEGER, INTENT(in) :: ngm         ! the number of g vectors
   REAL(DP), INTENT(in) :: g(3,ngm)   ! the coordinates of G vectors
   REAL(DP), INTENT(in) :: ecut       ! the cut-off energy
   INTEGER, INTENT(out) :: ngq        ! the number of q+G vectors inside the "ecut sphere"
   INTEGER, INTENT(out) :: igq(npwqx) ! the correspondence q+G <-> G
   REAL(DP), INTENT(out) :: gq(npwqx) ! the moduli of q+G
   !
   INTEGER :: ng   ! counter on   G vectors
   INTEGER :: nk   ! counter on k+G vectors
   REAL(DP) :: qq  ! |k+G|^2
   REAL(DP) :: q2x ! upper bound for |G|
   !
   ! ... first we count the number of q+G vectors inside the cut-off sphere
   !
   q2x = ( SQRT( SUM(q(:)**2) ) + SQRT( ecut ) )**2
   !
   ngq = 0
   igq(:) = 0
   gq (:) = 0.0_dp
   !
   DO ng = 1, ngm
      qq = SUM( ( q(:) + g(:,ng) )**2 )
      IF ( qq <= eps8 ) qq = 0.d0
      !
      ! ... here if |q+G|^2 <= Ecut
      !
      IF ( qq <= ecut ) THEN
         ngq = ngq + 1
         IF ( ngq > npwqx ) &
            CALL errore( 'gq_sort', 'array gq out-of-bounds', 1 )
         !
         gq(ngq) = qq
         !
         ! set the initial value of index array
         igq(ngq) = ng
      ELSE
         ! if |G| > |q| + SQRT( Ecut )  stop search and order vectors
         IF ( SUM( g(:,ng)**2 ) > ( q2x + eps8 ) ) EXIT
      ENDIF
   ENDDO
   !
   IF ( ng > ngm ) &
      CALL infomsg( 'gq_sort', 'unexpected exit from do-loop')
   !
   ! ... order vector gq keeping initial position in index
   !
   CALL hpsort_eps( ngq, gq, igq, eps8 )
   !
   ! ... now order true |q+G|
   !
   DO nk = 1, ngq
      gq(nk) = SUM( (q(:) + g(:,igq(nk)) )**2 )
   ENDDO
   !
END SUBROUTINE
!
!----------------------------------------------------------------------------
!SUBROUTINE gq_l2gmap( ngm, ig_l2g, ngk, igk, igk_l2g )
!  !----------------------------------------------------------------------------
!  !
!  ! ... This subroutine maps local G+k index to the global G vector index
!  ! ... the mapping is used to collect wavefunctions subsets distributed
!  ! ... across processors.
!  ! ... Written by Carlo Cavazzoni
!  !
!  IMPLICIT NONE
!  !
!  ! ... Here the dummy variables
!  !
!  INTEGER, INTENT(IN)  :: ngm, ngk, igk(ngk), ig_l2g(ngm)
!  INTEGER, INTENT(OUT) :: igk_l2g(ngk)
!  INTEGER              :: ig
!  !
!  ! ... input: mapping between local and global G vector index
!  !
!  DO ig = 1, ngk
!     !
!     igk_l2g(ig) = ig_l2g(igk(ig))
!     !
!  END DO
!  !
!END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE gq_l2gmap_kdip( npw_g, ngk_g, ngk, igk_l2g, igk_l2g_kdip )
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine maps local G+k index to the global G vector index
  ! ... the mapping is used to collect wavefunctions subsets distributed
  ! ... across processors.
  ! ... This map is used to obtained the G+k grids related to each kpt
  !
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp,       ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... Here the dummy variables
  !
  INTEGER, INTENT(IN)  :: npw_g, ngk_g, ngk
  INTEGER, INTENT(IN)  :: igk_l2g(ngk)
  INTEGER, INTENT(OUT) :: igk_l2g_kdip(ngk)
  !
  INTEGER, ALLOCATABLE :: igwk_(:), itmp(:), igwk_lup(:)
  INTEGER              :: ig, ig_, ngg
  !
  !
  ALLOCATE( itmp( npw_g ) )
  ALLOCATE( igwk_( ngk_g ) )
  !
  itmp(:)  = 0
  !
  !
  DO ig = 1, ngk
     !
     itmp(igk_l2g(ig)) = igk_l2g(ig)
     !
  END DO
  !
  CALL mp_sum( itmp, intra_bgrp_comm )
  !
  ngg = 0
  DO ig = 1, npw_g
     !
     IF ( itmp(ig) == ig ) THEN
        !
        ngg = ngg + 1
        !
        igwk_(ngg) = ig
        !
     END IF
     !
  END DO
  !
  IF ( ngg /= ngk_g ) CALL errore( 'gk_l2gmap_kdip', 'unexpected dimension in ngg', 1 )
  !
  ALLOCATE( igwk_lup( npw_g ) )
  !
!$omp parallel private(ig_, ig)
!$omp workshare
  igwk_lup = 0
!$omp end workshare
!$omp do
  DO ig_ = 1, ngk_g
     igwk_lup(igwk_(ig_)) = ig_
  ENDDO
!$omp end do
!$omp do
  DO ig = 1, ngk
     igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
  ENDDO
!$omp end do
!$omp end parallel
  !
  DEALLOCATE( igwk_lup )
  !
  DEALLOCATE( itmp, igwk_ )
  !
  RETURN
  !
END SUBROUTINE
