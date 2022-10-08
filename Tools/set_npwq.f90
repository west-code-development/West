!
! Copyright (C) 2015-2022 M. Govoni
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
  USE kinds,                  ONLY : DP
  USE westcom,                ONLY : npwq,npwq_g,npwqx,ngq,ngq_g,igq_q,l_use_ecutrho,fftdriver,dfft_io
  USE mp,                     ONLY : mp_max, mp_sum
  USE mp_global,              ONLY : intra_bgrp_comm
  USE gvect,                  ONLY : ig_l2g,ngm,ngmx,g
  USE gvecw,                  ONLY : gcutw
  USE pwcom,                  ONLY : npwx
  USE control_flags,          ONLY : gamma_only
  USE types_bz_grid,          ONLY : q_grid
  USE fft_base,               ONLY : dffts
  USE klist,                  ONLY : ngk
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER, EXTERNAL :: n_plane_waves
  REAL(DP), ALLOCATABLE :: gq2kin(:)
  INTEGER :: iq, ig
  !
  IF ( gamma_only ) THEN
     !
     IF( l_use_ecutrho ) THEN
        npwq      = ngm
        npwqx     = ngmx
        fftdriver = 'Rho'
        dfft_io   = dffts
     ELSE
        npwq      = ngk(1)
        npwqx     = npwx
        fftdriver = 'Wave'
        CALL set_dfft_ecutwf( dfft_io )
     ENDIF
     npwq_g=MAXVAL(ig_l2g(1:npwq))
     CALL mp_max(npwq_g,intra_bgrp_comm)
     !
  ELSE
     !
     IF( l_use_ecutrho ) CALL errore("set_npwq", "Rho grid not implemented with q-points",1)
     fftdriver = 'Wave'
     !
     CALL set_dfft_ecutwf( dfft_io )
     !
     npwqx = n_plane_waves( gcutw, q_grid%np, q_grid%p_cart, g, ngm )
     !
     ALLOCATE( gq2kin(npwqx) )
     ALLOCATE( ngq(q_grid%np) )
     ALLOCATE( igq_q(npwqx,q_grid%np) )
     igq_q(:,:) = 0
     DO iq = 1, q_grid%np
        CALL gq_sort( q_grid%p_cart(:,iq), ngm, g, gcutw, ngq(iq), igq_q(:,iq), gq2kin )
     ENDDO
     !
     DEALLOCATE(gq2kin)
     !
     ! ... compute the global number of q+G vectors for each q-point
     !
     ALLOCATE( ngq_g(q_grid%np) )
     !
     ngq_g = 0
     ngq_g(:) = ngq(:)
     CALL mp_sum( ngq_g, intra_bgrp_comm )
     !
     ! ... compute the maximum G vector index among all q+G in processors
     !
     npwq_g = 0
     DO iq = 1, q_grid%np
        DO ig = 1, ngq(iq)
           npwq_g = MAX( npwq_g, ig_l2g(igq_q(ig,iq)) )
        ENDDO
     ENDDO
     !
     CALL mp_max( npwq_g, intra_bgrp_comm )
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
   USE kinds,                 ONLY : DP
   USE constants,             ONLY : eps8
   USE westcom,               ONLY : npwqx
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
   gq (:) = 0._DP
   !
   DO ng = 1, ngm
      qq = SUM( ( q(:) + g(:,ng) )**2 )
      IF ( qq <= eps8 ) qq = 0._DP
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
!-----------------------------------------------------------------------
SUBROUTINE gq_l2gmap_kdip( npw_g, ngk_g, ngk, igk_l2g, igk_l2g_kdip )
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine maps local G+k index to the global G vector index
  ! ... the mapping is used to collect wavefunctions subsets distributed
  ! ... across processors.
  ! ... This map is used to obtained the G+k grids related to each kpt
  !
  USE mp_bands,               ONLY : intra_bgrp_comm
  USE mp,                     ONLY : mp_sum
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
!$OMP PARALLEL PRIVATE(ig_, ig)
!$OMP WORKSHARE
  igwk_lup = 0
!$OMP END WORKSHARE
!$OMP DO
  DO ig_ = 1, ngk_g
     igwk_lup(igwk_(ig_)) = ig_
  ENDDO
!$OMP END DO
!$OMP DO
  DO ig = 1, ngk
     igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
  !
  DEALLOCATE( igwk_lup )
  !
  DEALLOCATE( itmp, igwk_ )
  !
END SUBROUTINE
 !
 SUBROUTINE set_dfft_ecutwf (dfft)
    !
    ! OUTPUT  : dfft = FFT descriptor with real space set accorfing to ecutwfc
    !
    USE kinds,                ONLY : DP
    USE klist,                ONLY : xk, nks
    USE cell_base,            ONLY : at,bg
    USE control_flags,        ONLY : gamma_only
    USE fft_types,            ONLY : fft_type_descriptor, fft_type_init
    USE gvecw,                ONLY : gcutw
    USE mp,                   ONLY : mp_max
    USE mp_bands,             ONLY : ntask_groups
    USE mp_global,            ONLY : intra_bgrp_comm,inter_pool_comm
    USE stick_base,           ONLY : sticks_map
    USE gvecs,                ONLY : gcutms
    USE realus,               ONLY : real_space
    USE symm_base,            ONLY : fft_fact
    USE mp_bands,             ONLY : nyfft,nproc_bgrp
    USE command_line_options, ONLY : nmany_
    !
    IMPLICIT NONE
    !
    !
    ! I/O
    !
    TYPE ( fft_type_descriptor ), INTENT(OUT) :: dfft ! customized fft descriptor
    !
    ! Workspace
    !
    REAL(DP) :: gkcut
    TYPE( sticks_map ) :: smap
    INTEGER :: ik
    !
    LOGICAL :: lpara
    lpara =  ( nproc_bgrp > 1 )
    !
    ! ... calculate gkcut = max |k+G|^2, in (2pi/a)^2 units
    !
    IF (nks == 0) THEN
       !
       ! if k-points are automatically generated (which happens later)
       ! use max(bg)/2 as an estimate of the largest k-point
       !
       gkcut = 0.5_DP * MAX ( &
          SQRT (SUM(bg (1:3, 1)**2) ), &
          SQRT (SUM(bg (1:3, 2)**2) ), &
          SQRT (SUM(bg (1:3, 3)**2) ) )
    ELSE
       gkcut = 0._DP
       DO ik = 1, nks
          gkcut = MAX (gkcut, SQRT ( SUM(xk (1:3, ik)**2) ) )
       ENDDO
    ENDIF
    !
    gkcut = (SQRT (gcutw) + gkcut)**2
    !
    ! ... find maximum value among all the processors
    !
    CALL mp_max (gkcut, inter_pool_comm )
    !
    dfft%has_task_groups = (ntask_groups >1) .AND. .NOT. real_space
    CALL fft_type_init( dfft, smap, "wave", gamma_only, lpara, intra_bgrp_comm, at, bg, gkcut, &
    & MAX(gcutms/gkcut/4.0_DP,1.0_DP), fft_fact=fft_fact,nyfft=nyfft,nmany=nmany_)
    !
END SUBROUTINE
