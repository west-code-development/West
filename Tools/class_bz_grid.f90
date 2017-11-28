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
! Matteo Gerosa
!
!-----------------------------------------------------------------------
MODULE class_bz_grid
!-----------------------------------------------------------------------
   !
   USE kinds,            ONLY : DP
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   TYPE, PUBLIC :: bz_grid
      !
      INTEGER :: np1 = 0
      INTEGER :: np2 = 0
      INTEGER :: np3 = 0
      INTEGER :: nps = 0
      REAL(DP), ALLOCATABLE :: xp_cryst(:,:)
      REAL(DP), ALLOCATABLE :: xp_cart(:,:)
      REAL(DP), ALLOCATABLE :: wp(:)
      LOGICAL, ALLOCATABLE :: l_gammap(:)
      !
      ! Used only for k+q/k-q grids
      INTEGER, ALLOCATABLE :: index_kq(:,:)
      INTEGER, ALLOCATABLE :: index_q(:,:)
      REAL(DP), ALLOCATABLE :: g0(:,:,:) 
      COMPLEX(DP), ALLOCATABLE :: phase(:)
      !
      CONTAINS
      !
      PROCEDURE :: init => k_grid_init
      PROCEDURE :: init_kq => kq_grid_init
      PROCEDURE :: init_q => q_grid_init
      PROCEDURE :: get_phase => get_phase
      !
   END TYPE bz_grid
   !
   CONTAINS
   !
   !
   SUBROUTINE k_grid_init( grid, grid_type )
      !
      USE cell_base,        ONLY : at, bg
      USE klist,            ONLY : xk, wk, nkstot
      USE start_k,          ONLY : nk1, nk2, nk3
      USE westcom,          ONLY : nq
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(OUT) :: grid
      CHARACTER(LEN=1), INTENT(IN) :: grid_type
      !
      ! Workspace
      !
      INTEGER :: iq, ik, iq1, iq2, iq3, ip
      REAL(DP) :: dq1, dq2, dq3
      REAL(DP) :: eps = 1.D-8
      !
      SELECT CASE( grid_type )
      CASE ( 'K', 'k')
         !
         grid%np1 = nk1
         grid%np2 = nk2
         grid%np3 = nk3
         !
         grid%nps = nkstot
         !
         ALLOCATE ( grid%xp_cryst(3,grid%nps) )
         ALLOCATE ( grid%xp_cart(3,grid%nps) )
         ALLOCATE ( grid%wp(grid%nps) )
         !
         DO ik = 1, nkstot
            grid%xp_cart(:,ik) = xk(:,ik)
         ENDDO
         ! 
         DO ik = 1, nkstot
            grid%xp_cryst(:,ik) = xk(:,ik)
         ENDDO
         CALL cryst_to_cart( grid%nps, grid%xp_cryst, at, -1 )
         !
         DO ik = 1, nkstot
            grid%wp(ik) = wk(ik)
         ENDDO
         !
      CASE ( 'Q', 'q')
         !
         grid%np1 = nq(1)
         grid%np2 = nq(2)
         grid%np3 = nq(3)
         !
         grid%nps = nq(1) * nq(2) * nq(3)
         !
         ALLOCATE ( grid%xp_cryst(3,grid%nps) )
         ALLOCATE ( grid%xp_cart(3,grid%nps) )
         ALLOCATE ( grid%wp(grid%nps) )
         !
         dq1 = 1._DP / DBLE(nq(1))
         dq2 = 1._DP / DBLE(nq(2))
         dq3 = 1._DP / DBLE(nq(3))
         iq = 0
         DO iq1 = 1, nq(1)
            DO iq2 = 1, nq(2)
               DO iq3 = 1, nq(3)
                  iq = iq + 1
                  grid%xp_cryst(1,iq) = DBLE( iq1 - 1 ) * dq1
                  grid%xp_cryst(2,iq) = DBLE( iq2 - 1 ) * dq2
                  grid%xp_cryst(3,iq) = DBLE( iq3 - 1 ) * dq3
               ENDDO
            ENDDO
         ENDDO
         !
         grid%xp_cart = grid%xp_cryst
         CALL cryst_to_cart( grid%nps, grid%xp_cart, bg, +1 )
         !
         grid%wp = 1._DP / DBLE(grid%nps)
         !
      END SELECT
      !
      ALLOCATE( grid%l_gammap(grid%nps) )
      grid%l_gammap(:)=.FALSE.
      DO ip = 1, grid%nps
         grid%l_gammap(ip) = ( ALL( ABS ( grid%xp_cryst(:,ip) ) .LT. eps ) )
      ENDDO
      !
   END SUBROUTINE
   !
   !
   SUBROUTINE kq_grid_init( kqgrid, kgrid, qgrid, sig )
      !
      USE cell_base,        ONLY : at, bg
      USE klist,            ONLY : xk, wk, nkstot
      USE pwcom,            ONLY : nspin
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(OUT) :: kqgrid
      CLASS(bz_grid), INTENT(IN) :: kgrid, qgrid
      INTEGER, INTENT(IN) :: sig
      !
      ! Workspace
      !
      INTEGER :: ik, iq, ikq, ipol, ig
      INTEGER :: nks, nqs, nkqs
      INTEGER, ALLOCATABLE :: new_ikq(:), temp_index_ikq(:)
      REAL(DP) :: xkq(3)
      REAL(DP), ALLOCATABLE :: temp_xkq(:,:), temp_wkq(:)
      REAL(DP) :: dxk(3)
      LOGICAL :: xk_not_found
      CHARACTER(LEN=256) :: message
      CHARACTER(LEN=1) :: csig
      REAL(DP) :: eps = 1.D-8
      !
      IF ( sig == +1 ) THEN
         csig = '+'
      ELSE
         csig = '-'
      ENDIF
      !
      nks = kgrid%nps
      nqs = qgrid%nps
      !
      ALLOCATE( temp_xkq(3,nks), temp_wkq(nks) )
      ALLOCATE( new_ikq(nks), temp_index_ikq(nks) )
      !
      temp_xkq = kgrid%xp_cryst
      temp_wkq = kgrid%wp
      !
      ALLOCATE( kqgrid%index_kq(nks,nqs) )
      ALLOCATE( kqgrid%g0(3,nks,nqs) )
      !
      nkqs = 0
      new_ikq = 0
      kqgrid%g0 = 0._DP
      !
      DO ik = 1, nks
         DO iq = 1, nqs
            !
            xkq(:) = kgrid%xp_cryst(:,ik) + DBLE(sig) * qgrid%xp_cryst(:,iq) 
            !
            xk_not_found = .TRUE.
            !
            DO ikq = 1, nks
               IF ( xk_not_found ) THEN
                  dxk(:) = xkq(:)-temp_xkq(:,ikq) - NINT( xkq(:)-temp_xkq(:,ikq) )
                  IF ( ALL( ABS( dxk ) < eps ) ) THEN
                     xk_not_found = .FALSE.
                     IF ( new_ikq(ikq) == 0 ) THEN
                        nkqs = nkqs + 1
                        temp_index_ikq(nkqs) = ikq
                        new_ikq(ikq) = nkqs
                     ENDIF
                     kqgrid%index_kq(ik,iq) = new_ikq(ikq)
                     kqgrid%g0(:,ik,iq) = xkq(:) - temp_xkq(:,ikq)
                     CALL cryst_to_cart(1, kqgrid%g0(1,ik,iq), bg, +1)
                  ENDIF
               ENDIF
            ENDDO ! ikq
            !
            IF ( xk_not_found ) THEN
               WRITE(*,*) ik, iq, nkqs
               WRITE(*,*) xkq(:)
               WRITE(message,'(a,a,a)')'k ',TRIM( csig ),' q does not belong to k-point grid '
               CALL errore( 'kq_grid_init', TRIM( message ), (ik-1) * nqs + iq )
            ENDIF
            !
         ENDDO ! iq
         !
      ENDDO ! ik
      !
      ALLOCATE( kqgrid%xp_cryst(3,nspin*nkqs) , kqgrid%xp_cart(3,nspin*nkqs), kqgrid%wp(nspin*nkqs) )
      !
      DO ik = 1, nkqs
         ikq = temp_index_ikq(ik)
         kqgrid%xp_cryst(:,ik) = temp_xkq(:,ikq)
         kqgrid%wp(ik) = temp_wkq(ikq)
      ENDDO
      !
      kqgrid%xp_cart = kqgrid%xp_cryst
      CALL cryst_to_cart(nkqs, kqgrid%xp_cart, bg, +1)
      !
      IF ( nspin == 2 ) THEN
         !
         DO ik = 1, nks/2
             DO iq =1, nqs
               kqgrid%index_kq(nks/2+ik,iq) = kqgrid%index_kq(ik,iq) + nkqs
             ENDDO
         ENDDO
         !
         DO ikq = 1, nkqs
            kqgrid%xp_cart(:,ikq + nkqs) = kqgrid%xp_cart(:,ikq)
            kqgrid%xp_cryst(:,ikq + nkqs) = kqgrid%xp_cryst(:,ikq)
         ENDDO
         !
         nkqs = 2*nkqs
      ENDIF
      !
      kqgrid%nps = nkqs
      !
      DEALLOCATE( new_ikq, temp_index_ikq )
      DEALLOCATE( temp_xkq, temp_wkq )
      !
   END SUBROUTINE
   !
   !
   SUBROUTINE q_grid_init( qgrid, kgrid, k1grid )
      !
      USE cell_base,        ONLY : at, bg
      USE klist,            ONLY : xk, wk, nkstot
      USE pwcom,            ONLY : nspin
      USE westcom,          ONLY : nq
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(OUT) :: qgrid
      CLASS(bz_grid), INTENT(IN) :: kgrid, k1grid
      !
      ! Workspace
      !
      INTEGER :: ik, ikk, iq, iqq, iq1, iq2, iq3, ipol, ig
      INTEGER :: nks, nks1, nqs, temp_nqs
      INTEGER, ALLOCATABLE :: new_iq(:), temp_index_iq(:)
      REAL(DP) :: xq(3)
      REAL(DP), ALLOCATABLE :: temp_xq(:,:), temp_wq(:)
      REAL(DP) :: dxq(3), dq1, dq2, dq3
      LOGICAL :: xq_not_found
      CHARACTER(LEN=256) :: message
      CHARACTER(LEN=1) :: csig
      REAL(DP) :: eps = 1.D-8
      !
      nks = kgrid%nps
      nks1 = k1grid%nps
      temp_nqs = nq(1)*nq(2)*nq(3)
      !
      ALLOCATE( temp_xq(3,temp_nqs), temp_wq(temp_nqs) )
      ALLOCATE( new_iq(temp_nqs), temp_index_iq(temp_nqs) )
      !
      dq1 = 1._DP / DBLE(nq(1))
      dq2 = 1._DP / DBLE(nq(2))
      dq3 = 1._DP / DBLE(nq(3))
      iq = 0
      DO iq1 = 1, nq(1)
         DO iq2 = 1, nq(2)
            DO iq3 = 1, nq(3)
               iq = iq + 1
               temp_xq(1,iq) = DBLE( iq1 - 1 ) * dq1
               temp_xq(2,iq) = DBLE( iq2 - 1 ) * dq2
               temp_xq(3,iq) = DBLE( iq3 - 1 ) * dq3
            ENDDO
         ENDDO
      ENDDO
      temp_wq = 1 / DBLE(nq(1)*nq(2)*nq(3))
      !
      ALLOCATE( qgrid%index_q(nks,nks) )
      ALLOCATE( qgrid%g0(3,nks,nks) )
      !
      nqs = 0
      new_iq = 0
      qgrid%g0 = 0._DP
      !
      DO ik = 1, nks
         !
         DO ikk = 1, nks1
            !
            xq(:) = kgrid%xp_cryst(:,ik) - k1grid%xp_cryst(:,ikk) 
            !
            xq_not_found = .TRUE.
            !
            DO iq = 1, temp_nqs
               IF ( xq_not_found ) THEN
                  dxq(:) = xq(:)-temp_xq(:,iq) - NINT( xq(:)-temp_xq(:,iq) )
                  IF ( ALL( ABS( dxq ) < eps ) ) THEN
                     xq_not_found = .FALSE.
                     IF ( new_iq(iq) == 0 ) THEN
                        nqs = nqs + 1
                        temp_index_iq(nqs) = iq
                        new_iq(iq) = nqs
                     ENDIF
                     qgrid%index_q(ik,ikk) = new_iq(iq)
                     qgrid%g0(:,ik,ikk) = xq(:) - temp_xq(:,iq)
                     CALL cryst_to_cart(1, qgrid%g0(1,ik,ikk), bg, +1)
                  ENDIF
               ENDIF
            ENDDO ! iq
            !
            IF ( xq_not_found ) THEN
               WRITE(*,*) ik, ikk, nqs
               WRITE(*,*) xq(:)
               WRITE(message,'(a)') " k - k' is not a commensurate with k-point grid "
               CALL errore( 'q_grid_init', TRIM( message ), (ik-1) * nqs + ikk )
            ENDIF
            !
         ENDDO ! ikk
         !
      ENDDO ! ik
      !
      ALLOCATE( qgrid%xp_cryst(3,nqs) , qgrid%xp_cart(3,nqs), qgrid%wp(nqs) )
      !
      DO iqq = 1, nqs
         iq = temp_index_iq(iqq)
         qgrid%xp_cryst(:,iqq) = temp_xq(:,iq)
         qgrid%wp(iqq) = temp_wq(iq)
      ENDDO
      !
      qgrid%xp_cart = qgrid%xp_cryst
      CALL cryst_to_cart(nqs, qgrid%xp_cart, bg, +1)
      !
      ALLOCATE( qgrid%l_gammap(nqs) )
      qgrid%l_gammap(:)=.FALSE.
      DO iq = 1, nqs
         qgrid%l_gammap(iq) = ( ALL( ABS ( qgrid%xp_cryst(:,iq) ) .LT. eps ) )
      ENDDO
      !
      qgrid%nps = nqs
      !
      DEALLOCATE( new_iq, temp_index_iq )
      DEALLOCATE( temp_xq, temp_wq )
      !
   END SUBROUTINE
   !
   SUBROUTINE get_phase( grid_aux, ik, ikk )
      !
      USE cell_base,        ONLY : bg
      USE gvecs,            ONLY : ngms, nls
      USE gvect,            ONLY : g
      USE fft_base,         ONLY : dffts
      USE fft_interfaces,   ONLY : invfft
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(INOUT) :: grid_aux
      INTEGER, INTENT(IN) :: ik, ikk
      !
      ! Workspace
      !
      INTEGER :: ipol, ig, ig0
      REAL(DP) :: eps = 1.D-8
      !
      ig0 = 0
      DO ig = 1, ngms
         IF ( ALL ( ABS( g(:,ig) - grid_aux%g0(:,ik,ikk) ) < eps ) ) THEN
            ig0 = ig
         ENDIF
      ENDDO
      !
      IF (.NOT. ALLOCATED(grid_aux%phase) ) ALLOCATE ( grid_aux%phase(dffts%nnr) )
      grid_aux%phase = (0._DP, 0._DP)
      IF ( ig0 > 0 ) THEN
         grid_aux%phase( nls(ig0) ) = (1._DP, 0._DP)
      ENDIF
      ! phase = exp(-iG_0*r)
      CALL invfft( 'Wave', grid_aux%phase, dffts )
      grid_aux%phase(1:dffts%nnr) = CONJG( grid_aux%phase(1:dffts%nnr) )
      !
   END SUBROUTINE
   ! 
END MODULE
