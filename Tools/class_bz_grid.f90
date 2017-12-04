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
      INTEGER :: ngrid(3) = (/ 1, 1, 1 /)       ! number of points in each direction
      INTEGER :: np = 1                         ! total number of points = ngrid(1) * ngrid(2) * ngrid(3) 
      INTEGER :: ns = 1                         ! total number of spin = nspin
      INTEGER :: nps = 1                        ! total number of points and spins = np * ns
      INTEGER,ALLOCATABLE :: ip(:)              ! given ips --> ip   
      INTEGER,ALLOCATABLE :: is(:)              ! given ips --> is   
      REAL(DP), ALLOCATABLE :: p_cryst(:,:)     ! coordinates of point p in crystal                [ 1:np  ] 
      REAL(DP), ALLOCATABLE :: p_cart(:,:)      ! coordinates of point p in cart ( tpiba units )   [ 1:np  ]
      REAL(DP), ALLOCATABLE :: weight(:)        ! weight of point p (sum of weights = nspin)       [ 1:nps ]
      LOGICAL, ALLOCATABLE :: l_pIsGamma(:)     ! .true. if point p = (0,0,0), else .false.        [ 1:np  ]
      !
!      ! Used only for k+q/k-q grids  
!      INTEGER, ALLOCATABLE :: index_kq(:,:)     ! given ik and iq  => index of kp = k+/-q (reported in 1BZ) 
!      INTEGER, ALLOCATABLE :: index_q(:,:)      ! given ik and ikp => index of q = k - kp (reported in 1BZ)
!      REAL(DP), ALLOCATABLE :: g0(:,:,:)        ! (3,ik,iq) or (3,ik,ikp) --> G0 that brings back the point in 1BZ 
!      COMPLEX(DP), ALLOCATABLE :: phase(:)      ! given ik and ikp ==> e^(iG0.r) 
      !
      CONTAINS
      !
      PROCEDURE :: init => k_or_q_grid_init
      PROCEDURE :: find => findp
      PROCEDURE :: add => addp
!      PROCEDURE :: init_kq => kq_grid_init
!      PROCEDURE :: init_q => q_grid_init
!      PROCEDURE :: get_phase => get_phase
      !
   END TYPE bz_grid
   !
   CONTAINS
   !
   !
   !
   SUBROUTINE k_or_q_grid_init( this, grid_type )
      !
      ! ... grid_type = [ "K", "Q" ]
      !
      USE cell_base,        ONLY : at, bg
      USE klist,            ONLY : xk, wk, nkstot
      USE start_k,          ONLY : nk1, nk2, nk3
      USE pwcom,            ONLY : nspin
      USE westcom,          ONLY : nq
      USE constants,        ONLY : eps8
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(INOUT) :: this
      CHARACTER(LEN=1), INTENT(IN) :: grid_type
      !
      ! Workspace
      !
      INTEGER :: ip, iq1, iq2, iq3, ips
      INTEGER :: i, j, k
      !
      SELECT CASE( grid_type )
      CASE ( 'K', 'k')
         !
         this%ngrid(1:3) = (/ nk1, nk2, nk3 /) 
         this%np = this%ngrid(1) * this%ngrid(2) * this%ngrid(3) 
         this%ns = nspin 
         this%nps = nkstot   !    =   np * ns  
         !
         ! generate p-vectors in cart 
         !
         ALLOCATE ( this%p_cart (3,this%np) )
         DO ip = 1, this%np
            this%p_cart(:,ip) = xk(:,ip)
         ENDDO
         !
         ! generate p-vectors in cryst  
         !
         ALLOCATE ( this%p_cryst  (3,this%np) )
         this%p_cryst(:,:) = this%p_cart(:,:) 
         CALL cryst_to_cart( this%nps, this%p_cryst, at, -1 )
         !
         ! set weights 
         !
         DO ips = 1, this%nps
           this%weight(ips) = wk(ips)
         ENDDO
         !
      CASE ( 'Q', 'q')
         !
         this%ngrid(1:3) = nq(1:3) 
         this%np = this%ngrid(1) * this%ngrid(2) * this%ngrid(3) 
         this%ns = 1
         this%nps = this%np 
         !
         ! generate p-vectors in cryst 
         !
         ALLOCATE ( this%p_cryst  (3,this%np) )
         ip = 0
         DO iq1 = 1, nq(1)
            DO iq2 = 1, nq(2)
               DO iq3 = 1, nq(3)
                  ip = ip + 1
                  this%p_cryst(1,ip) = DBLE( iq1 - 1 ) / DBLE( this%ngrid(1) ) 
                  this%p_cryst(2,ip) = DBLE( iq2 - 1 ) / DBLE( this%ngrid(2) )
                  this%p_cryst(3,ip) = DBLE( iq3 - 1 ) / DBLE( this%ngrid(3) )
               ENDDO
            ENDDO
         ENDDO
         !
         ! generate p-vectors in cart
         !
         ALLOCATE ( this%p_cart (3,this%np) )
         this%p_cart(:,:) = this%p_cryst(:,:)
         CALL cryst_to_cart( this%nps, this%p_cart, bg, +1 )
         !
         ! set weights 
         !
         ALLOCATE ( this%weight (this%nps)   )
         !
         this%weight = 1._DP / DBLE(this%np)
         !
      END SELECT
      !
      ALLOCATE ( this%ip( this%nps ) ) 
      ALLOCATE ( this%is( this%nps ) ) 
      !
      ! generate map ips --> ip and is
      !
      this%ip = 0 
      this%is = 0 
      k = 0 
      DO i = 1, this%ns
         DO j = 1, this%np
            k = k+1 
            this%ip(k) = j
            this%is(k) = i
         ENDDO
      ENDDO
      !
      ALLOCATE( this%l_pIsGamma(this%np) )
      this%l_pIsGamma(:)=.FALSE.
      DO ip = 1, this%np
         this%l_pIsGamma(ip) = ( ALL( ABS ( this%p_cryst(:,ip) ) .LT. eps8 ) )
      ENDDO
      !
   END SUBROUTINE
   !
   !
   FUNCTION findp(this,p,unit_type) RESULT(ip)
      !
      ! ... ip is the index of p (unit_type = [ "cryst", "cart"])
      ! ... if on exit ip == 0 --> p is not commensurate with this grid
      !
      USE constants, ONLY : eps8 
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      CLASS(bz_grid),INTENT(IN) :: this 
      REAL(DP),INTENT(IN) :: p(3)
      CHARACTER(LEN=*),INTENT(IN) :: unit_type
      INTEGER :: ip
      !
      ! Workspace
      ! 
      INTEGER :: i
      ! 
      ip = 0 
      SELECT CASE( unit_type ) 
      CASE("cryst")
         DO i = 1, this%np
            IF( ( ALL( ABS ( p(:) - this%p_cryst(:,i) ) .LT. eps8 ) ) ) THEN 
               ip = i 
               EXIT
            ENDIF
         ENDDO
      CASE("cart")
         DO i = 1, this%np
            IF( ( ALL( ABS ( p(:) - this%p_cart(:,i) ) .LT. eps8 ) ) ) THEN  
               ip = i 
               EXIT
            ENDIF
         ENDDO
      CASE DEFAULT
         CALL errore( 1, "unit_type not supported, supported only cryst or cart" )  
      END SELECT
      !
   END FUNCTION
   !
   !
   SUBROUTINE addp( this, pin1, pin2, pout, g0, unit_type )
      !
      ! ... out : pout and g0 
      ! ... pout = pin1 + pin2 - g0   ( g0 makes sure that pout is in 1BZ ) 
      ! ... unit_type determines the units of pin1, pin2 and pout, g0  
      !
      USE cell_base,        ONLY : at, bg
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(IN) :: this
      REAL(DP), INTENT(IN) :: pin1(3), pin2(3) 
      REAL(DP), INTENT(OUT) :: pout(3)
      REAL(DP), INTENT(OUT) :: g0(3)
      CHARACTER(LEN=*),INTENT(IN) :: unit_type
      !
      ! Workspace
      !
      REAL(DP) :: ptemp(3)  
      !
      SELECT CASE(unit_type) 
      CASE("cryst","cart")
      CASE DEFAULT
         CALL errore( 1, "unit_type not supported, supported only cryst or cart" )  
      END SELECT 
      !   
      ptemp = pin1 + pin2
      IF( unit_type == "cart" ) CALL cryst_to_cart( 1, ptemp, at, -1 )
      !
      ! ptemp is now in cryst  
      !
      g0 = NINT( ptemp ) ! in cryst 
      pout = ptemp - g0  ! in cryst 
      !
      IF( unit_type == "cart" ) THEN 
         CALL cryst_to_cart( 1, pout , bg, 1 )
         CALL cryst_to_cart( 1, g0   , bg, 1 )
      ENDIF  
      !
   END SUBROUTINE
!     !
!     IF ( sig == +1 ) THEN
!        csig = '+'
!     ELSE
!        csig = '-'
!     ENDIF
!     !
!     nks = kgrid%nps
!     nqs = qgrid%nps
!     !
!     ALLOCATE( temp_xkq(3,nks), temp_wkq(nks) )
!     ALLOCATE( new_ikq(nks), temp_index_ikq(nks) )
!     !
!     temp_xkq = kgrid%xp_cryst
!     temp_wkq = kgrid%wp
!     !
!     ALLOCATE( kqgrid%index_kq(nks,nqs) )
!     ALLOCATE( kqgrid%g0(3,nks,nqs) )
!     !
!     nkqs = 0
!     new_ikq = 0
!     kqgrid%g0 = 0._DP
!     !
!     DO ik = 1, nks
!        DO iq = 1, nqs
!           !
!           xkq(:) = kgrid%xp_cryst(:,ik) + DBLE(sig) * qgrid%xp_cryst(:,iq) 
!           !
!           xk_not_found = .TRUE.
!           !
!           DO ikq = 1, nks
!              IF ( xk_not_found ) THEN
!                 dxk(:) = xkq(:)-temp_xkq(:,ikq) - NINT( xkq(:)-temp_xkq(:,ikq) )
!                 IF ( ALL( ABS( dxk ) < eps ) ) THEN
!                    xk_not_found = .FALSE.
!                    IF ( new_ikq(ikq) == 0 ) THEN
!                       nkqs = nkqs + 1
!                       temp_index_ikq(nkqs) = ikq
!                       new_ikq(ikq) = nkqs
!                    ENDIF
!                    kqgrid%index_kq(ik,iq) = new_ikq(ikq)
!                    kqgrid%g0(:,ik,iq) = xkq(:) - temp_xkq(:,ikq)
!                    CALL cryst_to_cart(1, kqgrid%g0(1,ik,iq), bg, +1)
!                 ENDIF
!              ENDIF
!           ENDDO ! ikq
!           !
!           IF ( xk_not_found ) THEN
!              WRITE(*,*) ik, iq, nkqs
!              WRITE(*,*) xkq(:)
!              WRITE(message,'(a,a,a)')'k ',TRIM( csig ),' q does not belong to k-point grid '
!              CALL errore( 'kq_grid_init', TRIM( message ), (ik-1) * nqs + iq )
!           ENDIF
!           !
!        ENDDO ! iq
!        !
!     ENDDO ! ik
!     !
!     ALLOCATE( kqgrid%xp_cryst(3,nspin*nkqs) , kqgrid%xp_cart(3,nspin*nkqs), kqgrid%wp(nspin*nkqs) )
!     !
!     DO ik = 1, nkqs
!        ikq = temp_index_ikq(ik)
!        kqgrid%xp_cryst(:,ik) = temp_xkq(:,ikq)
!        kqgrid%wp(ik) = temp_wkq(ikq)
!     ENDDO
!     !
!     kqgrid%xp_cart = kqgrid%xp_cryst
!     CALL cryst_to_cart(nkqs, kqgrid%xp_cart, bg, +1)
!     !
!     IF ( nspin == 2 ) THEN
!        !
!        DO ik = 1, nks/2
!            DO iq =1, nqs
!              kqgrid%index_kq(nks/2+ik,iq) = kqgrid%index_kq(ik,iq) + nkqs
!            ENDDO
!        ENDDO
!        !
!        DO ikq = 1, nkqs
!           kqgrid%xp_cart(:,ikq + nkqs) = kqgrid%xp_cart(:,ikq)
!           kqgrid%xp_cryst(:,ikq + nkqs) = kqgrid%xp_cryst(:,ikq)
!        ENDDO
!        !
!        nkqs = 2*nkqs
!     ENDIF
!     !
!     kqgrid%nps = nkqs
!     !
!     DEALLOCATE( new_ikq, temp_index_ikq )
!     DEALLOCATE( temp_xkq, temp_wkq )
!      !
!   END SUBROUTINE
   !
   !
!  SUBROUTINE q_grid_init( qgrid, kgrid, k1grid )
!     !
!     USE cell_base,        ONLY : at, bg
!     USE klist,            ONLY : xk, wk, nkstot
!     USE pwcom,            ONLY : nspin
!     USE westcom,          ONLY : nq
!     !
!     IMPLICIT NONE
!     !
!     ! I/O
!     !
!     CLASS(bz_grid), INTENT(OUT) :: qgrid
!     CLASS(bz_grid), INTENT(IN) :: kgrid, k1grid
!     !
!     ! Workspace
!     !
!     INTEGER :: ik, ikk, iq, iqq, iq1, iq2, iq3, ipol, ig
!     INTEGER :: nks, nks1, nqs, temp_nqs
!     INTEGER, ALLOCATABLE :: new_iq(:), temp_index_iq(:)
!     REAL(DP) :: xq(3)
!     REAL(DP), ALLOCATABLE :: temp_xq(:,:), temp_wq(:)
!     REAL(DP) :: dxq(3), dq1, dq2, dq3
!     LOGICAL :: xq_not_found
!     CHARACTER(LEN=256) :: message
!     CHARACTER(LEN=1) :: csig
!     REAL(DP) :: eps = 1.D-8
!     !
!     nks = kgrid%nps
!     nks1 = k1grid%nps
!     temp_nqs = nq(1)*nq(2)*nq(3)
!     !
!     ALLOCATE( temp_xq(3,temp_nqs), temp_wq(temp_nqs) )
!     ALLOCATE( new_iq(temp_nqs), temp_index_iq(temp_nqs) )
!     !
!     dq1 = 1._DP / DBLE(nq(1))
!     dq2 = 1._DP / DBLE(nq(2))
!     dq3 = 1._DP / DBLE(nq(3))
!     iq = 0
!     DO iq1 = 1, nq(1)
!        DO iq2 = 1, nq(2)
!           DO iq3 = 1, nq(3)
!              iq = iq + 1
!              temp_xq(1,iq) = DBLE( iq1 - 1 ) * dq1
!              temp_xq(2,iq) = DBLE( iq2 - 1 ) * dq2
!              temp_xq(3,iq) = DBLE( iq3 - 1 ) * dq3
!           ENDDO
!        ENDDO
!     ENDDO
!     temp_wq = 1 / DBLE(nq(1)*nq(2)*nq(3))
!     !
!     ALLOCATE( qgrid%index_q(nks,nks) )
!     ALLOCATE( qgrid%g0(3,nks,nks) )
!     !
!     nqs = 0
!     new_iq = 0
!     qgrid%g0 = 0._DP
!     !
!     DO ik = 1, nks
!        !
!        DO ikk = 1, nks1
!           !
!           xq(:) = kgrid%xp_cryst(:,ik) - k1grid%xp_cryst(:,ikk) 
!           !
!           xq_not_found = .TRUE.
!           !
!           DO iq = 1, temp_nqs
!              IF ( xq_not_found ) THEN
!                 dxq(:) = xq(:)-temp_xq(:,iq) - NINT( xq(:)-temp_xq(:,iq) )
!                 IF ( ALL( ABS( dxq ) < eps ) ) THEN
!                    xq_not_found = .FALSE.
!                    IF ( new_iq(iq) == 0 ) THEN
!                       nqs = nqs + 1
!                       temp_index_iq(nqs) = iq
!                       new_iq(iq) = nqs
!                    ENDIF
!                    qgrid%index_q(ik,ikk) = new_iq(iq)
!                    qgrid%g0(:,ik,ikk) = xq(:) - temp_xq(:,iq)
!                    CALL cryst_to_cart(1, qgrid%g0(1,ik,ikk), bg, +1)
!                 ENDIF
!              ENDIF
!           ENDDO ! iq
!           !
!           IF ( xq_not_found ) THEN
!              WRITE(*,*) ik, ikk, nqs
!              WRITE(*,*) xq(:)
!              WRITE(message,'(a)') " k - k' is not a commensurate with k-point grid "
!              CALL errore( 'q_grid_init', TRIM( message ), (ik-1) * nqs + ikk )
!           ENDIF
!           !
!        ENDDO ! ikk
!        !
!     ENDDO ! ik
!     !
!     ALLOCATE( qgrid%xp_cryst(3,nqs) , qgrid%xp_cart(3,nqs), qgrid%wp(nqs) )
!     !
!     DO iqq = 1, nqs
!        iq = temp_index_iq(iqq)
!        qgrid%xp_cryst(:,iqq) = temp_xq(:,iq)
!        qgrid%wp(iqq) = temp_wq(iq)
!     ENDDO
!     !
!     qgrid%xp_cart = qgrid%xp_cryst
!     CALL cryst_to_cart(nqs, qgrid%xp_cart, bg, +1)
!     !
!     ALLOCATE( qgrid%l_gammap(nqs) )
!     qgrid%l_gammap(:)=.FALSE.
!     DO iq = 1, nqs
!        qgrid%l_gammap(iq) = ( ALL( ABS ( qgrid%xp_cryst(:,iq) ) .LT. eps ) )
!     ENDDO
!     !
!     qgrid%nps = nqs
!     !
!     DEALLOCATE( new_iq, temp_index_iq )
!     DEALLOCATE( temp_xq, temp_wq )
!     !
!  END SUBROUTINE
!  !
!  SUBROUTINE get_phase( grid_aux, ik, ikk )
!     !
!     USE cell_base,        ONLY : bg
!     USE gvecs,            ONLY : ngms, nls
!     USE gvect,            ONLY : g
!     USE fft_base,         ONLY : dffts
!     USE fft_interfaces,   ONLY : invfft
!     !
!     IMPLICIT NONE
!     !
!     ! I/O
!     !
!     CLASS(bz_grid), INTENT(INOUT) :: grid_aux
!     INTEGER, INTENT(IN) :: ik, ikk
!     !
!     ! Workspace
!     !
!     INTEGER :: ipol, ig, ig0
!     REAL(DP) :: eps = 1.D-8
!     !
!     ig0 = 0
!     DO ig = 1, ngms
!        IF ( ALL ( ABS( g(:,ig) - grid_aux%g0(:,ik,ikk) ) < eps ) ) THEN
!           ig0 = ig
!        ENDIF
!     ENDDO
!     !
!     IF (.NOT. ALLOCATED(grid_aux%phase) ) ALLOCATE ( grid_aux%phase(dffts%nnr) )
!     grid_aux%phase = (0._DP, 0._DP)
!     IF ( ig0 > 0 ) THEN
!        grid_aux%phase( nls(ig0) ) = (1._DP, 0._DP)
!     ENDIF
!     ! phase = exp(-iG_0*r)
!     CALL invfft( 'Wave', grid_aux%phase, dffts )
!     grid_aux%phase(1:dffts%nnr) = CONJG( grid_aux%phase(1:dffts%nnr) )
!     !
!  END SUBROUTINE
   ! 
END MODULE
