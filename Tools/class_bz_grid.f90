!
! Copyright (C) 2015-2021 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `LICENSE'
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
   USE kinds,               ONLY : DP
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   TYPE, PUBLIC :: bz_grid
      !
      INTEGER :: ngrid(3) = (/ 1, 1, 1 /)       ! number of points in each direction
      INTEGER :: np = 1                         ! total number of points
      INTEGER :: ns = 1                         ! total number of spin = nspin_lsda
      INTEGER :: nps = 1                        ! total number of points and spins = np * ns
      INTEGER, ALLOCATABLE :: ip(:)             ! given ips --> ip
      INTEGER, ALLOCATABLE :: is(:)             ! given ips --> is
      REAL(DP), ALLOCATABLE :: p_cryst(:,:)     ! coordinates of point p in crystal                [ 1:np  ]
      REAL(DP), ALLOCATABLE :: p_cart(:,:)      ! coordinates of point p in cart ( tpiba units )   [ 1:np  ]
      REAL(DP), ALLOCATABLE :: weight(:)        ! weight of point p (sum of weights = nspin)       [ 1:nps ]
      LOGICAL, ALLOCATABLE :: l_pIsGamma(:)     ! .true. if point p = (0,0,0), else .false.        [ 1:np  ]
      !
      CONTAINS
      !
      PROCEDURE :: init => k_or_q_grid_init
      PROCEDURE :: find => findp
      PROCEDURE :: ipis2ips => from_ip_and_is_to_ips
      !
   END TYPE bz_grid
   !
   CONTAINS
   !
   SUBROUTINE k_or_q_grid_init( this, grid_type )
      !
      ! ... grid_type = [ 'K', 'Q' ]
      !
      USE cell_base,        ONLY : at, bg
      USE klist,            ONLY : xk, wk, nkstot
      USE start_k,          ONLY : nk1, nk2, nk3
      USE noncollin_module, ONLY : nspin_lsda
      USE control_flags,    ONLY : gamma_only
      USE constants,        ONLY : eps8
      USE westcom,          ONLY : qlist
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(INOUT) :: this
      CHARACTER(1), INTENT(IN) :: grid_type
      !
      ! Workspace
      !
      INTEGER :: ip, iq1, iq2, iq3, ips
      INTEGER :: i, j, k
      INTEGER :: iqlist
      !
      SELECT CASE( grid_type )
      CASE ( 'K', 'k')
         !
         ! This is a workaround to prevent ngrid(:) to be set to (/ 0, 0, 0 /) in the gamma_only case (espresso default)
         IF ( .NOT. gamma_only ) this%ngrid(1:3) = (/ nk1, nk2, nk3 /)
         this%np = this%ngrid(1) * this%ngrid(2) * this%ngrid(3)
         this%ns = nspin_lsda ! = 1 if nspin = 1 (unpolarized) or nspin = 4 (noncollinear)
         !                      = 2 if nspin = 2 (collinear)
         this%nps = nkstot    ! = np * ns
         IF ( this%nps /= this%np * this%ns ) CALL errore( 'types_bz_grid', 'nps /= np * ns', 1 )
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
         CALL cryst_to_cart( this%np, this%p_cryst, at, -1 )
         !
         ! set weights
         !
         ALLOCATE ( this%weight (this%nps) )
         DO ips = 1, this%nps
           this%weight(ips) = wk(ips)
         ENDDO
         !
      CASE ( 'Q', 'q')
         !
         IF ( .NOT. gamma_only ) this%ngrid(1:3) = (/ nk1, nk2, nk3 /)
         this%np = SIZE(qlist)
         this%ns = 1
         this%nps = this%np
         !
         ! generate p-vectors in cryst
         !
         ALLOCATE ( this%p_cryst  (3,this%np) )
         iqlist = 0
         ip = 0
         DO iq1 = 1, this%ngrid(1)
            DO iq2 = 1, this%ngrid(2)
               DO iq3 = 1, this%ngrid(3)
                  ip = ip + 1
                  IF ( ANY(qlist(:) == ip) ) THEN
                     iqlist = iqlist + 1
                     this%p_cryst(1,iqlist) = REAL( iq1 - 1, KIND=DP ) / REAL( this%ngrid(1), KIND=DP )
                     this%p_cryst(2,iqlist) = REAL( iq2 - 1, KIND=DP ) / REAL( this%ngrid(2), KIND=DP )
                     this%p_cryst(3,iqlist) = REAL( iq3 - 1, KIND=DP ) / REAL( this%ngrid(3), KIND=DP )
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         !
         ! generate p-vectors in cart
         !
         ALLOCATE ( this%p_cart (3,this%np) )
         this%p_cart(:,:) = this%p_cryst(:,:)
         CALL cryst_to_cart( this%np, this%p_cart, bg, +1 )
         !
         ! set weights
         !
         ALLOCATE ( this%weight (this%np)   )
         !
         this%weight = 1._DP / REAL(this%np,KIND=DP)
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
      this%l_pIsGamma(:) = .FALSE.
      DO ip = 1, this%np
         this%l_pIsGamma(ip) = ( ALL( ABS ( this%p_cryst(:,ip) ) .LT. eps8 ) )
      ENDDO
      !
   END SUBROUTINE
   !
   SUBROUTINE findp( this, p, unit_type, ip, g0 )
      !
      ! ... ip is the index of p (unit_type = [ 'cryst', 'cart'])
      ! ... if on exit ip == 0 --> p is not commensurate with this grid
      ! ... g0 relates p to an equivalent vector inside the 1BZ
      !
      USE constants,        ONLY : eps8
      USE cell_base,        ONLY : at, bg
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(IN) :: this
      REAL(DP), INTENT(IN) :: p(3)
      CHARACTER(*), INTENT(IN) :: unit_type
      INTEGER, INTENT(OUT) :: ip
      REAL(DP), INTENT(OUT) :: g0(3)
      !
      ! Workspace
      !
      INTEGER :: i
      REAL(DP) :: deltap(3)
      !
      SELECT CASE(unit_type)
      CASE('cryst','cart')
      CASE DEFAULT
         CALL errore( 'types_bz_grid', 'unit_type not supported, supported only cryst or cart', 1 )
      END SELECT
      !
      ! The search must be performed in crystalline coordinates
      !
      IF ( unit_type == 'cart' ) CALL cryst_to_cart( 1, p, at, -1 )
      !
      ip = 0
      DO i = 1, this%np
         deltap(:) = p(:) - this%p_cryst(:,i) - NINT( p(:) - this%p_cryst(:,i) )
         IF ( ALL ( ABS ( deltap ) .LT. eps8 ) ) THEN
            g0(:) = p(:) - this%p_cryst(:,i)
            ip = i
            EXIT
         ENDIF
      ENDDO
      !
      ! Tranform g0 back to cartesian coordinates if needed
      !
      IF ( unit_type == 'cart' ) CALL cryst_to_cart( 1, g0, bg, 1 )
      !
   END SUBROUTINE
   !
   FUNCTION from_ip_and_is_to_ips(this,ip,is) RESULT(ips)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CLASS(bz_grid), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: ip, is
      INTEGER :: ips
      !
      ips = ip + (is-1) * this%np
      !
   END FUNCTION
   !
END MODULE
