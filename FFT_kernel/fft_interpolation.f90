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
! He Ma, Marco Govoni
!
MODULE fourier_interpolation
 ! -----------------------------------------------------------------
 !
 IMPLICIT NONE
 !
 CONTAINS
 !
 SUBROUTINE get_G2R_mapping (n1, n2, n3, n, nl, ndim)
    ! 
    ! INPUT  : n1, n2, n3 = dimension of arbitrary R grid 
    !          n          = actual number of PW
    !          ndim       - dimensions of nl, 1: (n), 2:(n,-n)
    ! OUTPUT : nl         = computed mapping from G to R space (i,e. from [1,n] to [1, n1*n2*n3] )
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : g
    USE cell_base,     ONLY : at, bg, tpiba
    !
    IMPLICIT NONE
    ! 
    ! I/O
    !
    INTEGER, INTENT(IN)  :: n1, n2, n3, n, ndim
    INTEGER, INTENT(OUT) :: nl(n,ndim)        ! mapping from [1, npw] to [1, n1*n2*n3]  
    !
    ! Workspace
    !
    INTEGER  :: ig, i_dim
    INTEGER  :: h, k, l, hmax, kmax, lmax, hidx, kidx, lidx
    !
    nl = 0
    !
    ! Maximum Miller indexes associated to the R grid 
    !
    hmax = (n1 - 1) / 2
    kmax = (n2 - 1) / 2
    lmax = (n3 - 1) / 2
    !
    DO i_dim = 1, ndim 
       DO ig = 1, n
          !
          ! Get the current Miller index
          !
          SELECT CASE (i_dim)
          CASE(1)  
             h = NINT(SUM(g(:,ig) * at(:,1)))
             k = NINT(SUM(g(:,ig) * at(:,2)))
             l = NINT(SUM(g(:,ig) * at(:,3)))
          CASE(2)
             h = NINT(SUM(-g(:,ig) * at(:,1)))
             k = NINT(SUM(-g(:,ig) * at(:,2)))
             l = NINT(SUM(-g(:,ig) * at(:,3)))
          CASE DEFAULT
             CALL errore("","Bad mapping given: only 1 .OR. 2",1) 
          END SELECT
          !
          IF( ABS(h) <= hmax .AND. ABS(k) <= kmax .AND. ABS(l) <= lmax ) THEN 
             !
             ! derive position of G vector in (n1, n2, n3) grid
             !
             IF ( h >= 0 ) THEN
                hidx = h + 1
             ELSE
                hidx = h + n1 + 1
             ENDIF
             !
             IF ( k >= 0 ) THEN
                kidx = k + 1
             ELSE
                kidx = k + n2 + 1
             ENDIF
             !
             IF ( l >= 0 ) THEN
                lidx = l + 1
             ELSE
                lidx = l + n3 + 1
             ENDIF
             !
             nl(ig,i_dim) = hidx + (kidx - 1) * n1 + (lidx - 1) * n1 * n2
             !
          ENDIF
          !
       ENDDO
    ENDDO
    !
 END SUBROUTINE
 !
 !
 SUBROUTINE single_fwfft_fromArbitraryRGrid (fr, n1, n2, n3, n, nx, ndim, nl, fg, igk)
   !
   ! Note that the interpolation needs to be called by all processors within one band group 
   !
   ! FWFFT : R ---> G
   !
   ! INPUT  : n1, n2, n3 = dimension of arbitrary R grid 
   !          n     = actual number of PW
   !          nx    = leading dimendion for fg
   !          ndim  = 1, 2
   !          fr    = ONE COMPLEX array containing ONE function in R space (note that the array is not distributed, i.e. dimension = n1*n2*n3 )
   !          nl    = pre-computed mapping from G to R space (i,e. from [1,n] to [1, n1*n2*n3] )
   ! OUTPUT : fg    = ONE COMPLEX array containing ONE functions in G space (note that the array is distributed ) 
   !
   USE kinds,         ONLY : DP
   USE fft_scalar,    ONLY : cfft3d
   USE mp_bands,      ONLY : intra_bgrp_comm, me_bgrp
   USE mp,            ONLY : mp_bcast
   !
   ! I/O
   !
   INTEGER,     INTENT(IN)       :: n1, n2, n3, n, nx, ndim
   INTEGER,     INTENT(IN)       :: nl(n,ndim)
   COMPLEX(DP), INTENT(IN)       :: fr(n1*n2*n3)
   COMPLEX(DP), INTENT(OUT)      :: fg(nx)
   INTEGER, INTENT(IN), OPTIONAL :: igk(n)
   !
   ! Workspace
   !
   INTEGER :: ig, idx
   !
   ! FFT is serial...
   !
   IF ( me_bgrp == 0 ) THEN
      !
      CALL cfft3d( fr, n1, n2, n3, n1, n2, n3, 1, -1)
      !
   ENDIF
   !
   ! ... result is broadcast
   !
   CALL mp_bcast( fr, 0, intra_bgrp_comm )
   !
   IF( PRESENT(igk) ) THEN
      DO ig=1, n
         !
         idx = nl(igk(ig),1)
         IF (idx > 0) fg(ig) = fr(nl(igk(ig),1))
         !
      ENDDO
   ELSE
      DO ig=1, n
         !
         idx = nl(ig,1)
         IF (idx > 0) fg(ig) = fr(nl(ig,1))
         !
      ENDDO
   ENDIF
   !
   DO ig = (n+1), nx
      !
      fg(ig) = (0.0_DP,0.0_DP)
      !
   ENDDO
   !
 ENDSUBROUTINE
 !
 !
 SUBROUTINE single_invfft_toArbitraryRGrid (fr, n1, n2, n3, n, nx, ndim, nl, fg, igk)
   !
   ! Note that the interpolation needs to be called by all processors within one band group 
   !
   ! FWFFT : R ---> G
   !
   ! INPUT  : n1, n2, n3 = dimension of arbitrary R grid 
   !          n     = actual number of PW
   !          nx    = leading dimendion for fg
   !          ndim  = 1,2
   !          fr    = ONE COMPLEX array containing ONE function in R space (note that the array is not distributed, i.e. dimension = n1*n2*n3 )
   !          nl    = pre-computed mapping from G to R space (i,e. from [1,n] to [1, n1*n2*n3] )
   ! OUTPUT : fg    = ONE COMPLEX array containing ONE functions in G space (note that the array is distributed ) 
   !
   USE kinds,         ONLY : DP
   USE fft_scalar,    ONLY : cfft3d
   USE mp_bands,      ONLY : intra_bgrp_comm, me_bgrp
   USE mp,            ONLY : mp_bcast
   !
   ! I/O
   !
   INTEGER,     INTENT(IN) :: n1, n2, n3, n, nx, ndim
   INTEGER,     INTENT(IN) :: nl(n,ndim)
   COMPLEX(DP), INTENT(OUT) :: fr(n1*n2*n3)
   COMPLEX(DP), INTENT(IN):: fg(nx)
   INTEGER, INTENT(IN), OPTIONAL :: igk(n)
   !
   ! Workspace
   !
   INTEGER :: ig, idx
   !
   fr = 0._DP
   !
   IF( PRESENT(igk) ) THEN
      DO ig=1, n
         !
         idx = nl(igk(ig),1)
         IF (idx > 0) THEN 
            fr(nl(igk(ig),1)) = fg(ig)
         ENDIF
         !
      ENDDO
   ELSE 
      DO ig=1, n
         !
         idx = nl(ig,1)
         IF (idx > 0) THEN 
            fr(nl(ig,1)) = fg(ig)
         ENDIF
         !
      ENDDO
      IF ( ndim == 2 ) THEN
         DO ig=1, n
            !
            idx = nl(ig,2)
            IF (idx > 0) THEN 
               fr(nl(ig,2)) = DCONJG(fg(ig))
            ENDIF
            !
         ENDDO
      ENDIF
   ENDIF
   !
   ! ... result is gathered
   !
   CALL mp_sum( fr, intra_bgrp_comm )
   !
   ! FFT is serial...
   !
   IF ( me_bgrp == 0 ) THEN
      !
      CALL cfft3d( fr, n1, n2, n3, n1, n2, n3, 1, 1)
      !
   ENDIF
   !
 ENDSUBROUTINE
 !
END MODULE
