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
! He Ma, Marco Govoni, Han Yang
!
MODULE fourier_interpolation
  ! -----------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
  !
  SUBROUTINE ft_interpolate(fr, n1, n2, n3, fg, npw, npwx)
    !
    ! Fourier transform a set of functions f from R space (fr) to G space (fg)
    ! fr can be defined on arbitrary R space grid (n1, n2, n3), as long as it is smoother than dffts
    ! fg is defined on the internal G vectors of Quantum Espresso
    !
    ! An maximum ecut for will first be determined, fg has size npwx but all G components
    ! with energy larger than ecut will be 0
    !
    ! Gamma-trick is assumed, fr is real, k points are not supported yet
    !
    ! This subroutine should be called simultaneously within one image (band group), fr only need
    ! to be allocated on root, fg will be distributed within image (band group)
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : g
    USE control_flags, ONLY : gamma_only
    USE cell_base,     ONLY : at,bg,tpiba
    USE mp_images,     ONLY : me_image
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n1, n2, n3, npw, npwx
    REAL(DP), INTENT(IN) :: fr(n1*n2*n3)
    COMPLEX(DP), INTENT(OUT) :: fg(npwx)
    !
    INTEGER :: ig, ng
    INTEGER :: h, k, l, hmax, kmax, lmax, hidx, kidx, lidx
    REAL(DP) :: ecut, ecut1, ecut2, ecut3
    REAL(DP) :: g2
    INTEGER, ALLOCATABLE :: nls(:) ! mapping from [1, npw] to [1, n1*n2*n3]
    COMPLEX(DP), ALLOCATABLE :: auxr(:)
    !
    IF(.NOT. gamma_only) CALL errore('ft_interpolate', 'only Gamma point is supported', 1)
    !
    ALLOCATE(nls(npw))
    nls = 0
    !
    ALLOCATE(auxr(n1*n2*n3))
    !
    ! find maximum ecut in G space to represent R space grid (n1, n2, n3)
    !
    hmax = (n1 - 1) / 2
    kmax = (n2 - 1) / 2
    lmax = (n3 - 1) / 2
    ecut1 = 0.5_DP * hmax**2 * SUM((bg(:,1)*tpiba)**2)
    ecut2 = 0.5_DP * kmax**2 * SUM((bg(:,2)*tpiba)**2)
    ecut3 = 0.5_DP * lmax**2 * SUM((bg(:,3)*tpiba)**2)
    ecut = MIN(ecut1, ecut2, ecut3)
    !
    ng = 0
    DO ig = 1,npw
       !
       g2 = SUM(g(:,ig)**2) * tpiba**2
       !
       IF(g2 < 2*ecut) THEN
          !
          ng = ng + 1
          IF(ng > npw) CALL errore('ft_interpolate', 'size of fg is not enough to support fr', 1)
          !
          ! Compute Miller index, derive position of G vector in (n1, n2, n3) grid
          !
          h = NINT(SUM(g(:,ig) * at(:,1)))
          k = NINT(SUM(g(:,ig) * at(:,2)))
          l = NINT(SUM(g(:,ig) * at(:,3)))
          !
          IF(h >= 0) THEN
             hidx = h + 1
          ELSE
             hidx = h + n1 + 1
          ENDIF
          !
          IF(k >= 0) THEN
             kidx = k + 1
          ELSE
             kidx = k + n2 + 1
          ENDIF
          !
          IF(l >= 0) THEN
             lidx = l + 1
          ELSE
             lidx = l + n3 + 1
          ENDIF
          !
          nls(ng) = hidx + (kidx - 1) * n1 + (lidx - 1) * n1 * n2
          !
       ENDIF
       !
    ENDDO
    !
    fg(:) = (0._DP,0._DP)
    auxr(:) = (0._DP,0._DP)
    !
    IF(me_image == 0) THEN
       auxr(:) = CMPLX(fr,KIND=DP)
    ENDIF
    !
    CALL ft_interpolate_single_at_gamma(auxr, n1, n2, n3, npw, npwx, nls, fg)
    !
    DEALLOCATE(auxr)
    DEALLOCATE(nls)
    !
  END SUBROUTINE
  !
  SUBROUTINE ft_interpolate_single_at_gamma(fr, n1, n2, n3, npw, npwx, nls, fg)
    !
    USE kinds,         ONLY : DP
    USE fft_scalar,    ONLY : cfft3d
    USE mp_images,     ONLY : intra_image_comm,me_image
    USE mp,            ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n1, n2, n3, npw, npwx
    INTEGER, INTENT(IN) :: nls(npw)
    COMPLEX(DP), INTENT(IN) :: fr(n1*n2*n3)
    COMPLEX(DP),  INTENT(OUT) :: fg(npwx)
    !
    INTEGER :: ig, idx
    !
    IF(me_image == 0) THEN
       CALL cfft3d(fr, n1, n2, n3, n1, n2, n3, 1, -1)
    ENDIF
    !
    CALL mp_bcast(fr, 0, intra_image_comm)
    !
    DO ig = 1,npw
       idx = nls(ig)
       IF(idx > 0) fg(ig) = fr(nls(ig))
    ENDDO
    !
    DO ig = npw+1,npwx
       fg(ig) = (0._DP,0._DP)
    ENDDO
    !
  END SUBROUTINE
  !
END MODULE
