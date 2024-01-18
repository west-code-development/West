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
!----------------------------------------------------------------------------
SUBROUTINE dump_r ( auxr, fname )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE westcom,               ONLY : westpp_format,westpp_r0, westpp_nr, westpp_rmax
  USE fft_base,              ONLY : dffts
  USE cubefile,              ONLY : write_wfc_cube_r
  USE fft_at_gamma,          ONLY : single_fwfft_gamma
  USE fft_at_k,              ONLY : single_fwfft_k
  USE control_flags,         ONLY : gamma_only
  USE gvect,                 ONLY : ngm
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP),INTENT(IN) :: auxr(dffts%nnr)
  CHARACTER(LEN=*),INTENT(IN) :: fname
  !
  ! Workspace
  !
  COMPLEX(DP),ALLOCATABLE :: auxg(:)
  COMPLEX(DP),ALLOCATABLE :: auxr_(:)
  !$acc declare device_resident(auxr_)
  !
  ! Workspace
  !
  INTEGER :: i
  LOGICAL :: lgate(5)
  !
  lgate=.FALSE.
  DO i = 1, 7
     IF( westpp_format(i:i) == 'c' .OR. westpp_format(i:i) == 'C' ) lgate(1) = .TRUE. ! generate fname.cube
     IF( westpp_format(i:i) == 'x' .OR. westpp_format(i:i) == 'X' ) lgate(2) = .TRUE. ! generate fname.plavx
     IF( westpp_format(i:i) == 'y' .OR. westpp_format(i:i) == 'Y' ) lgate(3) = .TRUE. ! generate fname.plavy
     IF( westpp_format(i:i) == 'z' .OR. westpp_format(i:i) == 'Z' ) lgate(4) = .TRUE. ! generate fname.plavz
     IF( westpp_format(i:i) == 's' .OR. westpp_format(i:i) == 'S' ) lgate(5) = .TRUE. ! generate fname.spavr
  ENDDO
  !
  !
  IF( lgate(1) ) THEN
     CALL write_wfc_cube_r ( dffts, TRIM(fname)//'.cube', auxr )
  ENDIF
  !
  IF( lgate(2) ) THEN
     CALL write_wfc_1D_r ( dffts, TRIM(fname)//'.plavx', auxr, 1 )
  ENDIF
  !
  IF( lgate(3) ) THEN
     CALL write_wfc_1D_r ( dffts, TRIM(fname)//'.plavy', auxr, 2 )
  ENDIF
  !
  IF( lgate(4) ) THEN
     CALL write_wfc_1D_r ( dffts, TRIM(fname)//'.plavz', auxr, 3 )
  ENDIF
  !
  IF( lgate(5) ) THEN
     !
     ALLOCATE(auxg(ngm))
     !$acc enter data create(auxg)
     ALLOCATE(auxr_(dffts%nnr))
     !
     !$acc kernels present(auxr_,auxr)
     auxr_(:) = CMPLX( auxr, 0._DP, KIND = DP)
     !$acc end kernels
     !
     IF( gamma_only ) THEN
        CALL single_fwfft_gamma(dffts,ngm,ngm,auxr_,auxg,'Rho')
     ELSE
        CALL single_fwfft_k(dffts,ngm,ngm,auxr_,auxg,'Rho')
     ENDIF
     !
     !$acc update host(auxg)
     CALL write_wfc_spav ( TRIM(fname)//'.spavr', auxg, westpp_r0, westpp_nr, westpp_rmax )
     !
     !$acc exit data delete(auxg)
     DEALLOCATE(auxg)
     DEALLOCATE(auxr_)
     !
  ENDIF
  !
END SUBROUTINE
