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
! Marco Govoni
!
!----------------------------------------------------------------------------
SUBROUTINE dump_r ( auxr, fname )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE westcom,               ONLY : westpp_format,westpp_r0, westpp_nr, westpp_rmax,npwq0,npwq0x,fftdriver, &
                                  & westpp_calculation
  USE fft_base,              ONLY : dffts
  USE cubefile,              ONLY : write_wfc_cube_r 
  USE fft_at_gamma,          ONLY : single_fwfft_gamma
  USE fft_at_k,              ONLY : single_fwfft_k
  USE pwcom,                 ONLY : npw, npwx
  USE control_flags,         ONLY : gamma_only 
  USE gvect,                 ONLY : ngm 
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP) :: auxr(dffts%nnr)
  CHARACTER(LEN=*),INTENT(IN) :: fname
  COMPLEX(DP),ALLOCATABLE :: auxg(:)
  COMPLEX(DP),ALLOCATABLE :: auxr_(:)
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
     CALL write_wfc_cube_r ( dffts, 2001, TRIM(fname)//".cube", auxr )    
  ENDIF
  !
  IF( lgate(2) ) THEN 
     CALL write_wfc_1D_r ( dffts, 2002, TRIM(fname)//".plavx", auxr, 1 )    
  ENDIF
  !
  IF( lgate(3) ) THEN 
     CALL write_wfc_1D_r ( dffts, 2003, TRIM(fname)//".plavy", auxr, 2 )    
  ENDIF
  !
  IF( lgate(4) ) THEN 
     CALL write_wfc_1D_r ( dffts, 2004, TRIM(fname)//".plavz", auxr, 3 )    
  ENDIF
  !
  IF( lgate(5) ) THEN 
     !
     ALLOCATE(auxg(ngm))
     ALLOCATE(auxr_(dffts%nnr))
     auxr_ = CMPLX( auxr, 0.d0, KIND = DP)
     IF( gamma_only ) THEN 
        CALL single_fwfft_gamma(dffts,ngm,ngm,auxr_,auxg,'Smooth')
     ELSE
        CALL single_fwfft_k(dffts,ngm,ngm,auxr_,auxg,'Smooth')
     ENDIF
     CALL write_wfc_spav ( 2005, TRIM(fname)//".spavr", auxg, westpp_r0, westpp_nr, westpp_rmax )
     DEALLOCATE(auxg)
     DEALLOCATE(auxr_)
     !
  ENDIF 
  !
END SUBROUTINE
