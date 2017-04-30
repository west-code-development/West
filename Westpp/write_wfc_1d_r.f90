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
! -------------------------------------------------------------------
SUBROUTINE write_wfc_1D_r ( dfft, iu, fname, wfc_distr, ipol)
  ! -------------------------------------------------------------------
  !
  USE pwcom,                 ONLY : npw,npwx
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : celldm, at, bg, omega
  USE ions_base,             ONLY : nat, tau, atm, ityp
  USE scatter_mod,           ONLY : gather_grid
  USE fft_types,             ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER,INTENT(IN) :: iu
  CHARACTER(LEN=*),INTENT(IN) :: fname
  REAL(DP),INTENT(IN) :: wfc_distr(dfft%nnr)
  INTEGER,INTENT(IN) :: ipol
  !
  ! Workspace
  ! 
  REAL(DP)         :: alat
  INTEGER          :: nr1, nr2, nr3, nr1x, nr2x, nr3x
  INTEGER :: ir,ntot
  REAL(DP) :: wfc_gat(dfft%nr1x*dfft%nr2x*dfft%nr3x)
  REAL(DP),ALLOCATABLE :: integral(:)
  REAL(DP) :: conversion_1, conversion_2
  !
  wfc_gat=0.0_DP
  CALL gather_grid(dfft,wfc_distr,wfc_gat)
  !
  alat = celldm(1)
  nr1 = dfft%nr1
  nr2 = dfft%nr2
  nr3 = dfft%nr3
  nr1x= dfft%nr1x
  nr2x= dfft%nr2x
  nr3x= dfft%nr3x
  !
  ! ipol
  !
  SELECT CASE(ipol)
  CASE(1)
     !
     ! X
     !
     ntot=nr1
     ALLOCATE(integral(ntot))
     CALL calc_int_x(wfc_gat,nr1,nr2,nr3,integral)
     IF( dfft%mype == dfft%root ) THEN 
        OPEN(UNIT=iu,FILE=TRIM(fname))
     ENDIF
     !
  CASE(2)
     !
     ! Y
     ! 
     ntot=nr2
     ALLOCATE(integral(ntot))
     CALL calc_int_y(wfc_gat,nr1,nr2,nr3,integral)
     IF( dfft%mype == dfft%root ) THEN 
        OPEN(UNIT=iu,FILE=TRIM(fname))
     ENDIF
     !
  CASE(3)
     !
     ! Z
     !
     ntot=nr3
     ALLOCATE(integral(ntot))
     CALL calc_int_z(wfc_gat,nr1,nr2,nr3,integral)
     IF( dfft%mype == dfft%root ) THEN 
        OPEN(UNIT=iu,FILE=TRIM(fname))
     ENDIF
     !
  CASE DEFAULT 
     !
     CALL errore('write_1D_r','ipol must be 1,2,3',ipol)
     !
  END SELECT
  !
  ! Write
  !
  conversion_1 = alat*at(ipol,ipol)/REAL(ntot,KIND=DP) 
  conversion_2 = 1._DP/REAL(nr1*nr2*nr3,KIND=DP) 
  !
  IF( dfft%mype == dfft%root ) THEN 
     WRITE(iu,'("#",a12,a13)') 'r[au]', 'f(r)'
     DO ir=1,ntot
        WRITE(iu,'(2es14.6)') REAL(ir-1,KIND=DP)*conversion_1, integral(ir)*conversion_2
     ENDDO
     WRITE(iu,'(2es14.6)') REAL(ntot,KIND=DP)*conversion_1, integral(1)*conversion_2
     CLOSE(iu)
  ENDIF
  !
  DEALLOCATE(integral)
  !
END SUBROUTINE
!
!
SUBROUTINE calc_int_x(func,nr1,nr2,nr3,integral) 
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: func(nr1,nr2,nr3)
  INTEGER :: nr1,nr2,nr3
  REAL(DP) :: integral(nr1)
  !
  INTEGER :: i1,i2,i3
  !
  integral = 0._DP
  DO i1=1,nr1
     DO i2=1,nr2
        DO i3=1,nr3
           !
           integral(i1) = integral(i1) + func(i1,i2,i3) 
           !
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
!
!
SUBROUTINE calc_int_y(func,nr1,nr2,nr3,integral) 
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: func(nr1,nr2,nr3)
  INTEGER :: nr1,nr2,nr3
  REAL(DP) :: integral(nr2)
  !
  INTEGER :: i1,i2,i3
  !
  integral = 0._DP
  DO i1=1,nr1
     DO i2=1,nr2
        DO i3=1,nr3
           !
           integral(i2) = integral(i2) + func(i1,i2,i3) 
           !
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
!
!
SUBROUTINE calc_int_z(func,nr1,nr2,nr3,integral) 
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: func(nr1,nr2,nr3)
  INTEGER :: nr1,nr2,nr3
  REAL(DP) :: integral(nr3)
  !
  INTEGER :: i1,i2,i3
  !
  integral = 0._DP
  DO i1=1,nr1
     DO i2=1,nr2
        DO i3=1,nr3
           !
           integral(i3) = integral(i3) + func(i1,i2,i3) 
           !
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE
