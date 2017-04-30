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
MODULE cubefile 
 ! -----------------------------------------------------------------
 !
 IMPLICIT NONE
 !
 CONTAINS
 ! 
 !-----------------------------------------------------------------
   SUBROUTINE write_wfc_cube_r ( dfft, iu, fname, wfc_distr )
   ! -----------------------------------------------------------------
   !
   USE kinds,                 ONLY : DP
   USE cell_base,             ONLY : celldm, at, bg
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
   ! 
   ! Workspace
   !
   REAL(DP)         :: alat
   INTEGER          :: i, nt, i1, i2, i3, at_num, ir
   INTEGER, EXTERNAL  :: atomic_number
   REAL(DP)    :: at_chrg, tpos(3), inpos(3)
   REAL(DP) :: wfc_gat(dfft%nr1x*dfft%nr2x*dfft%nr3x)
   !
   wfc_gat=0.0_DP
   CALL gather_grid(dfft,wfc_distr,wfc_gat)
   !
   !      WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
   !      TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
   !      THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
   ! 
   !      LINE   FORMAT      CONTENTS
   !      ===============================================================
   !       1     A           TITLE
   !       2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
   !       3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
   !       4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
   !       #ATOMS LINES OF ATOM COORDINATES:
   !       ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
   !       REST: 6E13.5      CUBE DATA
   ! 
   !      ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
   !
   alat = celldm(1)
   !
   IF( dfft%mype == dfft%root ) THEN
      !
      OPEN(UNIT=iu,FILE=TRIM(ADJUSTL(fname)))
      !
      WRITE(iu,*) 'Cubfile created from WEST calculation'
      WRITE(iu,*) 'Comment'
      !                        origin is forced to (0.0,0.0,0.0)
      WRITE(iu,'(I5,3F12.6)') nat, 0.0_DP, 0.0_DP, 0.0_DP
      WRITE(iu,'(I5,3F12.6)') dfft%nr1, (alat*at(i,1)/DBLE(dfft%nr1),i=1,3)
      WRITE(iu,'(I5,3F12.6)') dfft%nr2, (alat*at(i,2)/DBLE(dfft%nr2),i=1,3)
      WRITE(iu,'(I5,3F12.6)') dfft%nr3, (alat*at(i,3)/DBLE(dfft%nr3),i=1,3)
      !
      DO i=1,nat
         !
         nt = ityp(i)
         ! find atomic number for this atom.
         at_num = atomic_number(TRIM(atm(nt)))
         at_chrg= DBLE(at_num)
         ! at_chrg could be alternatively set to valence charge
         ! positions are in cartesian coordinates and a.u.
         !
         ! wrap coordinates back into cell.
         tpos = MATMUL( TRANSPOSE(bg), tau(:,i) )
         tpos = tpos - NINT(tpos - 0.5_DP)
         inpos = alat * MATMUL( at, tpos )
         WRITE(iu,'(I5,5F12.6)') at_num, at_chrg, inpos
         !
      ENDDO
      !
      CALL actual_write_cube(iu,dfft%nr1,dfft%nr2,dfft%nr3,wfc_gat)
      !
      CLOSE(iu)
      !
   ENDIF
   !
 END SUBROUTINE
 ! 
 !-----------------------------------------------------------------
   SUBROUTINE read_wfc_cube_r ( dfft, iu, fname, wfc_distr )
   ! -----------------------------------------------------------------
   !
   USE kinds,                 ONLY : DP
   USE cell_base,             ONLY : celldm, at, bg
   USE ions_base,             ONLY : nat, tau, atm, ityp
   USE scatter_mod,           ONLY : scatter_grid
   USE fft_types,             ONLY : fft_type_descriptor
   !
   IMPLICIT NONE
   !
   ! I/O 
   !
   TYPE(fft_type_descriptor), INTENT(IN) :: dfft
   INTEGER,INTENT(IN) :: iu
   CHARACTER(LEN=*),INTENT(IN) :: fname
   REAL(DP),INTENT(OUT) :: wfc_distr(dfft%nnr)
   ! 
   ! Workspace
   !
   REAL(DP)         :: alat
   INTEGER          :: i, nt, i1, i2, i3, at_num, ir
   INTEGER, EXTERNAL  :: atomic_number
   REAL(DP)    :: at_chrg, tpos(3), inpos(3)
   REAL(DP) :: wfc_gat(dfft%nr1x*dfft%nr2x*dfft%nr3x)
   !
   !
   !      WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
   !      TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
   !      THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
   ! 
   !      LINE   FORMAT      CONTENTS
   !      ===============================================================
   !       1     A           TITLE
   !       2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
   !       3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
   !       4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
   !       #ATOMS LINES OF ATOM COORDINATES:
   !       ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
   !       REST: 6E13.5      CUBE DATA
   ! 
   !      ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
   !
   wfc_gat = 0._DP
   wfc_distr=0.0_DP
   !
   alat = celldm(1)
   !
   IF( dfft%mype == dfft%root ) THEN
      !
      OPEN(UNIT=iu,FILE=TRIM(ADJUSTL(fname)))
      !
      READ(iu,*) 
      READ(iu,*) 
      READ(iu,*) 
      READ(iu,*) 
      READ(iu,*) 
      READ(iu,*) 
      DO i=1,nat
         READ(iu,*) 
      ENDDO
      !
      CALL actual_read_cube(iu,dfft%nr1,dfft%nr2,dfft%nr3,wfc_gat)
      !
      CLOSE(iu)
      !
   ENDIF
   !
   CALL scatter_grid(dfft,wfc_gat,wfc_distr)
   !
 END SUBROUTINE
 !
 !
 SUBROUTINE actual_write_cube(iu,nr1,nr2,nr3,func) 
   !
   USE kinds,  ONLY : DP
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER,INTENT(IN) :: iu,nr1,nr2,nr3
   REAL(DP),INTENT(IN) :: func(nr1,nr2,nr3)
   !
   ! Workspace
   !
   INTEGER :: i1,i2,i3
   !
   DO i1=1,nr1
      DO i2=1,nr2
         !WRITE(iu,'(6E13.5)') (func(i1,i2,i3),i3=1,nr3)
         DO i3=1,nr3
            WRITE(iu,'(1E13.5)') func(i1,i2,i3)
         ENDDO
      ENDDO
   ENDDO
   !
 END SUBROUTINE
 !
 !
 SUBROUTINE actual_read_cube(iu,nr1,nr2,nr3,func) 
   !
   USE kinds,  ONLY : DP
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER,INTENT(IN) :: iu,nr1,nr2,nr3
   REAL(DP),INTENT(OUT) :: func(nr1,nr2,nr3)
   !
   ! Workspace
   !
   INTEGER :: i1,i2,i3
   !
   DO i1=1,nr1
      DO i2=1,nr2
         !READ(iu,'(6E13.5)') (func(i1,i2,i3),i3=1,nr3)
         DO i3=1,nr3
            READ(iu,*) func(i1,i2,i3)
         ENDDO
      ENDDO
   ENDDO
   !
 END SUBROUTINE
 !
END MODULE
