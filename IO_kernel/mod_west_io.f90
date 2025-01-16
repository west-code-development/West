!
! Copyright (C) 2015-2025 M. Govoni
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
MODULE west_io
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! Base64 was changed to binary in order to improve I/O performance.
  !
  ! A simple format is used here:
  ! A header consisting of HD_LENGTH=32 integers, followed by raw data.
  ! Currently only 3 integers are used in the header, storing:
  ! (1) the version identifier of the format,
  ! (2) the endianness (0 for GE, 1 for LE),
  ! (3) the length of the raw data (number of COMPLEX DP entries),
  ! (4-32) not used yet.
  !
  INTEGER, PARAMETER :: HD_LENGTH = 32
  INTEGER, PARAMETER :: HD_VERSION = 210405
  INTEGER, PARAMETER :: HD_ID_VERSION = 1
  INTEGER, PARAMETER :: HD_ID_LITTLE_ENDIAN = 2
  INTEGER, PARAMETER :: HD_ID_DIMENSION = 3
  !
  INTERFACE serial_data_write
    MODULE PROCEDURE serial_d1_data_write, serial_d2_data_write, serial_d3_data_write, &
                     serial_d4_data_write, serial_d5_data_write, serial_d6_data_write, &
                     serial_z1_data_write, serial_z2_data_write, serial_z3_data_write, &
                     serial_z4_data_write, serial_z5_data_write, serial_z6_data_write
  END INTERFACE
  !
  INTERFACE serial_data_read
    MODULE PROCEDURE serial_d1_data_read, serial_d2_data_read, serial_d3_data_read, &
                     serial_d4_data_read, serial_d5_data_read, serial_d6_data_read, &
                     serial_z1_data_read, serial_z2_data_read, serial_z3_data_read, &
                     serial_z4_data_read, serial_z5_data_read, serial_z6_data_read
  END INTERFACE
  !
  CONTAINS
    !
    ! REMOVE FILE IF IT EXISTS
    !
    SUBROUTINE remove_if_present(fname)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(*),INTENT(in) :: fname
      !
      ! Workspace
      !
      LOGICAL :: l_aux
      INTEGER :: iun
      !
      INQUIRE(FILE=TRIM(ADJUSTL(fname)),EXIST=l_aux)
      !
      IF(l_aux) THEN
        OPEN(NEWUNIT=iun,FILE=fname,STATUS='OLD')
        CLOSE(UNIT=iun,STATUS='DELETE')
      ENDIF
      !
    END SUBROUTINE
    !
    ! CHECK IF FILE IS PRESENT
    !
    LOGICAL FUNCTION file_is_present(lproc,fname,suffix)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      CHARACTER(*),INTENT(IN) :: suffix
      !
      ! Workspace
      !
      LOGICAL :: l_aux
      !
      IF(.NOT. lproc) RETURN
      !
      INQUIRE(FILE=TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(suffix)),EXIST=l_aux)
      file_is_present = l_aux
      !
    END FUNCTION
    !
    ! *************
    ! *** WRITE ***
    ! *************
    !
    SUBROUTINE serial_d1_data_write(lproc,fname,dummy,d1)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1
      REAL(DP),INTENT(IN) :: dummy(d1)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d2_data_write(lproc,fname,dummy,d1,d2)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2
      REAL(DP),INTENT(IN) :: dummy(d1,d2)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d3_data_write(lproc,fname,dummy,d1,d2,d3)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3
      REAL(DP),INTENT(IN) :: dummy(d1,d2,d3)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d4_data_write(lproc,fname,dummy,d1,d2,d3,d4)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4
      REAL(DP),INTENT(IN) :: dummy(d1,d2,d3,d4)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d5_data_write(lproc,fname,dummy,d1,d2,d3,d4,d5)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5
      REAL(DP),INTENT(IN) :: dummy(d1,d2,d3,d4,d5)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d6_data_write(lproc,fname,dummy,d1,d2,d3,d4,d5,d6)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5,d6
      REAL(DP),INTENT(IN) :: dummy(d1,d2,d3,d4,d5,d6)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z1_data_write(lproc,fname,dummy,d1)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1
      COMPLEX(DP),INTENT(IN) :: dummy(d1)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z2_data_write(lproc,fname,dummy,d1,d2)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2
      COMPLEX(DP),INTENT(IN) :: dummy(d1,d2)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z3_data_write(lproc,fname,dummy,d1,d2,d3)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3
      COMPLEX(DP),INTENT(IN) :: dummy(d1,d2,d3)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z4_data_write(lproc,fname,dummy,d1,d2,d3,d4)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4
      COMPLEX(DP),INTENT(IN) :: dummy(d1,d2,d3,d4)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z5_data_write(lproc,fname,dummy,d1,d2,d3,d4,d5)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5
      COMPLEX(DP),INTENT(IN) :: dummy(d1,d2,d3,d4,d5)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z6_data_write(lproc,fname,dummy,d1,d2,d3,d4,d5,d6)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5,d6
      COMPLEX(DP),INTENT(IN) :: dummy(d1,d2,d3,d4,d5,d6)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! ************
    ! *** READ ***
    ! ************
    !
    SUBROUTINE serial_d1_data_read(lproc,fname,dummy,d1)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1
      REAL(DP),INTENT(OUT) :: dummy(d1)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d2_data_read(lproc,fname,dummy,d1,d2)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2
      REAL(DP),INTENT(OUT) :: dummy(d1,d2)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d3_data_read(lproc,fname,dummy,d1,d2,d3)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3
      REAL(DP),INTENT(OUT) :: dummy(d1,d2,d3)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d4_data_read(lproc,fname,dummy,d1,d2,d3,d4)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4
      REAL(DP),INTENT(OUT) :: dummy(d1,d2,d3,d4)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d5_data_read(lproc,fname,dummy,d1,d2,d3,d4,d5)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5
      REAL(DP),INTENT(OUT) :: dummy(d1,d2,d3,d4,d5)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_d6_data_read(lproc,fname,dummy,d1,d2,d3,d4,d5,d6)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5,d6
      REAL(DP),INTENT(OUT) :: dummy(d1,d2,d3,d4,d5,d6)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z1_data_read(lproc,fname,dummy,d1)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1
      COMPLEX(DP),INTENT(OUT) :: dummy(d1)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z2_data_read(lproc,fname,dummy,d1,d2)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2
      COMPLEX(DP),INTENT(OUT) :: dummy(d1,d2)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z3_data_read(lproc,fname,dummy,d1,d2,d3)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3
      COMPLEX(DP),INTENT(OUT) :: dummy(d1,d2,d3)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z4_data_read(lproc,fname,dummy,d1,d2,d3,d4)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4
      COMPLEX(DP),INTENT(OUT) :: dummy(d1,d2,d3,d4)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z5_data_read(lproc,fname,dummy,d1,d2,d3,d4,d5)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5
      COMPLEX(DP),INTENT(OUT) :: dummy(d1,d2,d3,d4,d5)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    SUBROUTINE serial_z6_data_read(lproc,fname,dummy,d1,d2,d3,d4,d5,d6)
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: d1,d2,d3,d4,d5,d6
      COMPLEX(DP),INTENT(OUT) :: dummy(d1,d2,d3,d4,d5,d6)
      !
      ! Workspace
      !
      CHARACTER(512) :: fullname
      INTEGER :: iun
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iun,FILE=TRIM(fullname),IOSTAT=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iun,IOSTAT=ierr) dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iun,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
END MODULE
