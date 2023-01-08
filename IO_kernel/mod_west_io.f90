!
! Copyright (C) 2015-2022 M. Govoni
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
  ! a header consisting of HD_LENGTH=32 integers, followed by raw data.
  ! Currently only 3 integers are used in the header, storing:
  ! (1) the version identifier of the format
  ! (2) the endianness (0 for GE, 1 for LE)
  ! (3) the length of the raw data (number of COMPLEX DP entries)
  ! (4-32) not used (yet).
  !
  INTEGER, PARAMETER :: HD_LENGTH = 32
  INTEGER, PARAMETER :: HD_VERSION = 210405
  INTEGER, PARAMETER :: HD_ID_VERSION = 1
  INTEGER, PARAMETER :: HD_ID_LITTLE_ENDIAN = 2
  INTEGER, PARAMETER :: HD_ID_DIMENSION = 3
  !
#if defined(__SX6)
#  define DIRECT_IO_FACTOR 1
#else
#  define DIRECT_IO_FACTOR 8
#endif
  !
  INTERFACE serial_data_write
    MODULE PROCEDURE serial_i0_data_write, serial_i1_data_write, serial_i2_data_write, serial_i3_data_write, serial_i4_data_write, &
                     serial_d0_data_write, serial_d1_data_write, serial_d2_data_write, serial_d3_data_write, serial_d4_data_write, &
                     serial_z0_data_write, serial_z1_data_write, serial_z2_data_write, serial_z3_data_write, serial_z4_data_write
  END INTERFACE
  !
  INTERFACE serial_data_read
    MODULE PROCEDURE serial_i0_data_read, serial_i1_data_read, serial_i2_data_read, serial_i3_data_read, serial_i4_data_read, &
                     serial_d0_data_read, serial_d1_data_read, serial_d2_data_read, serial_d3_data_read, serial_d4_data_read, &
                     serial_z0_data_read, serial_z1_data_read, serial_z2_data_read, serial_z3_data_read, serial_z4_data_read
  END INTERFACE
  !
  ! Note:
  !
  ! Parallel I/O using MPI I/O may be found in the mod_mpiio module (mpiio.f90).
  ! Parallel I/O using MPI I/O has been removed from this file.
  ! See the git history of this file if needed.
  !
  CONTAINS
    !
    ! REMOVE FILE IF IT EXISTS
    !
    SUBROUTINE remove_if_present(fname)
      IMPLICIT NONE
      CHARACTER(*),INTENT(in) :: fname
      LOGICAL :: l_aux
      INTEGER :: iunit
      !
      INQUIRE(FILE=TRIM(ADJUSTL(fname)),EXIST=l_aux)
      !
      IF(l_aux) THEN
        OPEN(NEWUNIT=iunit,FILE=fname,STATUS='OLD')
        CLOSE(UNIT=iunit,STATUS='DELETE')
      ENDIF
      !
    END SUBROUTINE
    !
    ! CHECK IF FILE IS PRESENT
    !
    LOGICAL FUNCTION file_is_present(lproc,fname,suffix)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      CHARACTER(*),INTENT(IN) :: suffix
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
    ! WRITE I0
    !
    SUBROUTINE serial_i0_data_write(lproc,fname,i0dummy)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: i0dummy
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr)
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,'("Shape: (",a14,")")') 'scalar'
      WRITE(iunit,'(i14)',IOSTAT=ierr) i0dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE I1
    !
    SUBROUTINE serial_i1_data_write(lproc,fname,i1dummy,n)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      INTEGER,INTENT(IN) :: i1dummy(n)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr)
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,'("Shape: (",i14,")")') n
      DO i = 1,n
         WRITE(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i1dummy(i)
         IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE I2
    !
    SUBROUTINE serial_i2_data_write(lproc,fname,i2dummy,n,m)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      INTEGER,INTENT(IN) :: i2dummy(n,m)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i,j
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr)
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,'("Shape: (",i14,",",i14,")")') n,m
      DO j = 1,m
         DO i = 1,n
            WRITE(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i2dummy(i,j)
            IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
         ENDDO
         WRITE(iunit,*) ''
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE I3
    !
    SUBROUTINE serial_i3_data_write(lproc,fname,i3dummy,n,m,l)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      INTEGER,INTENT(IN) :: i3dummy(n,m,l)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i,j,k
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr)
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,'("Shape: (",i14,",",i14,",",i14,")")') n,m,l
      DO k = 1,l
         DO j = 1,m
            DO i = 1,n
               WRITE(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i3dummy(i,j,k)
               IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
            ENDDO
            WRITE(iunit,*) ''
         ENDDO
         WRITE(iunit,*) ''
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE I4
    !
    SUBROUTINE serial_i4_data_write(lproc,fname,i4dummy,n,m,l,q)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l,q
      INTEGER,INTENT(IN) :: i4dummy(n,m,l,q)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i,j,k,p
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr)
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,'("Shape: (",i14,",",i14,",",i14,",",i14,")")') n,m,l,q
      DO p = 1,q
         DO k = 1,l
            DO j = 1,m
               DO i = 1,n
                  WRITE(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i4dummy(i,j,k,p)
                  IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
               ENDDO
               WRITE(iunit,*) ''
            ENDDO
            WRITE(iunit,*) ''
         ENDDO
         WRITE(iunit,*) ''
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D0
    !
    SUBROUTINE serial_d0_data_write(lproc,fname,d0dummy)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      REAL(DP),INTENT(IN) :: d0dummy
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) d0dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D1
    !
    SUBROUTINE serial_d1_data_write(lproc,fname,d1dummy,n)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      REAL(DP),INTENT(IN) :: d1dummy(n)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) d1dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D2
    !
    SUBROUTINE serial_d2_data_write(lproc,fname,d2dummy,n,m)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      REAL(DP),INTENT(IN) :: d2dummy(n,m)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) d2dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D3
    !
    SUBROUTINE serial_d3_data_write(lproc,fname,d3dummy,n,m,l)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      REAL(DP),INTENT(IN) :: d3dummy(n,m,l)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) d3dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D4
    !
    SUBROUTINE serial_d4_data_write(lproc,fname,d4dummy,n,m,l,q)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l,q
      REAL(DP),INTENT(IN) :: d4dummy(n,m,l,q)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) d4dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z0
    !
    SUBROUTINE serial_z0_data_write(lproc,fname,z0dummy)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(IN) :: z0dummy
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) z0dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z1
    !
    SUBROUTINE serial_z1_data_write(lproc,fname,z1dummy,n)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      COMPLEX(DP),INTENT(IN) :: z1dummy(n)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) z1dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z2
    !
    SUBROUTINE serial_z2_data_write(lproc,fname,z2dummy,n,m)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      COMPLEX(DP),INTENT(IN) :: z2dummy(n,m)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) z2dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z3
    !
    SUBROUTINE serial_z3_data_write(lproc,fname,z3dummy,n,m,l)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      COMPLEX(DP),INTENT(IN) :: z3dummy(n,m,l)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) z3dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z4
    !
    SUBROUTINE serial_z4_data_write(lproc,fname,z4dummy,n,m,l,q)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l,q
      COMPLEX(DP),INTENT(IN) :: z4dummy(n,m,l,q)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(iunit,IOSTAT=ierr) z4dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! ************
    ! *** READ ***
    ! ************
    !
    ! READ I0
    !
    SUBROUTINE serial_i0_data_read(lproc,fname,i0dummy)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(OUT) :: i0dummy
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,*)
      READ(iunit,'(i14)',IOSTAT=ierr) i0dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ I1
    !
    SUBROUTINE serial_i1_data_read(lproc,fname,i1dummy,n)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      INTEGER,INTENT(OUT) :: i1dummy(n)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,*)
      DO i = 1,n
         READ(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i1dummy(i)
         IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ I2
    !
    SUBROUTINE serial_i2_data_read(lproc,fname,i2dummy,n,m)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      INTEGER,INTENT(OUT) :: i2dummy(n,m)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i,j
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,*)
      DO j = 1,m
         DO i = 1,n
            READ(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i2dummy(i,j)
            IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
         ENDDO
         READ(iunit,*)
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ I3
    !
    SUBROUTINE serial_i3_data_read(lproc,fname,i3dummy,n,m,l)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      INTEGER,INTENT(OUT) :: i3dummy(n,m,l)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i,j,k
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,*)
      DO k = 1,l
         DO j = 1,m
            DO i = 1,n
               READ(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i3dummy(i,j,k)
               IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
            ENDDO
            READ(iunit,*)
         ENDDO
         READ(iunit,*)
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ I4
    !
    SUBROUTINE serial_i4_data_read(lproc,fname,i4dummy,n,m,l,q)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l,q
      INTEGER,INTENT(OUT) :: i4dummy(n,m,l,q)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      INTEGER :: i,j,k,p
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,*)
      DO p = 1,q
         DO k = 1,l
            DO j = 1,m
               DO i = 1,n
                  READ(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i4dummy(i,j,k,p)
                  IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
               ENDDO
               READ(iunit,*)
            ENDDO
            READ(iunit,*)
         ENDDO
         READ(iunit,*)
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ D0
    !
    SUBROUTINE serial_d0_data_read(lproc,fname,d0dummy)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      REAL(DP),INTENT(OUT) :: d0dummy
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) d0dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ D1
    !
    SUBROUTINE serial_d1_data_read(lproc,fname,d1dummy,n)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      REAL(DP),INTENT(OUT) :: d1dummy(n)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) d1dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ D2
    !
    SUBROUTINE serial_d2_data_read(lproc,fname,d2dummy,n,m)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      REAL(DP),INTENT(OUT) :: d2dummy(n,m)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) d2dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ D3
    !
    SUBROUTINE serial_d3_data_read(lproc,fname,d3dummy,n,m,l)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      REAL(DP),INTENT(OUT) :: d3dummy(n,m,l)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) d3dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ D4
    !
    SUBROUTINE serial_d4_data_read(lproc,fname,d4dummy,n,m,l,q)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l,q
      REAL(DP),INTENT(OUT) :: d4dummy(n,m,l,q)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) d4dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ Z0
    !
    SUBROUTINE serial_z0_data_read(lproc,fname,z0dummy)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(OUT) :: z0dummy
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) z0dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ Z1
    !
    SUBROUTINE serial_z1_data_read(lproc,fname,z1dummy,n)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      COMPLEX(DP),INTENT(OUT) :: z1dummy(n)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) z1dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ Z2
    !
    SUBROUTINE serial_z2_data_read(lproc,fname,z2dummy,n,m)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      COMPLEX(DP),INTENT(OUT) :: z2dummy(n,m)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) z2dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ Z3
    !
    SUBROUTINE serial_z3_data_read(lproc,fname,z3dummy,n,m,l)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      COMPLEX(DP),INTENT(OUT) :: z3dummy(n,m,l)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) z3dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! READ Z4
    !
    SUBROUTINE serial_z4_data_read(lproc,fname,z4dummy,n,m,l,q)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l,q
      COMPLEX(DP),INTENT(OUT) :: z4dummy(n,m,l,q)
      !
      CHARACTER(512) :: fullname
      INTEGER :: iunit
      INTEGER :: ierr
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = TRIM(ADJUSTL(fname))//'.dat'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr,FORM='unformatted',ACCESS='stream',STATUS='old')
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      READ(iunit,IOSTAT=ierr) z4dummy
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot read file: '//TRIM(fullname),ABS(ierr))
      !
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
    ! ###############
    ! ## OUTPUT #####
    ! ###############
    !
    ! OUTPUT
    !
    SUBROUTINE serial_table_output(lproc,fname,d2dummy,nrow,ncol,header)
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nrow,ncol
      CHARACTER(*),INTENT(IN) :: header(ncol)
      REAL(DP),INTENT(IN) :: d2dummy(nrow,ncol)
      !
      CHARACTER(512) :: fullname
      CHARACTER(512) :: format_string
      INTEGER :: iunit
      INTEGER :: ierr,i,j
      REAL(DP) :: help(ncol)
      !
      IF(.NOT. lproc) RETURN
      !
      fullname = 'o-'//TRIM(ADJUSTL(fname))//'.tab'
      ierr = 0
      !
      OPEN(NEWUNIT=iunit,FILE=TRIM(fullname),IOSTAT=ierr)
      IF(ierr /= 0) CALL errore('WEST/IO','Cannot open file: '//TRIM(fullname),ABS(ierr))
      !
      WRITE(format_string,*) '("#",',ncol,'(a16))'
      WRITE(iunit,TRIM(format_string)) header
      WRITE(format_string,*) '(" ",',ncol,'(f16.6))'
      DO i = 1,nrow
         DO j = 1,ncol
            help(j) = d2dummy(i,j)
         ENDDO
         WRITE(iunit,TRIM(format_string),IOSTAT=ierr) help(1:ncol)
         IF(ierr /= 0) CALL errore('WEST/IO','Cannot write file: '//TRIM(fullname),ABS(ierr))
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
    END SUBROUTINE
    !
END MODULE
