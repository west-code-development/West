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
!-----------------------------------------------------------------------
MODULE west_io
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
#if defined(__SX6)
#  define DIRECT_IO_FACTOR 1
#else
#  define DIRECT_IO_FACTOR 8 
#endif
  !
  INTERFACE serial_data_write
    MODULE PROCEDURE serial_i0_data_write, serial_i1_data_write, serial_i2_data_write, serial_i3_data_write, & 
                     serial_d0_data_write, serial_d1_data_write, serial_d2_data_write, serial_d3_data_write, &
                     serial_z0_data_write, serial_z1_data_write, serial_z2_data_write, serial_z3_data_write
  END INTERFACE
  !
  INTERFACE serial_data_read
    MODULE PROCEDURE serial_i0_data_read, serial_i1_data_read, serial_i2_data_read, serial_i3_data_read, & 
                     serial_d0_data_read, serial_d1_data_read, serial_d2_data_read, serial_d3_data_read, &
                     serial_z0_data_read, serial_z1_data_read, serial_z2_data_read, serial_z3_data_read
  END INTERFACE
  !
  INTERFACE parallel_data_write
    MODULE PROCEDURE parallel_i1_data_write, parallel_d1_data_write, parallel_z1_data_write
  END INTERFACE
  !
  INTERFACE parallel_data_read
    MODULE PROCEDURE parallel_i1_data_read, parallel_d1_data_read, parallel_z1_data_read
  END INTERFACE
  !
  INTERFACE parallel_irrdata_write
    MODULE PROCEDURE parallel_i1_irrdata_write, parallel_d1_irrdata_write, parallel_z1_irrdata_write
  END INTERFACE
  !
  INTERFACE parallel_irrdata_read
    MODULE PROCEDURE parallel_i1_irrdata_read, parallel_d1_irrdata_read, parallel_z1_irrdata_read
  END INTERFACE
  !
  CONTAINS
    !
    ! CHECK IF FILE IS PRESENT
    !
    LOGICAL FUNCTION file_is_present(lproc,fname,suffix)
      LOGICAL,INTENT(IN) :: lproc
      CHARACTER(*),INTENT(IN) :: fname
      CHARACTER(*),INTENT(IN) :: suffix
      LOGICAL :: l_aux
      !
      IF(.NOT.lproc) RETURN
      !
      INQUIRE(FILE=TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(suffix)),EXIST=l_aux)
      file_is_present=l_aux
      !
    END FUNCTION
    !
    ! *************
    ! *************
    ! *** WRITE ***
    ! *************
    ! *************
    !
    ! WRITE I0
    !
    SUBROUTINE serial_i0_data_write(lproc,iunit,fname,i0dummy)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: i0dummy
      !
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      WRITE(iunit,'("Shape: (",a14,")")') 'scalar'
      WRITE(iunit,'(i14)',IOSTAT=ierr) i0dummy
      CLOSE(iunit,IOSTAT=ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE I1
    !
    SUBROUTINE serial_i1_data_write(lproc,iunit,fname,i1dummy,n)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      INTEGER,INTENT(IN) :: i1dummy(n)
      !
      INTEGER :: ierr
      INTEGER :: i
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      WRITE(iunit,'("Shape: (",i14,")")') n
      DO i=1,n
         WRITE(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i1dummy(i)
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE I2
    !
    SUBROUTINE serial_i2_data_write(lproc,iunit,fname,i2dummy,n,m)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      INTEGER,INTENT(IN) :: i2dummy(n,m)
      !
      INTEGER :: ierr
      INTEGER :: i,j
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      WRITE(iunit,'("Shape: (",i14,",",i14,")")') n,m
      DO j=1,m
         DO i=1,n
            WRITE(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i2dummy(i,j)
         ENDDO
         WRITE(iunit,*) ''
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE I3
    !
    SUBROUTINE serial_i3_data_write(lproc,iunit,fname,i3dummy,n,m,l)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      INTEGER,INTENT(IN) :: i3dummy(n,m,l)
      !
      INTEGER :: ierr
      INTEGER :: i,j,k
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      WRITE(iunit,'("Shape: (",i14,",",i14,",",i14,")")') n,m,l
      DO k=1,l
         DO j=1,m
            DO i=1,n
               WRITE(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i3dummy(i,j,k)
            ENDDO
            WRITE(iunit,*) ''
         ENDDO
         WRITE(iunit,*) ''
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D0
    !
    SUBROUTINE serial_d0_data_write(lproc,iunit,fname,d0dummy)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      REAL(DP),INTENT(IN) :: d0dummy
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(1, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) d0dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D1
    !
    SUBROUTINE serial_d1_data_write(lproc,iunit,fname,d1dummy,n)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      REAL(DP),INTENT(IN) :: d1dummy(n)
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(n, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) d1dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D2
    !
    SUBROUTINE serial_d2_data_write(lproc,iunit,fname,d2dummy,n,m)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      REAL(DP),INTENT(IN) :: d2dummy(n,m)
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(n*m, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) d2dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE D3
    !
    SUBROUTINE serial_d3_data_write(lproc,iunit,fname,d3dummy,n,m,l)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      REAL(DP),INTENT(IN) :: d3dummy(n,m,l)
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(n*m*l, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) d3dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z0 
    !
    SUBROUTINE serial_z0_data_write(lproc,iunit,fname,z0dummy)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(IN),TARGET :: z0dummy
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) z0dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z1
    !
    SUBROUTINE serial_z1_data_write(lproc,iunit,fname,z1dummy,n)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      COMPLEX(DP),INTENT(IN),TARGET :: z1dummy(n)
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2*n, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) z1dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z2 
    !
    SUBROUTINE serial_z2_data_write(lproc,iunit,fname,z2dummy,n,m)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      COMPLEX(DP),INTENT(IN),TARGET :: z2dummy(n,m)
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2*n*m, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) z2dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE Z3 
    !
    SUBROUTINE serial_z3_data_write(lproc,iunit,fname,z3dummy,n,m,l)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      COMPLEX(DP),INTENT(IN),TARGET :: z3dummy(n,m,l)
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2*n*m*l, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr,FORM='unformatted',&
           & STATUS='unknown',ACCESS='direct',RECL=unf_recl)
      WRITE(iunit,REC=1,IOSTAT=ierr) z3dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! ************
    ! ************
    ! *** READ ***
    ! ************
    ! ************
    !
    ! READ I0
    !
    SUBROUTINE serial_i0_data_read(lproc,iunit,fname,i0dummy)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(OUT) :: i0dummy
      !
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      READ(iunit,*) 
      READ(iunit,'(i14)',IOSTAT=ierr) i0dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! READ I1
    !
    SUBROUTINE serial_i1_data_read(lproc,iunit,fname,i1dummy,n)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      INTEGER,INTENT(OUT) :: i1dummy(n)
      !
      INTEGER :: ierr
      INTEGER :: i
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      READ(iunit,*) 
      DO i=1,n
         READ(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i1dummy(i)
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! READ I2
    !
    SUBROUTINE serial_i2_data_read(lproc,iunit,fname,i2dummy,n,m)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      INTEGER,INTENT(OUT) :: i2dummy(n,m)
      !
      INTEGER :: ierr
      INTEGER :: i,j
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      READ(iunit,*) 
      DO j=1,m
         DO i=1,n
            READ(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i2dummy(i,j)
         ENDDO
         READ(iunit,*) 
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! READ I3
    !
    SUBROUTINE serial_i3_data_read(lproc,iunit,fname,i3dummy,n,m,l)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      INTEGER,INTENT(OUT) :: i3dummy(n,m,l)
      !
      INTEGER :: ierr
      INTEGER :: i,j,k
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".dat",IOSTAT=ierr)
      READ(iunit,*) 
      DO k=1,l
         DO j=1,m
            DO i=1,n
               READ(iunit,'(i14)',ADVANCE='NO',IOSTAT=ierr) i3dummy(i,j,k)
            ENDDO
            READ(iunit,*) 
         ENDDO
         READ(iunit,*) 
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.dat',ierr)
      !
    END SUBROUTINE
    !
    ! READ D0
    !
    SUBROUTINE serial_d0_data_read(lproc,iunit,fname,d0dummy)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      REAL(DP),INTENT(OUT) :: d0dummy
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(1, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) d0dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ D1
    !
    SUBROUTINE serial_d1_data_read(lproc,iunit,fname,d1dummy,n)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      REAL(DP),INTENT(OUT) :: d1dummy(n)
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(n, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) d1dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ D2
    !
    SUBROUTINE serial_d2_data_read(lproc,iunit,fname,d2dummy,n,m)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      REAL(DP),INTENT(OUT) :: d2dummy(n,m)
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(n*m, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) d2dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ D3
    !
    SUBROUTINE serial_d3_data_read(lproc,iunit,fname,d3dummy,n,m,l)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      REAL(DP),INTENT(OUT) :: d3dummy(n,m,l)
      !
      INTEGER :: ierr
      INTEGER(8) :: unf_recl
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(n*m*l, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) d3dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ Z0
    !
    SUBROUTINE serial_z0_data_read(lproc,iunit,fname,z0dummy)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(OUT),TARGET :: z0dummy
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) z0dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ Z1
    !
    SUBROUTINE serial_z1_data_read(lproc,iunit,fname,z1dummy,n)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n
      COMPLEX(DP),INTENT(OUT),TARGET :: z1dummy(n)
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2*n, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) z1dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ Z2
    !
    SUBROUTINE serial_z2_data_read(lproc,iunit,fname,z2dummy,n,m)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m
      COMPLEX(DP),INTENT(OUT),TARGET :: z2dummy(n,m)
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2*n*m, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) z2dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ Z3
    !
    SUBROUTINE serial_z3_data_read(lproc,iunit,fname,z3dummy,n,m,l)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: n,m,l
      COMPLEX(DP),INTENT(OUT),TARGET :: z3dummy(n,m,l)
      !
      INTEGER(8) :: unf_recl
      INTEGER :: ierr
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      unf_recl = DIRECT_IO_FACTOR * INT(2*n*m*l, KIND=KIND(unf_recl))
      OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(fname))//".raw",IOSTAT=ierr, &
         & FORM='unformatted', STATUS='unknown', ACCESS='direct', RECL=unf_recl)
      READ(iunit,REC=1,IOSTAT=ierr) z3dummy
      CLOSE(iunit,IOSTAT=ierr)
      ! 
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! ###############
    ! ###############
    ! ## PARALLEL ###
    ! ###############
    ! ###############
    !
    ! WRITE I1, PARALLEL
    !
    SUBROUTINE parallel_i1_data_write(fname,i1dummy,nloc,offset,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: offset
      INTEGER,INTENT(IN) :: i1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id
      INTEGER :: sizeofdatum
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_TYPE_SIZE(MPI_INTEGER,sizeofdatum,ierr)
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,file_id,ierr)
      bytes_displacement = offset * sizeofdatum 
      CALL MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_WRITE(file_id,i1dummy,nloc,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE_I1, PARALLEL AND WITH MAP
    !
    SUBROUTINE parallel_i1_irrdata_write(fname,i1dummy,nloc,map,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: map(nloc)
      INTEGER,INTENT(IN) :: i1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id,filetype
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,file_id,ierr)
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nloc,1,map,MPI_INTEGER,filetype,ierr)
      CALL MPI_TYPE_COMMIT(filetype, ierr)
      bytes_displacement = 0
      call MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_INTEGER,filetype,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_WRITE_ALL(file_id,i1dummy,nloc,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !       
    END SUBROUTINE
!   !
!   ! WRITE I2, PARALLEL ! 
!   !
!   SUBROUTINE parallel_i2_data_write(fname,i2dummy,nloc,mloc,offset1,offset2,comm)
!     USE parallel_include
!     CHARACTER(*),INTENT(IN) :: fname
!     INTEGER,INTENT(IN) :: nloc,mloc
!     INTEGER,INTENT(IN) :: offset1,offset2
!     INTEGER,INTENT(IN) :: i2dummy(nloc,mloc)
!     INTEGER,INTENT(IN) :: comm ! communicator 
!     !
!     INTEGER :: ierr,file_id
!     INTEGER :: sizeofdatum
!     INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
!     !
!     ierr = 0
!     !
!     CALL MPI_TYPE_SIZE(MPI_INTEGER,sizeofdatum,ierr)
!     CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,file_id,ierr)
!     bytes_displacement = nloc * offset * sizeofdatum 
!     CALL MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierr)
!     CALL MPI_FILE_WRITE(file_id,i2dummy,nloc,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
!     CALL MPI_FILE_CLOSE(file_id,ierr)
!     !
!     IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
!     !
!   END SUBROUTINE
    !
    ! WRITE D1, PARALLEL
    !
    SUBROUTINE parallel_d1_data_write(fname,d1dummy,nloc,offset,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: offset
      REAL(DP),INTENT(IN) :: d1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id
      INTEGER :: sizeofdatum
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,sizeofdatum,ierr)
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,file_id,ierr)
      bytes_displacement = offset * sizeofdatum 
      CALL MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_WRITE(file_id,d1dummy,nloc,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE_D1, PARALLEL AND WITH MAP
    !
    SUBROUTINE parallel_d1_irrdata_write(fname,d1dummy,nloc,map,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: map(nloc)
      REAL(DP),INTENT(IN) :: d1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id,filetype
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,file_id,ierr)
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nloc,1,map,MPI_DOUBLE_PRECISION,filetype,ierr)
      CALL MPI_TYPE_COMMIT(filetype, ierr)
      bytes_displacement = 0
      call MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_DOUBLE_PRECISION,filetype,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_WRITE_ALL(file_id,d1dummy,nloc,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !       
    END SUBROUTINE
    !
    ! WRITE Z1, PARALLEL
    !
    SUBROUTINE parallel_z1_data_write(fname,z1dummy,nloc,offset,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: offset
      COMPLEX(DP),INTENT(IN) :: z1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id
      INTEGER :: sizeofdatum
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,sizeofdatum,ierr)
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,file_id,ierr)
      bytes_displacement = offset * sizeofdatum 
      CALL MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_WRITE(file_id,z1dummy,nloc,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! WRITE_Z1, PARALLEL AND WITH MAP
    !
    SUBROUTINE parallel_z1_irrdata_write(fname,z1dummy,nloc,map,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: map(nloc)
      COMPLEX(DP),INTENT(IN) :: z1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id,filetype
      INTEGER :: sizeofdatum
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,file_id,ierr)
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nloc,1,map,MPI_DOUBLE_COMPLEX,filetype,ierr)
      CALL MPI_TYPE_COMMIT(filetype, ierr)
      bytes_displacement = 0
      call MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_DOUBLE_COMPLEX,filetype,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_WRITE_ALL(file_id,z1dummy,nloc,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !       
    END SUBROUTINE
    !
    ! READ I1, PARALLEL
    !
    SUBROUTINE parallel_i1_data_read(fname,i1dummy,nloc,offset,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: offset
      INTEGER,INTENT(OUT) :: i1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id
      INTEGER :: sizeofdatum
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_TYPE_SIZE(MPI_INTEGER,sizeofdatum,ierr)
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_RDONLY,MPI_INFO_NULL,file_id,ierr)
      bytes_displacement = offset * sizeofdatum 
      CALL MPI_FILE_READ_AT(file_id,bytes_displacement,i1dummy,nloc,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ_I1, PARALLEL AND WITH MAP
    !
    SUBROUTINE parallel_i1_irrdata_read(fname,i1dummy,nloc,map,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: map(nloc)
      INTEGER,INTENT(OUT) :: i1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id,filetype
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_RDONLY,MPI_INFO_NULL,file_id,ierr)
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nloc,1,map,MPI_INTEGER,filetype,ierr)
      CALL MPI_TYPE_COMMIT(filetype, ierr)
      bytes_displacement = 0
      call MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_INTEGER,filetype,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_READ_ALL(file_id,i1dummy,nloc,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !       
    END SUBROUTINE
    !
    ! READ D1, PARALLEL
    !
    SUBROUTINE parallel_d1_data_read(fname,d1dummy,nloc,offset,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: offset
      REAL(DP),INTENT(OUT) :: d1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id
      INTEGER :: sizeofdatum
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,sizeofdatum,ierr)
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_RDONLY,MPI_INFO_NULL,file_id,ierr)
      bytes_displacement = offset * sizeofdatum 
      CALL MPI_FILE_READ_AT(file_id,bytes_displacement,d1dummy,nloc,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ_D1, PARALLEL AND WITH MAP
    !
    SUBROUTINE parallel_d1_irrdata_read(fname,d1dummy,nloc,map,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: map(nloc)
      REAL(DP),INTENT(OUT) :: d1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id,filetype
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_RDONLY,MPI_INFO_NULL,file_id,ierr)
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nloc,1,map,MPI_DOUBLE_PRECISION,filetype,ierr)
      CALL MPI_TYPE_COMMIT(filetype, ierr)
      bytes_displacement = 0
      call MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_DOUBLE_PRECISION,filetype,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_READ_ALL(file_id,d1dummy,nloc,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !       
    END SUBROUTINE
    !
    ! READ Z1, PARALLEL
    !
    SUBROUTINE parallel_z1_data_read(fname,z1dummy,nloc,offset,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: offset
      COMPLEX(DP),INTENT(OUT) :: z1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id
      INTEGER :: sizeofdatum
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,sizeofdatum,ierr)
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_RDONLY,MPI_INFO_NULL,file_id,ierr)
      bytes_displacement = offset * sizeofdatum 
      CALL MPI_FILE_READ_AT(file_id,bytes_displacement,z1dummy,nloc,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !
    END SUBROUTINE
    !
    ! READ_Z1, PARALLEL AND WITH MAP
    !
    SUBROUTINE parallel_z1_irrdata_read(fname,z1dummy,nloc,map,comm)
      USE parallel_include
      CHARACTER(*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nloc
      INTEGER,INTENT(IN) :: map(nloc)
      COMPLEX(DP),INTENT(OUT) :: z1dummy(nloc)
      INTEGER,INTENT(IN) :: comm ! communicator 
      !
      INTEGER :: ierr,file_id,filetype
      INTEGER(KIND=MPI_OFFSET_KIND) :: bytes_displacement
      !
      ierr = 0
      !
      CALL MPI_FILE_OPEN(comm,TRIM(ADJUSTL(fname))//".raw",MPI_MODE_RDONLY,MPI_INFO_NULL,file_id,ierr)
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nloc,1,map,MPI_DOUBLE_COMPLEX,filetype,ierr)
      CALL MPI_TYPE_COMMIT(filetype, ierr)
      bytes_displacement = 0
      call MPI_FILE_SET_VIEW(file_id,bytes_displacement,MPI_DOUBLE_COMPLEX,filetype,'native',MPI_INFO_NULL,ierr)
      CALL MPI_FILE_READ_ALL(file_id,z1dummy,nloc,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(file_id,ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot RD F:'//TRIM(ADJUSTL(fname))//'.raw',ierr)
      !       
    END SUBROUTINE
    !
    ! ###############
    ! ###############
    ! ## OUTPUT #####
    ! ###############
    ! ###############
    !
    ! OUTPUT
    !
    SUBROUTINE serial_table_output(lproc,iunit,fname,d2dummy,nrow,ncol,header)
      LOGICAL,INTENT(IN) :: lproc
      INTEGER,INTENT(IN) :: iunit 
      CHARACTER(LEN=*),INTENT(IN) :: fname
      INTEGER,INTENT(IN) :: nrow,ncol
      CHARACTER(LEN=*),INTENT(IN) :: header(ncol)
      REAL(DP),INTENT(IN) :: d2dummy(nrow,ncol)
      !
      INTEGER :: ierr,i,j
      CHARACTER(LEN=128) :: format_string
      REAL(DP) :: help(ncol)
      !
      IF(.NOT.lproc) RETURN
      !
      ierr = 0
      !
      OPEN(UNIT=iunit,FILE="o-"//TRIM(ADJUSTL(fname))//".tab",IOSTAT=ierr)
      WRITE(format_string,*) '("#",',ncol,'(a16))'
      WRITE(iunit,TRIM(format_string)) header
      WRITE(format_string,*) '(" ",',ncol,'(f16.6))'
      DO i=1,nrow
         DO j=1,ncol
            help(j) = d2dummy(i,j)
         ENDDO
         WRITE(iunit,TRIM(format_string),IOSTAT=ierr) help(1:ncol)
      ENDDO
      CLOSE(iunit,IOSTAT=ierr)
      !
      IF(ierr/=0) CALL errore('WEST/IO', 'Cannot WR F:'//"o-"//TRIM(ADJUSTL(fname))//'.tab',ierr)
      !
    END SUBROUTINE
    !
    !
END MODULE
