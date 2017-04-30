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
MODULE io_push
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP 
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTERFACE io_push_value
     MODULE PROCEDURE io_push_d0, io_push_i0, io_push_l0, io_push_c1
  END INTERFACE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_bar()
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      WRITE(stdout,'(5x,"--------------------------------------------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_title(labelin)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      CHARACTER(LEN=*) :: labelin
      !
      ! Workspace
      !
      CHARACTER(LEN=68) :: labelout
      !
      WRITE(labelout,"(a68)") labelin
      !
      WRITE(stdout,'(5x,"  ")')
      CALL io_push_bar
      WRITE(stdout,'(5x,a68)') ADJUSTL(labelout)
      CALL io_push_bar
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_d0(l_desc_in,d0,numa)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: l_desc_in
      INTEGER,INTENT(IN) :: numa
      REAL(DP),INTENT(IN) :: d0
      CHARACTER(LEN=numa) :: l_desc_out
      !
      WRITE(l_desc_out,'(a)') ADJUSTL(TRIM(l_desc_in))
      !
      WRITE(stdout,'(5x,a," = ",f14.6)') l_desc_out, d0
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_es0(l_desc_in,es0,numa)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: l_desc_in
      INTEGER,INTENT(IN) :: numa
      REAL(DP),INTENT(IN) :: es0
      CHARACTER(LEN=numa) :: l_desc_out
      !
      WRITE(l_desc_out,'(a)') ADJUSTL(TRIM(l_desc_in))
      !
      WRITE(stdout,'(5x,a," = ",es14.6)') l_desc_out, es0
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_i0(l_desc_in,i0,numa)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: l_desc_in
      INTEGER,INTENT(IN) :: numa
      INTEGER,INTENT(IN) :: i0
      CHARACTER(LEN=numa) :: l_desc_out
      !
      WRITE(l_desc_out,'(a)') ADJUSTL(TRIM(l_desc_in))
      !
      WRITE(stdout,'(5x,a," = ",i14)') l_desc_out, i0
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_l0(l_desc_in,l0,numa)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: l_desc_in
      INTEGER,INTENT(IN) :: numa
      LOGICAL,INTENT(IN) :: l0
      CHARACTER(LEN=numa) :: l_desc_out
      !
      WRITE(l_desc_out,'(a)') ADJUSTL(TRIM(l_desc_in))
      !
      WRITE(stdout,'(5x,a," = ",l14)') l_desc_out, l0
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_c256(l_desc_in,c256,numa)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: l_desc_in
      INTEGER,INTENT(IN) :: numa
      CHARACTER(LEN=256),INTENT(IN) :: c256
      CHARACTER(LEN=numa) :: l_desc_out
      !
      WRITE(l_desc_out,'(a)') ADJUSTL(TRIM(l_desc_in))
      !
      WRITE(stdout,'(5x,a," = ",a)') l_desc_out, TRIM(c256)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_push_c1(l_desc_in,c1,numa)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: l_desc_in
      INTEGER,INTENT(IN) :: numa
      CHARACTER(LEN=*),INTENT(IN) :: c1
      CHARACTER(LEN=numa) :: l_desc_out
      CHARACTER(LEN=14) :: cout
      !
      WRITE(l_desc_out,'(a)') ADJUSTL(TRIM(l_desc_in))
      WRITE(cout,'(a14)') ADJUSTR(TRIM(c1))
      !
      WRITE(stdout,'(5x,a," = ",a14)') l_desc_out, cout
      !
    END SUBROUTINE
    !
END MODULE
