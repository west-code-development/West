!
! Copyright (C) 2015-2022 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `LICENSE'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Marco Govoni
!
!-----------------------------------------------------------------------
MODULE sort_tools
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP,i8b
  !
  IMPLICIT NONE
  !
  INTERFACE swap
     MODULE PROCEDURE swap_i4, swap_i8
  END INTERFACE

  INTERFACE downheap
     MODULE PROCEDURE downheap_i4, downheap_i8
  END INTERFACE

  INTERFACE heapsort
     MODULE PROCEDURE heapsort_i4, heapsort_i8
  END INTERFACE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE heapsort_i8(length,array,idx)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: length
      INTEGER(i8b), INTENT(INOUT) :: array(length)
      INTEGER, INTENT(OUT) :: idx(length)
      !
      ! Workspace
      !
      INTEGER :: top
      INTEGER :: i
      !
      DO i = 1,length
         idx(i) = i
      ENDDO
      !
      top = length/2
      !
      DO i = top,1,-1
         CALL downheap(length,array,idx,i,length)
      ENDDO
      !
      i = length
      !
      DO WHILE(i > 1)
         CALL swap(length,array,1,i)
         CALL swap(length,idx,1,i)
         !
         i = i-1
         !
         CALL downheap(length,array,idx,1,i)
      ENDDO
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE downheap_i8(length,a,b,top,bottom)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: length
      INTEGER(i8b), INTENT(INOUT) :: a(length)
      INTEGER, INTENT(INOUT) :: b(length)
      INTEGER, INTENT(IN) :: top
      INTEGER, INTENT(IN) :: bottom
      !
      ! Workspace
      !
      INTEGER :: v
      INTEGER :: w
      !
      v = top
      w = 2*v
      !
      DO WHILE(w <= bottom)
         IF(w+1 <= bottom) then
            IF(a(w+1) > a(w)) THEN
               w = w+1
            END IF
         ENDIF
         !
         IF(a(v) >= a(w)) THEN
            RETURN
         ELSE
            CALL swap(length,a,v,w)
            CALL swap(length,b,v,w)
            !
            v = w
            w = 2*v
         ENDIF
      ENDDO
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE heapsort_i4(length,array,idx)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: length
      INTEGER, INTENT(INOUT) :: array(length)
      INTEGER, INTENT(OUT) :: idx(length)
      !
      ! Workspace
      !
      INTEGER :: top
      INTEGER :: i
      !
      DO i = 1,length
         idx(i) = i
      ENDDO
      !
      top = length/2
      !
      DO i = top,1,-1
         CALL downheap(length,array,idx,i,length)
      ENDDO
      !
      i = length
      !
      DO WHILE(i > 1)
         CALL swap(length,array,1,i)
         CALL swap(length,idx,1,i)
         !
         i = i-1
         !
         CALL downheap(length,array,idx,1,i)
      ENDDO
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE downheap_i4(length,a,b,top,bottom)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: length
      INTEGER, INTENT(INOUT) :: a(length)
      INTEGER, INTENT(INOUT) :: b(length)
      INTEGER, INTENT(IN) :: top
      INTEGER, INTENT(IN) :: bottom
      !
      ! Workspace
      !
      INTEGER :: v
      INTEGER :: w
      !
      v = top
      w = 2*v
      !
      DO WHILE(w <= bottom)
         IF(w+1 <= bottom) then
            IF(a(w+1) > a(w)) THEN
               w = w+1
            END IF
         ENDIF
         !
         IF(a(v) >= a(w)) THEN
            RETURN
         ELSE
            CALL swap(length,a,v,w)
            CALL swap(length,b,v,w)
            !
            v = w
            w = 2*v
         ENDIF
      ENDDO
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE swap_i8(length,array,i,j)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: length
      INTEGER(i8b), INTENT(INOUT) :: array(length)
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN) :: j
      !
      ! Workspace
      !
      INTEGER(i8b) :: tmp
      !
      tmp = array(i)
      array(i) = array(j)
      array(j) = tmp
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE swap_i4(length,array,i,j)
    !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: length
      INTEGER, INTENT(INOUT) :: array(length)
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN) :: j
      !
      ! Workspace
      !
      INTEGER :: tmp
      !
      tmp = array(i)
      array(i) = array(j)
      array(j) = tmp
      !
    END SUBROUTINE
    !
END MODULE
