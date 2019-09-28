MODULE conversions 
   USE kinds, ONLY : DP
   IMPLICIT NONE 
   CONTAINS
    !
    FUNCTION ltoa(l) RESULT(res)
       CHARACTER(:),ALLOCATABLE :: res
       LOGICAL,INTENT(IN) :: l
       CHARACTER(4) :: t="true"
       CHARACTER(5) :: f="false"
       IF( l ) THEN 
          res = t
       ELSE 
          res =f 
       ENDIF
    END FUNCTION
    !
    FUNCTION itoa(i) RESULT(res)
       CHARACTER(:),ALLOCATABLE :: res
       INTEGER,INTENT(IN) :: i
       CHARACTER(RANGE(i)+2) :: tmp
       WRITE(tmp,'(I0)') i
       res = TRIM(tmp)
    END FUNCTION
    !
    FUNCTION dtoa(d) RESULT(res)
       CHARACTER(:),ALLOCATABLE :: res
       REAL(DP),INTENT(IN) :: d
       CHARACTER(14) :: tmp
       WRITE(tmp,'(ES14.6)') d
       res = TRIM(tmp)
    END FUNCTION
END MODULE
