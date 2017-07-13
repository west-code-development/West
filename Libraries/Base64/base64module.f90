!
! Copyright (C) 2015-2017 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file: 
! Huihuo Zheng, Marco Govoni
!
!-------------------------------------------------------------------
module base64_module
  !-------------------------------------------------------------------
  !
  USE, INTRINSIC :: ISO_C_Binding,   ONLY  : C_DOUBLE, C_DOUBLE_COMPLEX, C_CHAR, C_SIGNED_CHAR, C_PTR, C_NULL_PTR, C_INT
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER :: ICHAR
  LOGICAL, PARAMETER :: lbigendian = ( ICHAR( TRANSFER(1,'a') ) == 0 )
  !
  INTERFACE
     !
     SUBROUTINE base64_init() BIND(C, NAME="b64init") 
       IMPORT
     END SUBROUTINE
     !
     SUBROUTINE base64_encode_double(from, n, to) BIND(C, NAME="encode_double")
       IMPORT 
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       REAL(C_DOUBLE), INTENT(IN) :: from(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: to(*)
     END SUBROUTINE
     !
     SUBROUTINE base64_decode_double(from, n, to) BIND(C, NAME="decode_double")
       IMPORT 
       INTEGER(C_int), VALUE, INTENT(IN) :: n
       REAL(C_DOUBLE), INTENT(IN) :: to(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: from(*)
     END SUBROUTINE
     !
     SUBROUTINE base64_encode_complex(from, n, to) BIND(C, NAME="encode_complex")
       IMPORT 
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: from(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: to(*)
     END SUBROUTINE
     !
     SUBROUTINE base64_decode_complex(from, n, to) BIND(C, NAME="decode_complex")
       IMPORT
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: to(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: from(*)
     END SUBROUTINE
     !
     SUBROUTINE base64_byteswap_complex(n, to) BIND(C, NAME="byteswap_complex") 
       IMPORT
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT) :: to(*)
     END SUBROUTINE  
     !   
     SUBROUTINE base64_byteswap_double(n, to) BIND(C, NAME="byteswap_double") 
       IMPORT 
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       REAL(C_DOUBLE), INTENT(INOUT) :: to(*)
     END SUBROUTINE
     !
  END INTERFACE
  !
  CONTAINS 
     !
     INTEGER FUNCTION lenbase64( nbytes )
        IMPLICIT NONE 
        INTEGER,INTENT(IN) :: nbytes
        lenbase64 = ( ( nbytes + 2 ) / 3 ) * 4
     END FUNCTION
     !
     LOGICAL FUNCTION islittleendian( )
        IMPLICIT NONE 
        islittleendian = (.NOT.lbigendian)
     END FUNCTION
  !
END MODULE
