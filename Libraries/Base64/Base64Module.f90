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
! Huihuo Zheng
!
!-------------------------------------------------------------------
module Base64_module
  !-------------------------------------------------------------------
  !
  USE, INTRINSIC :: ISO_C_Binding,   ONLY  : C_DOUBLE, C_DOUBLE_COMPLEX, C_CHAR, C_SIGNED_CHAR, C_PTR, C_NULL_PTR, C_INT
  !
  IMPLICIT NONE
  !
  TYPE Base64_type
     TYPE(C_PTR) :: object = C_NULL_PTR
  END TYPE Base64_type
  !
  INTERFACE
     !
     FUNCTION C_Base64__new() RESULT(this) BIND(C, NAME="Base64__new") 
       IMPORT
       TYPE(C_PTR) :: this
     END FUNCTION C_Base64__new
     !
     FUNCTION C_Base64__delete(this) BIND (C, NAME="Base64__delete") 
       IMPORT
       TYPE(C_PTR), VALUE :: this
     END FUNCTION C_Base64__delete
     !
     SUBROUTINE C_Base64__encode_double(this, from, n, to) BIND(C, NAME="Base64__encode_double")
       IMPORT 
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       REAL(C_DOUBLE), INTENT(IN) :: from(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: to(*)
       TYPE(C_PTR), VALUE :: this
     END SUBROUTINE C_Base64__encode_double
     !
     SUBROUTINE C_Base64__decode_double(this, from, n, to) BIND(C, NAME="Base64__decode_double")
       IMPORT 
       INTEGER(C_int), VALUE, INTENT(IN) :: n
       REAL(C_DOUBLE), INTENT(IN) :: to(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: from(*)
       TYPE(C_PTR), VALUE :: this
     END SUBROUTINE C_Base64__decode_double
     !
     SUBROUTINE C_Base64__encode_complex(this, from, n, to) BIND(C, NAME="Base64__encode_complex")
       IMPORT 
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: from(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: to(*)
       TYPE(C_PTR), VALUE :: this
     END SUBROUTINE C_Base64__encode_complex
     !
     SUBROUTINE C_Base64__decode_complex(this, from, n, to) BIND(C, NAME="Base64__decode_complex")
       IMPORT
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: to(*)
       CHARACTER(C_CHAR), INTENT(OUT) :: from(*)
       TYPE(C_PTR), VALUE :: this
     END SUBROUTINE C_Base64__decode_complex
     !
     SUBROUTINE C_Base64__byteswap_complex(this, n, to) BIND(C, NAME="Base64__byteswap_complex") 
       IMPORT
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT) :: to(*)
       TYPE(C_PTR), VALUE :: this
     END SUBROUTINE C_Base64__byteswap_complex  
     !   
     SUBROUTINE C_Base64__byteswap_double(this, n, to) BIND(C, NAME="Base64__byteswap_double") 
       IMPORT 
       INTEGER(C_INT), VALUE, INTENT(IN) :: n
       REAL(C_DOUBLE), INTENT(INOUT) :: to(*)
       TYPE(C_PTR), VALUE :: this
     END SUBROUTINE C_Base64__byteswap_double
     !
  END INTERFACE
  !
END MODULE
