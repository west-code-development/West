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
! Marco Govoni
!
! -------------------------------------------------------------------
MODULE function3d
 ! -----------------------------------------------------------------
 !
 IMPLICIT NONE
 !
 INTERFACE write_function3d
    MODULE PROCEDURE write_function3d_real !, write_function3d_complex
 END INTERFACE
 !
 INTERFACE read_function3d
    MODULE PROCEDURE read_function3d_real !, read_function3d_complex
 END INTERFACE
 !
 CONTAINS
 !
 !-----------------------------------------------------------------
   SUBROUTINE write_function3d_real ( fname, f_r, dfft )
   ! -----------------------------------------------------------------
   !
   USE kinds,                       ONLY : DP
   USE cell_base,                   ONLY : celldm, at
   USE control_flags,               ONLY : gamma_only
   USE mp_bands,                    ONLY : me_bgrp
   USE scatter_mod,                 ONLY : gather_grid
   USE fft_types,                   ONLY : fft_type_descriptor
   USE forpy_mod,  ONLY: call_py, call_py_noret, import_py, module_py
   USE forpy_mod,  ONLY: tuple, tuple_create
   USE forpy_mod,  ONLY: dict, dict_create
   USE forpy_mod,  ONLY: list, list_create
   USE forpy_mod,  ONLY: object, cast
   USE forpy_mod,  ONLY: exception_matches, KeyError, err_clear, err_print
   USE conversions, ONLY : ltoa, itoa, dtoa
   USE base64_module
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   CHARACTER(LEN=*),INTENT(IN) :: fname
   TYPE(fft_type_descriptor), INTENT(IN) :: dfft
   REAL(DP),INTENT(IN) :: f_r(dfft%nnr)
   !
   ! Workspace
   !
   CHARACTER(LEN=:),ALLOCATABLE :: charbase64
   REAL(DP),ALLOCATABLE :: f_r_gathered(:), f_r_gathered_nopadded(:)
   INTEGER :: ndim, nbytes, nlen
   TYPE(tuple) :: args
   TYPE(dict) :: kwargs
   TYPE(module_py) :: pymod
   TYPE(object) :: return_obj
   INTEGER :: return_int
   INTEGER :: IERR
   !
   ! Gather the function
   !
   ALLOCATE(f_r_gathered(dfft%nr1x*dfft%nr2x*dfft%nr3x)); f_r_gathered = 0._DP
   CALL gather_grid(dfft,f_r,f_r_gathered)
   !
   IF( me_bgrp == 0 ) THEN
      !
      ALLOCATE(f_r_gathered_nopadded(dfft%nr1*dfft%nr2*dfft%nr3)); f_r_gathered_nopadded = 0._DP
      CALL remove_padding_real(dfft,f_r_gathered,f_r_gathered_nopadded)
      !
      ! Encode
      !
      ndim = dfft%nr1*dfft%nr2*dfft%nr3
      nbytes = SIZEOF(f_r_gathered_nopadded(1)) * ndim
      nlen = lenbase64(nbytes)
      ALLOCATE(CHARACTER(LEN=nlen) :: charbase64)
      IF (.NOT. islittleendian()) CALL base64_byteswap_double(nbytes,f_r_gathered_nopadded(1:ndim))
      CALL base64_encode_double(f_r_gathered_nopadded(1:ndim), ndim, charbase64)
      DEALLOCATE(f_r_gathered_nopadded)
      !
      IERR = import_py(pymod, "west_function3d")
      !
      IERR = tuple_create(args, 1)
      IERR = args%setitem(0, TRIM(ADJUSTL(fname)) )
      IERR = dict_create(kwargs)
      IERR = kwargs%setitem("name","delta_v")
      IERR = kwargs%setitem("domain",'{'// &
      & '"a":['//dtoa(celldm(1)*at(1,1))//","//dtoa(celldm(1)*at(2,1))//","//dtoa(celldm(1)*at(3,1))//'],'// &
      & '"b":['//dtoa(celldm(1)*at(1,2))//","//dtoa(celldm(1)*at(2,2))//","//dtoa(celldm(1)*at(3,2))//'],'// &
      & '"c":['//dtoa(celldm(1)*at(1,3))//","//dtoa(celldm(1)*at(2,3))//","//dtoa(celldm(1)*at(3,3))//'] }')
      IERR = kwargs%setitem("grid","["//itoa(dfft%nr1)//","//itoa(dfft%nr2)//","//itoa(dfft%nr3)//"]" )
      IERR = kwargs%setitem("grid_function",charbase64)
      IERR = kwargs%setitem("dtype","double")
      !
      IERR = call_py_noret(pymod, "base64_to_function3D", args, kwargs)
      !
      CALL kwargs%destroy
      CALL args%destroy
      CALL pymod%destroy
      !
   ENDIF
   !
   DEALLOCATE(f_r_gathered)
   !
 END SUBROUTINE
 !
 !-----------------------------------------------------------------
   SUBROUTINE read_function3d_real ( fname, f_r, dfft )
   ! -----------------------------------------------------------------
   !
   USE kinds,                       ONLY : DP
   USE cell_base,                   ONLY : celldm, at
   USE control_flags,               ONLY : gamma_only
   USE mp_bands,                    ONLY : me_bgrp
   USE scatter_mod,                 ONLY : scatter_grid
   USE fft_types,                   ONLY : fft_type_descriptor
   USE forpy_mod,  ONLY: call_py, call_py_noret, import_py, module_py
   USE forpy_mod,  ONLY: tuple, tuple_create
   USE forpy_mod,  ONLY: dict, dict_create
   USE forpy_mod,  ONLY: list, list_create
   USE forpy_mod,  ONLY: object, cast
   USE forpy_mod,  ONLY: exception_matches, KeyError, err_clear, err_print
   USE base64_module
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   CHARACTER(LEN=*),INTENT(IN) :: fname
   TYPE(fft_type_descriptor), INTENT(IN) :: dfft
   REAL(DP),INTENT(OUT) :: f_r(dfft%nnr)
   !
   ! Workspace
   !
   CHARACTER(LEN=:),ALLOCATABLE :: charbase64
   REAL(DP),ALLOCATABLE :: f_r_gathered(:), f_r_gathered_nopadded(:)
   INTEGER :: ndim, nbytes, nlen
   TYPE(tuple) :: args
   TYPE(dict) :: kwargs, return_dict
   TYPE(module_py) :: pymod
   TYPE(object) :: return_obj
   INTEGER :: return_int
   INTEGER :: IERR
   !
   ALLOCATE(f_r_gathered(dfft%nr1x*dfft%nr2x*dfft%nr3x)); f_r_gathered = 0._DP
   !
   IF( me_bgrp == 0 ) THEN
      !
      ! Decode
      !
      !
      IERR = import_py(pymod, "west_function3d")
      !
      IERR = tuple_create(args, 1)
      IERR = args%setitem(0, TRIM(ADJUSTL(fname)) )
      IERR = dict_create(kwargs)
      !
      IERR = call_py(return_obj,pymod, "function3D_to_base64", args, kwargs)
      !
      IERR = cast(return_dict, return_obj)
      !
      ndim = dfft%nr1*dfft%nr2*dfft%nr3
      nbytes = SIZEOF(f_r_gathered_nopadded(1)) * ndim
      nlen = lenbase64(nbytes)
      ALLOCATE(CHARACTER(LEN=nlen) :: charbase64)
      IERR = return_dict%getitem(charbase64, "grid_function")
      ALLOCATE(f_r_gathered_nopadded(dfft%nr1*dfft%nr2*dfft%nr3)); f_r_gathered_nopadded = 0._DP
      CALL base64_decode_double(charbase64(1:nlen), ndim, f_r_gathered_nopadded(1:ndim))
      IF (.NOT. islittleendian()) CALL base64_byteswap_double(nbytes,f_r_gathered_nopadded(1:ndim))
      !
      CALL kwargs%destroy
      CALL args%destroy
      CALL return_obj%destroy
      CALL return_dict%destroy
      CALL pymod%destroy
      !
      CALL add_padding_real(dfft,f_r_gathered_nopadded,f_r_gathered)
      DEALLOCATE(f_r_gathered_nopadded)
      !
   ENDIF
   !
   CALL scatter_grid(dfft,f_r_gathered,f_r)
   !
   DEALLOCATE(f_r_gathered)
   !
 END SUBROUTINE
 !
 !
 SUBROUTINE add_padding_real(dfft,f_r_gathered_nopadded,f_r_gathered)
   USE kinds, ONLY :DP
   USE fft_types,                   ONLY : fft_type_descriptor
   IMPLICIT NONE
   TYPE(fft_type_descriptor), INTENT(IN) :: dfft
   REAL(DP),INTENT(IN) :: f_r_gathered_nopadded(dfft%nr1*dfft%nr2*dfft%nr3)
   REAL(DP),INTENT(OUT) :: f_r_gathered(dfft%nr1x*dfft%nr2x*dfft%nr3x)
   INTEGER :: i,j,k,ir_notpadded,ir_padded
   IF( dfft%nr1 == dfft%nr1x .AND. dfft%nr2 == dfft%nr2x .AND. dfft%nr3 == dfft%nr3x) THEN
      f_r_gathered = f_r_gathered_nopadded
   ELSE
      f_r_gathered = 0._DP
      DO k = 1, dfft%nr3
         DO j = 1, dfft%nr2
            DO i = 1, dfft%nr1
               ir_notpadded = (i-1)*dfft%nr1 *dfft%nr2  + (j-1)*dfft%nr2  + k
               ir_padded    = (i-1)*dfft%nr1x*dfft%nr2x + (j-1)*dfft%nr2x + k
               f_r_gathered(ir_padded) = f_r_gathered_nopadded(ir_notpadded)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
 END SUBROUTINE
 !
 SUBROUTINE remove_padding_real(dfft,f_r_gathered,f_r_gathered_nopadded)
   USE kinds, ONLY :DP
   USE fft_types,                   ONLY : fft_type_descriptor
   IMPLICIT NONE
   TYPE(fft_type_descriptor), INTENT(IN) :: dfft
   REAL(DP),INTENT(IN) :: f_r_gathered(dfft%nr1x*dfft%nr2x*dfft%nr3x)
   REAL(DP),INTENT(OUT) :: f_r_gathered_nopadded(dfft%nr1*dfft%nr2*dfft%nr3)
   INTEGER :: i,j,k,ir_notpadded,ir_padded
   IF( dfft%nr1 == dfft%nr1x .AND. dfft%nr2 == dfft%nr2x .AND. dfft%nr3 == dfft%nr3x) THEN
      f_r_gathered_nopadded = f_r_gathered
   ELSE
      f_r_gathered_nopadded = 0._DP
      DO k = 1, dfft%nr3
         DO j = 1, dfft%nr2
            DO i = 1, dfft%nr1
               ir_notpadded = (i-1)*dfft%nr1 *dfft%nr2  + (j-1)*dfft%nr2  + k
               ir_padded    = (i-1)*dfft%nr1x*dfft%nr2x + (j-1)*dfft%nr2x + k
               f_r_gathered_nopadded(ir_notpadded) = f_r_gathered(ir_padded)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
 END SUBROUTINE
 !
END MODULE
