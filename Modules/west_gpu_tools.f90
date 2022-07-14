!
! Copyright (C) 2020-2021
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
!-----------------------------------------------------------------------
MODULE west_gpu_tools
   !-----------------------------------------------------------------------
   !
   USE kinds,         ONLY : DP
   USE ISO_C_BINDING, ONLY : C_INT
#if defined(__CUDA)
   !
   IMPLICIT NONE
   !
   INTERFACE memcpy_H2D
      MODULE PROCEDURE memcpy_H2D_c16
      MODULE PROCEDURE memcpy_H2D_r8
   END INTERFACE
   !
   INTERFACE memcpy_D2H
      MODULE PROCEDURE memcpy_D2H_c16
      MODULE PROCEDURE memcpy_D2H_r8
   END INTERFACE
   !
   INTERFACE acc_memcpy_H2D
      SUBROUTINE acc_memcpy_H2D_c16(dev,src,bytes) BIND(C,NAME='acc_memcpy_to_device')
      USE kinds,         ONLY : DP
      USE ISO_C_BINDING, ONLY : C_INT
      COMPLEX(DP) :: dev(*)
      COMPLEX(DP) :: src(*)
      INTEGER(C_INT), VALUE :: bytes
      END SUBROUTINE
      !
      SUBROUTINE acc_memcpy_H2D_r8(dev,src,bytes) BIND(C,NAME='acc_memcpy_to_device')
      USE kinds,         ONLY : DP
      USE ISO_C_BINDING, ONLY : C_INT
      REAL(DP) :: dev(*)
      REAL(DP) :: src(*)
      INTEGER(C_INT), VALUE :: bytes
      END SUBROUTINE
   END INTERFACE
   !
   INTERFACE acc_memcpy_D2H
      SUBROUTINE acc_memcpy_D2H_c16(dest,dev,bytes) BIND(C,NAME='acc_memcpy_from_device')
      USE kinds,         ONLY : DP
      USE ISO_C_BINDING, ONLY : C_INT
      COMPLEX(DP) :: dest(*)
      COMPLEX(DP) :: dev(*)
      INTEGER(C_INT), VALUE :: bytes
      END SUBROUTINE
      !
      SUBROUTINE acc_memcpy_D2H_r8(dest,dev,bytes) BIND(C,NAME='acc_memcpy_from_device')
      USE kinds,         ONLY : DP
      USE ISO_C_BINDING, ONLY : C_INT
      REAL(DP) :: dest(*)
      REAL(DP) :: dev(*)
      INTEGER(C_INT), VALUE :: bytes
      END SUBROUTINE
   END INTERFACE
   !
   CONTAINS
   !
   !-----------------------------------------------------------------------
   SUBROUTINE memcpy_H2D_c16(dev,src,length)
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   COMPLEX(DP), INTENT(INOUT) :: dev(*)
   COMPLEX(DP), INTENT(IN) :: src(*)
   INTEGER, INTENT(IN) :: length
   !
   ! Workspace
   !
   INTEGER(C_INT) :: bytes
   COMPLEX(DP) :: dummy
   !
   bytes = INT(length,KIND=C_INT)*SIZEOF(dummy)
   !
   !$acc data present(dev)
   !$acc host_data use_device(dev)
   CALL acc_memcpy_H2D(dev,src,bytes)
   !$acc end host_data
   !$acc end data
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE memcpy_H2D_r8(dev,src,length)
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   REAL(DP), INTENT(INOUT) :: dev(*)
   REAL(DP), INTENT(IN) :: src(*)
   INTEGER, INTENT(IN) :: length
   !
   ! Workspace
   !
   INTEGER(C_INT) :: bytes
   REAL(DP) :: dummy
   !
   bytes = INT(length,KIND=C_INT)*SIZEOF(dummy)
   !
   !$acc data present(dev)
   !$acc host_data use_device(dev)
   CALL acc_memcpy_H2D(dev,src,bytes)
   !$acc end host_data
   !$acc end data
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE memcpy_D2H_c16(dest,dev,length)
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   COMPLEX(DP), INTENT(OUT) :: dest(*)
   COMPLEX(DP), INTENT(INOUT) :: dev(*)
   INTEGER, INTENT(IN) :: length
   !
   ! Workspace
   !
   INTEGER(C_INT) :: bytes
   COMPLEX(DP) :: dummy
   !
   bytes = INT(length,KIND=C_INT)*SIZEOF(dummy)
   !
   !$acc data present(dev)
   !$acc host_data use_device(dev)
   CALL acc_memcpy_D2H(dest,dev,bytes)
   !$acc end host_data
   !$acc end data
   !
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE memcpy_D2H_r8(dest,dev,length)
   !-----------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   REAL(DP), INTENT(OUT) :: dest(*)
   REAL(DP), INTENT(INOUT) :: dev(*)
   INTEGER, INTENT(IN) :: length
   !
   ! Workspace
   !
   INTEGER(C_INT) :: bytes
   REAL(DP) :: dummy
   !
   bytes = INT(length,KIND=C_INT)*SIZEOF(dummy)
   !
   !$acc data present(dev)
   !$acc host_data use_device(dev)
   CALL acc_memcpy_D2H(dest,dev,bytes)
   !$acc end host_data
   !$acc end data
   !
   END SUBROUTINE
#endif
END MODULE
