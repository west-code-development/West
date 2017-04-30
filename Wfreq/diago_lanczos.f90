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
SUBROUTINE diago_lanczos( bnorm, diago, subdiago, braket, ndim ) 
  !-----------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE westcom,             ONLY : n_lanczos
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP),INTENT(IN) :: bnorm
  REAL(DP),INTENT(INOUT) :: diago(n_lanczos)
  REAL(DP),INTENT(INOUT) :: subdiago(n_lanczos-1)
  REAL(DP),INTENT(INOUT) :: braket(ndim,n_lanczos)
  INTEGER,INTENT(IN) :: ndim
  !
  ! Workspace
  !
  INTEGER :: info,il,ip
  REAL(DP),ALLOCATABLE :: z(:,:)
  REAL(DP),ALLOCATABLE :: work(:)
  REAL(DP),ALLOCATABLE :: braket_tilde(:,:)
  !
  ALLOCATE( z(n_lanczos,n_lanczos), work(MAX(1,2*n_lanczos-2)), braket_tilde(ndim,n_lanczos) )
  !
  ! Diagonalize
  !
  CALL DSTEQR ('I', n_lanczos, diago, subdiago, z, n_lanczos, work, info)
  !
  ! braket_tilde(ip,il) = brak(ip,:) * z(:,il)
  !  
  CALL DGEMM('N','N',ndim,n_lanczos,n_lanczos,bnorm,braket,ndim,z,n_lanczos,0._DP,braket_tilde,ndim)
  ! 
  ! braket_tilde(ip,il) = braket_tilde(ip,il) * z(1,il)
  !
  DO il=1,n_lanczos
     DO ip = 1, ndim
        braket(ip,il) = braket_tilde(ip,il) * z(1,il)
     ENDDO
  ENDDO 
  !
  DEALLOCATE( z, work, braket_tilde )
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE diago_lanczos_complex( bnorm, diago, subdiago, braket, ndim ) 
  !-----------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE westcom,             ONLY : n_lanczos
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP),INTENT(IN) :: bnorm
  REAL(DP),INTENT(INOUT) :: diago(n_lanczos)
  REAL(DP),INTENT(INOUT) :: subdiago(n_lanczos-1)
  COMPLEX(DP),INTENT(INOUT) :: braket(ndim,n_lanczos)
  INTEGER,INTENT(IN) :: ndim
  !
  ! Workspace
  !
  INTEGER :: info,il,ip
  COMPLEX(DP) :: znorm, zzero
  COMPLEX(DP),ALLOCATABLE :: z(:,:)
  REAL(DP),ALLOCATABLE :: work(:)
  COMPLEX(DP),ALLOCATABLE :: braket_tilde(:,:)
  !
  znorm = CMPLX( bnorm, 0._DP, KIND=DP )
  zzero = ( 0._DP, 0._DP ) 
  !
  ALLOCATE( z(n_lanczos,n_lanczos), work(MAX(1,2*n_lanczos-2)), braket_tilde(ndim,n_lanczos) )
  !
  ! Diagonalize
  !
  CALL ZSTEQR ('I', n_lanczos, diago, subdiago, z, n_lanczos, work, info)
  !
  ! braket_tilde(ip,il) = brak(ip,:) * z(:,il)
  !  
  CALL ZGEMM('N','N',ndim,n_lanczos,n_lanczos,znorm,braket,ndim,z,n_lanczos,zzero,braket_tilde,ndim)
  ! 
  ! braket_tilde(ip,il) = braket_tilde(ip,il) * z(1,il)
  !
  DO il=1,n_lanczos
     DO ip = 1, ndim
        braket(ip,il) = braket_tilde(ip,il) * z(1,il)
     ENDDO
  ENDDO 
  !
  DEALLOCATE( z, work, braket_tilde )
  !
END SUBROUTINE
