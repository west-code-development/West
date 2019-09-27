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
!-------------------------------------------------------------------------
SUBROUTINE apply_operator (m,dvg,dng,tr2,iq)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE mp,                    ONLY : mp_barrier
  USE westcom,               ONLY : npwqx,npwq
  USE mp_world,              ONLY : world_comm
  USE types_coulomb,         ONLY : pot3D
  USE dfpt_module,           ONLY : dfpt
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN), OPTIONAL :: iq
  REAL(DP),INTENT(IN), OPTIONAL :: tr2
  INTEGER, INTENT(IN) :: m
  COMPLEX(DP), INTENT(IN) :: dvg(npwqx,m)
  COMPLEX(DP), INTENT(OUT) :: dng(npwqx,m)
  !
  ! Workspace
  !
  INTEGER :: ipert, ig
  INTEGER :: iq_
  !
  REAL(DP) :: tr2_
  !
  COMPLEX(DP), ALLOCATABLE ::aux_g(:,:)
  !
  LOGICAL :: l_outsource
  !
  CALL mp_barrier( world_comm )
  !
  l_outsource = .FALSE.
  !
  IF(PRESENT(iq)) THEN 
     iq_ = iq
  ELSE
     iq_ = 1
  ENDIF
  IF(PRESENT(tr2)) THEN 
     tr2_ = tr2
  ELSE
     tr2_ = 1.d-8
  ENDIF
  !
  dng=0._DP
  !
  ALLOCATE( aux_g(npwqx,m) ); aux_g=0._DP
  !
  DO CONCURRENT (ipert = 1:m, ig = 1:npwq)
     aux_g(ig,ipert) = dvg(ig,ipert) * pot3D%sqvc(ig) ! perturbation acts only on body 
  ENDDO
  !
  IF( l_outsource ) THEN
     CALL calc_outsourced(m,aux_g,dng)
  ELSE 
     CALL dfpt(m,aux_g,dng,tr2_,iq_)
  ENDIF 
  !
  DO CONCURRENT (ipert = 1:m, ig = 1:npwq) 
     dng(ig,ipert) = dng(ig,ipert) * pot3D%sqvc(ig) ! perturbation acts only on body  
  ENDDO
  !
  CALL mp_barrier( world_comm )
  !
END SUBROUTINE
!
!
!-------------------------------------------------------------------------
SUBROUTINE calc_outsourced (m,dvg,dng)
  !-----------------------------------------------------------------------
  !
  USE kinds,   ONLY : DP
  USE westcom, ONLY : npwqx      
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m
  COMPLEX(DP), INTENT(IN) :: dvg(npwqx,m)
  COMPLEX(DP), INTENT(OUT) :: dng(npwqx,m)
  !
  dng = dvg ! placeholder
  !
END SUBROUTINE 
