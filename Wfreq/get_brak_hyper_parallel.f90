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
SUBROUTINE get_brak_hyper_parallel(dvpsi,NRHS,NLSTEPS,x,brak,idistr)
  !-----------------------------------------------------------------------
  !
  ! ... brak = < dvpsi | x >
  !
  USE kinds,                ONLY : DP 
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE random_numbers,       ONLY : randy
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left
  USE io_global,            ONLY : stdout, ionode
  USE gvect,                ONLY : g,nl,gstart,ngm_g,ig_l2g
  USE cell_base,            ONLY : tpiba2
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE constants,            ONLY : tpi,fpi
  USE pwcom,                ONLY : npw,npwx
!  USE control_flags,        ONLY : gamma_only 
  USE noncollin_module,     ONLY : noncolin,npol 
  USE class_idistribute,    ONLY : idistribute
  !
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute),INTENT(IN) :: idistr
  COMPLEX(DP) :: dvpsi(npwx*npol,idistr%nlocx)
  INTEGER,INTENT(IN) :: NRHS,NLSTEPS
  COMPLEX(DP),INTENT(IN) :: x(npwx*npol,NRHS,NLSTEPS)
  REAL(DP),INTENT(OUT) :: brak(idistr%nglob,NLSTEPS,NRHS)
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i3,ip,ig,i1_glob,i2_glob,il
  INTEGER :: nblock_i
  INTEGER,ALLOCATABLE :: perm(:)
  INTEGER,ALLOCATABLE :: lblock(:)
  INTEGER :: icycl,aux_int
  !
  REAL(DP), EXTERNAL :: DDOT
!  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock( "brak" )
  !
  ALLOCATE( perm(0:nimage-1), lblock(0:nimage-1) )
  ! 
  ! -------------------
  ! BEGIN 
  ! brak = < dvpsi, x > 
  ! -------------------
  !
  brak = 0._DP
  !
  ! Initialize to zero
  !
  lblock=0 
  perm=0
  !
  lblock(my_image_id)=idistr%nloc
  perm(my_image_id)=my_image_id
  !
  CALL mp_sum(lblock,inter_image_comm)
  CALL mp_sum(perm,inter_image_comm)
  !
  DO icycl=0,nimage-1
     !
     nblock_i=lblock(perm(my_image_id))
     !
     DO i1=1,nblock_i
        !
        i1_glob = idistr%l2g( i1, perm(my_image_id) )
        !
!        IF(gamma_only) THEN
           DO il=1,NLSTEPS
              DO i2=1,NRHS
                 brak(i1_glob,il,i2) = 2.0_DP * DDOT(2*npw,dvpsi(1,i1),1,x(1,i2,il),1)
                 IF(gstart==2) brak(i1_glob,il,i2) = brak(i1_glob,il,i2) - REAL(dvpsi(1,i1),KIND=DP) * REAL(x(1,i2,il),KIND=DP)
              ENDDO
           ENDDO
!        ELSE
!           DO il=1,NLSTEPS
!              DO i2=1,NRHS
!                 brak(i1_glob,il,i2) = ZDOTC(npw,dvpsi(1,i1),1,x(1,i2,il),1)
!              ENDDO
!           ENDDO
!           IF(noncolin) THEN
!              DO il=1,NLSTEPS
!                 DO i2=1,NRHS
!                    brak(i1_glob,il,i2) = brak(i1_glob,il,i2) + ZDOTC(npw,dvpsi(1+npwx,i1),1,x(1+npwx,i2,il),1)
!                 ENDDO
!              ENDDO
!           ENDIF
!        ENDIF
     ENDDO
     !
     ! Update perm array according to a circular_shift_left
     !
     aux_int = perm(0)
     !
     DO i1=0,nimage-2 
        perm(i1) = perm(i1+1)
     ENDDO
     perm(nimage-1)=aux_int
     !
     ! Cycle the dvpsi array 
     ! 
     CALL mp_circular_shift_left( dvpsi, icycl, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( perm, lblock ) 
  !
  ! Syncronize brak
  !
  CALL mp_sum(brak,intra_bgrp_comm) 
  !
  CALL stop_clock( "brak" )
  !
END SUBROUTINE 
!
!
!-----------------------------------------------------------------------
SUBROUTINE get_brak_hyper_parallel_complex(dvpsi,NRHS,NLSTEPS,x,brak,idistr)
  !-----------------------------------------------------------------------
  !
  ! ... brak = < dvpsi | x >
  !
  USE kinds,                ONLY : DP 
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE random_numbers,       ONLY : randy
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left
  USE io_global,            ONLY : stdout, ionode
  USE gvect,                ONLY : g,nl,gstart,ngm_g,ig_l2g
  USE cell_base,            ONLY : tpiba2
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE constants,            ONLY : tpi,fpi
  USE pwcom,                ONLY : npw,npwx
!  USE control_flags,        ONLY : gamma_only 
  USE noncollin_module,     ONLY : noncolin,npol 
  USE class_idistribute,    ONLY : idistribute
  !
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  TYPE(idistribute),INTENT(IN) :: idistr
  COMPLEX(DP) :: dvpsi(npwx*npol,idistr%nlocx)
  INTEGER,INTENT(IN) :: NRHS,NLSTEPS
  COMPLEX(DP),INTENT(IN) :: x(npwx*npol,NRHS,NLSTEPS)
  COMPLEX(DP),INTENT(OUT) :: brak(idistr%nglob,NLSTEPS,NRHS)
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i3,ip,ig,i1_glob,i2_glob,il
  INTEGER :: nblock_i
  INTEGER,ALLOCATABLE :: perm(:)
  INTEGER,ALLOCATABLE :: lblock(:)
  INTEGER :: icycl,aux_int
  !
!  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock( "brak" )
  !
  ALLOCATE( perm(0:nimage-1), lblock(0:nimage-1) )
  ! 
  ! -------------------
  ! BEGIN 
  ! brak = < dvpsi, x > 
  ! -------------------
  !
  brak = 0._DP
  !
  ! Initialize to zero
  !
  lblock=0 
  perm=0
  !
  lblock(my_image_id)=idistr%nloc
  perm(my_image_id)=my_image_id
  !
  CALL mp_sum(lblock,inter_image_comm)
  CALL mp_sum(perm,inter_image_comm)
  !
  DO icycl=0,nimage-1
     !
     nblock_i=lblock(perm(my_image_id))
     !
     DO i1=1,nblock_i
        !
        i1_glob = idistr%l2g( i1, perm(my_image_id) )
        !
!        IF(gamma_only) THEN
!           DO il=1,NLSTEPS
!              DO i2=1,NRHS
!                 brak(i1_glob,il,i2) = 2.0_DP * DDOT(2*npw,dvpsi(1,i1),1,x(1,i2,il),1)
!                 IF(gstart==2) brak(i1_glob,il,i2) = brak(i1_glob,il,i2) - REAL(dvpsi(1,i1),KIND=DP) * REAL(x(1,i2,il),KIND=DP)
!              ENDDO
!           ENDDO
!        ELSE
           DO il=1,NLSTEPS
              DO i2=1,NRHS
                 brak(i1_glob,il,i2) = ZDOTC(npw,dvpsi(1,i1),1,x(1,i2,il),1)
              ENDDO
           ENDDO
           IF(noncolin) THEN
              DO il=1,NLSTEPS
                 DO i2=1,NRHS
                    brak(i1_glob,il,i2) = brak(i1_glob,il,i2) + ZDOTC(npw,dvpsi(1+npwx,i1),1,x(1+npwx,i2,il),1)
                 ENDDO
              ENDDO
           ENDIF
!        ENDIF
     ENDDO
     !
     ! Update perm array according to a circular_shift_left
     !
     aux_int = perm(0)
     !
     DO i1=0,nimage-2 
        perm(i1) = perm(i1+1)
     ENDDO
     perm(nimage-1)=aux_int
     !
     ! Cycle the dvpsi array 
     ! 
     CALL mp_circular_shift_left( dvpsi, icycl, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( perm, lblock ) 
  !
  ! Syncronize brak
  !
  CALL mp_sum(brak,intra_bgrp_comm) 
  !
  CALL stop_clock( "brak" )
  !
END SUBROUTINE 
