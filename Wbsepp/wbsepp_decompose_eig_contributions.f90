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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
SUBROUTINE wbsepp_decompose_eig_contributions ( )
  !
  ! ... This pp reads eig-values and -vectors from davidson diago
  ! ... and project them on KS empty eig-vectors
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi
  USE io_global,            ONLY : stdout
  USE distribution_center,  ONLY : pert,aband
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE westcom,              ONLY : nbnd_occ,ev, &
                                   dvg_exc, d0psi, n_plep_read_from_file
  !wbsecom combined into westcom
  !USE wbsecom,              ONLY : dvg_exc, d0psi, n_plep_read_from_file
  USE pwcom,                ONLY : nks,npwx,npw
  USE plep_db,              ONLY : plep_db_read
  USE mp_world,             ONLY : mpime
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_image_comm, intra_bgrp_comm
  USE wavefunctions_module, ONLY : evc
  USE control_flags,        ONLY : gamma_only
  USE wvfct,                ONLY : nbnd
  USE gvect,                ONLY : gstart
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER    :: nvec,nbndx_occ,nbndx_emp
  INTEGER    :: iocc, iemp, ipol
  INTEGER    :: il1,ig1
  REAL(DP)   :: temp, anormr
  COMPLEX(DP):: anorm
  REAL(DP), ALLOCATABLE :: os_strength(:,:)
  REAL(DP), ALLOCATABLE :: projection_matrix(:,:,:)
  REAL(DP), EXTERNAL    :: DDOT
  INTEGER :: num_w, i
  REAL(DP):: w_start, w_end, delta_w, epsil, w_i
  REAL(DP), ALLOCATABLE :: delta_func(:)
  !
  ! ... INITIALIZATION
  !
  nvec  = n_plep_read_from_file
  nbndx_occ = MAXVAL(nbnd_occ)
  nbndx_emp = nbnd - nbndx_occ
  !
  ! CHECK IF KS EMPTY STATES WERE COMPUTED
  !
  IF (nbndx_emp < 1) THEN
     !
     CALL errore('chidiago', 'This eigenvectors decompositions needs KS empty stats, rerun pwscf with nbnd>nbnd_occ',1)
     !
  ENDIF
  !
  CALL start_clock( 'eig_decompose' )
  !
  ! ... DISTRIBUTE nvec
  !
  pert = idistribute()
  CALL pert%init(nvec,'i','nvec',.TRUE.)
  !
  ! ... DISTRIBUTE nband
  !
  aband = idistribute()
  CALL aband%init(nbndx_occ,'b','band_paralel',.TRUE.)
  !
  CALL wbse_memory_report() ! Before allocating I report the memory required.
  !
  ! READ EIGENVALUES AND VECTORS FROM OUTPUT
  !
  CALL plep_db_read( n_plep_read_from_file )
  !
  ! CHECK IF KS EMPTY STATES WERE COMPUTED
  !
  IF (nbnd <= nbndx_occ) THEN
     !
     CALL errore('eig_decompose', 'This eigenvectors decompositions needs empty stats',1)
     !
  ENDIF
  !
  ALLOCATE (d0psi(npwx,nbndx_occ,nks,3))
  !
  CALL solve_e_psi()
  !
  ALLOCATE (projection_matrix(nbndx_occ, nbndx_emp, nvec))
  ALLOCATE (os_strength(3, nvec))
  !
  DO il1 = 1, pert%nloc
     !
     ig1 = pert%l2g(il1)
     !
     IF( ig1 < 1 .OR. ig1 > nvec ) CYCLE
     !
     DO iocc = 1, nbndx_occ
        !
        DO iemp = 1, nbndx_emp
           !
           IF (gamma_only) THEN
              !
              anorm = 0.0_DP
              !
              anorm = 2.0_DP * DDOT(2*npw,dvg_exc(1,iocc,1,il1),1,evc(1,nbndx_occ+iemp),1)
              !
              IF (gstart==2) THEN
                 !
                 anorm = anorm - DBLE(dvg_exc(1,iocc,1,il1)) * DBLE(evc(1,nbndx_occ+iemp))
                 !
              ENDIF
              !
           ELSE
              !
              anorm = 0.0_DP
              !
              anorm  = DDOT(2*npw,dvg_exc(1,iocc,1,il1),1,evc(1,nbndx_occ+iemp),1)
              !
           ENDIF
           !
           CALL mp_sum(anorm,intra_bgrp_comm)
           !
           projection_matrix(iocc,iemp,ig1) = SQRT(DBLE(anorm)*DBLE(anorm) + AIMAG(anorm)*AIMAG(anorm))
           !
        ENDDO
        !
     ENDDO
     !
     DO ipol = 1, 3
        !
        anormr = 0.0_DP
        !
        DO iocc = 1, nbndx_occ
           !
           IF (gamma_only) THEN
              !
              anormr = anormr + 2.0_DP*DDOT(2*npw,dvg_exc(1,iocc,1,il1),1,d0psi(1,iocc,1,ipol),1 )
              !
              IF (gstart==2) THEN
                 !
                 anormr = anormr - DBLE(dvg_exc(1,iocc,1,il1)) * DBLE(d0psi(1,iocc,1,ipol))
                 !
              ENDIF
              !
           ELSE
              !
              anormr  = anormr + DDOT(2*npw,dvg_exc(1,iocc,1,il1),1,d0psi(1,iocc,1,ipol),1)
              !
           ENDIF
           !
        ENDDO
        !
        CALL mp_sum(anormr,intra_bgrp_comm)
        !
        !os_strength(ipol, ig1) =  anormr*anormr*2.0d0/(pi*3)
        !
        ! see the formular in Eq.26, JCTC, 12, 3969 (2016)
        os_strength(ipol, ig1) =  anormr*anormr*2.0d0/3.0d0
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(projection_matrix, inter_image_comm)
  CALL mp_sum(os_strength,       inter_image_comm)
  !
  ! ... Print out results
  !
  WRITE(stdout, "( /,5x,' *----------*   THE PRINCIPLE PROJECTED COMPONENTS   *----------*')")
  !
  DO ig1 = 1, nvec
     !
     WRITE(stdout, &
         & "(   /, 5x, ' #     Exciton : | ', i8,' |','   ','Excitation energy : | ', f12.6)") &
         ig1, ev(ig1)
     !
     WRITE(stdout, "( /,5x,' OSCILLATOR STRENGTH : ')")
     WRITE(stdout, "( /,15x,' XX ', 12x, 'YY', 12x, 'ZZ' 12x, 'TOTAL' )")
     WRITE(stdout, "( /,15x, f12.6, 5x, f12.6, 5x, f12.6, 5x, f12.6)") &
          os_strength(1, ig1)*ev(ig1), os_strength(2, ig1)*ev(ig1), os_strength(3, ig1)*ev(ig1), &
          (os_strength(1, ig1)+ os_strength(2, ig1)+ os_strength(3, ig1))*ev(ig1)
     !
     WRITE(stdout, &
         & "( /,  10x,'     OCCUPIED STATES', 5x, 'UNOCCUPIED STATES', 5x, 'PROJECTION' )")
     DO iocc = 1, nbndx_occ
        !
        DO iemp = 1, nbndx_emp
           !
           temp = projection_matrix(iocc,iemp,ig1)
           !
           IF (temp >= 0.1) THEN
              !
              WRITE(stdout, "(16x, i8, 16x, i8, 14x, f12.6, 2x)" ) iocc, iemp, temp
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  WRITE(stdout, "( /,5x, ' *----------*   DETAIL OF PROJECTED COMPONENTS   *----------*' )")
  !
  DO ig1 = 1, nvec
     !
     WRITE(stdout, &
         & "(   /, 5x, ' #     Exciton : | ', i8,' |','   ','Excitation energy : | ', f12.6)") &
         ig1, ev(ig1)
     !
     WRITE(stdout, &
         & "( /,  10x,'     OCCUPIED STATES', 5x, 'UNOCCUPIED STATES', 5x, 'PROJECTION' )")
     !
     DO iocc = 1, nbndx_occ
        !
        DO iemp = 1, nbndx_emp
           !
           WRITE(stdout, "(16x, i8, 16x, i8, 14x, f12.6, 2x)" ) iocc, iemp, projection_matrix(iocc,iemp,ig1)
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  w_start = 0
  w_end   = 1.0
  delta_w = 0.001
  epsil   = 0.01
  num_w   = INT((w_end - w_start)/delta_w)
  ALLOCATE (delta_func(num_w))
  w_i     = 0.0
  delta_func(:) = 0.0
  DO ig1 = 1, nvec
     DO i = 1, num_w
        w_i = w_start + delta_w*(i-1)
        IF (w_i <= w_end) THEN
           anormr = (os_strength(1, ig1)+ os_strength(2, ig1)+ os_strength(3, ig1))*ev(ig1) !w_i
           delta_func(i) = delta_func(i) + anormr*(epsil/((w_i-ev(ig1))*(w_i-ev(ig1)) + epsil*epsil))*(1.0/pi)
        ENDIF
     ENDDO
  ENDDO
  !
  DO i = 1, num_w
     WRITE(stdout,*) w_start + delta_w*(i-1), delta_func(i)
  ENDDO
  !
  DEALLOCATE(delta_func)
  !
  DEALLOCATE( projection_matrix )
  DEALLOCATE( os_strength )
  DEALLOCATE( ev )
  DEALLOCATE( dvg_exc )
  DEALLOCATE( d0psi )
  !
  CALL stop_clock( 'eig_decompose')
  !
  RETURN
  !
END SUBROUTINE
