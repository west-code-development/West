!
! Copyright (C) 2015-2025 M. Govoni
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
SUBROUTINE do_setup
  !-----------------------------------------------------------------------
  !
  USE json_module,            ONLY : json_file
  USE pwcom,                  ONLY : ngk,nbnd,nkstot,nspin,nelec,nelup,neldw,isk
  USE fixed_occ,              ONLY : f_inp
  USE kinds,                  ONLY : DP
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : intra_bgrp_comm,nproc_bgrp,me_bgrp
  USE io_global,              ONLY : stdout
  USE lsda_mod,               ONLY : lsda
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin,npol,lspinorb
  USE cell_base,              ONLY : omega,celldm,at,bg,tpiba
  USE fft_base,               ONLY : dfftp,dffts
  USE gvect,                  ONLY : ngm,ecutrho
  USE gvecw,                  ONLY : ecutwfc
  USE io_push,                ONLY : io_push_title,io_push_value,io_push_bar,io_push_es0
  USE westcom,                ONLY : logfile
  USE mp_world,               ONLY : mpime,root
  USE types_bz_grid,          ONLY : k_grid,q_grid
  !
  IMPLICIT NONE
  !
  TYPE(json_file) :: json
  INTEGER :: iunit
  INTEGER :: ik, iq, iks, spin, ip
  INTEGER, ALLOCATABLE :: ngm_i(:), npw_i(:)
  REAL(DP) :: alat
  CHARACTER(LEN=6) :: cik, ciq, cip
  !
  CALL start_clock('do_setup')
  !
  ! INIT PW
  !
  CALL init_pw_arrays(nbnd)
  !
  CALL set_dirs()
  !
  ! INIT K, Q GRIDS
  !
  CALL k_grid%init('K')
  !
  CALL q_grid%init('Q')
  !
  IF ( ANY ( (q_grid%ngrid(:) - k_grid%ngrid(:)) /= 0   ) ) THEN
     CALL errore('do_setup','q-point grid must be the same as k-point grid',1)
  ENDIF
  !
  IF( mpime == root ) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
  !
  IF ( lsda ) THEN
     IF ( INT( nelup ) == 0 .AND. INT( neldw ) == 0 ) THEN
        DO iks = 1, k_grid%nps
           spin = isk(iks)
           !
           SELECT CASE(spin)
           CASE(1)
              nelup = SUM( f_inp(:,1) )
           CASE(2)
              neldw = SUM( f_inp(:,2) )
           END SELECT
           !
        ENDDO
     ENDIF
     IF ( INT( nelup ) == 0 .AND. INT( neldw ) == 0 ) THEN
        CALL errore( 'do_setup','nelup = 0 and neldw = 0 ',1)
     ENDIF
  ENDIF
  !
  ! SYSTEM OVERVIEW
  !
  ALLOCATE( npw_i(0:nproc_bgrp-1), ngm_i(0:nproc_bgrp-1) )
  npw_i = 0
  ngm_i = 0
  npw_i(me_bgrp) = ngk(1)
  ngm_i(me_bgrp) = ngm
  CALL mp_sum( npw_i, intra_bgrp_comm )
  CALL mp_sum( ngm_i, intra_bgrp_comm )
  IF( mpime == root ) THEN
     DO ip = 0, nproc_bgrp-1
        WRITE(cip,'(i6)') ip+1
        CALL json%add('system.basis.npw.proc('//TRIM(ADJUSTL(cip))//')',npw_i(ip))
        CALL json%add('system.basis.ngm.proc('//TRIM(ADJUSTL(cip))//')',ngm_i(ip))
        CALL json%add('system.basis.npw.min',MINVAL(npw_i(:)))
        CALL json%add('system.basis.npw.max',MAXVAL(npw_i(:)))
        CALL json%add('system.basis.npw.sum',SUM(npw_i(:)))
        CALL json%add('system.basis.ngm.min',MINVAL(ngm_i(:)))
        CALL json%add('system.basis.ngm.max',MAXVAL(ngm_i(:)))
        CALL json%add('system.basis.ngm.sum',SUM(ngm_i(:)))
     ENDDO
  ENDIF
  DEALLOCATE( npw_i, ngm_i )
  !
  CALL io_push_title('System Overview')
  CALL io_push_value('gamma_only',gamma_only,20)
  IF( mpime == root ) CALL json%add('system.basis.gamma_only',gamma_only)
  CALL io_push_value('ecutwfc [Ry]',ecutwfc,20)
  IF( mpime == root ) CALL json%add('system.basis.ecutwfc:ry',ecutwfc)
  CALL io_push_value('ecutrho [Ry]',ecutrho,20)
  IF( mpime == root ) CALL json%add('system.basis.ecutrho:ry',ecutrho)
  CALL io_push_es0('omega [au^3]',omega,20)
  IF( mpime == root ) CALL json%add('system.cell.units','a.u.')
  IF( mpime == root ) CALL json%add('system.cell.omega',omega)
  CALL io_push_value('nbnd',nbnd,20)
  IF( mpime == root ) CALL json%add('system.electron.nbnd',nbnd)
  CALL io_push_value('nkstot',nkstot,20)
  IF( mpime == root ) CALL json%add('system.electron.nkstot',nkstot)
  CALL io_push_value('nspin',nspin,20)
  IF( mpime == root ) CALL json%add('system.electron.nspin',nspin)
  CALL io_push_value('nelec',nelec,20)
  IF( mpime == root ) CALL json%add('system.electron.nelec',nelec)
  IF(nspin == 2) THEN
     CALL io_push_value('nelup',nelup,20)
     IF( mpime == root ) CALL json%add('system.electron.nelup',nelup)
     CALL io_push_value('neldw',neldw,20)
     IF( mpime == root ) CALL json%add('system.electron.neldw',neldw)
  ENDIF
  CALL io_push_value('npol',npol,20)
  IF( mpime == root ) CALL json%add('system.electron.npol',npol)
  CALL io_push_value('lsda',lsda,20)
  IF( mpime == root ) CALL json%add('system.electron.lsda',lsda)
  CALL io_push_value('noncolin',noncolin,20)
  IF( mpime == root ) CALL json%add('system.electron.noncolin',noncolin)
  CALL io_push_value('lspinorb',lspinorb,20)
  IF( mpime == root ) CALL json%add('system.electron.lspinorb',lspinorb)
  CALL io_push_bar
  !
  alat = celldm(1)
  !
  WRITE( stdout, '(/5x,"3DFFT grid")')
  WRITE( stdout, '( 8x,"s : (",i4,",",i4,",",i4,")")') dffts%nr1, dffts%nr2, dffts%nr3
  WRITE( stdout, '( 8x,"p : (",i4,",",i4,",",i4,")")') dfftp%nr1, dfftp%nr2, dfftp%nr3
  WRITE( stdout, '(/5x,"Direct Lattice Cell [a.u.]")')
  WRITE( stdout, '( 8x,"a1 = (",3f14.7,")")') alat*at(1:3,1)
  WRITE( stdout, '( 8x,"a2 = (",3f14.7,")")') alat*at(1:3,2)
  WRITE( stdout, '( 8x,"a3 = (",3f14.7,")")') alat*at(1:3,3)
  WRITE( stdout, '(/5x,"Reciprocal Lattice Cell [a.u.]")')
  WRITE( stdout, '( 8x,"b1 = (",3f14.7,")")') tpiba*bg(1:3,1)
  WRITE( stdout, '( 8x,"b2 = (",3f14.7,")")') tpiba*bg(1:3,2)
  WRITE( stdout, '( 8x,"b3 = (",3f14.7,")")') tpiba*bg(1:3,3)
  WRITE( stdout, *)
  IF( mpime == root ) THEN
     CALL json%add('system.3dfft.s',(/ dffts%nr1, dffts%nr2, dffts%nr3 /) )
     CALL json%add('system.3dfft.p',(/ dfftp%nr1, dfftp%nr2, dfftp%nr3 /) )
     CALL json%add('system.cell.a1',alat*at(1:3,1))
     CALL json%add('system.cell.a2',alat*at(1:3,2))
     CALL json%add('system.cell.a3',alat*at(1:3,3))
     CALL json%add('system.cell.b1',tpiba*bg(1:3,1))
     CALL json%add('system.cell.b2',tpiba*bg(1:3,2))
     CALL json%add('system.cell.b3',tpiba*bg(1:3,3))
     CALL json%add('system.cell.alat',alat)
     CALL json%add('system.cell.tpiba',tpiba)
  ENDIF
  !
  WRITE( stdout, '(/5x,"Brillouin Zone sampling [cryst. coord.]")')
  WRITE( stdout, * )
  DO ik = 1, k_grid%np
     WRITE( cik, '(i6)') ik
     WRITE( stdout, '(8x,"k(",i6.6,") = (",3f14.7,")")') ik, k_grid%p_cryst(1:3,ik)
     IF( mpime == root ) THEN
        CALL json%add('system.bzsamp.k('//TRIM(ADJUSTL(cik))//').id',ik)
        CALL json%add('system.bzsamp.k('//TRIM(ADJUSTL(cik))//').crystcoord',k_grid%p_cryst(1:3,ik))
     ENDIF
  ENDDO
  !
  ! q-point grid
  !
  IF (.NOT. gamma_only ) THEN
     WRITE( stdout, * )
     DO iq = 1, q_grid%np
        WRITE( ciq, '(i6)') iq
        WRITE( stdout, '(8x,"q(",i6.6,") = (",3f14.7,")")') iq, q_grid%p_cryst(1:3,iq)
        IF( mpime == root ) THEN
           CALL json%add('system.bzsamp.q('//TRIM(ADJUSTL(ciq))//').id',iq)
           CALL json%add('system.bzsamp.q('//TRIM(ADJUSTL(ciq))//').crystcoord',q_grid%p_cryst(1:3,iq))
        ENDIF
     ENDDO
  ENDIF
  !
  IF( mpime == root ) THEN
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     CALL json%destroy()
  ENDIF
  !
  CALL stop_clock('do_setup')
  !
END SUBROUTINE
