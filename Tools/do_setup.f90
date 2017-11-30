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
!-----------------------------------------------------------------------
SUBROUTINE do_setup
  !-----------------------------------------------------------------------
  !
  USE json_module,            ONLY : json_file
  USE pwcom,                  ONLY : npw,nbnd,nkstot,xk,wk,nspin,nelec,nelup,neldw,et,wg,&
                                   & lspinorb,domag,lsda,isk,nks,two_fermi_energies,ngk
  USE fixed_occ,              ONLY : tfixed_occ,f_inp
  USE kinds,                  ONLY : DP
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : intra_bgrp_comm,npool,nbgrp,nproc_bgrp,me_bgrp
  USE mp_pools,               ONLY : intra_pool_comm, inter_pool_comm, &
                                     my_pool_id, nproc_pool, kunit
  USE io_global,              ONLY : stdout
  USE lsda_mod,               ONLY : current_spin,lsda
  USE constants,              ONLY : rytoev
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin,npol
  USE cell_base,              ONLY : omega,celldm,at
  USE fft_base,               ONLY : dfftp,dffts
  USE gvecs,                  ONLY : ngms_g, ngms
  USE gvect,                  ONLY : ngm_g, ngm, ecutrho
  USE gvecw,                  ONLY : ecutwfc
  USE io_push
  USE westcom,                ONLY : logfile
  USE mp_world,               ONLY : mpime, root
  USE class_bz_grid,          ONLY : bz_grid
  USE types_bz_grid,          ONLY : k_grid, q_grid, kmq_grid, kpq_grid
  !
  IMPLICIT NONE
  !
  TYPE(json_file) :: json
  INTEGER :: iunit
  INTEGER :: auxi,ib
  INTEGER :: ipol,ik,iq,npwx_g, nkbl, nkl, nkr, iks, ike, spin, ip
  INTEGER,ALLOCATABLE :: ngm_i(:), npw_i(:) 
  INTEGER, ALLOCATABLE :: ngk_g(:)
!  REAL(DP) :: xkg(3)
  REAL(DP) :: alat
  CHARACTER(LEN=6) :: cik, ciq, cip
  !
  CALL start_clock('do_setup')
  !
  ! INIT PW
  !
  CALL init_pw_arrays(nbnd)
  CALL set_iks_l2g()
  !
  CALL set_dirs()
  !
  ! INIT K GRID
  !
  k_grid = bz_grid()
  CALL k_grid%init('K')
  !
  IF ( .NOT. gamma_only ) THEN
     !
     ! INIT Q GRID
     !
     ! initialize q-point grid
     q_grid = bz_grid()
     CALL q_grid%init('Q')
     IF ( ANY ( q_grid%ngrid(:) - k_grid%ngrid(:) /= 0   ) ) THEN
        CALL errore( 'do_setup','q-point grid must be the same as k-point grid ',1)
     ENDIF
     ! initialize (k-q) grid
     !kmq_grid = bz_grid()
     !CALL kmq_grid%init_kq( k_grid, q_grid, -1 )
     !! initialize (k-q) grid
     !kpq_grid = bz_grid()
     !CALL kpq_grid%init_kq( k_grid, q_grid, +1 )
  ENDIF
  !
  IF( mpime == root ) THEN 
     CALL json%initialize()
     CALL json%load_file(filename=TRIM(logfile))
  ENDIF
  !
  IF ( lsda ) THEN
     IF ( INT( nelup ) == 0 .AND. INT( neldw ) == 0 ) THEN
     !IF ( .NOT. two_fermi_energies ) THEN
        DO iks = 1, nks
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
  npw_i(me_bgrp) = npw
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
! IF ( gamma_only ) THEN
!    auxi = npw
!    CALL mp_sum(auxi,intra_bgrp_comm)
!    CALL io_push_value('glob. #G',auxi,20)
!    IF( mpime == root ) CALL json%add('system.basis.globg',auxi)
! ELSE
!    ALLOCATE( ngk_g(nkstot) )
!    !npool = nproc_image / nproc_pool
!    nkbl = nkstot / kunit
!    nkl = kunit * ( nkbl / npool )
!    nkr = ( nkstot - nkl * npool ) / kunit
!    IF ( my_pool_id < nkr ) nkl = nkl + kunit
!    iks = nkl*my_pool_id + 1
!    IF ( my_pool_id >= nkr ) iks = iks + nkr*kunit
!    ike = iks + nkl - 1
!    ngk_g = 0
!    ngk_g(iks:ike) = ngk(1:nks)
!    CALL mp_sum( ngk_g, inter_pool_comm )
!    CALL mp_sum( ngk_g, intra_pool_comm )
!    ngk_g = ngk_g / nbgrp
!    npwx_g = MAXVAL( ngk_g(1:nkstot) )
!    CALL io_push_value('glob. #PW',npwx_g,20)
!    IF( mpime == root ) CALL json%add('system.basis.globpw',npwx_g)
!    DEALLOCATE( ngk_g )
! ENDIF
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
  CALL io_push_value('domag',domag,20)
  IF( mpime == root ) CALL json%add('system.electron.domag',domag)
  CALL io_push_bar
  !
  alat = celldm(1)
  !
  WRITE( stdout, '(/5x,"sFFT : (",i4,",",i4,",",i4,")")') dffts%nr1, dffts%nr2, dffts%nr3
  WRITE( stdout, '(/5x,"pFFT : (",i4,",",i4,",",i4,")")') dfftp%nr1, dfftp%nr2, dfftp%nr3
  WRITE( stdout, '(/5x,"Cell [a.u.]          = ",3f14.6)') alat*at(1,1:3)
  WRITE( stdout, '( 5x,"                     = ",3f14.6)') alat*at(2,1:3)
  WRITE( stdout, '( 5x,"                     = ",3f14.6)') alat*at(3,1:3)
  WRITE( stdout, '( 5x," ")')
  IF( mpime == root ) THEN 
     CALL json%add('system.basis.sFFT',(/ dffts%nr1, dffts%nr2, dffts%nr3 /) )
     CALL json%add('system.basis.pFFT',(/ dfftp%nr1, dfftp%nr2, dfftp%nr3 /) )
     CALL json%add('system.cell.a1',alat*at(1:3,1))
     CALL json%add('system.cell.a2',alat*at(1:3,2))
     CALL json%add('system.cell.a3',alat*at(1:3,3))
     CALL json%add('system.cell.alat',alat)
  ENDIF
  !
  WRITE( stdout, '(5x,"number of ks points = ",i6)') k_grid%nps
  IF( mpime == root ) CALL json%add('system.kpt.nkstot',k_grid%nps)
  WRITE( stdout, '(23x,"cart. coord. in units 2pi/alat")')
  DO iks = 1, k_grid%nps
     ik = k_grid%ip(iks)
     WRITE( cik, '(i6)') ik
     WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') iks, &
          (k_grid%p_cart(ipol,ik) , ipol = 1, 3) , k_grid%weight(iks)
     IF( mpime == root ) CALL json%add('system.kpt.k('//TRIM(ADJUSTL(cik))//').cartcoord:tpiba',k_grid%p_cart(1:3,ik))
     IF( mpime == root ) CALL json%add('system.kpt.k('//TRIM(ADJUSTL(cik))//').weight',k_grid%weight(iks))
  ENDDO
  WRITE( stdout, '(/23x,"cryst. coord.")')
  DO ik = 1, k_grid%nps
     ik = k_grid%ip(iks)
     WRITE( cik, '(i6)') ik
     WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
          ik, (k_grid%p_cryst(ipol,ik) , ipol = 1, 3) , k_grid%weight(iks)
     IF( mpime == root ) CALL json%add('system.kpt.k('//TRIM(ADJUSTL(cik))//').crystcoord',k_grid%p_cryst(1:3,ik))
  ENDDO
  !
  ! q-point grid
  !
  IF (.NOT. gamma_only ) THEN
     WRITE( stdout, * )
     WRITE( stdout, '(5x,"number of q points = ",i6)') q_grid%np
     IF( mpime == root ) CALL json%add('system.qpt.nqtot',q_grid%np)
     WRITE( stdout, '(23x,"cart. coord. in units 2pi/alat")')
     DO iq = 1, q_grid%np
        WRITE( ciq, '(i6)') iq
        WRITE( stdout, '(8x,"q(",i5,") = (",3f12.7,")")') iq, &
             (q_grid%p_cart(ipol, iq) , ipol = 1, 3)
        IF( mpime == root ) CALL json%add('system.qpt.q('//TRIM(ADJUSTL(ciq))//').cartcoord:tpiba',q_grid%p_cart(1:3,iq))
     ENDDO
     WRITE( stdout, '(/23x,"cryst. coord.")')
     DO iq = 1, q_grid%np
        WRITE( ciq, '(i6)') iq
        WRITE( stdout, '(8x,"q(",i5,") = (",3f12.7,")")') &
             iq, (q_grid%p_cryst(ipol,iq) , ipol = 1, 3)
        IF( mpime == root ) CALL json%add('system.qpt.q('//TRIM(ADJUSTL(ciq))//').crystcoord',q_grid%p_cryst(1:3,iq))
     ENDDO
  ENDIF
  !
  !
  IF( mpime == root ) THEN
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print_file( iunit )
     CLOSE( iunit )
     CALL json%destroy()
  ENDIF 
  !
  !
  CALL stop_clock('do_setup')
  !
END SUBROUTINE 
