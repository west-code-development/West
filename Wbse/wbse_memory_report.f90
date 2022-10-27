!
! Copyright (C) 2015-2021 M. Govoni
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
SUBROUTINE wbse_init_memory_report()
  !----------------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE fft_base,            ONLY : dfftp
  USE mp_world,            ONLY : mpime,root
  USE westcom,             ONLY : nbnd_occ,logfile
  USE json_module,         ONLY : json_file
  !
  IMPLICIT NONE
  !
  TYPE(json_file) :: json
  INTEGER :: iunit
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  REAL(DP) :: mem_tot, mem_partial
  !
  CALL pw_memory_report()
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load( filename=TRIM(logfile) )
     CALL json%add( 'memory.units', 'Mb' )
     !
  ENDIF
  !
  mem_tot = 0.0_DP
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WBSE_INIT global")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = (1.0_DP/Mb)*complex_size*dfftp%nnr*nbnd_occ(1)
  WRITE(stdout,'(5x,"[MEM] evc_in_r                ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, dfftp%nnr, nbnd_occ(1)
  IF( mpime == root ) CALL json%add( 'memory.evc_in_r', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Total estimate          ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,*)
  !
  IF( mpime == root ) THEN
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     CALL json%destroy()
     !
  ENDIF
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_memory_report()
  !----------------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE wvfct,               ONLY : npwx
  USE control_flags,       ONLY : gamma_only
  USE mp_world,            ONLY : mpime,root
  USE westcom,             ONLY : nbnd_occ,n_pdep_basis,npwqx,logfile
  USE distribution_center, ONLY : pert
  USE noncollin_module,    ONLY : npol
  USE json_module,         ONLY : json_file
  !
  IMPLICIT NONE
  !
  TYPE(json_file) :: json
  INTEGER :: iunit
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  REAL(DP) :: mem_tot, mem_partial
  !
  CALL pw_memory_report()
  !
  IF( mpime == root ) THEN
     !
     CALL json%initialize()
     CALL json%load( filename=TRIM(logfile) )
     CALL json%add( 'memory.units', 'Mb' )
     !
  ENDIF
  !
  mem_tot = 0._DP
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WBSE global")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwqx*nbnd_occ(1)*pert%nlocx
  WRITE(stdout,'(5x,"[MEM] dvg_exc                 ",f10.2," Mb", 5x,"(",i7,",",i5,",",i5,")")') &
     mem_partial, npwqx, nbnd_occ(1), pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.dvg_exc', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwqx*nbnd_occ(1)*pert%nlocx
  WRITE(stdout,'(5x,"[MEM] dng_exc                 ",f10.2," Mb", 5x,"(",i7,","i5,",",i5,")")') &
     mem_partial, npwqx, nbnd_occ(1), pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.dng_exc', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  IF( gamma_only ) THEN
     mem_partial = (1.0_DP/Mb)*real_size*n_pdep_basis*pert%nlocx
  ELSE
     mem_partial = (1.0_DP/Mb)*complex_size*n_pdep_basis*pert%nlocx
  ENDIF
  WRITE(stdout,'(5x,"[MEM] hr_distr                ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, n_pdep_basis, pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.hr_distr', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  IF( gamma_only ) THEN
     mem_partial = (1.0_DP/Mb)*real_size*n_pdep_basis*pert%nlocx
  ELSE
     mem_partial = (1.0_DP/Mb)*complex_size*n_pdep_basis*pert%nlocx
  ENDIF
  WRITE(stdout,'(5x,"[MEM] vr_distr                ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, n_pdep_basis, pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.vr_distr', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Total estimate          ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] ")')
  !
  mem_tot = 0._DP
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WBSE temporary")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*npol*nbnd_occ(1)
  WRITE(stdout,'(5x,"[MEM] dvpsi                   ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, nbnd_occ(1)
  IF( mpime == root ) CALL json%add( 'memory.dvpsi', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*npol*nbnd_occ(1)
  WRITE(stdout,'(5x,"[MEM] dpsi                    ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, nbnd_occ(1)
  IF( mpime == root ) CALL json%add( 'memory.dpsi', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Total estimate          ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,*)
  !
  IF( mpime == root ) THEN
     !
     OPEN( NEWUNIT=iunit,FILE=TRIM(logfile) )
     CALL json%print( iunit )
     CLOSE( iunit )
     CALL json%destroy()
     !
  ENDIF
  !
END SUBROUTINE
