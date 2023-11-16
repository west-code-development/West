!
! Copyright (C) 2015-2023 M. Govoni
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
SUBROUTINE wstat_memory_report()
  !----------------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE wvfct,               ONLY : npwx,nbnd
  USE control_flags,       ONLY : gamma_only
  USE mp_global,           ONLY : nbgrp
  USE mp_world,            ONLY : mpime,root
  USE westcom,             ONLY : nbndval0x,n_pdep_basis,npwqx,logfile
  USE distribution_center, ONLY : pert
  USE noncollin_module,    ONLY : npol
  USE json_module,         ONLY : json_file
  !
  IMPLICIT NONE
  !
  TYPE(json_file) :: json
  INTEGER :: iunit
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  INTEGER :: nbndloc
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
  nbndloc = nbndval0x/nbgrp
  !
  mem_tot = 0.0_DP
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WSTAT global")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwqx*pert%nlocx
  WRITE(stdout,'(5x,"[MEM] dvg                     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwqx, pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.dvg', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwqx*pert%nlocx
  WRITE(stdout,'(5x,"[MEM] dng                     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwqx, pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.dng', mem_partial )
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
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Total estimate          ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM]")')
  !
  mem_tot = 0.0_DP
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WSTAT temporary")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*npol*nbndloc
  WRITE(stdout,'(5x,"[MEM] dvpsi                   ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, nbndloc
  IF( mpime == root ) CALL json%add( 'memory.dvpsi', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*npol*nbndloc
  WRITE(stdout,'(5x,"[MEM] dpsi                    ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, nbndloc
  IF( mpime == root ) CALL json%add( 'memory.dpsi', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  IF( .NOT. gamma_only ) THEN
     mem_partial = (1.0_DP/Mb)*complex_size*nbnd*npwx
     WRITE(stdout,'(5x,"[MEM] evckmq                  ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
        mem_partial, npwx, nbnd
     IF( mpime == root ) CALL json%add( 'memory.evckmq', mem_partial )
     mem_tot = mem_tot + mem_partial
  ENDIF
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*npol*nbndloc*4
  WRITE(stdout,'(5x,"[MEM] Sternheimer workspace   ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, nbndloc*4
  IF( mpime == root ) CALL json%add( 'memory.sternheimer', mem_partial )
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
