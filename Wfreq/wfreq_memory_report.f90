!
! Copyright (C) 2015-2024 M. Govoni
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
SUBROUTINE wfreq_memory_report()
  !----------------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE wvfct,               ONLY : npwx,nbnd
  USE uspp,                ONLY : nkb
  USE control_flags,       ONLY : gamma_only
  USE mp_world,            ONLY : mpime,root
  USE westcom,             ONLY : npwqx,l_macropol,l_skip_nl_part_of_hcomr,n_lanczos,logfile
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
  mem_tot = 0.0_DP
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WFREQ temporary")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*npol*pert%nlocx
  WRITE(stdout,'(5x,"[MEM] dvpsi                   ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.dvpsi', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwqx*pert%nlocx
  WRITE(stdout,'(5x,"[MEM] pertg_all               ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwqx, pert%nlocx
  IF( mpime == root ) CALL json%add( 'memory.pertg_all', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  IF( .NOT. gamma_only ) THEN
     mem_partial = (1.0_DP/Mb)*complex_size*nbnd*npwx
     WRITE(stdout,'(5x,"[MEM] evckpq                  ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
        mem_partial, npwx, nbnd
     IF( mpime == root ) CALL json%add( 'memory.evckpq', mem_partial )
     mem_tot = mem_tot + mem_partial
  ENDIF
  !
  IF( l_macropol .AND. .NOT. l_skip_nl_part_of_hcomr .AND. nkb > 0 ) THEN
     mem_partial = (1.0_DP/Mb)*complex_size*nkb*npwx*2
     WRITE(stdout,'(5x,"[MEM] Hr commutator workspace ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
        mem_partial, npwx, nkb*2
     IF( mpime == root ) CALL json%add( 'memory.hcomr', mem_partial )
     mem_tot = mem_tot + mem_partial
  ENDIF
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*npol*pert%nlocx*n_lanczos
  WRITE(stdout,'(5x,"[MEM] Lanczos workspace       ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, pert%nlocx*n_lanczos
  IF( mpime == root ) CALL json%add( 'memory.lanczos', mem_partial )
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
