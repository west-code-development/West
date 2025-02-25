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
SUBROUTINE pw_memory_report()
  !----------------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE wvfct,               ONLY : npwx,nbnd
  USE fft_base,            ONLY : dfftp
  USE gvect,               ONLY : ngl
  USE gvecs,               ONLY : ngms
  USE uspp,                ONLY : nkb
  USE mp_world,            ONLY : mpime,root
  USE westcom,             ONLY : logfile
  USE json_module,         ONLY : json_file
  !
  IMPLICIT NONE
  !
  TYPE(json_file) :: json
  INTEGER :: iunit
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  REAL(DP) :: mem_tot, mem_partial
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
  WRITE(stdout,'(/,5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: QE")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = (1.0_DP/Mb)*complex_size*nbnd*npwx
  WRITE(stdout,'(5x,"[MEM] Kohn-Sham Wavefunctions ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx, nbnd
  IF( mpime == root ) CALL json%add( 'memory.evc', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*complex_size*nkb*npwx
  WRITE(stdout,'(5x,"[MEM] NL pseudopotentials     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx, nkb
  IF( mpime == root ) CALL json%add( 'memory.nlpp', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*complex_size*dfftp%nnr
  WRITE(stdout,'(5x,"[MEM] Each V/rho on FFT grid  ",f10.2," Mb", 5x,"(",i7,")")') &
     mem_partial, dfftp%nnr
  IF( mpime == root ) CALL json%add( 'memory.rhor', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*real_size*ngms
  WRITE(stdout,'(5x,"[MEM] Each G-vector array     ",f10.2," Mb", 5x,"(",i7,")")') &
     mem_partial, ngms
  IF( mpime == root ) CALL json%add( 'memory.rhog', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = (1.0_DP/Mb)*real_size*ngl
  WRITE(stdout,'(5x,"[MEM] G-vector shells         ",f10.2," Mb", 5x,"(",i7,")")') &
     mem_partial, ngl
  IF( mpime == root ) CALL json%add( 'memory.gshells', mem_partial )
  mem_tot = mem_tot + mem_partial
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Total estimate          ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM]")')
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
