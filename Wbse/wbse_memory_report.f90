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
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_memory_report()
  !----------------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : stdout
  USE wvfct,               ONLY : npwx
  USE control_flags,       ONLY : gamma_only
  USE mp_global,           ONLY : nbgrp
  USE mp_world,            ONLY : mpime,root
  USE pwcom,               ONLY : nks
  USE westcom,             ONLY : l_bse,l_hybrid_tddft,l_lanczos,nbnd_occ,n_trunc_bands,&
                                & n_pdep_basis,npwqx,logfile
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
  nbndloc = (nbnd_occ(1)-n_trunc_bands-1)/nbgrp+1
  !
  IF( .NOT. l_lanczos ) THEN
     mem_tot = 0._DP
     WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
     WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WBSE global")')
     WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
     WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
     WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
     !
     mem_partial = (1.0_DP/Mb)*complex_size*npwqx*nbndloc*pert%nlocx
     WRITE(stdout,'(5x,"[MEM] dvg_exc                 ",f10.2," Mb", 5x,"(",i7,",",i5,",",i5,")")') &
        mem_partial, npwqx, nbndloc, pert%nlocx
     IF( mpime == root ) CALL json%add( 'memory.dvg_exc', mem_partial )
     mem_tot = mem_tot + mem_partial
     !
     mem_partial = (1.0_DP/Mb)*complex_size*npwqx*nbndloc*pert%nlocx
     WRITE(stdout,'(5x,"[MEM] dng_exc                 ",f10.2," Mb", 5x,"(",i7,","i5,",",i5,")")') &
        mem_partial, npwqx, nbndloc, pert%nlocx
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
  ENDIF
  !
  mem_tot = 0._DP
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WBSE temporary")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  IF( .NOT. l_lanczos ) THEN
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
  ELSE
     mem_partial = (1.0_DP/Mb)*complex_size*npwx*nbndloc*nks*3
     WRITE(stdout,'(5x,"[MEM] d0psi                   ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
        mem_partial, npwx, nbndloc*nks*3
     IF( mpime == root ) CALL json%add( 'memory.d0psi', mem_partial )
     mem_tot = mem_tot + mem_partial
     !
     mem_partial = (1.0_DP/Mb)*complex_size*npwx*nbndloc
     WRITE(stdout,'(5x,"[MEM] evc1                    ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
        mem_partial, npwx, nbndloc
     IF( mpime == root ) CALL json%add( 'memory.evc1', mem_partial )
     mem_tot = mem_tot + mem_partial
     !
     mem_partial = (1.0_DP/Mb)*complex_size*npwx*nbndloc*2
     WRITE(stdout,'(5x,"[MEM] Lanczos workspace       ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
        mem_partial, npwx, nbndloc*2
     IF( mpime == root ) CALL json%add( 'memory.lanczos', mem_partial )
     mem_tot = mem_tot + mem_partial
  ENDIF
  !
  mem_partial = (1.0_DP/Mb)*complex_size*npwx*nbndloc*2
  WRITE(stdout,'(5x,"[MEM] Liouville workspace     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx, nbndloc*2
  IF( mpime == root ) CALL json%add( 'memory.liouville', mem_partial )
  mem_tot = mem_tot + mem_partial
  !
  IF( l_bse .OR. l_hybrid_tddft ) THEN
     mem_partial = (1.0_DP/Mb)*complex_size*npwx*(nbnd_occ(1)-n_trunc_bands)*2
     WRITE(stdout,'(5x,"[MEM] kernel                  ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
        mem_partial, npwx, (nbnd_occ(1)-n_trunc_bands)*2
     IF( mpime == root ) CALL json%add( 'memory.kernel', mem_partial )
     mem_tot = mem_tot + mem_partial
  ENDIF
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
