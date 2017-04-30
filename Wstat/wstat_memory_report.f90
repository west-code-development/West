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
SUBROUTINE wstat_memory_report()
  !----------------------------------------------------------------------------
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout
  USE wvfct,          ONLY : npwx, nbnd, nbndx
  USE basis,          ONLY : natomwfc
  USE fft_base,       ONLY : dfftp,dffts
  USE gvect,          ONLY : ngl, ngm, ngm_g
  USE gvecs,          ONLY : ngms_g, ngms
  USE uspp,           ONLY : nkb
  USE control_flags,  ONLY : isolve, nmix, gamma_only, lscf
  USE mp_global,      ONLY : np_ortho
  USE westcom,        ONLY : nbnd_occ,n_pdep_basis,npwq0x
  USE distribution_center,  ONLY : pert
  USE noncollin_module,      ONLY : noncolin,npol
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  REAL(DP) :: mem_tot, mem_partial
  !
  WRITE(stdout,'(/,5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: QE")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  ! the conversions to double prevent integer overflow in very large run
  !
  mem_tot = 0._DP
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = DBLE(complex_size*nbnd*npwx)/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] Kohn-Sham Wavefunctions ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx,nbnd
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = DBLE(complex_size*nkb*npwx)/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] NL pseudopotentials     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx, nkb
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = DBLE(complex_size*dfftp%nnr)/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] Each V/rho on FFT grid  ",f10.2," Mb", 5x,"(",i7,")")') &
     mem_partial, dfftp%nnr
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = DBLE(real_size*ngms)/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] Each G-vector array     ",f10.2," Mb", 5x,"(",i7,")")') &
     mem_partial, ngms
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = DBLE(real_size*ngl)/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] G-vector shells         ",f10.2," Mb", 5x,"(",i7,")")') &
     mem_partial, ngl
  mem_tot = mem_tot + mem_partial
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE( stdout, '(5x,"[MEM] TOT                     ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] ")')
  !
  !
  !
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WSTAT global")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  mem_tot = 0._DP
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = DBLE(complex_size*npwq0x*pert%nlocx)/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] dvg                     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwq0x, pert%nlocx
  mem_tot = mem_tot + mem_partial
  !
  mem_partial = DBLE(complex_size*npwq0x*pert%nlocx)/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] dng                     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwq0x, pert%nlocx
  mem_tot = mem_tot + mem_partial
  !
  IF( gamma_only ) THEN  
     mem_partial = DBLE(real_size*n_pdep_basis*pert%nlocx)/DBLE(Mb)
  ELSE
     mem_partial = DBLE(complex_size*n_pdep_basis*pert%nlocx)/DBLE(Mb)
  ENDIF
  WRITE( stdout, '(5x,"[MEM] hr_distr                ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, n_pdep_basis, pert%nlocx
  mem_tot = mem_tot + mem_partial 
  !
  IF( gamma_only ) THEN  
     mem_partial = DBLE(real_size*n_pdep_basis*pert%nlocx)/DBLE(Mb)
  ELSE
     mem_partial = DBLE(complex_size*n_pdep_basis*pert%nlocx)/DBLE(Mb)
  ENDIF
  WRITE( stdout, '(5x,"[MEM] vr_distr                ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, n_pdep_basis, pert%nlocx
  mem_tot = mem_tot + mem_partial 
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE( stdout, '(5x,"[MEM] TOT                     ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] ")')
  !
  !
  mem_tot = 0._DP
  !
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] **Memory** analysis: WSTAT temporary")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE(stdout,'(5x,"[MEM] Allocated arrays      ",5x,"est. size (Mb)", 5x,"dimensions")')
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
  mem_partial = DBLE(complex_size*npwx*npol*nbnd_occ(1))/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] dvpsi                   ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, nbnd_occ(1)
  mem_tot = mem_tot + mem_partial 
  !
  mem_partial = DBLE(complex_size*npwx*npol*nbnd_occ(1))/DBLE(Mb)
  WRITE( stdout, '(5x,"[MEM] dpsi                    ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
     mem_partial, npwx*npol, nbnd_occ(1)
  mem_tot = mem_tot + mem_partial
  !
!  mem_partial = DBLE(complex_size*dffts%nnr)/DBLE(Mb)
!  WRITE( stdout, '(5x,"[MEM] aux_r                   ",f10.2," Mb", 5x,"(",i7,")")') &
!     mem_partial, dffts%nnr
!  mem_tot1 = mem_tot1 + mem_partial 
!  !
!  mem_partial = DBLE(complex_size*npwx)/DBLE(Mb)
!  WRITE( stdout, '(5x,"[MEM] aux_g                   ",f10.2," Mb", 5x,"(",i7,")")') &
!     mem_partial, npwx
!  mem_tot1 = mem_tot1 + mem_partial 
!  !
!  IF(.NOT.gamma_only) THEN
!     mem_partial = DBLE(complex_size*npwx)/DBLE(Mb)
!     WRITE( stdout, '(5x,"[MEM] dpsic                   ",f10.2," Mb", 5x,"(",i7,")")') &
!        mem_partial, dffts%nnr
!     mem_tot1 = mem_tot1 + mem_partial 
!  ENDIF
!  !
!  mem_partial = DBLE(complex_size*npwx*pert%nlocx)/DBLE(Mb)
!  WRITE( stdout, '(5x,"[MEM] dhg                     ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
!     mem_partial, npwx, pert%nlocx
!  mem_tot2 = mem_tot2 + mem_partial
  !
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  WRITE( stdout, '(5x,"[MEM] Total estimate          ",f10.2," Mb", 5x)') mem_tot
  WRITE(stdout,'(5x,"[MEM] ----------------------------------------------------------")')
  !
END SUBROUTINE
