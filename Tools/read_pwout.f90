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
! Huihuo Zheng, Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE read_pwout() ! to be sync'd with PW/src/read_file.f90
  !----------------------------------------------------------------------------
  !
  ! Wrapper routine, for compatibility
  !
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, tmp_dir, wfc_dir
  USE io_global,            ONLY : stdout, ionode
  USE buffers,              ONLY : open_buffer, close_buffer
  USE wvfct,                ONLY : nbnd, npwx, npw
  USE noncollin_module,     ONLY : npol
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE uspp,                 ONLY : becsum
  USE scf,                  ONLY : rho
  USE realus,               ONLY : betapointlist, &
                                   init_realspace_vars,real_space
  USE dfunct,               ONLY : newd
  USE ldaU,                 ONLY : lda_plus_u, U_projection
  USE pw_restart,           ONLY : pw_readfile
  USE control_flags,        ONLY : io_level
  USE klist,                ONLY : init_igk
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
  USE mp_images,            ONLY : my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_bcast
  USE wavefunctions_module, ONLY : evc
  USE gvect,                ONLY : ngm_g, ecutrho
  USE gvecs,                ONLY : ngms_g, dual
  USE gvecw,                ONLY : ecutwfc
  USE fft_base,             ONLY : dfftp
  USE fft_base,             ONLY : dffts
  USE wvfct,                ONLY : npwx
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE 
  INTEGER :: ierr
  LOGICAL :: exst
  CHARACTER( 256 )  :: dirname
  !
  !
  ierr = 0 
  !
  ! ... Read the contents of the xml data file
  !
  IF ( ionode ) WRITE( stdout, '(/,5x,A,/,5x,A)') &
     'Reading data from directory:', TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  CALL start_clock('read_xml')
  CALL read_xml_file ( )
  CALL stop_clock('read_xml')
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  CALL start_clock('read_wave')
  !
  wfc_dir = tmp_dir
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  IF( my_image_id == 0 ) THEN 
     CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  ENDIF 
  !
  ! ... Allocate and compute k+G indices and number of plane waves
  ! ... FIXME: should be read from file, not re-computed
  !
  CALL init_igk ( npwx, ngm, g, gcutw ) 
  !
  IF( my_image_id == 0 ) THEN
     CALL pw_readfile( 'wave', ierr )
  ENDIF
  ! pw
  CALL mp_bcast( ecutwfc,    0, inter_image_comm )
  CALL mp_bcast( ecutrho,    0, inter_image_comm )
  CALL mp_bcast( dual,       0, inter_image_comm )
  CALL mp_bcast( gamma_only, 0, inter_image_comm )
  CALL mp_bcast( dfftp%nr1,  0, inter_image_comm )
  CALL mp_bcast( dfftp%nr2,  0, inter_image_comm )
  CALL mp_bcast( dfftp%nr3,  0, inter_image_comm )
  CALL mp_bcast( ngm_g,      0, inter_image_comm )
  CALL mp_bcast( npw,        0, inter_image_comm )
  CALL mp_bcast( dffts%nr1,  0, inter_image_comm )
  CALL mp_bcast( dffts%nr2,  0, inter_image_comm )
  CALL mp_bcast( dffts%nr3,  0, inter_image_comm )
  CALL mp_bcast( ngms_g,     0, inter_image_comm )
  ! wfc 
  CALL mp_bcast(evc(:,:),0,inter_image_comm)
  !
  CALL stop_clock('read_wave')
  !
  ! ... Assorted initialization: pseudopotentials, PAW
  ! ... Not sure which ones (if any) should be done here
  !
  CALL init_us_1()
  !
  IF (lda_plus_u .AND. (U_projection == 'pseudo')) CALL init_q_aeps()
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  !
  IF ( real_space ) THEN
    CALL betapointlist()
    CALL init_realspace_vars()
    IF( ionode ) WRITE(stdout,'(5x,"Real space initialisation completed")')
  ENDIF
  CALL newd()
  !
  IF( my_image_id == 0 ) CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE
