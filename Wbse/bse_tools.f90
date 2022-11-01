!
! Copyright (C) 2015-2022 M. Govoni
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
SUBROUTINE read_bse_pots_g2g( rhog, fixed_band_i, fixed_band_j, ispin, single_only)
  !
  ! ... this routine writes the pot-density in xml format
  !     seems the path is wbse_init_save_dir
  !
  USE kinds,          ONLY : DP
  USE pdep_io,        ONLY : pdep_read_G_and_distribute
  USE fft_at_gamma,   ONLY : single_invfft_gamma
  USE fft_at_k,       ONLY : single_invfft_k
  USE westcom,        ONLY : wbse_init_save_dir
  USE pwcom,          ONLY : npwx
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rhog(npwx)
  INTEGER, INTENT(IN) :: fixed_band_i, fixed_band_j, ispin
  LOGICAL, INTENT(IN) :: single_only
  !
  ! local variables
  !
  INTEGER             :: band_i, band_j
  CHARACTER(LEN=6)    :: my_labeli
  CHARACTER(LEN=6)    :: my_labelj
  CHARACTER(LEN=6)    :: my_spin
  CHARACTER(LEN=256)  :: file_base
  !
  IF (fixed_band_j < fixed_band_i) THEN
     !
     ! p_31 -> p_13
     !
     band_i = fixed_band_j
     !
     band_j = fixed_band_i
     !
  ELSE
     !
     band_i = fixed_band_i
     !
     band_j = fixed_band_j
     !
  ENDIF
  !
  WRITE (my_labeli,'(i6.6)') band_i
  WRITE (my_labelj,'(i6.6)') band_j
  WRITE (my_spin,  '(i1)') ispin
  !
  rhog = (0.0_DP, 0.0_DP)
  !
  file_base = TRIM(wbse_init_save_dir)//'/E'//TRIM(ADJUSTL(my_labeli))//&
           '_'//TRIM(ADJUSTL(my_labelj))//'_'//TRIM(ADJUSTL(my_spin))//'.dat'
  CALL pdep_read_G_and_distribute(file_base, rhog(:))
  !
END SUBROUTINE
!
SUBROUTINE read_bse_pots_g2r( rho_all, fixed_band_i, fixed_band_j, ispin, single_only)
  !
  ! ... this routine writes the pot-density in xml format
  !     call by wbse_bse_kernel -> wbse.x  but using wbse_init.x data
  !     folder: wbse_init_save_dir
  !
  USE kinds,          ONLY : DP
  USE control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY : dfftp, dffts
  USE pdep_io,        ONLY : pdep_read_G_and_distribute
  USE gvect,          ONLY : ngm, ngmx
  USE fft_at_gamma,   ONLY : single_invfft_gamma
  USE fft_at_k,       ONLY : single_invfft_k
  USE westcom,        ONLY : wbse_init_save_dir
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rho_all(dfftp%nnr)
  INTEGER, INTENT(IN) :: fixed_band_i, fixed_band_j, ispin
  LOGICAL, INTENT(IN) :: single_only
  !
  ! local variables
  !
  INTEGER             :: band_i, band_j
  REAL(DP)            :: rhoaux1(dfftp%nnr)
  REAL(DP)            :: rhoaux2(dfftp%nnr)
  COMPLEX(DP), ALLOCATABLE :: aux_g(:), aux_r(:)
  CHARACTER(LEN=6)    :: my_labeli
  CHARACTER(LEN=6)    :: my_labelj
  CHARACTER(LEN=6)    :: my_spin
  CHARACTER(LEN=256)  :: file_base
  !
  ALLOCATE (aux_g(ngmx), aux_r(dfftp%nnr))
  !
  rhoaux1(:) = 0.0_DP
  rhoaux2(:) = 0.0_DP
  !
  IF (fixed_band_j < fixed_band_i) THEN
     !
     ! p_31 -> p_13
     !
     band_i = fixed_band_j
     !
     band_j = fixed_band_i
     !
  ELSE
     !
     band_i = fixed_band_i
     !
     band_j = fixed_band_j
     !
  ENDIF
  !
  WRITE (my_labeli,'(i6.6)') band_i
  WRITE (my_labelj,'(i6.6)') band_j
  WRITE (my_spin,  '(i1)') ispin
  !
  aux_g = (0.0_DP, 0.0_DP)
  !
  file_base = TRIM(wbse_init_save_dir)//'/E'//TRIM(ADJUSTL(my_labeli))//&
              '_'//TRIM(ADJUSTL(my_labelj))//'_'//TRIM(ADJUSTL(my_spin))//'.dat'
  CALL pdep_read_G_and_distribute(file_base, aux_g(:))
  !
  ! G -> R
  !
  aux_r = (0.0_DP, 0.0_DP)
  !
  IF (gamma_only) THEN
     !
     CALL single_invfft_gamma(dffts,ngm,ngmx,aux_g,aux_r,'Rho')
     !
  ELSE
     !
     CALL single_invfft_k(dffts,ngm,ngmx,aux_g,aux_r,'Rho') ! no igk
     !
  ENDIF
  !
  rhoaux1(:) = REAL(aux_r(:),KIND=DP)
  !
  IF (.NOT. single_only) THEN
     !
     IF (fixed_band_j < (fixed_band_i+1)) THEN
        !
        ! p_21 -> p_12
        !
        band_i = fixed_band_j
        !
        band_j = fixed_band_i+1
        !
     ELSE
        !
        band_i = fixed_band_i+1
        !
        band_j = fixed_band_j
        !
     ENDIF
     !
     WRITE (my_labeli,'(i6.6)') band_i
     WRITE (my_labelj,'(i6.6)') band_j
     WRITE (my_spin,  '(i1)') ispin
     !
     aux_g = (0.0_DP, 0.0_DP)
     !
     file_base = TRIM(wbse_init_save_dir)//'/E'//TRIM(ADJUSTL(my_labeli))//&
                '_'//TRIM(ADJUSTL(my_labelj))//'_'//TRIM(ADJUSTL(my_spin))//'.dat'
     CALL pdep_read_G_and_distribute(file_base,aux_g(:))
     !
     ! G -> R
     !
     aux_r = (0.0_DP, 0.0_DP)
     !
     IF (gamma_only) THEN
        !
        CALL single_invfft_gamma(dffts,ngm,ngmx,aux_g,aux_r,'Rho')
        !
     ELSE
        !
        CALL single_invfft_k(dffts,ngm,ngmx,aux_g,aux_r,'Rho') ! no igk
        !
     ENDIF
     !
     rhoaux2(:) = REAL(aux_r(:),KIND=DP)
     !
  ENDIF
  !
  ! total : rho(:) = [rho_ij, rho_i+1,j]
  !
  rho_all(:) = (0.0_DP, 0.0_DP)
  !
  rho_all(:) = CMPLX (rhoaux1(:), rhoaux2(:), KIND=DP)
  !
  DEALLOCATE (aux_g, aux_r)
  !
END SUBROUTINE
!
SUBROUTINE write_umatrix_and_omatrix (oumat_dim,ispin,umatrix,omatrix)
  !
  ! ...   This subroutine writes orbs to unit emptyunitc0
  !
  USE kinds,          ONLY : DP
  USE iotk_module
  USE mp_world,       ONLY : world_comm
  USE io_global,      ONLY : stdout
  USE mp,             ONLY : mp_barrier
  USE mp_global,      ONLY : my_pool_id,my_bgrp_id,me_bgrp,root_bgrp
  USE westcom,        ONLY : wbse_init_save_dir
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,           INTENT(IN)   :: oumat_dim,ispin
  REAL(DP),          INTENT(IN)   :: omatrix(oumat_dim,oumat_dim)
  COMPLEX(DP),       INTENT(IN)   :: umatrix(oumat_dim,oumat_dim)
  !
  ! Workspace
  !
  INTEGER :: ierr, iunout
  CHARACTER(LEN=256)  :: filename
  !character (len=:), allocatable :: filename
  CHARACTER(LEN=4)    :: my_spin
  !
  ! BARRIER
  !
  CALL mp_barrier(world_comm)
  !
  WRITE(my_spin,'(i1)') ispin
  !
  filename =  TRIM(wbse_init_save_dir) // '/u_matrix.wan.occ.' // TRIM(my_spin) // '.dat'
  !
  WRITE(stdout,'(5X,"Writing Omatrix & Umatrix to ", A256 )') TRIM(filename)
  !
  IF(my_pool_id /= 0) RETURN
  IF(my_bgrp_id /= 0) RETURN
  !
  ! Resume all components
  !
  IF (me_bgrp==root_bgrp) THEN
     !
     !
     ! ... open XML descriptor
     !
     CALL iotk_free_unit(iunout,ierr)
     CALL iotk_open_write( iunout, FILE = TRIM(filename), BINARY = .False., IERR = ierr )
     !
     CALL iotk_write_begin(iunout,'OUMATRIX_SIZE')
     CALL iotk_write_dat(iunout,'oumat_dim',oumat_dim)
     CALL iotk_write_end(iunout,'OUMATRIX_SIZE')
     !
     CALL iotk_write_begin(iunout, 'UMATRIX_ELE')
     CALL iotk_write_dat(iunout, 'umat_ele', umatrix(:,:))
     CALL iotk_write_end(iunout, 'UMATRIX_ELE')
     !
     CALL iotk_write_begin(iunout, 'OMATRIX_ELE')
     CALL iotk_write_dat(iunout, 'omat_ele', omatrix(:,:))
     CALL iotk_write_end(iunout, 'OMATRIX_ELE')
     !
     CALL iotk_close_write(iunout)
     !
  ENDIF
  !
  ! BARRIER
  !
  CALL mp_barrier(world_comm)
  !
END SUBROUTINE
!
SUBROUTINE read_umatrix_and_omatrix (oumat_dim,ispin,umatrix,omatrix)
  !
  USE kinds,          ONLY : DP
  USE iotk_module
  USE io_global,      ONLY : stdout,ionode
  USE mp_world,       ONLY : world_comm,root
  USE mp,             ONLY : mp_bcast,mp_barrier
  USE mp_global,      ONLY : intra_image_comm
  USE westcom,        ONLY : wbse_init_save_dir
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,           INTENT(IN)   :: oumat_dim,ispin
  REAL(DP),          INTENT(INOUT):: omatrix(oumat_dim,oumat_dim)
  COMPLEX(DP),       INTENT(INOUT):: umatrix(oumat_dim,oumat_dim)
  !
  ! Workspace
  !
  INTEGER :: ierr,iunout
  INTEGER :: oumat_dim_tmp
  CHARACTER(LEN=256)  :: filename
  CHARACTER(LEN=4)    :: my_spin
  REAL(DP),    ALLOCATABLE :: omatrix_tmp(:,:)
  COMPLEX(DP), ALLOCATABLE :: umatrix_tmp(:,:)
  !
  !
  ! BARRIER
  !
  CALL mp_barrier(world_comm)
  !
  WRITE(my_spin,'(i1)') ispin
  filename = TRIM(wbse_init_save_dir)//'/u_matrix.wan.occ.'//TRIM(my_spin)//'.dat'
  !
  WRITE(stdout,'(/,5X,"Reading Omatrix & Umatrix from ", A256 )') filename
  !
  ierr = 0
  IF ( ionode ) THEN
     CALL iotk_free_unit( iunout, ierr )
     CALL iotk_open_read( iunout, FILE = TRIM( filename ), BINARY = .False., IERR = ierr )
  ENDIF
  !
  CALL mp_bcast( ierr, root, world_comm )
  !
  IF ( ierr /=0 ) CALL errore( 'write_umatrix_and_ovl_matrix', 'cannot open file for reading', ierr )
  !
  IF ( ionode ) THEN
     !
     CALL iotk_scan_begin(iunout, 'OUMATRIX_SIZE' )
     CALL iotk_scan_dat(  iunout, 'oumat_dim', oumat_dim_tmp)
     CALL iotk_scan_end(  iunout, 'OUMATRIX_SIZE' )
     !
  ENDIF
  !
  CALL mp_bcast( oumat_dim_tmp, 0, intra_image_comm )
  !
  ALLOCATE(umatrix_tmp(oumat_dim_tmp,oumat_dim_tmp))
  ALLOCATE(omatrix_tmp(oumat_dim_tmp,oumat_dim_tmp))
  !
  IF ( ionode ) THEN
     !
     CALL iotk_scan_begin(iunout, 'UMATRIX_ELE')
     CALL iotk_scan_dat(  iunout, 'umat_ele', umatrix_tmp(:,:))
     CALL iotk_scan_end(  iunout, 'UMATRIX_ELE')
     !
  ENDIF
  !
  IF ( ionode ) THEN
     !
     CALL iotk_scan_begin(iunout, 'OMATRIX_ELE')
     CALL iotk_scan_dat(  iunout, 'omat_ele', omatrix_tmp(:,:))
     CALL iotk_scan_end(  iunout, 'OMATRIX_ELE')
     !
  ENDIF
  !
  IF ( ionode ) CALL iotk_close_read( iunout )
  !
  CALL mp_bcast( umatrix_tmp, 0 , intra_image_comm )
  CALL mp_bcast( omatrix_tmp, 0 , intra_image_comm )
  !
  umatrix(:,:) = 0.0
  omatrix(:,:) = 0.0
  umatrix(1:oumat_dim_tmp,1:oumat_dim_tmp) = umatrix_tmp(1:oumat_dim_tmp,1:oumat_dim_tmp)
  omatrix(1:oumat_dim_tmp,1:oumat_dim_tmp) = omatrix_tmp(1:oumat_dim_tmp,1:oumat_dim_tmp)
  !
  DEALLOCATE(umatrix_tmp)
  DEALLOCATE(omatrix_tmp)
  !
END SUBROUTINE
!
SUBROUTINE read_pwscf_wannier_orbs ( ne, npw, c_emp, filename )
  !
  ! ... This subroutine reads wannier orbital from unit emptyunit
  !
  USE kinds,          ONLY : DP
  USE iotk_module
  USE io_global,      ONLY : stdout
  USE gvect,          ONLY : ig_l2g
  USE mp_wave,        ONLY : splitwf
  USE mp,             ONLY : mp_get,mp_size,mp_rank,mp_bcast,mp_max
  USE mp_global,      ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm, &
                           & my_pool_id,my_bgrp_id,inter_bgrp_comm,inter_pool_comm
  USE klist,          ONLY : ngk, igk_k
  USE wvfct,          ONLY : npwx
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT (IN)   :: npw,ne
  COMPLEX(DP), INTENT(INOUT) :: c_emp(npw,ne)
  CHARACTER(LEN=*),INTENT(IN):: filename
  !
  INTEGER :: ierr, iun, ig, i
  INTEGER :: ngw_l, ngw_g
  !
  CHARACTER(LEN=256)       :: fileempty
  COMPLEX(DP), ALLOCATABLE :: ctmp(:)
  INTEGER, ALLOCATABLE     :: igk_l2g(:)
  !
  ! ... Subroutine Body
  !
  ALLOCATE ( igk_l2g( npwx ) )
  !
  ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
  !
  igk_l2g = 0
  DO ig = 1, ngk(1)
     igk_l2g(ig) = ig_l2g(igk_k(ig,1))
  ENDDO
  !
  fileempty = TRIM(filename)
  !
  WRITE(stdout,'(/,5X,"Reading pwscf wannier orbs from ", A8 )') fileempty
  !
  ngw_l    = npw!MAXVAL(igk_l2g(:))
  ngw_g    = MAXVAL(igk_l2g(:))
  CALL mp_max(ngw_g,intra_bgrp_comm)
  !
  ALLOCATE( ctmp(ngw_g) )
  !
  IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
     !
     ! ONLY ROOT W/IN BGRP READS
     !
     IF (me_bgrp==root_bgrp) THEN
        !
        CALL iotk_free_unit( iun, ierr )
        CALL iotk_open_read( iun, FILE = TRIM(fileempty), BINARY = .TRUE., IERR = ierr)
        !
     ENDIF
     !
  ENDIF
  !
  CALL mp_bcast(ierr,0,inter_bgrp_comm)
  !
  IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
     !
     ! ONLY ROOT W/IN BGRP READS
     !
     IF (me_bgrp==root_bgrp) THEN
        !
        CALL iotk_scan_begin( iun, 'WANWFC_GSPACE' )
        !
     ENDIF
     !
  ENDIF
  !
  DO i = 1, ne
     !
     IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
        !
        ! ONLY ROOT W/IN BGRP READS
        !
        IF (me_bgrp==root_bgrp) THEN
           !
           CALL iotk_scan_dat( iun, &
                             'wfc' // iotk_index(i), ctmp(:))
           !
        ENDIF
        !
        CALL splitwf( c_emp(:,i), ctmp, ngw_l, igk_l2g, &
                         me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm)
        !
     ENDIF
     !
  ENDDO
  !
  IF (my_pool_id==0.AND.my_bgrp_id==0) THEN
     !
     ! ONLY ROOT W/IN BGRP READS
     !
     IF (me_bgrp==root_bgrp) THEN
        !
        CALL iotk_scan_end( iun, 'WANWFC_GSPACE' )
        CALL iotk_close_read( iun )
        !
     ENDIF
     !
  ENDIF
  !
  DEALLOCATE(ctmp, igk_l2g)
  !
  CALL mp_bcast(c_emp,0,inter_bgrp_comm)
  CALL mp_bcast(c_emp,0,inter_pool_comm)
  !
END SUBROUTINE
