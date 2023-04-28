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
MODULE wbse_io
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
  !
  SUBROUTINE read_bse_pots_g(rhog,fixed_band_i,fixed_band_j,ispin)
    !
    USE kinds,          ONLY : DP
    USE pdep_io,        ONLY : pdep_read_G_and_distribute
    USE westcom,        ONLY : wbse_init_save_dir,l_bse,l_reduce_io,tau_is_read,tau_all,n_tau,&
                             & n_trunc_bands
    USE pwcom,          ONLY : npwx
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    COMPLEX(DP), INTENT(OUT) :: rhog(npwx)
    INTEGER, INTENT(IN) :: fixed_band_i,fixed_band_j,ispin
    !
    ! Workspace
    !
    INTEGER :: band_i,band_j,iread
    CHARACTER :: my_spin
    CHARACTER(LEN=6) :: my_labeli,my_labelj
    CHARACTER(LEN=256) :: fname
    !
    band_i = MIN(fixed_band_i,fixed_band_j)
    band_j = MAX(fixed_band_i,fixed_band_j)
    !
    IF(l_reduce_io) THEN
       !
       iread = tau_is_read(band_i,band_j,ispin)
       IF(iread > 0) THEN
          rhog(:) = tau_all(:,iread)
          RETURN
       ENDIF
       !
    ENDIF
    !
    WRITE(my_labeli,'(i6.6)') band_i+n_trunc_bands
    WRITE(my_labelj,'(i6.6)') band_j+n_trunc_bands
    WRITE(my_spin,'(i1)') ispin
    !
    IF(l_bse) THEN
       fname = TRIM(wbse_init_save_dir)//'/int_W'//my_labeli//'_'//my_labelj//'_'//my_spin//'.dat'
    ELSE
       fname = TRIM(wbse_init_save_dir)//'/int_v'//my_labeli//'_'//my_labelj//'_'//my_spin//'.dat'
    ENDIF
    CALL pdep_read_G_and_distribute(fname,rhog)
    !
    IF(l_reduce_io) THEN
       !
       n_tau = n_tau+1
       tau_is_read(band_i,band_j,ispin) = n_tau
       tau_all(:,n_tau) = rhog
       !
    ENDIF
    !
  END SUBROUTINE
  !
  SUBROUTINE write_umatrix_and_omatrix(oumat_dim,ispin,umatrix,omatrix)
    !
    USE kinds,          ONLY : DP,i8b
    USE mp_world,       ONLY : world_comm
    USE io_global,      ONLY : stdout
    USE mp,             ONLY : mp_barrier
    USE mp_world,       ONLY : mpime,root
    USE westcom,        ONLY : wbse_init_save_dir
    USE west_io,        ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,HD_ID_DIMENSION
    USE base64_module,  ONLY : islittleendian
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER, INTENT(IN) :: oumat_dim,ispin
    REAL(DP), INTENT(IN) :: omatrix(oumat_dim,oumat_dim)
    COMPLEX(DP), INTENT(IN) :: umatrix(oumat_dim,oumat_dim)
    !
    ! Workspace
    !
    INTEGER :: iun
    CHARACTER(LEN=256) :: fname
    CHARACTER :: my_spin
    INTEGER :: header(HD_LENGTH)
    INTEGER(i8b) :: offset
    !
    WRITE(my_spin,'(i1)') ispin
    fname = TRIM(wbse_init_save_dir)//'/o_and_u.'//my_spin//'.dat'
    !
    WRITE(stdout,'(/,5X,"Writing overlap and rotation matrices to ",A)') TRIM(fname)
    !
    ! Resume all components
    !
    IF(mpime == root) THEN
       !
       header = 0
       header(HD_ID_VERSION) = HD_VERSION
       header(HD_ID_DIMENSION) = oumat_dim
       IF(islittleendian()) THEN
          header(HD_ID_LITTLE_ENDIAN) = 1
       ENDIF
       !
       OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
       offset = 1
       WRITE(iun,POS=offset) header
       offset = offset+SIZEOF(header)
       WRITE(iun,POS=offset) umatrix(1:oumat_dim,1:oumat_dim)
       offset = offset+SIZEOF(umatrix)
       WRITE(iun,POS=offset) omatrix(1:oumat_dim,1:oumat_dim)
       CLOSE(iun)
       !
    ENDIF
    !
    ! BARRIER
    !
    CALL mp_barrier(world_comm)
    !
  END SUBROUTINE
  !
  SUBROUTINE read_umatrix_and_omatrix(oumat_dim,ispin,umatrix,omatrix)
    !
    USE kinds,          ONLY : DP,i8b
    USE io_global,      ONLY : stdout,ionode
    USE mp_world,       ONLY : world_comm
    USE mp,             ONLY : mp_bcast,mp_barrier
    USE mp_global,      ONLY : intra_image_comm
    USE westcom,        ONLY : wbse_init_save_dir
    USE west_io,        ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,HD_ID_DIMENSION
    USE base64_module,  ONLY : islittleendian
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    INTEGER, INTENT(IN) :: oumat_dim,ispin
    REAL(DP), INTENT(OUT) :: omatrix(oumat_dim,oumat_dim)
    COMPLEX(DP), INTENT(OUT) :: umatrix(oumat_dim,oumat_dim)
    !
    ! Workspace
    !
    INTEGER :: ierr,iun
    INTEGER :: oumat_dim_tmp
    CHARACTER(LEN=256) :: fname
    CHARACTER :: my_spin
    INTEGER :: header(HD_LENGTH)
    INTEGER(i8b) :: offset
    REAL(DP), ALLOCATABLE :: omatrix_tmp(:,:)
    COMPLEX(DP), ALLOCATABLE :: umatrix_tmp(:,:)
    !
    ! BARRIER
    !
    CALL mp_barrier(world_comm)
    !
    WRITE(my_spin,'(i1)') ispin
    fname = TRIM(wbse_init_save_dir)//'/o_and_u.'//my_spin//'.dat'
    !
    WRITE(stdout,'(/,5X,"Reading overlap and rotation matrices from ",A)') TRIM(fname)
    !
    IF(ionode) THEN
       !
       OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr /= 0) THEN
          CALL errore('read_umatrix_and_omatrix','Cannot read file: '//TRIM(fname),1)
       ENDIF
       !
       offset = 1
       READ(iun,POS=offset) header
       IF(HD_VERSION /= header(HD_ID_VERSION)) THEN
          CALL errore('read_umatrix_and_omatrix','Unknown file format: '//TRIM(fname),1)
       ENDIF
       IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
          .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) THEN
          CALL errore('read_umatrix_and_omatrix','Endianness mismatch: '//TRIM(fname),1)
       ENDIF
       oumat_dim_tmp = header(HD_ID_DIMENSION)
       !
    ENDIF
    !
    CALL mp_bcast(oumat_dim_tmp,0,intra_image_comm)
    !
    ALLOCATE(umatrix_tmp(oumat_dim_tmp,oumat_dim_tmp))
    ALLOCATE(omatrix_tmp(oumat_dim_tmp,oumat_dim_tmp))
    !
    IF(ionode) THEN
       !
       offset = offset+SIZEOF(header)
       READ(iun,POS=offset) umatrix_tmp(1:oumat_dim_tmp,1:oumat_dim_tmp)
       offset = offset+SIZEOF(umatrix_tmp)
       READ(iun,POS=offset) omatrix_tmp(1:oumat_dim_tmp,1:oumat_dim_tmp)
       CLOSE(iun)
       !
    ENDIF
    !
    CALL mp_bcast(umatrix_tmp,0,intra_image_comm)
    CALL mp_bcast(omatrix_tmp,0,intra_image_comm)
    !
    umatrix(:,:) = 0._DP
    omatrix(:,:) = 0._DP
    umatrix(1:oumat_dim_tmp,1:oumat_dim_tmp) = umatrix_tmp(1:oumat_dim_tmp,1:oumat_dim_tmp)
    omatrix(1:oumat_dim_tmp,1:oumat_dim_tmp) = omatrix_tmp(1:oumat_dim_tmp,1:oumat_dim_tmp)
    !
    DEALLOCATE(umatrix_tmp)
    DEALLOCATE(omatrix_tmp)
    !
  END SUBROUTINE
  !
END MODULE
