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
!----------------------------------------------------------------------------
SUBROUTINE do_eigenpot2 ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE io_push,               ONLY : io_push_title
  USE westcom,               ONLY : westpp_sign,westpp_range,westpp_save_dir,fftdriver,&
                                  & npwq,npwqx,dvg,ngq,igq_q,westpp_n_pdep_eigen_to_use
  USE fft_base,              ONLY : dffts
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k
  USE distribution_center,   ONLY : pert
  USE class_idistribute,     ONLY : idistribute
  USE control_flags,         ONLY : gamma_only
  USE pdep_db,               ONLY : pdep_db_read
  USE types_bz_grid,         ONLY : q_grid
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : psic=>psic_d
  USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu
#else
  USE wavefunctions,         ONLY : psic
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ir, local_j, global_j, iq, dffts_nnr
  REAL(DP),ALLOCATABLE :: auxr(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  CHARACTER(LEN=6) :: labeli
  CHARACTER(LEN=5) :: labelq
  LOGICAL :: l_print_pdep_read
  !
  pert = idistribute()
  CALL pert%init(westpp_n_pdep_eigen_to_use,'i','npdep',.TRUE.)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  ALLOCATE(auxr(dffts%nnr))
  !$acc enter data create(auxr)
  !
  dffts_nnr = dffts%nnr
  !
  CALL io_push_title('(E)igenpotentials')
  !
  CALL start_bar_type( barra, 'westpp', pert%nloc*q_grid%np)
  !
  DO iq = 1, q_grid%np
     !
     IF (gamma_only) THEN
        CALL pdep_db_read(westpp_n_pdep_eigen_to_use)
     ELSE
        IF (iq==1) THEN
           l_print_pdep_read = .TRUE.
        ELSE
           l_print_pdep_read = .FALSE.
        ENDIF
        CALL pdep_db_read(westpp_n_pdep_eigen_to_use,iq,l_print_pdep_read)
        npwq = ngq(iq)
     ENDIF
     !
     !$acc enter data copyin(dvg)
     !
     DO local_j=1,pert%nloc
        !
        ! local -> global
        !
        global_j = pert%l2g(local_j)
        IF( global_j < westpp_range(1) .OR. global_j > westpp_range(2) ) CYCLE
        !
        IF( gamma_only ) THEN
           !$acc host_data use_device(dvg)
           CALL single_invfft_gamma(dffts,npwq,npwqx,dvg(:,local_j),psic,TRIM(fftdriver))
           !$acc end host_data
        ELSE
           !$acc host_data use_device(dvg)
           CALL single_invfft_k(dffts,npwq,npwqx,dvg(:,local_j),psic,'Wave',igq_q(:,iq))
           !$acc end host_data
        ENDIF
        IF( westpp_sign ) THEN
           !$acc parallel loop present(auxr)
           DO ir = 1, dffts_nnr
              auxr(ir) = REAL( psic(ir), KIND=DP) *  ABS( REAL( psic(ir), KIND=DP) )
           ENDDO
           !$acc end parallel loop
        ELSE
           !$acc parallel loop present(auxr)
           DO ir = 1, dffts_nnr
              auxr(ir) = REAL( psic(ir), KIND=DP) *  REAL( psic(ir), KIND=DP)
           ENDDO
           !$acc end parallel loop
        ENDIF
        !
        WRITE(labeli,'(i6.6)') global_j
        WRITE(labelq,'(i5.5)') iq
        fname = TRIM( westpp_save_dir ) // '/eigQ'//labelq//'I'//labeli
        !$acc update host(auxr)
        CALL dump_r( auxr, fname)
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO
     !
     !$acc exit data delete(dvg)
     !
  ENDDO
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  !$acc exit data delete(auxr)
  DEALLOCATE( auxr )
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
