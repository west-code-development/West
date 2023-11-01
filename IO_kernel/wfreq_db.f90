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
MODULE wfreq_db
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    ! *****************************
    ! WFREQ WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE wfreq_db_write( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_barrier,mp_sum
      USE mp_global,            ONLY : inter_pool_comm
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_save_dir,qp_bands,n_bands,wfreq_calculation,n_spectralf,logfile,&
                                     & sigma_exx,sigma_vxcl,sigma_vxcnl,sigma_hf,sigma_z,sigma_eqplin,&
                                     & sigma_eqpsec,sigma_sc_eks,sigma_sc_eqplin,sigma_sc_eqpsec,sigma_diff,&
                                     & sigma_freq,sigma_spectralf,l_enable_off_diagonal,pijmap,n_pairs,&
                                     & sigma_vxcl_full,sigma_vxcnl_full,sigma_exx_full,sigma_hf_full,&
                                     & sigma_sc_eks_full,sigma_sc_eqplin_full,sigma_corr_full,occupation
      USE pwcom,                ONLY : et,nspin
      USE io_push,              ONLY : io_push_bar
      USE json_module,          ONLY : json_file
      USE constants,            ONLY : rytoev
      USE types_bz_grid,        ONLY : k_grid
      USE distribution_center,  ONLY : kpt_pool
      !
      IMPLICIT NONE
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: iks,iks_g,ib,ipair
      CHARACTER(LEN=6) :: my_label_k,my_label_b
      CHARACTER(LEN=10) :: label
      !
      TYPE(json_file) :: json
      INTEGER :: iun,i
      REAL(DP),ALLOCATABLE :: eks(:),occ(:,:)
      LOGICAL :: l_generate_plot,l_optics
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('wfreq_db')
      time_spent(1) = get_clock('wfreq_db')
      !
      ! get global copy of occupation
      !
      ALLOCATE(occ(n_bands,k_grid%nps))
      !
      occ(:,:) = 0._DP
      DO iks = 1,kpt_pool%nloc
         iks_g = kpt_pool%l2g(iks)
         DO ib = 1,n_bands
            occ(ib,iks_g) = occupation(qp_bands(ib),iks)
         ENDDO
      ENDDO
      !
      CALL mp_sum(occ,inter_pool_comm)
      !
      IF(nspin == 1) occ(:,:) = 2._DP*occ
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(logfile))
         !
         l_generate_plot = .FALSE.
         l_optics = .FALSE.
         DO i = 1,9
            IF(wfreq_calculation(i:i) == 'P') l_generate_plot = .TRUE.
            IF(wfreq_calculation(i:i) == 'O') l_optics = .TRUE.
         ENDDO
         !
         CALL json%add('output.Q.bandmap',qp_bands)
         !
         IF(ALLOCATED(pijmap)) THEN
            DO ipair = 1,n_pairs
               WRITE(label,'(i10)') ipair
               CALL json%add('output.Q.indexmap('//label//')',(/pijmap(1,ipair),pijmap(2,ipair)/))
            ENDDO
         ELSE
            DO ib = 1,n_bands
               WRITE(label,'(i10)') ib
               CALL json%add('output.Q.indexmap('//label//')',(/ib,ib/))
            ENDDO
         ENDIF
         !
         IF(l_generate_plot) CALL json%add('output.P.freqlist',sigma_freq*rytoev)
         !
         ALLOCATE(eks(n_bands))
         !
         DO iks = 1,k_grid%nps
            !
            DO ib = 1,n_bands
               eks(ib) = et(qp_bands(ib),iks)
            ENDDO
            !
            WRITE(my_label_k,'(i6.6)') iks
            !
            CALL json%add('output.Q.K'//my_label_k//'.eks',eks*rytoev)
            CALL json%add('output.Q.K'//my_label_k//'.z',sigma_z(:,iks))
            CALL json%add('output.Q.K'//my_label_k//'.eqpLin',sigma_eqplin(:,iks)*rytoev)
            CALL json%add('output.Q.K'//my_label_k//'.eqpSec',sigma_eqpsec(:,iks)*rytoev)
            CALL json%add('output.Q.K'//my_label_k//'.sigma_diff',sigma_diff(:,iks)*rytoev)
            CALL json%add('output.Q.K'//my_label_k//'.occupation',occ(:,iks))
            IF(.NOT. l_enable_off_diagonal) THEN
               CALL json%add('output.Q.K'//my_label_k//'.sigmax',sigma_exx(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.vxcl',sigma_vxcl(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.vxcnl',sigma_vxcnl(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.hf',sigma_hf(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eks.re',&
               & REAL(sigma_sc_eks(:,iks)*rytoev,KIND=DP))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eks.im',&
               & AIMAG(sigma_sc_eks(:,iks)*rytoev))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpLin.re',&
               & REAL(sigma_sc_eqplin(:,iks)*rytoev,KIND=DP))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpLin.im',&
               & AIMAG(sigma_sc_eqplin(:,iks)*rytoev))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpSec.re',&
               & REAL(sigma_sc_eqpsec(:,iks)*rytoev,KIND=DP))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpSec.im',&
               & AIMAG(sigma_sc_eqpsec(:,iks)*rytoev))
            ELSE
               CALL json%add('output.Q.K'//my_label_k//'.sigmax',sigma_exx_full(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.vxcl',sigma_vxcl_full(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.vxcnl',sigma_vxcnl_full(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.hf',sigma_hf_full(:,iks)*rytoev)
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eks.re',&
               & REAL(sigma_sc_eks_full(:,iks)*rytoev,KIND=DP))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eks.im',&
               & AIMAG(sigma_sc_eks_full(:,iks)*rytoev))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpLin.re',&
               & REAL(sigma_sc_eqplin_full(:,iks)*rytoev,KIND=DP))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpLin.im',&
               & AIMAG(sigma_sc_eqplin_full(:,iks)*rytoev))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpSec.re',&
               & REAL(sigma_corr_full(:,iks)*rytoev,KIND=DP))
               CALL json%add('output.Q.K'//my_label_k//'.sigmac_eqpSec.im',&
               & AIMAG(sigma_corr_full(:,iks)*rytoev))
            ENDIF
            !
            IF(l_generate_plot) THEN
               DO ib = 1,n_bands
                  WRITE(my_label_b,'(i6.6)') qp_bands(ib)
                  CALL json%add('output.P.K'//my_label_k//'.B'//my_label_b//'.sigmac.re',&
                  & REAL(sigma_spectralf(:,ib,iks),KIND=DP)*rytoev)
                  CALL json%add('output.P.K'//my_label_k//'.B'//my_label_b//'.sigmac.im',&
                  & AIMAG(sigma_spectralf(:,ib,iks))*rytoev)
               ENDDO
            ENDIF
            !
            IF(l_optics) THEN
               CALL json%add('output.O',"optics.json")
            ENDIF
            !
         ENDDO
         !
         DEALLOCATE(eks)
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(logfile))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      DEALLOCATE(occ)
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      time_spent(2) = get_clock('wfreq_db')
      CALL stop_clock('wfreq_db')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"SAVE written in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wfreq_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE qdet_db_write_eri(eri_w,eri_vc,eri_w_full)
    !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_save_dir,logfile,n_pairs
      USE pwcom,                ONLY : nspin
      USE io_push,              ONLY : io_push_bar
      USE json_module,          ONLY : json_file
      USE constants,            ONLY : rytoev
      !
      IMPLICIT NONE
      !
      COMPLEX(DP),INTENT(IN):: eri_w(n_pairs,n_pairs,nspin,nspin)
      REAL(DP),INTENT(IN),OPTIONAL:: eri_vc(n_pairs,n_pairs,nspin,nspin)
      COMPLEX(DP),INTENT(IN),OPTIONAL:: eri_w_full(n_pairs,n_pairs,nspin,nspin)
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: iks,jks
      CHARACTER(LEN=6) :: my_label_ik,my_label_jk,my_label_ipair
      !
      TYPE(json_file) :: json
      INTEGER :: iun,ipair
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('qdet_db')
      time_spent(1) = get_clock('qdet_db')
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(logfile))
         !
         DO iks = 1,nspin
            DO jks = 1,nspin
               !
               WRITE(my_label_ik,'(i6.6)') iks
               WRITE(my_label_jk,'(i6.6)') jks
               !
               DO ipair = 1,n_pairs
                  !
                  WRITE(my_label_ipair,'(i6.6)') ipair
                  !
                  IF(PRESENT(eri_vc)) THEN
                     CALL json%add('qdet.eri_vc.K'//my_label_ik//'.K'//my_label_jk//'.pair'//&
                     & my_label_ipair,eri_vc(:,ipair,jks,iks)*rytoev)
                  ENDIF
                  !
                  IF(PRESENT(eri_w_full)) THEN
                     CALL json%add('qdet.eri_w_full.K'//my_label_ik//'.K'//my_label_jk//'.pair'//&
                     & my_label_ipair,REAL(eri_w_full(:,ipair,jks,iks),KIND=DP)*rytoev)
                  ENDIF
                  !
                  CALL json%add('qdet.eri_w.K'//my_label_ik//'.K'//my_label_jk//'.pair'//&
                  & my_label_ipair,REAL(eri_w(:,ipair,jks,iks),KIND=DP)*rytoev)
                  !
               ENDDO
               !
            ENDDO
         ENDDO
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(logfile))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      time_spent(2) = get_clock('qdet_db')
      CALL stop_clock('qdet_db')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"SAVE written in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x "In location : ",a)') TRIM(wfreq_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE qdet_db_write_h1e(h1e)
    !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_save_dir,logfile,n_pairs
      USE pwcom,                ONLY : nspin
      USE io_push,              ONLY : io_push_bar
      USE json_module,          ONLY : json_file
      USE constants,            ONLY : rytoev
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN):: h1e(n_pairs,nspin)
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: iks
      CHARACTER(LEN=6) :: my_label_ik
      !
      TYPE(json_file) :: json
      INTEGER :: iun
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('qdet_db')
      time_spent(1) = get_clock('qdet_db')
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(logfile))
         !
         DO iks = 1,nspin
            WRITE(my_label_ik,'(i6.6)') iks
            CALL json%add('qdet.h1e.K'//my_label_ik,h1e(:,iks)*rytoev)
         ENDDO
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(logfile))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2) = get_clock('qdet_db')
      CALL stop_clock('qdet_db')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"SAVE written in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wfreq_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE qdet_db_write_overlap(overlap)
    !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : wfreq_save_dir,logfile,n_bands
      USE io_push,              ONLY : io_push_bar
      USE json_module,          ONLY : json_file
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN):: overlap(n_bands,n_bands)
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      TYPE(json_file) :: json
      INTEGER :: iun
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('qdet_db')
      time_spent(1) = get_clock('qdet_db')
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(logfile))
         !
         CALL json%add('qdet.overlap_ab',RESHAPE(overlap,[n_bands*n_bands]))
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(logfile))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2) = get_clock('qdet_db')
      CALL stop_clock('qdet_db')
      !
      WRITE(stdout,*)
      CALL io_push_bar()
      WRITE(stdout,'(5x,"SAVE written in ",a)') TRIM(human_readable_time(time_spent(2)-time_spent(1)))
      WRITE(stdout,'(5x,"In location : ",a)') TRIM(wfreq_save_dir)
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
END MODULE
