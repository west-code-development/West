!
! Copyright (C) 2015-2017 M. Govoni 
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
  !
  CONTAINS
    !
    !
    ! *****************************
    ! WFREQ WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE wfreq_db_write( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp_global,            ONLY : my_image_id
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : wfreq_save_dir,iks_l2g,qp_bandrange,wfreq_calculation,n_spectralf, &
                                     & sigma_exx,sigma_vxcl,sigma_vxcnl,sigma_hf,sigma_z,sigma_eqplin,sigma_eqpsec,sigma_sc_eks,&
                                     & sigma_sc_eqplin,sigma_sc_eqpsec,sigma_diff,sigma_freq,sigma_spectralf
      USE pwcom,                ONLY : npw,nbnd,nkstot,nspin,nelec,nelup,neldw,lspinorb,domag,lsda,nks,et
      USE io_push,              ONLY : io_push_bar
      USE json_module,          ONLY : json_file 
      USE constants,            ONLY : rytoev
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=512)    :: fname
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: iunout,global_j,local_j
      INTEGER :: ierr, iks, ib
      CHARACTER(LEN=6) :: my_label_k, my_label_b
      !
      TYPE(json_file) :: json 
      INTEGER :: iunit, i
      INTEGER,ALLOCATABLE :: ilist(:) 
      LOGICAL :: l_generate_plot, l_optics
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('wfreq_db')
      time_spent(1)=get_clock('wfreq_db')
      !
      IF ( mpime == root ) THEN
         !
         CALL json%initialize()
         !
         CALL add_intput_parameters_to_json_file( 3,(/1,2,3/) , json )
         !
         l_generate_plot = .FALSE.
         l_optics = .FALSE.
         DO i = 1,8
            IF( wfreq_calculation(i:i) == 'P' ) l_generate_plot = .TRUE.
            IF( wfreq_calculation(i:i) == 'O' ) l_optics = .TRUE.
         ENDDO
         !
         ALLOCATE(ilist(qp_bandrange(1):qp_bandrange(2)))
         DO ib = qp_bandrange(1),qp_bandrange(2)
            ilist(ib) = ib
         ENDDO
         CALL json%add('output.Q.bandmap',ilist(qp_bandrange(1):qp_bandrange(2)))
         DEALLOCATE(ilist)
         IF( l_generate_plot ) CALL json%add('output.P.freqlist',sigma_freq(1:n_spectralf)*rytoev)
         !
         DO iks = 1, nks 
            !
            WRITE(my_label_k,'(i6.6)') iks_l2g(iks)
            !
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigmax', sigma_exx(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.vxcl', sigma_vxcl(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.vxcnl', sigma_vxcnl(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.hf', sigma_hf(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.z', sigma_z(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.eks', et(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.eqpLin', sigma_eqplin(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.eqpSec', sigma_eqpsec(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigmac_eks.re', &
            & DBLE(sigma_sc_eks(qp_bandrange(1):qp_bandrange(2),iks)*rytoev))
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigmac_eks.im', &
            & AIMAG(sigma_sc_eks(qp_bandrange(1):qp_bandrange(2),iks)*rytoev))
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigmac_eqpLin.re', &
            & DBLE(sigma_sc_eqplin(qp_bandrange(1):qp_bandrange(2),iks)*rytoev))
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigmac_eqpLin.im', &
            & AIMAG(sigma_sc_eqplin(qp_bandrange(1):qp_bandrange(2),iks)*rytoev))
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigmac_eqpSec.re', &
            & DBLE(sigma_sc_eqpsec(qp_bandrange(1):qp_bandrange(2),iks)*rytoev))
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigmac_eqpSec.im', &
            & AIMAG(sigma_sc_eqpsec(qp_bandrange(1):qp_bandrange(2),iks)*rytoev))
            CALL json%add('output.Q.K'//TRIM(my_label_k)//'.sigma_diff', sigma_diff(qp_bandrange(1):qp_bandrange(2),iks)*rytoev)
            !
            IF( l_generate_plot ) THEN 
               DO ib = qp_bandrange(1), qp_bandrange(2)
                  WRITE(my_label_b,'(i6.6)') ib 
                  CALL json%add('output.P.K'//TRIM(my_label_k)//'.B'//TRIM(my_label_b)//'.sigmac.re',&
                  &DBLE(sigma_spectralf(1:n_spectralf,ib,iks))*rytoev)
                  CALL json%add('output.P.K'//TRIM(my_label_k)//'.B'//TRIM(my_label_b)//'.sigmac.im',&
                  &AIMAG(sigma_spectralf(1:n_spectralf,ib,iks))*rytoev)
               ENDDO
            ENDIF
            !
            IF( l_optics) THEN
               CALL json%add('output.O',"optics.json") 
            ENDIF 
            !
         ENDDO
         !
         OPEN( NEWUNIT=iunit, FILE=TRIM( wfreq_save_dir ) // '/' // TRIM('wfreq-save.json') )
         CALL json%print_file( iunit )
         CLOSE( iunit )
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
      time_spent(2)=get_clock('wfreq_db')
      CALL stop_clock('wfreq_db')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'Database written in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wfreq_save_dir )  
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
END MODULE
