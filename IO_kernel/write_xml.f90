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
MODULE write_xml
  !----------------------------------------------------------------------------
  !
  USE iotk_module
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_xml_dump( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier,mp_max
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout 
      USE io_push,              ONLY : io_push_bar
      USE westcom,              ONLY : wfreq_save_dir, wstat_save_dir, l_is_wstat_converged, ev
      USE westcom,              ONLY : qe_prefix, west_prefix, outdir 
      USE westcom,              ONLY : wstat_calculation, n_pdep_eigen, n_pdep_times, n_pdep_maxiter, n_dfpt_maxiter, &
                                     & n_pdep_read_from_file, trev_pdep, trev_pdep_rel, tr2_dfpt, l_minimize_exx_if_active, & 
                                     & l_kinetic_only, l_use_ecutrho
      USE control_flags,        ONLY : gamma_only
      USE gvecw,                ONLY : ecutwfc
      USE cell_base,            ONLY : omega
      USE pwcom,                ONLY : npw,nbnd,nkstot,nspin,nelec,nelup,neldw,lspinorb,domag,lsda
      USE noncollin_module,     ONLY : noncolin,npol
      USE fft_base,             ONLY : dfftp,dffts
      USE west_version,         ONLY : west_version_number
      USE mp_global,            ONLY : nimage,npool,nbgrp,nproc_image,nproc_pool,nproc_bgrp 
      USE mp_world,             ONLY : nproc 
      !
      IMPLICIT NONE
      !
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2), tot_walltime
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: iunout, ierr, nth 
      CHARACTER(LEN=9)  :: cdate, ctime
      CHARACTER(iotk_attlenx) :: attr
      !
#if defined(__OPENMP)
      INTEGER, EXTERNAL :: omp_get_max_threads
#endif
      !
      CALL date_and_tim( cdate, ctime )
#if defined(__OPENMP)
      nth = omp_get_max_threads()
#else
      nth = 0
#endif
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('wstat_xml')
      time_spent(1)=get_clock('wstat_xml')
      !
      tot_walltime = get_clock('WSTAT')
      CALL mp_max( tot_walltime, world_comm ) 
      !
      ! 1)  CREATE THE INPUT FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         !
         CALL iotk_open_write( iunout, FILE = TRIM( wstat_save_dir ) // '/' // TRIM("wstat.xml") , &
         & ROOT="wstat",BINARY=.FALSE.,SKIP_HEAD=.TRUE.,IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      CALL errore( 'wstat', 'cannot open wstat.xml for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         ! INPUT
         !
         CALL iotk_write_begin( iunout, "input" )
         !
         CALL iotk_write_begin( iunout, "input_west" )
         CALL iotk_write_dat( iunout, "qe_prefix"        , qe_prefix )
         CALL iotk_write_dat( iunout, "west_prefix"        , west_prefix )
         CALL iotk_write_dat( iunout, "outdir"        , outdir )
         CALL iotk_write_end( iunout, "input_west" )
         !
         CALL iotk_write_begin( iunout, "wstat_control" )
         CALL iotk_write_dat( iunout, "wstat_calculation"        , wstat_calculation)
         CALL iotk_write_dat( iunout, "n_pdep_eigen"             , n_pdep_eigen)
         CALL iotk_write_dat( iunout, "n_pdep_times"             , n_pdep_times)
         CALL iotk_write_dat( iunout, "n_pdep_maxiter"           , n_pdep_maxiter)
         CALL iotk_write_dat( iunout, "n_dfpt_maxiter"           , n_dfpt_maxiter)
         CALL iotk_write_dat( iunout, "n_pdep_read_from_file"    , n_pdep_read_from_file)
         CALL iotk_write_dat( iunout, "trev_pdep"                , trev_pdep)
         CALL iotk_write_dat( iunout, "trev_pdep_rel"            , trev_pdep_rel)
         CALL iotk_write_dat( iunout, "tr2_dfpt"                 , tr2_dfpt)
         CALL iotk_write_dat( iunout, "l_minimize_exx_if_active" , l_minimize_exx_if_active)
         CALL iotk_write_dat( iunout, "l_kinetic_only"           , l_kinetic_only)
         CALL iotk_write_dat( iunout, "l_use_ecutrho"            , l_use_ecutrho)
         CALL iotk_write_end( iunout, "wstat_control" )
         !
         CALL iotk_write_end( iunout, "input"  )
         !
         ! OUTPUT
         ! 
         CALL iotk_write_begin( iunout, "output" )
         !
         CALL iotk_write_begin( iunout, "pdepEigenval" )
         IF( ALLOCATED(ev) ) THEN 
           CALL iotk_write_attr(attr,"coll","list",FIRST=.TRUE.) 
           CALL iotk_write_dat( iunout, "ev" , ev(1:n_pdep_eigen), ATTR=attr)
         ENDIF
         CALL iotk_write_dat( iunout, "l_is_wstat_converged", .true. )
         CALL iotk_write_end( iunout, "pdepEigenval" )
         !
         CALL iotk_write_end( iunout, "output"  )
         !
         ! EXECUTION
         !
         CALL iotk_write_begin( iunout, "execution" )
         !
         CALL iotk_write_begin( iunout, "general" )
         CALL iotk_write_dat( iunout, "version"                  , west_version_number )
         CALL iotk_write_dat( iunout, "date"                     , cdate )
         CALL iotk_write_dat( iunout, "start_time"                     , ctime )
         CALL iotk_write_end( iunout, "general"  )
         !
         CALL iotk_write_begin( iunout, "parallel" )
         CALL iotk_write_dat( iunout, "nranks"                   , nproc )
         CALL iotk_write_dat( iunout, "nimage"                   , nimage )
         CALL iotk_write_dat( iunout, "npool"                    , npool )
         CALL iotk_write_dat( iunout, "nbgrp"                    , nbgrp )
         CALL iotk_write_dat( iunout, "nz"                       , nproc_bgrp )
         CALL iotk_write_dat( iunout, "nthreads"                 , nth )
         CALL iotk_write_end( iunout, "parallel"  )
         !
         CALL iotk_write_begin( iunout, "status" )
         CALL iotk_write_dat( iunout, "completed", .true. )
         CALL iotk_write_end( iunout, "status"  )
         !
         CALL iotk_write_begin( iunout, "timings" )
         CALL iotk_write_dat( iunout, "tot_time", tot_walltime )
         CALL iotk_write_end( iunout, "timings"  )
         !
         CALL iotk_write_end( iunout, "execution"  )
         !
         ! SYSINFO
         !
         CALL iotk_write_begin( iunout, "sysinfo" )
         !
         CALL iotk_write_dat( iunout, "gamma_only"               , gamma_only )
         CALL iotk_write_dat( iunout, "ecutwfc"                  , ecutwfc )
         CALL iotk_write_dat( iunout, "omega"                    , omega )
         CALL iotk_write_dat( iunout, "nbnd"                     , nbnd )
         CALL iotk_write_dat( iunout, "nkstot"                   , nkstot )
         CALL iotk_write_dat( iunout, "nspin"                    , nspin )
         CALL iotk_write_dat( iunout, "nelec"                    , nelec )
         CALL iotk_write_dat( iunout, "npol"                     , npol )
         CALL iotk_write_dat( iunout, "lsda"                     , lsda )
         CALL iotk_write_dat( iunout, "noncolin"                 , noncolin )
         CALL iotk_write_dat( iunout, "lspinorb"                 , lspinorb )
         CALL iotk_write_dat( iunout, "domag"                    , domag )
         CALL iotk_write_attr(attr,"coll","list",FIRST=.TRUE.) 
         CALL iotk_write_dat( iunout, "fft_grid_s"           , (/ dffts%nr1, dffts%nr2, dffts%nr3 /), ATTR=attr )
         CALL iotk_write_attr(attr,"coll","list",FIRST=.TRUE.) 
         CALL iotk_write_dat( iunout, "fft_grid_p"           , (/ dfftp%nr1, dfftp%nr2, dfftp%nr3 /), ATTR=attr )
         !
         CALL iotk_write_end( iunout, "sysinfo"  )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2)=get_clock('wstat_xml')
      CALL stop_clock('wstat_xml')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'wstat.xml file written in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wstat_save_dir )  
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wfreq_xml_dump( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier,mp_max
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout 
      USE io_push,              ONLY : io_push_bar
      USE westcom,              ONLY : l_is_wstat_converged, ev, iks_l2g
      USE westcom,              ONLY : wfreq_save_dir, qe_prefix, west_prefix, outdir 
      USE westcom,              ONLY : wstat_calculation, n_pdep_eigen, n_pdep_times, n_pdep_maxiter, n_dfpt_maxiter, &
                                     & n_pdep_read_from_file, trev_pdep, trev_pdep_rel, tr2_dfpt, l_minimize_exx_if_active, & 
                                     & l_kinetic_only, l_use_ecutrho
      USE westcom,              ONLY : wfreq_calculation, n_pdep_eigen_to_use, qp_bandrange, macropol_calculation, n_lanczos, &
                                     & n_imfreq, n_refreq, ecut_imfreq, ecut_refreq, wfreq_eta, n_secant_maxiter, trev_secant, &
                                     & l_enable_lanczos, l_enable_gwetot, div_kind_hf, o_restart_time, ecut_spectralf, &
                                     & n_spectralf
      USE westcom,              ONLY : sigma_exx,sigma_vxcl,sigma_vxcnl,sigma_hf,sigma_z,sigma_eqplin,sigma_eqpsec,sigma_sc_eks,&
                                     & sigma_sc_eqplin,sigma_sc_eqpsec,sigma_diff,sigma_freq,sigma_spectralf
      USE control_flags,        ONLY : gamma_only
      USE gvecw,                ONLY : ecutwfc
      USE cell_base,            ONLY : omega
      USE pwcom,                ONLY : npw,nbnd,nkstot,nspin,nelec,nelup,neldw,lspinorb,domag,lsda,nks,et
      USE noncollin_module,     ONLY : noncolin,npol
      USE fft_base,             ONLY : dfftp,dffts
      USE west_version,         ONLY : west_version_number
      USE mp_global,            ONLY : nimage,npool,nbgrp,nproc_image,nproc_pool,nproc_bgrp 
      USE mp_world,             ONLY : nproc 
      USE constants,            ONLY : rytoev
      !
      IMPLICIT NONE
      !
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2), tot_walltime
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: iunout, ierr, nth, iks, ibnd, i
      CHARACTER(LEN=9)  :: cdate, ctime
      CHARACTER(iotk_attlenx) :: attr
      CHARACTER(LEN=5) :: myglobk,myglobb
      LOGICAL :: l_generate_plot
      !
#if defined(__OPENMP)
      INTEGER, EXTERNAL :: omp_get_max_threads
#endif
      !
      CALL date_and_tim( cdate, ctime )
#if defined(__OPENMP)
      nth = omp_get_max_threads()
#else
      nth = 0
#endif
      !
      ! MPI BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! TIMING
      !
      CALL start_clock('wfreq_xml')
      time_spent(1)=get_clock('wfreq_xml')
      !
      tot_walltime = get_clock('WFREQ')
      CALL mp_max( tot_walltime, world_comm ) 
      !
      ! 1)  CREATE THE INPUT FILE
      !
      IF ( mpime == root ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_free_unit( iunout, ierr )
         !
         CALL iotk_open_write( iunout, FILE = TRIM( wfreq_save_dir ) // '/' // TRIM("wfreq.xml") , &
         & ROOT="wfreq",BINARY=.FALSE.,SKIP_HEAD=.TRUE.,IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, root, world_comm )
      CALL errore( 'wfreq', 'cannot open wfreq.xml for writing', ierr )
      !
      IF ( mpime == root ) THEN  
         !
         ! INPUT
         !
         CALL iotk_write_begin( iunout, "input" )
         !
         CALL iotk_write_begin( iunout, "input_west" )
         CALL iotk_write_dat( iunout, "qe_prefix"        , qe_prefix )
         CALL iotk_write_dat( iunout, "west_prefix"        , west_prefix )
         CALL iotk_write_dat( iunout, "outdir"        , outdir )
         CALL iotk_write_end( iunout, "input_west" )
         !
         CALL iotk_write_begin( iunout, "wstat_control" )
         CALL iotk_write_dat( iunout, "wstat_calculation"        , wstat_calculation)
         CALL iotk_write_dat( iunout, "n_pdep_eigen"             , n_pdep_eigen)
         CALL iotk_write_dat( iunout, "n_pdep_times"             , n_pdep_times)
         CALL iotk_write_dat( iunout, "n_pdep_maxiter"           , n_pdep_maxiter)
         CALL iotk_write_dat( iunout, "n_dfpt_maxiter"           , n_dfpt_maxiter)
         CALL iotk_write_dat( iunout, "n_pdep_read_from_file"    , n_pdep_read_from_file)
         CALL iotk_write_dat( iunout, "trev_pdep"                , trev_pdep)
         CALL iotk_write_dat( iunout, "trev_pdep_rel"            , trev_pdep_rel)
         CALL iotk_write_dat( iunout, "tr2_dfpt"                 , tr2_dfpt)
         CALL iotk_write_dat( iunout, "l_minimize_exx_if_active" , l_minimize_exx_if_active)
         CALL iotk_write_dat( iunout, "l_kinetic_only"           , l_kinetic_only)
         CALL iotk_write_dat( iunout, "l_use_ecutrho"            , l_use_ecutrho)
         CALL iotk_write_end( iunout, "wstat_control" )
         !
         CALL iotk_write_begin( iunout, "wfreq_control" )
         CALL iotk_write_dat( iunout, "wfreq_calculation"        , wfreq_calculation)
         CALL iotk_write_dat( iunout, "n_pdep_eigen_to_use"      , n_pdep_eigen_to_use)
         CALL iotk_write_dat( iunout, "qp_bandrange"             , qp_bandrange)
         CALL iotk_write_dat( iunout, "macropol_calculation"     , macropol_calculation)
         CALL iotk_write_dat( iunout, "n_lanczos"                , n_lanczos)
         CALL iotk_write_dat( iunout, "n_imfreq"                 , n_imfreq)
         CALL iotk_write_dat( iunout, "n_refreq"                 , n_refreq)
         CALL iotk_write_dat( iunout, "ecut_imfreq"              , ecut_imfreq)
         CALL iotk_write_dat( iunout, "ecut_refreq"              , ecut_refreq)
         CALL iotk_write_dat( iunout, "wfreq_eta"                , wfreq_eta)
         CALL iotk_write_dat( iunout, "n_secant_maxiter"         , n_secant_maxiter)
         CALL iotk_write_dat( iunout, "trev_secant"              , trev_secant)
         CALL iotk_write_dat( iunout, "l_enable_lanczos"         , l_enable_lanczos)
         CALL iotk_write_dat( iunout, "l_enable_gwetot"          , l_enable_gwetot)
         CALL iotk_write_dat( iunout, "div_kind_hf"              , div_kind_hf)
         CALL iotk_write_dat( iunout, "o_restart_time"           , o_restart_time)
         CALL iotk_write_dat( iunout, "ecut_spectralf"           , ecut_spectralf)
         CALL iotk_write_dat( iunout, "n_spectralf"              , n_spectralf)
         CALL iotk_write_end( iunout, "wfreq_control" )
         !
         CALL iotk_write_end( iunout, "input"  )
         !
         ! OUTPUT
         ! 
         CALL iotk_write_begin( iunout, "output" )
         !
         CALL iotk_write_attr(attr,"coll","kblist",FIRST=.TRUE.) 
         CALL iotk_write_begin( iunout, "gw", ATTR=attr )
         !
         l_generate_plot = .FALSE.
         DO i = 1,8
            IF( wfreq_calculation(i:i) == 'P' ) l_generate_plot = .TRUE.
         ENDDO
         !
         DO iks = 1, nks
            !
            !WRITE(myglobk,'(i5.5)') iks_l2g(iks)
            DO ibnd = qp_bandrange(1), qp_bandrange(2)
               !
               CALL iotk_write_attr(attr,"k",iks,FIRST=.TRUE.)
               CALL iotk_write_attr(attr,"b",ibnd)
               CALL iotk_write_begin( iunout, "state", ATTR=attr )
               !
               CALL iotk_write_dat( iunout, "sigmax"          , sigma_exx       (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "vxcl"            , sigma_vxcl      (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "vxcnl"           , sigma_vxcnl     (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "hf"              , sigma_hf        (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "z"               , sigma_z         (ibnd,iks) )
               CALL iotk_write_dat( iunout, "eks"             , et              (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "eqpLin"          , sigma_eqplin    (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "eqpSec"          , sigma_eqpsec    (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "sigmac_eks"      , sigma_sc_eks    (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "sigmac_eqpLin"   , sigma_sc_eqplin (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "sigmac_eqpSec"   , sigma_sc_eqpsec (ibnd,iks) * rytoev )
               CALL iotk_write_dat( iunout, "ediff"           , sigma_diff      (ibnd,iks) * rytoev )
               !
               IF( l_generate_plot ) THEN 
                  CALL iotk_write_attr(attr,"coll","list",FIRST=.TRUE.)
                  CALL iotk_write_dat( iunout, "sigmac" , sigma_spectralf(1:n_spectralf,ibnd,iks) * rytoev, ATTR=attr)
               ENDIF
               !
               CALL iotk_write_end( iunout, "state")
               !
            ENDDO
            !
         ENDDO
         !
         CALL iotk_write_end( iunout, "gw" )
         !
         IF( l_generate_plot ) THEN 
            !
            CALL iotk_write_attr(attr,"coll","list",FIRST=.TRUE.) 
            CALL iotk_write_dat( iunout, "sigmac_energylist", sigma_freq(1:n_spectralf) * rytoev, ATTR=attr )
            ! 
         ENDIF  
         !
         CALL iotk_write_end( iunout, "output"  )
         !
         ! EXECUTION
         !
         CALL iotk_write_begin( iunout, "execution" )
         !
         CALL iotk_write_begin( iunout, "general" )
         CALL iotk_write_dat( iunout, "version"                  , west_version_number )
         CALL iotk_write_dat( iunout, "date"                     , cdate )
         CALL iotk_write_dat( iunout, "start_time"                     , ctime )
         CALL iotk_write_end( iunout, "general"  )
         !
         CALL iotk_write_begin( iunout, "parallel" )
         CALL iotk_write_dat( iunout, "nranks"                   , nproc )
         CALL iotk_write_dat( iunout, "nimage"                   , nimage )
         CALL iotk_write_dat( iunout, "npool"                    , npool )
         CALL iotk_write_dat( iunout, "nbgrp"                    , nbgrp )
         CALL iotk_write_dat( iunout, "nz"                       , nproc_bgrp )
         CALL iotk_write_dat( iunout, "nthreads"                 , nth )
         CALL iotk_write_end( iunout, "parallel"  )
         !
         CALL iotk_write_begin( iunout, "status" )
         CALL iotk_write_dat( iunout, "completed", .true. )
         CALL iotk_write_end( iunout, "status"  )
         !
         CALL iotk_write_begin( iunout, "timings" )
         CALL iotk_write_dat( iunout, "tot_time", tot_walltime )
         CALL iotk_write_end( iunout, "timings"  )
         !
         CALL iotk_write_end( iunout, "execution"  )
         !
         ! SYSINFO
         !
         CALL iotk_write_begin( iunout, "sysinfo" )
         !
         CALL iotk_write_dat( iunout, "gamma_only"               , gamma_only )
         CALL iotk_write_dat( iunout, "ecutwfc"                  , ecutwfc )
         CALL iotk_write_dat( iunout, "omega"                    , omega )
         CALL iotk_write_dat( iunout, "nbnd"                     , nbnd )
         CALL iotk_write_dat( iunout, "nkstot"                   , nkstot )
         CALL iotk_write_dat( iunout, "nspin"                    , nspin )
         CALL iotk_write_dat( iunout, "nelec"                    , nelec )
         CALL iotk_write_dat( iunout, "npol"                     , npol )
         CALL iotk_write_dat( iunout, "lsda"                     , lsda )
         CALL iotk_write_dat( iunout, "noncolin"                 , noncolin )
         CALL iotk_write_dat( iunout, "lspinorb"                 , lspinorb )
         CALL iotk_write_dat( iunout, "domag"                    , domag )
         CALL iotk_write_attr(attr,"coll","list",FIRST=.TRUE.) 
         CALL iotk_write_dat( iunout, "fft_grid_s"           , (/ dffts%nr1, dffts%nr2, dffts%nr3 /), ATTR=attr )
         CALL iotk_write_attr(attr,"coll","list",FIRST=.TRUE.) 
         CALL iotk_write_dat( iunout, "fft_grid_p"           , (/ dfftp%nr1, dfftp%nr2, dfftp%nr3 /), ATTR=attr )
         !
         CALL iotk_write_end( iunout, "sysinfo"  )
         !
         ! ... close XML descriptor
         !
         CALL iotk_close_write( iunout )
         !
      END IF
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2)=get_clock('wfreq_xml')
      CALL stop_clock('wfreq_xml')
      !
      WRITE(stdout,'(  5x," ")')
      CALL io_push_bar()
      WRITE(stdout, "(5x, 'wfreq.xml file written in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
      WRITE(stdout, "(5x, 'In location : ',a)") TRIM( wfreq_save_dir )  
      CALL io_push_bar()
      !
    END SUBROUTINE
    !
    !
END MODULE
