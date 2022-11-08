!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE do_spectrum()
  !---------------------------------------------------------------------
  !
  ! This routine is modified from the TDDFT code.
  ! Pls do not distribute it to avoid the copyright problem.
  ! TODO: write this code in Python to make it become ours.
  !
  USE kinds,               ONLY : DP
  USE constants,           ONLY : pi,rytoev,evtonm,rytonm
  USE io_global,           ONLY : stdout,ionode
  USE environment,         ONLY : environment_start,environment_end
  USE mp_global,           ONLY : mp_startup,mp_global_end
  USE mp_world,            ONLY : world_comm
  USE mp,                  ONLY : mp_bcast,mp_barrier
  USE westcom,             ONLY : qe_prefix,itermax,itermax0,extrapolation,start,end,increment,&
                                & epsil,ipol,sym_op,verbosity,units,spin_channel
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  ! Constants
  !
  !  integer,  parameter :: dp=kind(0.d0)
  !  real(dp), parameter :: pi=3.14159265d0
  !  real(dp), parameter :: ry=13.6056981d0
  !  real(dp), parameter :: ev2nm=1239.84172
  !  real(dp), parameter :: ry2nm=91.1266519
  !
  ! User controlled variables
  !
  REAL(DP) :: omega(3)
  CHARACTER(LEN=256):: filename
  CHARACTER(LEN=256):: filename_plot
  !
  ! General use variables & counters
  !
  INTEGER :: n_ipol, i,j, info, ip, ip2, counter
  INTEGER :: iun1, iun2
  REAL(DP) :: norm0(3), average(3), av_amplitude(3), alpha_temp(3), scale, wl, degspin, f_sum
  COMPLEX(DP) :: omeg_c
  REAL(DP), ALLOCATABLE :: beta_store(:,:), gamma_store(:,:)
  COMPLEX(DP), ALLOCATABLE :: zeta_store(:,:,:)
  COMPLEX(DP) :: green(3,3) ! susceptibility chi
  COMPLEX(DP), ALLOCATABLE :: a(:), b(:), c(:), r(:,:)
  LOGICAL :: skip
  !
  ! For perceived color analysis
  !
  REAL(DP), PARAMETER :: vis_start    = 0.116829041_DP
  REAL(DP), PARAMETER :: vis_start_wl = 780._DP
  REAL(DP), PARAMETER :: vis_end      = 0.239806979_DP
  REAL(DP), PARAMETER :: vis_end_wl   = 380._DP
  LOGICAL  :: do_perceived
  REAL(DP) :: perceived_red   = 0._DP
  REAL(DP) :: perceived_green = 0._DP
  REAL(DP) :: perceived_blue  = 0._DP
  REAL(DP) :: perceived_renorm
  INTEGER  :: perceived_itermax,perceived_iter
  REAL(DP), ALLOCATABLE :: perceived_intensity(:), perceived_evaluated(:)
  !
  ! Subroutines etc.
  !
  COMPLEX(DP), EXTERNAL :: zdotc
  !
  ! User controlled variable initialisation
  !
  CALL mp_startup()
  !
  IF(ionode) THEN
     !
     IF(itermax0 < 151 .AND. TRIM(extrapolation) /= 'no') THEN
        WRITE(stdout,'(5x,"Itermax0 is less than 150, no extrapolation scheme can be used")')
        extrapolation = 'no'
     ENDIF
     !
     IF(ipol < 4) THEN
        n_ipol = 1
     ELSE
        n_ipol = 3
        ipol = 1
     ENDIF
     !
     ! Polarization symmetry
     !
     IF(.NOT. sym_op == 0) THEN
        !
        IF(sym_op == 1) THEN
           WRITE(stdout,'(5x,"All polarization axes will be considered to be equal")')
           n_ipol = 3
           ipol = 1
        ELSE
           CALL errore('do_spectrum','Unsupported symmetry operation',1)
        ENDIF
        !
     ENDIF
     !
     ! Terminator scheme
     !
     IF(TRIM(extrapolation) == 'no') THEN
        itermax = itermax0
     ENDIF
     !
     ! Check the units (Ry, eV, nm)
     !
     IF(units < 0 .OR. units > 2) CALL errore('do_spectrum','Unsupported unit system',1)
     !
     IF(units /= 0 .AND. verbosity > 4) THEN
        verbosity = 4
     ENDIF
     !
     ! Initialisation of coefficients
     !
     ALLOCATE(beta_store(n_ipol,itermax))
     ALLOCATE(gamma_store(n_ipol,itermax))
     ALLOCATE(zeta_store(n_ipol,n_ipol,itermax))
     !
     ALLOCATE(a(itermax))
     ALLOCATE(b(itermax-1))
     ALLOCATE(c(itermax-1))
     ALLOCATE(r(n_ipol,itermax))
     !
     a(:) = (0._DP,0._DP)
     b(:) = (0._DP,0._DP)
     c(:) = (0._DP,0._DP)
     r(:,:) = (0._DP,0._DP)
     !
     ! Read beta, gamma, and zeta coefficients
     !
     CALL read_b_g_z_file()
     !
     ! Optional: use an extrapolation scheme
     !
     CALL extrapolate()
     !
     ! Spectrum calculation
     !
     WRITE(stdout,'(/5x,"Data ready, starting to calculate observables...")')
     WRITE(stdout,'(/5x,"Broadening = ",f15.8," Ry")') epsil
     !
     filename = TRIM(qe_prefix) // '.plot_chi.dat'
     filename_plot = TRIM(qe_prefix) // '.abs_spectrum.dat'
     !
     WRITE(stdout,'(/5x,"Output file name: ",a)') TRIM(filename)
     WRITE(stdout,'(/5x,"Output file name: ",a)') TRIM(filename_plot)
     !
     WRITE(stdout,'(/,5x,"chi_i_j: dipole polarizability tensor in units of e^2*a_0^2/energy")')
     !
     IF(n_ipol == 3) THEN
        WRITE(stdout,'(/,5x,"S: oscillator strength in units of 1/energy")')
        WRITE(stdout,'(/,5x,"S(\hbar \omega) = 2m/( 3 \pi e^2 \hbar) \omega sum_j chi_j_j")')
        WRITE(stdout,'(/,5x,"S(\hbar \omega) satisfies the f-sum rule: \int_0^\infty dE S(E) = N_el ")')
     ELSE
        WRITE(stdout,'(/,5x,"Insufficent info for S")')
     ENDIF
     !
     ! Units
     !
     IF(units == 0) THEN
        WRITE(stdout,'(/,5x,"Functions are reported in \hbar.\omega Energy unit is (Ry)")')
     ELSEIF(units == 1) THEN
        WRITE(stdout,'(/,5x,"Functions are reported in \hbar.\omega Energy unit is (eV)")')
     ELSEIF(units == 2) THEN
        WRITE(stdout,'(/,5x,"Functions are reported in (nm), Energy unit is (eV)")')
     ENDIF
     !
     ! The static dipole polarizability / static charge-density susceptibility
     !
     WRITE(stdout,'(/,5x,"Static dipole polarizability tensor:")')
     !
     CALL calc_chi(0._DP,epsil,green(:,:))
     !
     DO ip = 1,n_ipol
        DO ip2 = 1,n_ipol
           IF(n_ipol == 3) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e15.8," + i",e15.8)') &
                           & ip2, ip, REAL(green(ip,ip2),KIND=DP), AIMAG(green(ip,ip2))
           IF(n_ipol == 1) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e15.8," + i",e15.8)') &
                           & ipol, ipol, REAL(green(ip,ip2),KIND=DP), AIMAG(green(ip,ip2))
        ENDDO
     ENDDO
     !
     ! Open the output file
     !
     OPEN(NEWUNIT=iun1,FILE=filename)
     OPEN(NEWUNIT=iun2,FILE=filename_plot)
     !
     !-----------------------------------------------------------------------------!
     !                          PERCEIVED COLOR ANALYSIS                           !
     !-----------------------------------------------------------------------------!
     !
     ! The perceived color analysis uses the perception fit from the following program:
     ! RGB VALUES FOR VISIBLE WAVELENGTHS by Dan Bruton (astro@tamu.edu)
     ! Let's see if the environment is suitable for perceived color analysis
     ! This is needed for optics, not for EELS.
     !
     do_perceived = .FALSE.
     !
     IF(units == 0 .AND. start < vis_start .AND. end > vis_end .AND. n_ipol == 3) THEN
        do_perceived = .TRUE.
        perceived_itermax = INT((vis_end-vis_start)/increment)
     ELSEIF(units == 2 .AND. start < vis_end_wl .AND. end > vis_start_wl .AND. n_ipol == 3) THEN
        do_perceived = .TRUE.
        perceived_itermax = INT((vis_start_wl-vis_end_wl)/increment)
     ENDIF
     !
     IF(do_perceived) THEN
        WRITE(stdout,'(/,5x,"Will attempt to calculate perceived color")')
        perceived_iter = 1
        ALLOCATE(perceived_intensity(perceived_itermax))
        ALLOCATE(perceived_evaluated(perceived_itermax))
        perceived_intensity(:) = 0._DP
        perceived_evaluated(:) = -1._DP
        perceived_renorm = -9999999._DP
     ENDIF
     !
     ! Header of the output plot file
     !
     IF(units == 0) THEN
        WRITE(iun1,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (Ry) Re(chi) (e^2*a_0^2/Ry) Im(chi) (e^2*a_0^2/Ry)")')
     ELSEIF(units == 1) THEN
        WRITE(iun1,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (eV) Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV)")')
     ELSEIF(units == 2) THEN
        WRITE(iun1,'("#Chi is reported as CHI_(i)_(j) wavelength (nm) Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV)")')
     ENDIF
     !
     IF(n_ipol == 3) THEN
        WRITE(iun2,'("# Trchi(w) satisfies the sum rule")')
     ENDIF
     !
     ! Start a loop on frequency
     !
     ! Units conversion and omega history
     !
     omega(1) = omega(2)
     omega(2) = omega(3)
     !
     IF(units == 0) THEN
        omega(3) = start
     ELSEIF(units == 1) THEN
        omega(3) = start/rytoev
     ELSEIF(units == 2) THEN
        omega(3) = rytonm/start
     ENDIF
     !
     !-------------------------------------------------------------!
     !                          FIRST STEP                         !
     !-------------------------------------------------------------!
     !
     ! In order to gain speed, we perform the first step seperately
     !
     IF(n_ipol == 3) THEN
        !
        ! Calculation of the susceptibility
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF(units == 1 .OR. units == 2) THEN
           green(:,:) = green(:,:)/rytoev
        ENDIF
        !
        DO ip = 1,n_ipol
           DO ip2 = 1,n_ipol
              WRITE(iun1,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
              & ip2, ip, start, REAL(green(ip,ip2),KIND=DP), AIMAG(green(ip,ip2))
           ENDDO
        ENDDO
        !
        alpha_temp(3) = omega(3) * AIMAG(green(1,1)+green(2,2)+green(3,3))/(pi*3._DP)
        !
        IF(units == 1 .OR. units == 2) THEN
           alpha_temp(3) = alpha_temp(3)*rytoev
        ENDIF
        !
        ! alpha is ready
        !
        WRITE(iun2,'(5x,"Trchi(w)=",2x,2(e15.8,2x))') start, alpha_temp(3)
        !
        ! This is for the f-sum rule
        !
        f_sum = increment*alpha_temp(3)/3._DP
        !
        start = start + increment
        !
     ENDIF
     !
     !--------------------------------------------------------------!
     !                       OMEGA LOOP                             !
     !--------------------------------------------------------------!
     !
     DO WHILE(start < end)
        !
        ! Units conversion and omega history
        !
        omega(1) = omega(2)
        omega(2) = omega(3)
        !
        IF(units == 0) THEN
           omega(3) = start
        ELSEIF(units == 1) THEN
           omega(3) = start/rytoev
        ELSEIF(units == 2) THEN
           omega(3) = rytonm/start
        ENDIF
        !
        ! Calculation of the susceptibility for a given frequency omega.
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF(units == 1 .OR. units == 2) THEN
           green(:,:) = green(:,:)/rytoev
        ENDIF
        !
        ! Writing of chi
        !
        DO ip = 1,n_ipol
           DO ip2 = 1,n_ipol
              IF(n_ipol == 3) WRITE(iun1,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
                              & ip2, ip, start, REAL(green(ip,ip2),KIND=DP), AIMAG(green(ip,ip2))
              IF(n_ipol == 1) WRITE(iun1,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
                              & ipol, ipol, start, REAL(green(ip,ip2),KIND=DP), AIMAG(green(ip,ip2))
           ENDDO
        ENDDO
        !
        IF(n_ipol == 3) THEN
           !
           ! Calculation of the absorption coefficient
           !
           alpha_temp(1) = alpha_temp(2)
           alpha_temp(2) = alpha_temp(3)
           alpha_temp(3) = omega(3) * AIMAG(green(1,1)+green(2,2)+green(3,3))/(pi*3._DP)
           !
           IF(units == 1 .OR. units == 2) THEN
              alpha_temp(3) = alpha_temp(3)*rytoev
           ENDIF
           !
           ! alpha is ready
           !
           WRITE(iun2,'(5x,"Trchi(w)=",2x,2(e15.8,2x))') start, alpha_temp(3)
           !
           IF(is_peak(omega(3),alpha_temp(3))) &
               WRITE(stdout,'(5x,"Possible peak at ",F15.8," Ry; Intensity=",E11.2)') omega(1),alpha_temp(1)
           !
           ! f-sum rule
           !
           f_sum = f_sum + integrator(increment,alpha_temp(3))
           !
           ! Perceived color analysis
           !
           IF(omega(3) < vis_end .AND. omega(3) > vis_start .AND. do_perceived) THEN
              !
              perceived_intensity(perceived_iter) = alpha_temp(3)
              perceived_evaluated(perceived_iter) = omega(3)
              perceived_iter = perceived_iter + 1
              !
              ! Renormalization to 1
              !
              IF(alpha_temp(3) > perceived_renorm) perceived_renorm = alpha_temp(3)
              !
           ENDIF
           !
        ENDIF
        !
        start = start + increment
        !
     ENDDO
     !
     !------------------------------------------------------------------!
     !                      END OF OMEGA LOOP                           !
     !------------------------------------------------------------------!
     !
     !------------------------------------------------------------------!
     !                          LAST STEP                               !
     !------------------------------------------------------------------!
     !
     ! In order to gain speed, we perform the last step seperately
     !
     IF(n_ipol == 3) THEN
        !
        ! Units conversion
        !
        IF(units == 0) THEN
           omega(3) = start
        ELSEIF(units == 1) THEN
           omega(3) = start/rytoev
        ELSEIF(units == 2) THEN
           omega(3) = rytonm/start
        ENDIF
        !
        ! Calculation of the susceptibility
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF(units == 1 .OR. units == 2) THEN
           green(:,:) = green(:,:)/rytoev
        ENDIF
        !
        DO ip = 1,n_ipol
           DO ip2 = 1,n_ipol
              WRITE(iun1,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
              & ip2, ip, start, REAL(green(ip,ip2),KIND=DP), AIMAG(green(ip,ip2))
           ENDDO
        ENDDO
        !
        ! Absorption coefficient
        !
        alpha_temp(3) = omega(3) * AIMAG(green(1,1)+green(2,2)+green(3,3))/(pi*3._DP)
        !
        IF(units == 1 .OR. units == 2) THEN
           alpha_temp(3) = alpha_temp(3)*rytoev
        ENDIF
        !
        ! alpha is ready
        !
        WRITE(iun2,'(5x,"Trchi(w)=",2x,2(e15.8,2x))') start, alpha_temp(3)
        !
        ! alpha is ready
        !
        f_sum = f_sum + increment*alpha_temp(3)/3._DP
        !
     ENDIF
     !
     CLOSE(iun1)
     !
     IF(n_ipol == 3) THEN
        WRITE(stdout,'(5x,"Integral of absorbtion coefficient ",F15.8)') f_sum
     ENDIF
     !
     !------------------------------------------------------------------------!
     !                      Perceived color analysis                          !
     !------------------------------------------------------------------------!
     !
     IF(ALLOCATED(perceived_intensity)) THEN
        !
        WRITE(stdout,'(5x,"Perceived color analysis is experimental")')
        perceived_intensity(:) = perceived_intensity(:)/perceived_renorm
        perceived_intensity(:) = 1._DP-perceived_intensity(:) !inverse spectrum
        !
        DO j = 1,INT(perceived_itermax/8)
           DO i = 1,perceived_itermax
              !
              wl = 91.1266519_DP/perceived_evaluated(i) !hc/hbar.omega=lambda (hbar.omega in rydberg units)
              !
              ! WARNING alpha_temp duty change: now contains R G and B
              !
              CALL wl_to_color(wl,alpha_temp(1),alpha_temp(2),alpha_temp(3))
              !
              ! Now the intensities
              ! First the degradation toward the end
              !
              IF(wl > 700._DP) THEN
                 scale = 0.3_DP+0.7_DP* (780._DP-wl)/(780._DP-700._DP)
              ELSEIF(wl < 420._DP) THEN
                 scale = 0.3_DP+0.7_DP*(wl-380._DP)/(420._DP-380._DP)
              ELSE
                 scale = 1._DP
              ENDIF
              !
              alpha_temp(:) = scale*alpha_temp(:)
              !
              ! Then the data from absorbtion spectrum
              !
              alpha_temp(:) = perceived_intensity(i)*alpha_temp(:)
              !
              ! The perceived color can also be calculated here
              !
              IF(j == 1) THEN
                 perceived_red = perceived_red+alpha_temp(1)
                 perceived_green = perceived_green+alpha_temp(2)
                 perceived_blue = perceived_blue+alpha_temp(3)
              ENDIF
              !
           ENDDO
        ENDDO
        !
        perceived_red = perceived_red/(1._DP*perceived_itermax)
        perceived_green = perceived_green/(1._DP*perceived_itermax)
        perceived_blue = perceived_blue/(1._DP*perceived_itermax)
        !
        WRITE(stdout,'(5x,"Perceived R G B ",3(F15.8,1X))') &
        & perceived_red,perceived_green,perceived_blue
        WRITE(stdout,'(5x,"Perceived R G B ",3(F15.8,1X))') &
        & perceived_red*255._DP,perceived_green*255._DP,perceived_blue*255._DP
        !
     ENDIF
     !
     !----------------------------------------------------------------------!
     !                   End of perceived color analysis                    !
     !----------------------------------------------------------------------!
     !
     ! f-sum rule
     !
     start = 0._DP
     f_sum = 0._DP
     !
     DO WHILE(start < end)
        f_sum = f_sum + integrator(increment,start)
        start = start + increment
     ENDDO
     !
     f_sum = f_sum + increment*start/3._DP
     !
     WRITE(stdout,'(5x,"Integral test:",F15.8,"  Actual: ",F15.8:)') f_sum, 0.5_DP*start*start
     !
     ! Deallocations
     !
     IF(ALLOCATED(perceived_intensity)) DEALLOCATE(perceived_intensity)
     IF(ALLOCATED(perceived_evaluated)) DEALLOCATE(perceived_evaluated)
     !
     IF(ALLOCATED(beta_store))  DEALLOCATE(beta_store)
     IF(ALLOCATED(gamma_store)) DEALLOCATE(gamma_store)
     IF(ALLOCATED(zeta_store))  DEALLOCATE(zeta_store)
     !
     DEALLOCATE(a)
     DEALLOCATE(b)
     DEALLOCATE(c)
     DEALLOCATE(r)
     !
  ENDIF
  !
  CALL mp_barrier(world_comm)
  CALL mp_global_end ()
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION is_peak(omeg,alpha)
    !---------------------------------------------------------------------
    !
    ! A simple algorithm for detecting peaks.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: omeg, alpha !x and y
    !
    ! Local variables
    !
    REAL(DP), SAVE :: omeg_save = 0._DP
    REAL(DP), SAVE :: thm1
    REAL(DP), SAVE :: h2m1
    REAL(DP), SAVE :: first_der_save = 9.E99_DP
    REAL(DP), SAVE :: alpha_save(3) = 0._DP
    INTEGER, SAVE :: current_iter = 0
    LOGICAL, SAVE :: trigger = .TRUE.
    REAL(DP) :: first_der, second_der
    !
    is_peak = .FALSE.
    !
    ! counter
    ! Rotate the variables
    !
    IF(current_iter < 3) THEN
        current_iter = current_iter + 1
        omeg_save = omeg
        alpha_save(current_iter) = alpha
        RETURN
    ELSE
        IF(current_iter == 3) THEN
           current_iter = current_iter + 1
           thm1 = omeg-omeg_save
           h2m1 = 1._DP/(thm1*thm1) !for second derivative
           thm1 = 0.5_DP/thm1       !for first derivative
        ENDIF
        alpha_save(1) = alpha_save(2) !t-h
        alpha_save(2) = alpha_save(3) !t
        alpha_save(3) = alpha         !t+h
    ENDIF
    !
    !The derivatives
    !
    first_der = (alpha_save(3)-alpha_save(1))*thm1
    second_der = (alpha_save(3)-2._DP*alpha_save(2)+alpha_save(1))*h2m1
    !
    IF(second_der < 0) THEN
       IF(trigger) THEN
          IF(ABS(first_der) < ABS(first_der_save)) THEN
             first_der_save = first_der
             RETURN
          ELSE
             is_peak = .TRUE.
             trigger = .FALSE.
             RETURN
          ENDIF
       ENDIF
    ELSE
       first_der_save = 9.E99_DP
       trigger = .TRUE.
    ENDIF
    !
  END FUNCTION
  !
  !-----------------------------------------------------------------------
  REAL(DP) FUNCTION integrator(dh,alpha)
    !---------------------------------------------------------------------
    !
    ! Calculates an integral every three points, using the Simpson's rule.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: dh, alpha !x and y
    !
    LOGICAL, SAVE :: flag = .TRUE.
    !
    ! COMPOSITE SIMPSON INTEGRATOR, (precision level ~ float)
    ! \int a b f(x) dx = ~ h/3 (f(a) + \sum_odd-n 2*f(a+n*h) + \sum_even-n 4*f(a+n*h) +f(b))
    !
    integrator = 0._DP
    !
    IF(flag) THEN
       !
       ! odd steps
       !
       integrator = (4._DP/3._DP)*dh*alpha
       flag = .FALSE.
       !
    ELSE
       !
       ! even steps
       !
       integrator = (2._DP/3._DP)*dh*alpha
       flag = .TRUE.
       !
    ENDIF
    !
  END FUNCTION
  !
  !-----------------------------------------------------------------------
  SUBROUTINE read_b_g_z_file()
    !---------------------------------------------------------------------
    !
    ! Reads beta, gamma, and zeta coefficients.
    !
    USE iotk_module
    !
    IMPLICIT NONE
    !
    INTEGER :: iun, ierr, is
    INTEGER :: nipol_input_tmp, nspin
    CHARACTER(LEN=3) :: ipol_label_tmp
    INTEGER :: nlan_tmp
    REAL(DP), ALLOCATABLE :: beta_store_tmp(:,:)
    REAL(DP), ALLOCATABLE :: gamma_store_tmp(:,:)
    COMPLEX(DP), ALLOCATABLE :: zeta_store_tmp(:,:,:)
    CHARACTER(LEN=4) :: my_ip
    CHARACTER(LEN=256) :: file_ip
    !
    IF(sym_op == 0) THEN
       !
       DO ip = 1,n_ipol
          !
          WRITE(my_ip,'(i1)')ip
          file_ip = 'summary.'//TRIM(my_ip)//'.xml'
          !
          CALL iotk_free_unit(iun, ierr)
          CALL iotk_open_read(iun, FILE=TRIM(file_ip), IERR=ierr)
          !
          CALL iotk_scan_begin(iun, 'SUMMARY')
          CALL iotk_scan_dat(iun, 'nspin', nspin)
          CALL iotk_scan_dat(iun, 'nipol_input', nipol_input_tmp)
          CALL iotk_scan_dat(iun, 'ipol_label', ipol_label_tmp)
          CALL iotk_scan_dat(iun, 'n_lanczos', nlan_tmp)
          CALL iotk_scan_end(iun, 'SUMMARY')
          !
          WRITE(stdout,*)
          WRITE(stdout,'(5x,a)') 'Reading beta, gamma, zeta of the polarzation ' &
          & //TRIM(ipol_label_tmp)//' from file '//TRIM(file_ip)
          !
          IF(nlan_tmp < itermax0) THEN
             CALL errore('read_b_g_z_file', 'Error in lriter_stop < itermax0, reduce itermax0',1)
          ENDIF
          !
          ALLOCATE(beta_store_tmp(nlan_tmp,nspin))
          ALLOCATE(gamma_store_tmp(nlan_tmp,nspin))
          !
          CALL iotk_scan_begin(iun, 'BETA_STORE')
          DO is = 1,nspin
             CALL iotk_scan_dat(iun, 'beta_store_k', beta_store_tmp(1:nlan_tmp,is))
          ENDDO
          CALL iotk_scan_end(iun, 'BETA_STORE')
          !
          CALL iotk_scan_begin(iun, 'GAMMA_STORE')
          DO is = 1,nspin
             CALL iotk_scan_dat(iun, 'gamma_store_k', gamma_store_tmp(1:nlan_tmp,is))
          ENDDO
          CALL iotk_scan_end(iun, 'GAMMA_STORE')
          !
          ALLOCATE(zeta_store_tmp(3,nlan_tmp,nspin))
          !
          CALL iotk_scan_begin(iun, 'ZETA_STORE')
          !
          DO ip2 = 1,3
             CALL iotk_scan_dat(iun, 'zeta_store_ipol_j', ip2)
             DO is = 1,nspin
                CALL iotk_scan_dat(iun, 'zeta_store_k', zeta_store_tmp(ip2,1:nlan_tmp,is))
             ENDDO
          ENDDO
          CALL iotk_scan_end(iun, 'ZETA_STORE')
          !
          CALL iotk_close_read(iun)
          !
          IF(nspin == 1) spin_channel = 1
          !
          norm0(ip)                    = beta_store_tmp(1, spin_channel)
          beta_store(ip,1:itermax0-1)  = beta_store_tmp(2:itermax0, spin_channel)
          gamma_store(ip,1:itermax0-1) = gamma_store_tmp(2:itermax0, spin_channel)
          beta_store(ip,itermax0)      = beta_store_tmp(itermax0, spin_channel)
          gamma_store(ip,itermax0)     = gamma_store_tmp(itermax0, spin_channel)
          !
          IF(n_ipol == 1) THEN
             zeta_store(1,1,1:itermax0) = zeta_store_tmp(ipol,1:itermax0,spin_channel)
          ELSE
             zeta_store(ip,:,1:itermax0) = zeta_store_tmp(:,1:itermax0,spin_channel)
          ENDIF
          !
          DEALLOCATE(beta_store_tmp)
          DEALLOCATE(gamma_store_tmp)
          DEALLOCATE(zeta_store_tmp)
          !
       ENDDO
       !
       IF(nspin == 1) degspin = 2
       IF(nspin == 2) degspin = 1
       !
       beta_store(:,itermax0+1:itermax)   = 0._DP
       gamma_store(:,itermax0+1:itermax)  = 0._DP
       zeta_store(:,:,itermax0+1:itermax) = (0._DP,0._DP)
       !
    ENDIF
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE extrapolate()
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    skip = .FALSE.
    !
    IF(TRIM(extrapolation) /= 'no') THEN
       !
       average = 0._DP
       av_amplitude = 0._DP
       !
       DO ip = 1,n_ipol
          !
          WRITE(stdout,'(/5x,"Polarization direction:",I1)') ip
          counter = 0
          !
          DO i = 151,itermax0
             IF(skip) THEN
                skip = .FALSE.
                CYCLE
             ENDIF
             !
             IF(MOD(i,2) == 1) THEN
                IF(i /= 151 .AND. ABS(beta_store(ip,i)-average(ip)/counter) > 2._DP) THEN
                   skip = .TRUE.
                ELSE
                   average(ip) = average(ip)+beta_store(ip,i)
                   av_amplitude(ip) = av_amplitude(ip)+beta_store(ip,i)
                   counter = counter+1
                ENDIF
             ELSE
                IF(i /= 151 .AND. ABS(beta_store(ip,i)-average(ip)/counter) > 2._DP) THEN
                   skip = .TRUE.
                ELSE
                   average(ip) = average(ip)+beta_store(ip,i)
                   av_amplitude(ip) = av_amplitude(ip)-beta_store(ip,i)
                   counter = counter+1
                ENDIF
             ENDIF
          ENDDO
          !
          average(ip) = average(ip)/counter
          av_amplitude(ip) = av_amplitude(ip)/counter
          !
          WRITE(stdout,'(5x,"Lanczos coefficients:")')
          WRITE(stdout,'(5x,"Average =",3F15.8)') average(ip)
          WRITE(stdout,'(5x,"Average oscillation amplitude =",F15.8)') av_amplitude(ip)
          !
       ENDDO
       !
       IF(TRIM(extrapolation) == 'constant') av_amplitude = 0
       !
       DO ip = 1,n_ipol
          DO i = itermax0,itermax
             IF(MOD(i,2) == 1) THEN
                beta_store(ip,i) = average(ip)+av_amplitude(ip)
                gamma_store(ip,i) = average(ip)+av_amplitude(ip)
             ELSE
                beta_store(ip,i) = average(ip)-av_amplitude(ip)
                gamma_store(ip,i) = average(ip)-av_amplitude(ip)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calc_chi(freq,broad,chi)
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: freq
    REAL(DP), INTENT(IN) :: broad
    COMPLEX(DP), INTENT(OUT) :: chi(:,:)
    !
    omeg_c = CMPLX(freq,broad,KIND=DP)
    !
    DO ip = 1,n_ipol
       !
       a(:) = omeg_c
       !
       DO i = 1,itermax-1
          b(i) = CMPLX(-beta_store(ip,i),KIND=DP)
          c(i) = CMPLX(-gamma_store(ip,i),KIND=DP)
       ENDDO
       !
       r(ip,:) = (0._DP,0._DP)
       r(ip,1) = (1._DP,0._DP)
       !
       ! |w_t|=(w-L) |1,0,0,...,0|
       !
       CALL zgtsv(itermax,1,b,a,c,r(ip,:),itermax,info)
       !
       IF(info /= 0) CALL errore('calc_chi','Unable to solve tridiagonal system',1)
       !
       ! p=-div.rho'
       ! p= chi . E
       ! Thus, chi = - <zeta|w_t>
       !
       ! Notice that brodening has a positive sign,
       ! thus the abs. coefficient is Im(tr(chi)) not -Im(Tr(chi)) as usual
       !
       DO ip2 = 1,n_ipol
           !
           chi(ip,ip2) = ZDOTC(itermax,zeta_store(ip,ip2,:),1,r(ip,:),1)
           !
           ! Multiplication with a norm
           !
           chi(ip,ip2) = chi(ip,ip2) * CMPLX(norm0(ip),KIND=DP)
           !
           ! The response charge density is defined as 2*evc0*q, see Eq. (43) in
           ! JCP 128, 154105 (2008).
           ! Therefore, the dipole is given by 2*degspin* zeta^T *
           ! (w-T^itermax)^-1 * e_1. See also Eq. (15) in that paper.
           ! Optics: The minus sign accounts for the negative electron charge
           ! (perturbation is -e E x, rather than E x)
           !
           chi(ip,ip2) = chi(ip,ip2) * CMPLX(-2._DP*degspin,KIND=DP)
           !
       ENDDO
       !
    ENDDO
    !
  END SUBROUTINE
  !
  !-----------------------------------------------------------------------
  SUBROUTINE wl_to_color(wavelength,red,green,blue)
    !---------------------------------------------------------------------
    !
    ! Gives the colour intensity of a given wavelength in terms of RGB (red, green and blue).
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: wavelength
    REAL(DP), INTENT(OUT) :: red,green,blue
    !
    IF(wavelength >= 380._DP .AND. wavelength <= 440._DP) THEN
       red = -1._DP*(wavelength-440._DP)/(440._DP-380._DP)
       green = 0._DP
       blue = 1._DP
    ENDIF
    !
    IF(wavelength >= 440._DP .AND. wavelength<=490._DP) THEN
       red = 0._DP
       green = (wavelength-440._DP)/(490._DP-440._DP)
       blue = 1._DP
    ENDIF
    !
    IF(wavelength >= 490._DP .AND. wavelength<=510._DP) THEN
       red = 0._DP
       green = 1._DP
       blue = -1._DP*(wavelength-510._DP)/(510._DP-490._DP)
    ENDIF
    !
    IF(wavelength >= 510._DP .AND. wavelength <= 580._DP) THEN
       red = (wavelength-510._DP)/(580._DP-510._DP)
       green = 1._DP
       blue = 0._DP
    ENDIF
    !
    IF(wavelength >= 580._DP .AND. wavelength <= 645._DP) THEN
       red = 1._DP
       green = -1._DP*(wavelength-645._DP)/(645._DP-580._DP)
       blue = 0._DP
    ENDIF
    !
    IF(wavelength >= 645._DP .AND. wavelength <= 780._DP) THEN
       red = 1._DP
       green = 0._DP
       blue = 0._DP
    ENDIF
    !
  END SUBROUTINE
  !
END SUBROUTINE
