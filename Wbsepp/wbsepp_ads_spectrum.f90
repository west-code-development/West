!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE wbsepp_ads_spectrum()
  !---------------------------------------------------------------------
  !
  ! This routine is modified from the TDDFT code.
  ! Pls do not distribute it to avoid the copyright problem.
  ! Todo: write this code in Python to make it become our
  !
  USE kinds,               ONLY : dp
  USE constants,           ONLY : pi,rytoev,evtonm,rytonm
  USE io_files,            ONLY : tmp_dir, nd_nmbr
  USE global_version,      ONLY : version_number
  USE io_global,           ONLY : stdout, ionode, ionode_id
  USE environment,         ONLY : environment_start,environment_end
  USE mp_global,           ONLY : mp_startup,mp_global_end, my_image_id
  USE mp_world,            ONLY : world_comm
  USE mp,                  ONLY : mp_bcast, mp_barrier
  USE westcom,             ONLY : qe_prefix
  USE westcom,           ONLY : itermax, itermax0, extrapolation, &
                                  start, end, increment, epsil, ipol, &
                                  sym_op, verbosity, units, spin_channel
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256), EXTERNAL :: trimcheck
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
  REAL(dp) :: omega(3)
  REAL(kind=dp) :: omegmax,delta_omeg
  CHARACTER(len=256):: filename
  CHARACTER(len=256):: filename_plot
  !
  ! General use variables & counters
  !
  INTEGER :: n_ipol, i,j, info, ip, ip2, counter, ios
  REAL(dp) :: norm0(3), volume, &
            & alat, q1, q2, q3, modulus_q, nelec, &
            & average(3), av_amplitude(3), &
            & alpha_temp(3), scale, wl, &
            & omeg, z1,z2, degspin, integration_function, start_save, f_sum
  !
  COMPLEX(kind=dp) :: omeg_c
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: beta_store, gamma_store
  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: zeta_store
  COMPLEX(kind=dp) :: green(3,3), &  ! susceptibility chi
                      eps(3,3),   &  ! dielectric function
                      epsm1(3,3)     ! inverse dielectric function
  COMPLEX(kind=dp), ALLOCATABLE :: a(:), b(:), c(:), r(:,:)
  LOGICAL :: skip, exst
  !
  ! For perceived color analysis
  !
  REAL(dp),PARAMETER :: vis_start    = 0.116829041, &
                      & vis_start_wl = 780
  REAL(dp),PARAMETER :: vis_end      = 0.239806979, &
                      & vis_end_wl   = 380
  REAL(dp) :: perceived_red   = 0.0d0, &
            & perceived_green = 0.0d0, &
            & perceived_blue  = 0.0d0
  REAL(dp) :: perceived_renorm
  INTEGER  :: perceived_itermax,perceived_iter
  REAL(dp), ALLOCATABLE :: perceived_intensity(:), &
                         & perceived_evaluated(:)
  LOGICAL  :: do_perceived
  !
  ! Subroutines etc.
  !
  COMPLEX(kind=dp), EXTERNAL :: zdotc
  !
  ! User controlled variable initialisation
  !
  CALL mp_startup ( )
  !
  IF (ionode) THEN
     WRITE (stdout,*) ' '
     WRITE (stdout,*) ' '
     WRITE(stdout,'( 10x, " Warning: Only a single CPU will be used ")')
  ENDIF
  !
  IF (ionode) THEN
     !
     IF (itermax0 < 151 .and. trim(extrapolation)/="no") THEN
        WRITE(stdout,'(5x, "Itermax0 is less than 150, no extrapolation scheme can be used ")')
        extrapolation="no"
     ENDIF
     !
     IF (ipol < 4) THEN
        n_ipol=1
     ELSE
        n_ipol=3
        ipol = 1
     ENDIF
     !
     ! Polarization symmetry
     !
     IF ( .not. sym_op == 0 ) THEN
        !
        IF (sym_op == 1) THEN
           WRITE(stdout,'(5x,"All polarization axes will be considered to be equal.")')
           n_ipol = 3
           ipol = 1
        ELSE
           CALL errore("lr_calculate_spectrum","Unsupported symmetry operation",1)
        ENDIF
        !
     ENDIF
     !
     ! Terminator scheme
     !
     IF (trim(extrapolation)=="no") THEN
        !
        itermax = itermax0
        !
     ENDIF
     !
     ! Check the units (Ry, eV, nm)
     !
     IF (units < 0 .or. units >2) CALL errore("lr_calculate_spectrum","Unsupported unit system",1)
     !
     IF ( units /= 0 .and. verbosity > 4) THEN
        WRITE(stdout,'(5x,"Such a high verbosity is not supported when &
                        & non-default units are used")')
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
     a(:) = (0.0d0,0.0d0)
     b(:) = (0.0d0,0.0d0)
     c(:) = (0.0d0,0.0d0)
     r(:,:) = (0.0d0,0.0d0)
     !
     ! Read beta, gamma, and zeta coefficients
     !
     CALL read_b_g_z_file_html()
     !
     ! Optional: use an extrapolation scheme
     !
     CALL extrapolate()
     !
     !  Spectrum calculation
     !
     WRITE (stdout,'(/5x,"Data ready, starting to calculate observables...")')
     WRITE (stdout,'(/5x,"Broadening = ",f15.8," Ry")') epsil
     !
     filename = trim(qe_prefix) // ".plot_chi.dat"
     filename_plot = trim(qe_prefix) // ".abs_spectrum.dat"
     !
     WRITE (stdout,'(/5x,"Output file name: ",A20)') filename
     WRITE (stdout,'(/5x,"Output file name: ",A20)') filename_plot
     !
     WRITE(stdout,'(/,5x,"chi_i_j: dipole polarizability tensor &
                                & in units of e^2*a_0^2/energy")')
     !
     IF (n_ipol == 3) THEN
        WRITE(stdout,'(/,5x,"S: oscillator strength in units of 1/energy")')
        WRITE(stdout,'(/,5x,"S(\hbar \omega) = 2m/( 3 \pi e^2 \hbar) &
                               & \omega sum_j chi_j_j")')
        WRITE(stdout,'(/,5x,"S(\hbar \omega) satisfies the f-sum rule: &
                               & \int_0^\infty dE S(E) = N_el ")')
     ELSE
        WRITE (stdout,'(/,5x,"Insufficent info for S")')
     ENDIF
     !
     ! Units
     !
     IF (units == 0) THEN      ! Ry
        WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega &
                              & Energy unit is (Ry)")')
     ELSEIF (units == 1) THEN  ! eV
        WRITE (stdout,'(/,5x,"Functions are reported in \hbar.\omega &
                              & Energy unit is (eV)")')
     ELSEIF (units == 2) THEN  ! nm
        WRITE (stdout,'(/,5x,"Functions are reported in (nm), &
                              & Energy unit is (eV) ")')
     ENDIF
     !
     ! The static dipole polarizability / static charge-density susceptibility
     !
     WRITE (stdout,'(/,5x,"Static dipole polarizability tensor:")')
     !
     CALL calc_chi(0.0d0,epsil,green(:,:))
     !
     DO ip=1,n_ipol
        DO ip2=1,n_ipol
           !
           IF (n_ipol == 3) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e15.8," + i",e15.8)') &
                                          & ip2, ip, dble(green(ip,ip2)), aimag(green(ip,ip2))
           !
           IF (n_ipol == 1) WRITE(stdout,'(5x,"chi_",i1,"_",i1,"=",2x,e15.8," + i",e15.8)') &
                                          & ipol, ipol, dble(green(ip,ip2)), aimag(green(ip,ip2))
           !
        ENDDO
     ENDDO
     !
     ! Open the output file
     !
     OPEN(17,file=filename,status="unknown")
     OPEN(18,file=filename_plot,status="unknown")
     !
     !-----------------------------------------------------------------------------!
     !                          PERCEIVED COLOR ANALYSIS                           !
     !-----------------------------------------------------------------------------!
     !
     ! The perceived color analysis uses the perception fit from the following program:
     ! RGB VALUES FOR VISIBLE WAVELENGTHS   by Dan Bruton (astro@tamu.edu)
     ! Let's see if the environment is suitable for perceived color analysis
     ! This is needed for optics, not for EELS.
     !
     do_perceived = .false.
     !
     IF (units == 0 .and. start<vis_start .and. end>vis_end .and. n_ipol == 3) THEN
        !
        WRITE (stdout,'(/,5x,"Will attempt to calculate perceived color")')
        do_perceived=.true.
        perceived_iter=1
        perceived_itermax=int((vis_end-vis_start)/increment)
        ALLOCATE(perceived_intensity(perceived_itermax))
        ALLOCATE(perceived_evaluated(perceived_itermax))
        perceived_intensity(:)=0.0d0
        perceived_evaluated(:)=-1.0d0
        perceived_renorm=-9999999.0d0
        !
     ELSEIF (units == 2 .and. start<vis_end_wl .and. end>vis_start_wl .and. n_ipol == 3) THEN
        !
        WRITE (stdout,'(/,5x,"Will attempt to calculate perceived color")')
        do_perceived=.true.
        perceived_iter=1
        perceived_itermax=int((vis_start_wl-vis_end_wl)/increment)
        ALLOCATE(perceived_intensity(perceived_itermax))
        ALLOCATE(perceived_evaluated(perceived_itermax))
        perceived_intensity(:)=0.0d0
        perceived_evaluated(:)=-1.0d0
        perceived_renorm=-9999999.0d0
        !
     ENDIF
     !
     ! Header of the output plot file
     !
     IF (units == 0) THEN
        WRITE (17,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (Ry) &
                      & Re(chi) (e^2*a_0^2/Ry) Im(chi) (e^2*a_0^2/Ry) ")')
     ELSEIF (units == 1) THEN
        WRITE (17,'("#Chi is reported as CHI_(i)_(j) \hbar \omega (eV) &
                      & Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV) ")')
     ELSEIF (units == 2) THEN
        WRITE (17,'("#Chi is reported as CHI_(i)_(j) wavelength (nm) &
                      & Re(chi) (e^2*a_0^2/eV) Im(chi) (e^2*a_0^2/eV) ")')
     ENDIF
     !
     IF (n_ipol == 3) THEN
        WRITE(18,'("# Trchi(w) satisfies the sum rule ")' )
     ENDIF
     !
     ! Start a loop on frequency
     !
     ! Units conversion and omega history
     !
     omega(1) = omega(2)
     omega(2) = omega(3)
     !
     IF (units == 0) THEN
        omega(3) = start
     ELSEIF (units == 1) THEN
        omega(3) = start/rytoev
     ELSEIF (units == 2) THEN
        omega(3) = rytonm/start
     ENDIF
     !
     !-------------------------------------------------------------!
     !                          FIRST STEP                         !
     !-------------------------------------------------------------!
     !
     ! In order to gain speed, we perform the first step seperately
     !
     IF (n_ipol == 3) THEN
        !
        ! Calculation of the susceptibility
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF (units == 1 .or. units == 2) THEN
           !
           green(:,:) = green(:,:)/rytoev
           !
        ENDIF
        !
        DO ip=1,n_ipol
           DO ip2=1,n_ipol
              !
              WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
              !
           ENDDO
        ENDDO
        !
        alpha_temp(3) = omega(3) * aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
        !
        IF (units == 1 .or. units == 2) THEN
           alpha_temp(3) = alpha_temp(3)*rytoev
        ENDIF
        !
        ! alpha is ready
        !
        WRITE(18,'(5x,"Trchi(w)=",2x,2(e15.8,2x))') start, alpha_temp(3)
        !
        ! This is for the f-sum rule
        !
        f_sum = increment*alpha_temp(3)/3.0d0
        !
        start = start + increment
        !
     ENDIF
     !
     !--------------------------------------------------------------!
     !                       OMEGA LOOP                             !
     !--------------------------------------------------------------!
     !
     DO WHILE (start < end)
        !
        ! Units conversion and omega history
        !
        omega(1) = omega(2)
        omega(2) = omega(3)
        !
        IF (units == 0) THEN
           omega(3) = start
        ELSEIF (units == 1) THEN
           omega(3) = start/rytoev
        ELSEIF (units == 2) THEN
           omega(3) = rytonm/start
        ENDIF
        !
        ! Calculation of the susceptibility for a given frequency omega.
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF (units == 1 .or. units == 2) THEN
           !
           green(:,:) = green(:,:)/rytoev
           !
        ENDIF
        !
        ! Writing of chi
        !
        DO ip=1,n_ipol
           !
           DO ip2=1,n_ipol
              !
              !eps(ip,ip2)=(1.d0,0.d0)-(32.d0*pi/omega)*green(ip,ip2)
              !
              IF (n_ipol == 3) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
                           & ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
              IF (n_ipol == 1) WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
                           & ipol, ipol, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
              !
           ENDDO
           !
        ENDDO
        !
        IF (n_ipol==3) THEN
           !
           ! Calculation of the absorption coefficient
           !
           alpha_temp(1) = alpha_temp(2)
           alpha_temp(2) = alpha_temp(3)
           alpha_temp(3) = omega(3) * aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
           !
           IF (units == 1 .or. units == 2) THEN
              alpha_temp(3) = alpha_temp(3)*rytoev
           ENDIF
           !
           ! alpha is ready
           !
           WRITE(18,'(5x,"Trchi(w)=",2x,2(e15.8,2x))') start, alpha_temp(3)
           !
           IF ( is_peak(omega(3),alpha_temp(3))) &
                WRITE(stdout,'(5x,"Possible peak at ",F15.8," Ry; &
                                  & Intensity=",E11.2)') omega(1),alpha_temp(1)
           !
           ! f-sum rule
           !
           f_sum = f_sum + integrator(increment,alpha_temp(3))
           !
           ! Perceived color analysis
           !
           IF ( omega(3)<vis_end .and. omega(3)>vis_start .and. do_perceived ) THEN
              !
              perceived_intensity(perceived_iter) = alpha_temp(3)
              perceived_evaluated(perceived_iter) = omega(3)
              perceived_iter = perceived_iter + 1
              !
              ! Renormalization to 1
              !
              IF (alpha_temp(3) > perceived_renorm) perceived_renorm = alpha_temp(3)
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
     IF (n_ipol == 3) THEN
        !
        ! Units conversion
        !
        IF (units == 0) THEN
           omega(3) = start
        ELSEIF (units == 1) THEN
           omega(3) = start/rytoev
        ELSEIF (units == 2) THEN
           omega(3) = rytonm/start
        ENDIF
        !
        ! Calculation of the susceptibility
        !
        CALL calc_chi(omega(3),epsil,green(:,:))
        !
        IF (units == 1 .or. units == 2) THEN
           green(:,:) = green(:,:)/rytoev
        ENDIF
        !
        DO ip=1,n_ipol
           DO ip2=1,n_ipol
              !
              WRITE(17,'(5x,"chi_",i1,"_",i1,"=",2x,3(e15.8,2x))') &
                  ip2, ip, start, dble(green(ip,ip2)), aimag(green(ip,ip2))
              !
           ENDDO
        ENDDO
        !
        ! Absorption coefficient
        !
        alpha_temp(3) = omega(3) * aimag(green(1,1)+green(2,2)+green(3,3))/(pi*3.d0)
        !
        IF (units == 1 .or. units == 2) THEN
           alpha_temp(3) = alpha_temp(3)*rytoev
        ENDIF
        !
        ! alpha is ready
        !
        WRITE(18,'(5x,"Trchi(w)=",2x,2(e15.8,2x))') start, alpha_temp(3)
        !
        ! alpha is ready
        !
        f_sum = f_sum + increment*alpha_temp(3)/3.0d0
        !
     ENDIF
     !
     CLOSE(17)
     !
     IF (n_ipol==3) THEN
        !
        ! S(w)=2m_e/(pi e^2 hbar)
        WRITE(stdout,'(5x,"Integral of absorbtion coefficient ",F15.8)') f_sum
        !
     ENDIF
     !
     !------------------------------------------------------------------------!
     !                      Perceived color analysis                          !
     !------------------------------------------------------------------------!
     !
     IF (allocated(perceived_intensity)) THEN
        !
        WRITE(stdout,'(5x,"Perceived color analysis is experimental")')
        perceived_intensity(:)=perceived_intensity(:)/perceived_renorm
        perceived_intensity(:)=1.0d0-perceived_intensity(:) !inverse spectrum
        !
        DO j=1, int(perceived_itermax/8)
           !
           DO i=1,perceived_itermax
              !
              wl=91.1266519/perceived_evaluated(i) !hc/hbar.omega=lambda (hbar.omega in rydberg units)
              !
              ! WARNING alpha_temp duty change: now contains R G and B
              !
              CALL wl_to_color(wl,alpha_temp(1),alpha_temp(2),alpha_temp(3))
              !
              ! Now the intensities
              ! First the degradation toward the end
              !
              IF (wl >700) THEN
                 scale=.3+.7* (780.-wl)/(780.-700.)
              ELSEIF (wl<420.) THEN
                 scale=.3+.7*(wl-380.)/(420.-380.)
              ELSE
                 scale=1.
              ENDIF
              !
              alpha_temp(:)=scale*alpha_temp(:)
              !
              ! Then the data from absorbtion spectrum
              !
              alpha_temp(:)=perceived_intensity(i)*alpha_temp(:)
              !
              ! The perceived color can also be calculated here
              !
              IF (j==1) THEN
                 perceived_red=perceived_red+alpha_temp(1)
                 perceived_green=perceived_green+alpha_temp(2)
                 perceived_blue=perceived_blue+alpha_temp(3)
              ENDIF
              !
              IF (alpha_temp(1)>1.0d0) PRINT *,alpha_temp(1)
              !
           ENDDO
           !
        ENDDO
        !
        perceived_red  = perceived_red  /(1.0d0*perceived_itermax)
        perceived_green= perceived_green/(1.0d0*perceived_itermax)
        perceived_blue = perceived_blue /(1.0d0*perceived_itermax)
        !
        WRITE(stdout,'(5x,"Perceived R G B ",3(F15.8,1X))') perceived_red,perceived_green,perceived_blue
        WRITE(stdout,'(5x,"Perceived R G B ",3(F15.8,1X))') perceived_red*255,perceived_green*255,perceived_blue*255
        !
     ENDIF
     !
     !----------------------------------------------------------------------!
     !                   End of perceived color analysis                    !
     !----------------------------------------------------------------------!
     !
     ! f-sum rule
     !
     start = 0.0
     f_sum = 0.0
     !
     DO WHILE (start<end)
        !
        f_sum = f_sum + integrator(increment,start)
        start = start + increment
        !
     ENDDO
     !
     f_sum = f_sum + increment*start/3.0d0
     !
     WRITE(stdout,'(5x,"Integral test:",F15.8,"  Actual: ",F15.8:)') &
                                             f_sum, 0.5*start*start
     !
     ! Deallocations
     !
     IF (allocated(perceived_intensity)) DEALLOCATE(perceived_intensity)
     IF (allocated(perceived_evaluated)) DEALLOCATE(perceived_evaluated)
     !
     IF (allocated(beta_store))  DEALLOCATE(beta_store)
     IF (allocated(gamma_store)) DEALLOCATE(gamma_store)
     IF (allocated(zeta_store))  DEALLOCATE(zeta_store)
     !
     DEALLOCATE(a)
     DEALLOCATE(b)
     DEALLOCATE(c)
     DEALLOCATE(r)
     !
  ENDIF
  !
  CALL mp_barrier (world_comm)
  CALL mp_global_end ()
  !
  RETURN
  !
CONTAINS

LOGICAL FUNCTION is_peak(omeg,alpha)
  !-------------------------------------------------------------------------
  !
  ! A simple algorithm for detecting peaks.
  ! Increments of omega between alpha steps should be constant
  ! omega must increase monothonically
  ! no checks performed!
  ! OBM 2010
  !
  IMPLICIT NONE
  !Input and output
  REAL(kind=dp),INTENT(in) :: omeg, alpha !x and y
  !
  ! Local variables
  !
  REAL(kind=dp),SAVE :: omeg_save = 0.0d0, &
                      & thm1, h2m1,&
                      & first_der_save=9.0d99
  REAL(kind=dp),SAVE :: alpha_save(3) = 0.0d0
  INTEGER, SAVE :: current_iter = 0
  LOGICAL, SAVE :: trigger=.true.
  REAL(kind=dp) :: first_der, second_der
  !
  is_peak = .false.
  ! counter
  ! Rotate the variables
  !
  IF (current_iter < 3) THEN
      current_iter = current_iter + 1
      omeg_save = omeg
      alpha_save(current_iter) = alpha
      RETURN
  ELSE
      IF (current_iter == 3) THEN
         current_iter = current_iter + 1
         thm1=(omeg-omeg_save)
         h2m1=1.0d0/(thm1*thm1) !for second derivative
         thm1=0.5d0/thm1        !for first derivative
         !thm1=0.083333333333333d0/thm1        !for first derivative
      ENDIF
      !alpha_save(1)=alpha_save(2) !t-2h
      !alpha_save(2)=alpha_save(3) !t-h
      !alpha_save(3)=alpha_save(4) !t
      !alpha_save(4)=alpha_save(5) !t+h
      !alpha_save(5)=alpha         !t+2h
      alpha_save(1)=alpha_save(2)  !t-h
      alpha_save(2)=alpha_save(3)  !t
      alpha_save(3)=alpha          !t+h
  ENDIF
  !
  !The derivatives
  first_der = (alpha_save(3)-alpha_save(1))*thm1
  second_der = (alpha_save(3)-2.0d0*alpha_save(2)+alpha_save(1))*h2m1
  ! second derivative corresponds to t, 3 steps before
  !first_der = (-alpha_save(5)+8.0d0*(alpha_save(4)-alpha_save(2))+alpha_save(1))*thm1
  !first derivative corresponds to t, 3 steps before
  !second_der = (alpha_save(4)-2.0d0*alpha_save(3)+alpha_save(2))*h2m1
  ! second derivative corresponds to t, 3 steps before
  !Decide
  !print *,"w",omeg-0.25d0/thm1,"f=",abs(first_der),"s=",second_der
  !print *,"w",omeg-0.5d0/thm1,"f=",abs(first_der),"s=",second_der
  !if (abs(first_der) < 1.0d-8 .and. second_der < 0 ) is_peak=.true.
  !
  IF (second_der < 0) THEN
     IF (trigger) THEN
        IF (abs(first_der) <abs(first_der_save)) THEN
           first_der_save = first_der
           RETURN
        ELSE
           is_peak=.true.
           trigger=.false.
           RETURN
        ENDIF
     ENDIF
  ELSE
     first_der_save=9.0d99
     trigger=.true.
  ENDIF
  !
  RETURN
  !
END FUNCTION is_peak

REAL(kind=dp) FUNCTION integrator(dh,alpha)
  !------------------------------------------------------------------------
  !
  ! This function calculates an integral every
  ! three points, using the Simpson's rule.
  !
  IMPLICIT NONE
  !Input and output
  REAL(kind=dp),INTENT(in) :: dh, alpha !x and y
  !internal
  LOGICAL,SAVE :: flag=.true.
  !
  ! COMPOSITE SIMPSON INTEGRATOR, (precision level ~ float)
  ! \int a b f(x) dx = ~ h/3 (f(a) + \sum_odd-n 2*f(a+n*h) + \sum_even-n 4*f(a+n*h) +f(b))
  !
  integrator = 0.0d0
  !
  IF (flag) THEN
     ! odd steps
     integrator = (4.0d0/3.0d0)*dh*alpha
     flag = .false.
  ELSE
     ! even steps
     integrator = (2.0d0/3.0d0)*dh*alpha
     flag = .true.
  ENDIF
  !
  RETURN
  !
END FUNCTION integrator

SUBROUTINE read_b_g_z_file_html()
  !------------------------------------------------------------------------
  !
  ! This subroutine reads the coefficients from the html file.
  !
  USE iotk_module
  !
  IMPLICIT NONE
  !
  INTEGER :: iun, ierr, is
  INTEGER :: nipol_input_tmp, nspin
  CHARACTER(LEN=3) :: ipol_label_tmp
  INTEGER :: n_lzstep_tmp
  REAL(DP), ALLOCATABLE :: beta_store_tmp(:,:)
  REAL(DP), ALLOCATABLE :: gamma_store_tmp(:,:)
  COMPLEX(DP), ALLOCATABLE :: zeta_store_tmp(:,:,:)
  CHARACTER(LEN=4)   :: my_ip
  CHARACTER(LEN=256) :: file_ip
  !
  IF (sym_op == 0) THEN
     !
     DO ip = 1, n_ipol
        !
        WRITE(my_ip,'(i1)')ip
        file_ip = 'summary.'//TRIM(my_ip)//'.xml'
        !
        CALL iotk_free_unit( iun, ierr )
        CALL iotk_open_read( iun, FILE = TRIM( file_ip ), IERR = ierr )
        !
        CALL iotk_scan_begin( iun, "SUMMARY")
        CALL iotk_scan_dat( iun, "nspin", nspin )
        CALL iotk_scan_dat( iun, "nipol_input", nipol_input_tmp )
        CALL iotk_scan_dat( iun, "ipol_label",  ipol_label_tmp )
        CALL iotk_scan_dat( iun, "n_lanczos",    n_lzstep_tmp )
        CALL iotk_scan_end( iun, "SUMMARY"  )
        !
        WRITE (stdout,*) ' '
        WRITE (stdout,*) ' Reading alpha beta zeta of the polarzation : ', ipol_label_tmp, &
                          ' from file : ', file_ip
        !
        IF (n_lzstep_tmp < itermax0) THEN
           CALL errore("read_b_g_z_file", "Error in lriter_stop < itermax0, reduce itermax0",1)
        ENDIF
        !
        ALLOCATE(beta_store_tmp (n_lzstep_tmp,nspin))
        ALLOCATE(gamma_store_tmp(n_lzstep_tmp,nspin))
        !
        CALL iotk_scan_begin( iun, "BETA_STORE")
        DO is = 1, nspin
           CALL iotk_scan_dat( iun, "beta_store_k", beta_store_tmp(1:n_lzstep_tmp,is) )
        ENDDO
        CALL iotk_scan_end( iun, "BETA_STORE"  )
        !
        CALL iotk_scan_begin( iun, "GAMMA_STORE")
        DO is = 1, nspin
           CALL iotk_scan_dat( iun, "gamma_store_k", gamma_store_tmp(1:n_lzstep_tmp,is) )
        ENDDO
        CALL iotk_scan_end( iun, "GAMMA_STORE"  )
        !
        ALLOCATE(zeta_store_tmp (3,n_lzstep_tmp,nspin))
        !
        CALL iotk_scan_begin( iun, "ZETA_STORE")
        DO ip2 = 1, 3
           CALL iotk_scan_dat( iun, "zeta_store_ipol_j", ip2 )
           DO is = 1, nspin
              CALL iotk_scan_dat( iun, "zeta_store_k", zeta_store_tmp(ip2,1:n_lzstep_tmp,is) )
           ENDDO
        ENDDO
        CALL iotk_scan_end( iun, "ZETA_STORE"  )
        !
        CALL iotk_close_read( iun )
        !
        IF (nspin == 1) spin_channel = 1
        !
        norm0(ip)                    =  beta_store_tmp(1, spin_channel)
        beta_store(ip,1:itermax0-1)  =  beta_store_tmp(2:itermax0, spin_channel)
        gamma_store(ip,1:itermax0-1) =  gamma_store_tmp(2:itermax0, spin_channel)
        beta_store(ip,itermax0)      =  beta_store_tmp(itermax0, spin_channel)
        gamma_store(ip,itermax0)     =  gamma_store_tmp(itermax0, spin_channel)
        !
        IF (n_ipol == 1) THEN
           !
           zeta_store(1,1,1:itermax0)=  zeta_store_tmp(ipol,1:itermax0,spin_channel)
           !
        ELSE
           !
           zeta_store(ip,:,1:itermax0) = zeta_store_tmp(:,1:itermax0,spin_channel)
           !
        ENDIF
        !
        DEALLOCATE(beta_store_tmp)
        DEALLOCATE(gamma_store_tmp)
        DEALLOCATE(zeta_store_tmp)
        !
     ENDDO
     !
     IF (nspin == 1) degspin  = 2
     IF (nspin == 2) degspin  = 1
     !
     beta_store(:,itermax0+1:itermax)   = 0.d0
     gamma_store(:,itermax0+1:itermax)  = 0.d0
     zeta_store(:,:,itermax0+1:itermax) = (0.d0,0.d0)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE read_b_g_z_file_html


SUBROUTINE extrapolate()
  !
  ! This subroutine applies the "extrapolation" scheme
  ! for extrapolating the reduced matrix.
  !
  IMPLICIT NONE
  !
  !  Terminatore
  !
  skip = .false.
  !
  IF (trim(extrapolation)/="no") THEN
   !
   average = 0.d0
   av_amplitude = 0.d0
   !
   DO ip=1,n_ipol
     !
     WRITE(stdout,'(/5x,"Polarization direction:",I1)') ip
     counter=0
     !
     DO i=151,itermax0
        !
        IF (skip .eqv. .true.) THEN
           skip=.false.
           CYCLE
        ENDIF
        !
        IF (mod(i,2)==1) THEN
           !
           IF ( i/=151 .and. abs( beta_store(ip,i)-average(ip)/counter ) > 2.d0 ) THEN
              !
              !if ( i.ne.151 .and. counter == 0) counter = 1
              skip=.true.
              !
           ELSE
              !
              average(ip)=average(ip)+beta_store(ip,i)
              av_amplitude(ip)=av_amplitude(ip)+beta_store(ip,i)
              counter=counter+1
              !
           ENDIF
           !
        ELSE
           !
           IF ( i/=151 .and. abs( beta_store(ip,i)-average(ip)/counter ) > 2.d0 ) THEN
              !
              skip=.true.
              !
           ELSE
              !
              average(ip)=average(ip)+beta_store(ip,i)
              av_amplitude(ip)=av_amplitude(ip)-beta_store(ip,i)
              counter=counter+1
              !
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
     average(ip)=average(ip)/counter
     av_amplitude(ip)=av_amplitude(ip)/counter
     !
     WRITE(stdout,'(5x,"Lanczos coefficients:")')
     WRITE(stdout,'(5x,"Average =",3F15.8)') average(ip)
     WRITE(stdout,'(5x,"Average oscillation amplitude =",F15.8)') av_amplitude(ip)
     !
   ENDDO
   !
   IF (trim(extrapolation)=="constant") av_amplitude=0
   !
   !
   DO ip=1,n_ipol
     !
     DO i=itermax0,itermax
        !
        IF (mod(i,2)==1) THEN
           !
           beta_store(ip,i)=average(ip)+av_amplitude(ip)
           gamma_store(ip,i)=average(ip)+av_amplitude(ip)
           !
        ELSE
           !
           beta_store(ip,i)=average(ip)-av_amplitude(ip)
           gamma_store(ip,i)=average(ip)-av_amplitude(ip)
           !
        ENDIF
        !
     ENDDO
     !
   ENDDO
   !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE extrapolate

SUBROUTINE calc_chi(freq,broad,chi)
  !-----------------------------------------------------------------------------
  !
  ! This subroutine Calculates the susceptibility.
  !
  IMPLICIT NONE
  !
  REAL(kind=dp), INTENT(in) :: freq
  REAL(kind=dp), INTENT(in) :: broad
  COMPLEX(kind=dp), INTENT(out) :: chi(:,:)
  !
  omeg_c = cmplx(freq,broad,dp)
  !
  DO ip =1, n_ipol
     !
     a(:) = omeg_c
     !
     DO i = 1,itermax-1
        !
        b(i) = cmplx(-beta_store(ip,i),0.0d0,dp)
        c(i) = cmplx(-gamma_store(ip,i),0.0d0,dp)
        !
     ENDDO
     !
     r(ip,:) = (0.0d0,0.0d0)
     r(ip,1) = (1.0d0,0.0d0)
     !
     ! |w_t|=(w-L) |1,0,0,...,0|
     !
     CALL zgtsv(itermax,1,b,a,c,r(ip,:),itermax,info)
     !
     IF (info /= 0) CALL errore("calc_chi", "Unable to solve tridiagonal system",1)
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
         chi(ip,ip2) = chi(ip,ip2) * cmplx(norm0(ip),0.0d0,dp)
         !
         ! The response charge density is defined as 2*evc0*q, see Eq. (43) in
         ! JCP 128, 154105 (2008).
         ! Therefore, the dipole is given by 2*degspin* zeta^T *
         ! (w-T^itermax)^-1 * e_1. See also Eq. (15) in that paper.
         ! Optics: The minus sign accounts for the negative electron charge
         ! (perturbation is -e E x, rather than E x)
         !
         chi(ip,ip2) = chi(ip,ip2) * cmplx(-2.d0*degspin, 0.d0, dp)
         !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE calc_chi

SUBROUTINE wl_to_color(wavelength,red,green,blue)
  !----------------------------------------------------------------------------
  !
  ! Gives the colour intensity of a given wavelength
  ! in terms of RGB (red, green and blue).
  !
  IMPLICIT NONE
  !
  REAL(kind=dp), INTENT(in) :: wavelength
  REAL(kind=dp), INTENT(out) :: red,green,blue
  !
  IF ((wavelength>=380.).and.(wavelength<=440.)) THEN
     red = -1.*(wavelength-440.)/(440.-380.)
     green = 0.
     blue = 1.
  ENDIF
  !
  IF ((wavelength>=440.).and.(wavelength<=490.)) THEN
     red = 0.
     green = (wavelength-440.)/(490.-440.)
     blue = 1.
  ENDIF
  !
  IF ((wavelength>=490.).and.(wavelength<=510.)) THEN
     red = 0.
     green = 1.
     blue = -1.*(wavelength-510.)/(510.-490.)
  ENDIF
  !
  IF ((wavelength>=510.).and.(wavelength<=580.)) THEN
     red = (wavelength-510.)/(580.-510.)
     green = 1.
     blue = 0.
  ENDIF
  !
  IF ((wavelength>=580.).and.(wavelength<=645.)) THEN
     red = 1.
     green = -1.*(wavelength-645.)/(645.-580.)
     blue = 0.
  ENDIF
  !
  IF ((wavelength>=645.).and.(wavelength<=780.)) THEN
     red = 1.
     green = 0.
     blue = 0.
  ENDIF
  !
  RETURN
  !
END SUBROUTINE wl_to_color

  !
END SUBROUTINE wbsepp_ads_spectrum

