module Tmatrix
    ! other modules used:
    use Basics
    use SpecialFunctions
    use ModulatedPulses
    use FourierTransforms
    use IntegrationModPulses
    use Parameters, only : k0, omega, modulated, dk, Nthetai_1, Nthetai_2, Nphii, pulse_type, thetai_1_max

    implicit none
    private

    ! Define public routines
    public :: calculalte_T_matrix_element, integrand, theta_curve
    
    ! Define working precision
    integer, parameter :: wp = dp


    ! a comment
contains

    ! Function to calculate element ka,kb point
    function calculalte_T_matrix_element(ka, kb, N, ns, ls, ms, time, bx, by, symmetry_factor) result(Tfinal)
        ! arguments by
        real(wp), dimension(3), intent(in) :: ka, kb                   ! outgoing momenta
        real(wp), intent(in) :: time                                   ! time 
        real(wp), intent(in) :: bx, by                                 ! impact parameter         
        integer, intent(in)  :: N                                      ! Number of states in superposition
        integer, intent(in), dimension(N)  :: ns, ls, ms       ! target state qunatum numbers
        !integer, intent(in)  :: Ntheta, Nphi                  ! number of quadrature points
        integer, intent(in)  :: symmetry_factor                ! whether symmetric or antisymmetric spatial wave function

        
        ! variables
        real(wp) :: kn, frontfactor, theta0
        integer  :: h
        real(wp) :: phi_points(Nphii)
        real(wp) :: theta_points(Nthetai_1), theta_points_mod(Nthetai_1, 3) ! arrays for theta_points
        ! since only three intervals needed for modulated pulses we just use the extra memeory. This
        ! also removes any trouble that might come with using allocation, which I suspect might also be time costly
        !
        integer :: Nintervals ! the number of theta intevals to integrate over (2 or 3)
        integer, dimension(3) :: Nthetas ! the number of quadrature points in each of the theta_i integration intervals
        integer :: Nthetai  ! takes the value of Nthetai_1 or Nthetai_2 depending onthe if statement in integration non-modulated Morimoto pulses

        ! specific for the case that pulse_type = modulPRA
        integer  :: bandnumber
        real(wp) :: d_k
        real(wp) :: theta_points2(Nthetai_2)

        ! output
        complex(wp) :: T, Tfinal
        T = cmplx(0._wp, 0._wp, wp) 

        ! the front factor
        ! first factor from plane waves, second a 4pi from Bethe integral and a 4pi from Bessel expansion used in Fourier transf.
        ! finally 1/sqrt(N) is from normalisation. This assumes equal population
        frontfactor = 1._wp/(2._wp*pi)**(9./2.) * (4._wp*pi)**2 *1._wp/sqrt(dble(N)) ! times (-i^l)

        ! integration points for phi
        call linspace(phi_points, 0._wp, 2._wp*pi, Nphii)
        call linspace(theta_points, 0._wp, thetai_1_max, Nthetai_1)
 
        ! loop over states in superposition
        do h=1, N
            ! calculate kn
            kn = sqrt(norm2(ka)**2 + norm2(kb)**2 + 1._wp/((ns(h))**2)) ! evaluating the energy delta function
            
            if (pulse_type .eq. 'modulatd') then
                ! If the type equals the very short(in k-space) pulse from Morimoto et al. 
                ! Here we can integrale theta only along a thin curve
                !                
                ! if pulse is modulated
                if (modulated .eqv. .true.) then 
                    !
                    call find_integration_limits(kn, Nthetas, Nintervals, theta_points_mod) ! creates the theta intervals and sets Nintervals
                    ! now we integrate. This modifies T 
                    call integrate_over_theta_phi_modulated(T, Nintervals, Nthetas, Nphii, theta_points_mod, &
                                                            & phi_points, kn, ka, kb, ns(h), ls(h), ms(h), &
                                                            & time, bx, by, symmetry_factor)
                    !
                    !
                ! if the pulse is unmodulated (simply a gaussian) but very thin! 
                else if (modulated .eqv. .false.) then
                    !
                    ! first find the nearest 
                    if (kn .ge. (k0 + 0.0002_wp)) then    ! if kn larger than k0. 
                        theta0 = theta_curve(kn, k0)
                        !call linspace(theta_points, max(0._wp, theta0-0.0005_wp), max(0.004_wp, theta0+0.0005_wp), Nthetai_2) ! uses max to avoid negative thetas
                        call linspace(theta_points, theta0-0.0012_wp, theta0+0.0012_wp, Nthetai_2)
                        Nthetai = Nthetai_2 ! many points needed
                    else 
                        call linspace(theta_points, 0._wp, 0.0035_wp, Nthetai_1)
                        Nthetai = Nthetai_1 ! few points needed
                    end if
                    !
                    ! now call integration
                    call integrate_over_theta_phi(T, Nthetai, Nphii, theta_points, phi_points, kn, ka, kb, ns(h), ls(h), ms(h), &
                                                                            & time, bx, by, symmetry_factor)
                    !
                end if
                !
            ! If now modulated but without interferece in kn's, modulated with sigma_trans = k0 * 1mrad
            else if (pulse_type .eq. 'modulPRA') then
                ! want to set up different theta_i integrals depending on kn, but only a single interval for each kn
                ! find band
                bandnumber = nint((kn - k0)/dk) ! nint round to nearest whole number
                d_k = kn - k0 - bandnumber*dk   ! distance from nearest band center

                if (d_k .le. -0.0001_wp) then       ! Smaller than central momentum
                    theta0 = theta_curve2(kn, N-1)
                    call linspace(theta_points2, theta0 - 0.0009_wp, theta0 + 0.0009_wp, Nthetai_2)
                    call integrate_over_theta_phi(T, Nthetai_2, Nphii, theta_points, phi_points, kn, ka, kb, ns(h), ls(h), ms(h), &
                    & time, bx, by, symmetry_factor)
                !

                else if (d_k .ge. 0.0001_wp) then   ! greater than central momentum
                    theta0 = theta_curve2(kn, N)
                    call linspace(theta_points2, theta0 - 0.0009_wp, theta0 + 0.0009_wp, Nthetai_2)
                    call integrate_over_theta_phi(T, Nthetai_2, Nphii, theta_points, phi_points, kn, ka, kb, ns(h), ls(h), ms(h), &
                    & time, bx, by, symmetry_factor)
                !
                else ! the case there kn is close to the band center momentum.
                    call integrate_over_theta_phi(T, Nthetai_1, Nphii, theta_points, phi_points, kn, ka, kb, ns(h), ls(h), ms(h), &
                    & time, bx, by, symmetry_factor)
                !
                !
                end if

            else 
                ! if the pulse is long in k-space such that there is no curve of central theta points to integrate around
                ! this is the case for the 100 as Shao and Starace pulse
                ! simply create a long theta range
                ! call linspace(theta_points, 0._wp, 0.004_wp, Nthetai_1)
                call integrate_over_theta_phi(T, Nthetai_1, Nphii, theta_points, phi_points, kn, ka, kb, ns(h), ls(h), ms(h), &
                                                                            & time, bx, by, symmetry_factor)
                !
            end if
            !
        end do ! end do over state in coherent superposition (h loop)

        ! finally we add the frontfactor
        Tfinal = frontfactor * T
    end function calculalte_T_matrix_element


    ! The integrand we are integrating over 
    function integrand(kn, ka, kb, n, l, m, theta, phi, bx, by, symmetry_factor) result(I)
        use Parameters, only : delay_time, pulse_to_use
        ! arguments
        real(wp), intent(in) :: kn, ka(3), kb(3), theta, phi, bx, by
        integer, intent(in)  :: n, l, m
        integer, intent(in)  :: symmetry_factor  ! Whether the spatial wave function is symmetric or not (=1 for symmetric, =-1 for anti symmetric) 
        
        ! variables
        real(wp) :: q(3), cos_theta_q, phi_q      ! concerning the target momentum
        real(wp) :: delta_squared_a, delta_squared_b, delta_factor, radial_of_q ! factors in contained in integrand (Hydrogen)
        complex(wp) :: sph_harm_factor, psi
        
        complex(wp) :: I ! the resulting integrand

        ! pulse factor
        !psi = gaussian_pulse(kn, theta)        ! single Gaussian
        !psi = energy_gaussian_pulse(kn, theta) ! single Gaussian in energy
        psi = pulse_to_use(kn, theta)           ! Modulated pulse

        ! change in momentum incoming electron
        ! calculated for both outgoing electrons due to symmetrization of fermions.
        delta_squared_a = (ka(1) - kn*sin(theta)*cos(phi))**2  +  (ka(2) - kn*sin(theta)*sin(phi))**2  +  (ka(3) - kn*cos(theta))**2 
        delta_squared_b = (kb(1) - kn*sin(theta)*cos(phi))**2  +  (kb(2) - kn*sin(theta)*sin(phi))**2  +  (kb(3) - kn*cos(theta))**2
        delta_factor    = 1._wp/delta_squared_a + symmetry_factor*1_wp/delta_squared_b

        ! calculate the target momentum q
        q = kn * [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)] - ka - kb
        cos_theta_q = q(3)/norm2(q)
        phi_q       = sign(1._wp, q(2)) * acos(possible_zero_divided_by_zero(q(1), norm2(q(:2))))

        ! The radial integral evaluated at |q|
        radial_of_q = fourier_hydrogen_orb(norm2(q), n, l) ! these integrals are evaluated and hardcoded in module "FourierTransforms"

        ! spherical harmonic
        sph_harm_factor = py_sph_harm(cos_theta_q, phi_q) ! sph_harm(cos_theta_q, phi_q, l, m)
        
        ! The integrand
        ! contains : psi(k_i), sine from spherical jacobian, radial fourier, 1/delta^2 from Bethe integral
        !            spherical harmonic from spherical fourier, phase from impact parameter
        !            and phase from time delay (z-impact parameter)
        I = psi * sin(theta) * radial_of_q * cmplx(delta_factor, 0._wp, wp) * sph_harm_factor &
                    & * exp(cmplx(0._wp, -kn*sin(theta)*(bx*cos(phi) + by*sin(phi)), wp)) !&
                    !& * exp(cmplx(0._wp, (kn)**2/2._wp*delay_time - kn*cos(theta)*k0*delay_time, wp))  
                    !& * exp(cmplx(0._wp, (kn*sin(theta))**2/2*delay_time, wp))
                    !& * exp(cmplx(0._wp, - kn*cos(theta)*k0*delay_time, wp))
                    
                    !& * exp(cmplx(0._wp, (kn)**2/2._wp*delay_time - kn*cos(theta)*k0*delay_time, wp))     
                    !& * exp(cmplx(0._wp, -kn*k0*cos(theta)*(-2._wp*pi*288._wp/7._wp)/4._wp * dble(m), wp))
                    !
                    !& * exp(cmplx(0._wp, (kn*sin(theta))**2/2*delay_time, wp))                        ! transvers delay
                    !& * exp(cmplx(0._wp, (kn)**2/2._wp*delay_time - kn*k0*cos(theta)*delay_time, wp)) ! full delay
        !
    end function integrand

    function theta_curve(kn, k0) result(theta0)
        ! function to calc. the central theta value to integrate around.
        ! will need to be updated for modulated pulses
        real(wp), intent(in) :: kn, k0
        real(wp) :: theta0

        theta0 = sqrt(kn**2 - k0**2)/kn
    end function theta_curve

    
end module Tmatrix