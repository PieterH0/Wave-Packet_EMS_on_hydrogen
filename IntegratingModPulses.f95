!
module IntegrationModPulses
    ! Module performs integrations over appropriate phi_i and theta_i intervals
    use Basics
    use Tmatrix, only : integrand
    use Parameters, only : k0, omega, dk, Nthetai_1, Nthetai_2
    use ModulatedPulses  
    
    implicit none
    private

    public find_integration_limits, integrate_over_theta_phi_modulated, integrate_over_theta_phi, set_up_kf_integrals, theta_curve2

    ! Define working precision
    integer, parameter :: wp = dp 


    contains
    subroutine integrate_over_theta_phi_modulated(T, Nintervals, Nthetas, Nphi, theta_points, phi_points, kn, ka, kb, &
                                & n, l, m, time, bx, by, symmetry_factor)
        
        ! arguments
        complex(wp), intent(inout) :: T     ! the T-matrix up to constant. What is returned.
        integer, intent(in) :: Nintervals   ! number of integration intervals over theta. Output from find_integration_limits
        integer, intent(in) :: Nphi         ! number of quadrature points phi integral
        integer, dimension(3), intent(in) :: Nthetas  ! number of quadrature points in each integration interval over theta_i
        real(wp), dimension(Nthetai_1, 3), intent(in) :: theta_points
        real(wp), dimension(Nphi), intent(in) :: phi_points
        integer, intent(in)  :: n, l, m     ! target state
        integer, intent(in)  :: symmetry_factor  ! Whether the spatial wave function is symmetric or not (=1 for symmetric, =-1 for anti symmetric) 
        real(wp), intent(in) :: kn, ka(3), kb(3) 
        real(wp), intent(in) :: time        ! time in phases coherence
        real(wp), intent(in) :: bx, by      ! impact parameters
        
        ! variables
        integer     :: i,j,k  ! for looping
        complex(wp) :: Tn ! T-matrix contribution of n'th target state
        real(wp)    :: dtheta, dphi ! integration steps 


        do i=1, Nintervals ! loop over theta intervals
            dtheta = theta_points(2,i) - theta_points(1,i)
            dphi   = phi_points(2) - phi_points(1)
            Tn     = cmplx(0._wp, 0._wp, wp)
            ! loops over theta and phi
            do j=1, Nthetas(i)
                do k=1, Nphi
                    !
                    Tn = Tn + integrand(kn, ka, kb, n, l, m, theta_points(j,i), phi_points(k), bx, by, symmetry_factor)
                    !
                end do
            end do
            ! Multiply remaining factors
            T = T + Tn * kn * dtheta * dphi * exp(cmplx(0._wp, -0.5_wp/(n**2) * time, wp))
        end do

    end subroutine integrate_over_theta_phi_modulated


    subroutine integrate_over_theta_phi(T, Ntheta, Nphi, theta_points, phi_points, kn, ka, kb, n, l, m, time, &
                                                                    & bx, by, symmetry_factor)
        ! integrating over theta and phi of the incoming momenta if the pulse is gaussian non-modulated. T
        ! non-modulated means find_integration_limits is not called beforehand. 
        ! arguments
        complex(wp), intent(inout) :: T     ! the T-matrix up to constant. What is returned.
        integer, intent(in) :: Ntheta, Nphi ! integration points ntheta for each interval)
        real(wp), dimension(Ntheta), intent(in) :: theta_points
        real(wp), dimension(Nphi), intent(in)   :: phi_points
        integer, intent(in)  :: n, l, m          ! target state
        integer, intent(in)  :: symmetry_factor  ! Whether the spatial wave function is symmetric or not (=1 for symmetric, =-1 for anti symmetric) 
        real(wp), intent(in) :: kn, ka(3), kb(3)
        real(wp), intent(in) :: time, bx, by ! delay time, impact parameter
        
        ! variables
        integer :: i,j    ! for looping
        complex(wp) :: Tn ! T-matrix contribution
        real(wp)    :: dtheta, dphi ! integration steps 


        dtheta = theta_points(2) - theta_points(1)
        dphi   = phi_points(2) - phi_points(1)

        Tn = cmplx(0._wp, 0._wp, wp)
        ! loops over theta and phi
        do i=1, Ntheta
            do j=1, Nphi
                !
                Tn = Tn + integrand(kn, ka, kb, n, l, m, theta_points(i), phi_points(j), bx, by, symmetry_factor)
                !
            end do
            ! correction to get proper trapeziodal (N = Nphi + 1) integration as end points are not necessarily zero.
            Tn = Tn - (integrand(kn, ka, kb, n, l, m, theta_points(i), phi_points(1), bx, by, symmetry_factor))/2._wp
            Tn = Tn - (integrand(kn, ka, kb, n, l, m, theta_points(i), phi_points(Nphi), bx, by, symmetry_factor))/2._wp
        end do
        ! Multiply remaining factors
        T = T + Tn * kn * dtheta * dphi * exp(cmplx(0._wp, -0.5_wp/(n**2) * time, wp))

    end subroutine integrate_over_theta_phi


    subroutine find_integration_limits(kn, Nthetas, Nintervals, theta_points)
        ! Subroutine to create the arrays of theta points to be evaluate in the theta integral over the 
        ! incoming momenta. 
        ! arguments
        real(wp), intent(in)   :: kn
        integer,  intent(out)  :: Nintervals  ! number of theta intervals to integrate over, either 2 or 3.
        real(wp), intent(out)  :: theta_points(Nthetai_1, 3)  
        integer,  intent(out)  :: Nthetas(3)   ! number of theta points in each integration interval, either [Nthetai_1, Nthetai_2, Nthetai_2] or [Nthetai_2, Nthetai_2, Nthetai_2]
        
        integer :: N      ! number of nearest band for kn
        real(wp), dimension(Nthetai_1) :: interval_points ! quadrature points in each interval of integration.
        real(wp) :: d_k   ! distance from nearest band center
        real(wp) :: theta0, theta1, theta2  ! Theta points to integrate around 

        N = nint((kn - k0)/dk) ! nint round to nearest whole number
        ! distance from nearest band center
        d_k = kn - k0 - N*dk

        if (d_k .le. -0.0001_wp) then       ! Smaller than central momentum
            Nintervals = 2
            Nthetas = [Nthetai_2, Nthetai_2, Nthetai_2]

            theta0 = theta_curve2(kn, N-1)
            theta1 = theta_curve2(kn, N-2)
            call linspace(interval_points, theta0 - 0.0009_wp, theta0 + 0.0009_wp, Nthetai_2)
            theta_points(:,1) = interval_points
            call linspace(interval_points, theta1 - 0.0009_wp, theta1 + 0.0009_wp, Nthetai_2)
            theta_points(:,2) = interval_points

            !
            ! last column is arbitrary 
        !
        else if (d_k .ge. 0.0001_wp) then   ! greater than central momentum
            Nintervals = 3
            Nthetas = [Nthetai_2, Nthetai_2, Nthetai_2]

            theta0 = theta_curve2(kn, N)
            theta1 = theta_curve2(kn, N-1)
            theta2 = theta_curve2(kn, N-2)

            call linspace(interval_points, theta0 - 0.0009_wp, theta0 + 0.0009_wp, Nthetai_2)
            theta_points(:,1) = interval_points
            call linspace(interval_points, theta1 - 0.0009_wp, theta1 + 0.0009_wp, Nthetai_2)
            theta_points(:,2) = interval_points
            call linspace(interval_points, theta2 - 0.0009_wp, theta2 + 0.0009_wp, Nthetai_2)
            theta_points(:,3) = interval_points
        !
        else                ! the case there kn is close to the band center momentum.
            Nintervals = 3
            Nthetas = [Nthetai_1, Nthetai_2, Nthetai_2]

            theta1 = theta_curve2(kn, N-1)
            theta2 = theta_curve2(kn, N-2)

            call linspace(interval_points, 0._wp, 0.004_wp, Nthetai_1)
            theta_points(:,1) = interval_points
            call linspace(interval_points, theta1 - 0.0009_wp, theta1 + 0.0009_wp, Nthetai_2)
            theta_points(:,2) = interval_points
            call linspace(interval_points, theta2 - 0.0009_wp, theta2 + 0.0009_wp, Nthetai_2)
            theta_points(:,3) = interval_points
        !
        end if

    end subroutine find_integration_limits

    function theta_curve2(kn, N) result(theta0)
        ! function to calc. the central theta value to integrate around
        ! this time for modulated pulses
        ! will need to be updated for modulated pulses
        real(wp), intent(in) :: kn
        integer, intent(in)  :: N 
        real(wp) :: theta0

        theta0 = sqrt(kn**2 - (k0+N*dk)**2)/kn
    end function theta_curve2


    subroutine set_up_kf_integrals(Nint, Nkpoints, min, max, kfpoints)
        ! Subroutine to set up the integration interval for a modulated electron pulse
        !
        ! Arguments:
        !   Nint : (int) number of intervals above zero
        !   Nkpoints : (int) number of k points in each integral
        !   min : (float) differnce between minimum point of integration and pulse band center (k0 + N dk)
        !   max : (float) differnce between maximum point of integration and pulse band center (k0 + N dk)

        ! Output:
        ! kfpoints : matrix (Nk)
        
        integer, intent(in)  :: Nint, Nkpoints
        real(wp), intent(in) :: min, max

        real(wp), dimension(Nkpoints), intent(out) :: kfpoints ! output matrix for quadrature points
        !real(wp) :: mink_f, maxk_f
        real(wp), allocatable :: kpoints_subinterval(:) ! kf points in a subinterval of the integral
        integer :: Nkf_in_subinterval 
        integer :: N, i, j

        Nkf_in_subinterval = Nkpoints/(2*Nint+1)
        allocate(kpoints_subinterval(Nkf_in_subinterval))


        dk = omega/k0
        ! loop over bands from -Nint to Nint
        ! do this by introducing i = 1 to 2Nint+1
        do i=1, 2*Nint+1
            N = -Nint - 1 + i ! band order
            j = i - 1
            !
            call linspace(kpoints_subinterval, sqrt((k0+(N-1)*dk - min)**2/2._wp - 1._wp/32._wp), &
                    & sqrt((k0 + (N)*dk + max)**2/2._wp - 1._wp/18._wp), Nkf_in_subinterval)
            ! fill the list of k-points
            kfpoints(j*Nkf_in_subinterval + 1 : i*Nkf_in_subinterval) = kpoints_subinterval
        end do

    end subroutine set_up_kf_integrals


end module IntegrationModPulses