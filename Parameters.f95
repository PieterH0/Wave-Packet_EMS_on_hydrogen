module Parameters
    use basics, only : dp, pi
    implicit none

    integer, parameter, private :: wp = dp

    !public delay_time

    real(wp) :: delay_time      ! a.u. Delay time for shifting the projectile pulse beam focus!
    ! --- Loading Parameter --- !
    !!! Default Arguments !!!
    real(wp) :: theta = pi/4._wp ! final polar angle of measurement (45 degr)
    real(wp) :: k0    =  27.11132398604629_wp  ! centra k-value

    ! phi_f integral
    integer  :: Nphi   = 100
    real(wp) :: minPhi = -0.02_wp
    real(wp) :: maxPhi =  0.02_wp

    ! k_f integral (depends on pulse simulated)
    integer  :: Nkpoints = 3500
    real(wp) :: mink = 19.1_wp
    real(wp) :: maxk = 19.2_wp
    logical  :: multiple = .false. ! whether or not to use multiple k_f integrals
    integer  :: Nint = 0           ! If multiple = .true. the number of side bands larger than zero, so 2*Nint + 1 in total.

    real(wp) :: time = 0.0    !(time parameter for target evolution t_c)

    integer  :: Nthetai_1 = 100 !(quadrature points internal theta integral for kn close to band center)
    integer  :: Nthetai_2 = 100 !(quadrature points internal theta integral for kn far from band center)
                                ! Nthetai_1 is one used for unmodulated pulses
    real(wp) :: thetai_1_max = 0.004 ! maximum theta value for integration Gaussian pulses
    integer  :: Nphii     =  30 !(quadrature points internal phi integral)

    ! impact parameter
    real(wp) :: bx = 0.0_wp 
    real(wp) :: by = 0.0_wp

    !!!!!!!!!      Variables to the projectile pulse       !!!!!!!!
    character(len=8) :: pulse_type = 'gaussian'
    logical  :: modulated  = .false.
    integer  :: g2 = 0
    real(wp) :: omega = 0._wp
    real(wp) :: dk = 0._wp                                   ! a.u. Default. True value calculated in main2
    real(wp) :: tbunch = 0._wp                               ! a.u. This is the default value. In the main2 the actual value is calculated
    real(wp) :: sigmatrans = 0.02711132398604629_wp          ! a.u. (1 mrad * k0)
    real(wp) :: sigmalong  = 0.01050490402250473_wp          ! a.u. (0.1 fs)
    real(wp) :: sigmatheta = 0.01_wp                         ! rad

    !!! Picking the correct pulse type !!! 
    ! need an interface containing function of same shape as the one we point to 
    abstract interface
        function some_function(x,y) result(z)
            implicit none
            real(8), intent(in) :: x,y
            complex(8) :: z
        end function some_function
        ! real(8) is to define with double precision. Somehow cannot get wp inside interface. 
    end interface

    ! defining the pointer 
    procedure(some_function), pointer :: pulse_to_use => null()

end module Parameters