! Module to contains all pulses to be integrated over. Both Gaussian and modulated pulses
!
module ModulatedPulses
use Basics
use Parameters, only : k0, dk, tbunch, sigmatrans, sigmalong, sigmatheta, g2
implicit none
private

! Define working precision 
integer, parameter :: wp = dp 

public :: trans_pulse, long_pulse, modulated_pulse, gaussian_pulse, energy_gaussian_pulse

    
contains
    ! Transversal Gaussian 
    function trans_pulse(ktrans) result(pulse)
        real(wp), intent(in)  :: ktrans
        real(wp) :: pulse
        real(wp) :: normfactor

        normfactor = 1._wp/(sqrt(2._wp * pi * sigmatrans**2))
        pulse = normfactor * exp(-ktrans**2/(4*sigmatrans**2))
    end function trans_pulse

    ! Longitudinal pulse (non-modulated)
    function long_pulse(klong) result(pulse)
        real(wp), intent(in)  :: klong
        real(wp) :: pulse
        real(wp) :: normfactor

        normfactor = 1._wp/((2._wp * pi * sigmalong**2)**(1./4.))

        pulse = normfactor *  exp(-(klong - k0)**2) * 1._wp/(2*pi*sigmalong**2)**(1./4.) * exp(-(klong)**2/(4._wp*sigmalong**2))
    end function long_pulse
  
    function modulated_pulse(kn, theta) result(final_pulse)
        ! arguments
        real(wp), intent(in) :: kn, theta ! length and polar angle of ki component 
        ! result
        complex(wp) :: final_pulse
        complex(wp) :: exp_factor
        integer :: N, Nmax
        
        ! pulse components
        real(wp)  :: klong, ktrans, normfactor
        klong  = kn * cos(theta)
        ktrans = kn * sin(theta)

        !Nmax = 6
        Nmax = nint((kn - k0)/dk) ! nint rounds to nearest whole number

        normfactor = 1._wp/((2._wp * pi * sigmalong**2)**(1./4.))

        exp_factor = exp(cmplx(0._wp, tbunch*k0*klong - tbunch*klong**2/2._wp, wp))

        final_pulse = normfactor * side_band(0, klong) * trans_pulse(ktrans) * exp_factor

        !final_pulse = normfactor * sum([(side_band(N, klong), N=Nmax-2, Nmax)]) * trans_pulse(ktrans) * exp_factor

    end function modulated_pulse

    function modulated_e_gauss(kn, theta) result(final_pulse)
        !!! NOT TESTED
        !!! BE CAREFUL!
        !arguments
        real(wp), intent(in) :: kn, theta ! length and polar angle of ki component 
        !result
        complex(wp) :: final_pulse

        real(wp) :: theta_gauss, frontfactor
        complex(wp) :: exp_factor
        
        theta_gauss = exp(-sin(theta)**2/(4._wp*sigmatheta**2))
        exp_factor  = exp(cmplx(0._wp, tbunch*k0*kn - tbunch*kn**2/2._wp, wp))

        final_pulse = cmplx(0._wp, 0._wp, wp) !frontfactor * sum([(side_band(N, klong), N=Nmax-2, Nmax)]) * theta_gauss

    end function modulated_e_gauss

    function gaussian_pulse(kn, theta) result(final_pulse)
        ! Gaussian pulse in that it is a single guassian in momentum. Not a gaussian in energy!
        ! (25-01-24) I removed the exponential factor which in all other calculations has been interted in the T-matrix in the T-matrix module. This 
        ! means that the exp_factor should not have been there meaning results with tp=/= 0 are wrong prior to today.

        ! arguments
        real(wp), intent(in) :: kn, theta ! length and polar angle of ki component 
        ! result
        complex(wp) :: final_pulse
        !complex(wp) :: exp_factor

        ! pulse components
        real(wp)  :: klong, ktrans, normfactor
        klong  = kn * cos(theta)
        ktrans = kn * sin(theta)



        normfactor  = 1._wp/((2._wp * pi * sigmalong**2)**(1./4.))
        !exp_factor  = exp(cmplx(0._wp, tp*k0*klong - tp*klong**2/2._wp, wp))
        final_pulse = cmplx(normfactor * exp(-(klong - k0)**2/(4._wp * sigmalong**2)) * trans_pulse(ktrans), 0._wp, wp) ! * exp_factor 

    end function gaussian_pulse

    function energy_gaussian_pulse(kn, theta) result(final_pulse)
        ! Gaussian pulse which is a gaussian in energy!
        ! arguments
        real(wp), intent(in) :: kn, theta ! length and polar angle of ki component 
        ! result
        complex(wp) :: final_pulse ! returned as complex to be consistent with 
        
        ! pulse components
        real(wp)  :: normfactor

        normfactor  = 1._wp/((8._wp * pi**3)**(1._wp/4._wp) * sqrt(sigmalong) * k0 * sigmatheta)
        final_pulse = cmplx(normfactor * exp(-(kn - k0)**2/(4._wp*sigmalong**2)) * &
                        & exp(-sin(theta)**2/(4._wp*sigmatheta**2)), 0._wp, wp)

    end function energy_gaussian_pulse

    ! Side band gaussian
    function side_band(N, klong) result(band)
        integer, intent(in) :: N
        real(wp) :: band
        real(wp) :: bandcenter
        real(wp), intent(in) :: klong

        bandcenter = k0 + real(N, 8) * dk

        band  = bessel_jn(N, real(g2, 8)) * exp(-(klong - bandcenter)**2/(4._wp * sigmalong**2))
    end function side_band

    
end module ModulatedPulses
