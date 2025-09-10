module FourierTransforms
    implicit none
    private

    ! Double floating point number
    integer, parameter :: dp = kind(1.0d0) 

    ! Now define the working precision of this module 
    integer, parameter :: wp = dp

    public :: fourier_hydrogen_orb, coherent_fourier_orb


contains
    ! We hard code the radial fourier transforms for the simple hydrogen states
    ! the integrals can be calculated analytically using Laplace transforms or looked up in B&J appendix A.
    function fourier_hydrogen_orb(q, n, l) result(orbital)
        ! quantum numbers
        integer, intent(in)  :: n, l
        ! momentum
        real(wp), intent(in) :: q

        ! probability amplitude
        real(wp) :: orbital


        if (n == 1 .and. l == 0) then
            orbital = 4._wp/(1._wp + q**2)**2
        else if (n == 2 .and. l == 0) then
            orbital = 1._wp/sqrt(2._wp) * (2._wp/(1._wp/4._wp + q**2)**2 - 1._wp/(1._wp/4._wp + q**2)**3)
        else if (n ==2 .and. l == 1) then
            orbital = 2._wp * q / sqrt(6._wp) * 1._wp/(1._wp/4._wp + q**2)**3
        else if (n == 3 .and. l == 1) then
            orbital = 864._wp/sqrt(6._wp) * q * (9._wp*q**2 - 1._wp)/(9._wp*q**2 + 1._wp)**4
        else if (n == 4 .and. l == 1) then
            orbital = 256._wp/sqrt(60._wp) * 4._wp*q/(16._wp*q**2 + 1._wp)**3 * &
                        (12._wp * (16._wp*q**2 - 1._wp)**2/(16._wp*q**2 + 1._wp)**2 - 2._wp)
        else
            write(*,*)"Orbital not yet known"
        end if

    end function fourier_hydrogen_orb

    ! hardcoded coherent superposition
    function coherent_fourier_orb(q, t, n1, n2, l1, l2) result(orbital)
        integer, intent(in)  :: n1, n2, l1, l2
        real(wp), intent(in) :: q, t
        real(wp) :: omega1, omega2, orbital1, orbital2

        complex(wp) :: orbital

        ! energies
        omega1 = 0.5_wp/n1**2
        omega2 = 0.5_wp/n2**2

        orbital1 = fourier_hydrogen_orb(q,n1,l1)
        orbital2 = fourier_hydrogen_orb(q,n2,l2)

        orbital = 1._wp/sqrt(2._wp) * orbital1  * exp(cmplx(0, -omega1*t, wp)) + &
                        1._wp/sqrt(2._wp) * orbital2 * exp(cmplx(0, -omega2*t, wp))
        
    end function coherent_fourier_orb


end module FourierTransforms