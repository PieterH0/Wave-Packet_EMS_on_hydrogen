module SpecialFunctions
! module containing special function such a spherical harmonics and Legendre polynomials 
use Basics
implicit none
private

! Define working precision 
integer, parameter :: wp = dp 

public :: sph_harm, assoc_legendre, py_sph_harm
    
contains
    ! Function to determine the spherical harmonics 
    function sph_harm(cos_theta, phi, l, m) result(sph_harmonic)
        integer     :: l,m 
        real(wp)    :: cos_theta, phi 
        complex(wp) :: sph_harmonic 
        real(wp)    :: sqrt_fac, plm
        integer     :: abs_m

        ! Calcualate legendre poly for positive m and check sign afterwards 
        abs_m = abs(m)
        sqrt_fac = sqrt((2._wp*l + 1._wp)/(4._wp*pi) * gamma(l-abs_m+1._wp)/gamma(l+abs_m+1._wp))
        plm = assoc_legendre(l, abs_m, cos_theta)

        if (m>=0) then 
            sph_harmonic = sqrt_fac * plm * exp(cmplx(0._wp, m*phi, wp))   
        else 
            sph_harmonic = (-1._wp)**abs_m * sqrt_fac * plm * exp(cmplx(0._wp, m*phi, wp))  
        end if

    end function sph_harm

    ! function to determine associated legendre polynomial from recursion
    function assoc_legendre(l, m, x) result(plm)
        integer, intent(in)  :: l, m
        real(wp), intent(in) :: x
        real(wp) :: plm      ! the result 
        real(wp) :: pmm      ! the m, m associoated legendre poly. Start of recursion
        real(wp) :: sqrt_fact, d_factorial_factor, pmm1, pmj
        integer  :: i, j     ! index for looping

        

        pmm = 1._wp ! the p00 polynomial
        if (m.gt.0) then
            sqrt_fact = sqrt((1._wp + x)*(1._wp - x))
            d_factorial_factor = 1._wp
            do i=1, m
                pmm = -pmm*d_factorial_factor * sqrt_fact
                d_factorial_factor = d_factorial_factor + 2
            end do
        end if

        ! now for l not equal m.
        if (l==m) then
            plm = pmm
        else
            ! if not first we calt pm(l = m+1)
            pmm1 = (2._wp * m + 1._wp) * x * pmm

            ! we check if this is the wanted poly.
            if (l == m+1) then
                plm = pmm1
            
            else
                ! if not we start recursion
                do j=m+2, l !j is the l value checked
                     pmj = ((2._wp * j - 1._wp) * x * pmm1 - (j + m - 1._wp) * pmm) / (j-m)
                     pmm = pmm1
                     pmm1= pmj
                enddo
                plm = pmj
            end if
        end if
    end function assoc_legendre

    function py_sph_harm(cos_theta, phi) result(sph_harmonic)
        ! Special for the py states, which is a superposition of m=-1 and the m=1 spherical harmonics
        real(wp)    :: cos_theta, phi ! arguments
        complex(wp) :: sph_harmonic   ! result



        sph_harmonic = sqrt(3._wp/(4._wp*pi)) * sqrt(1._wp - cos_theta**2) * sin(phi)
    end function 

end module SpecialFunctions