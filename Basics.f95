module Basics
    !
    ! This module contains basic definitions ans function to be used in other modules
    ! 
    implicit none
    private

    ! Different kinds of floating point numbers
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0) 

    ! Some constants 
    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp   

    ! Define the working precision for this module! 
    integer, parameter :: wp = dp

    ! Public constants and kinds 
    public sp, dp, pi   

    ! public subroutines
    public linspace, possible_zero_divided_by_zero, save_array


    
contains
    ! Just a standard linspace subroutine, nothing fancy 
    ! shamelessly stolen from Mads 
    subroutine linspace(arr, start, end, N)
        integer, intent(in) :: N 
        real(wp), intent(in)  :: start, end 
        real(wp), intent(out) :: arr(N)
        real(wp) :: step 
        integer  :: i 

        step = (end - start) / (N-1)

        do i = 1, N
            arr(i) = start + (i-1)*step 
        end do 
    end subroutine linspace

    function possible_zero_divided_by_zero(x, y) result(z)
        real(wp), intent(in) :: x, y
        real(wp) :: z
        if (abs(x) == 0._wp .and. y == 0._wp) then
            z = 1._wp
        else if (y == 0._wp) then
            stop "Division by zero!!!" 
        else
            z = x/y
        end if
    end function possible_zero_divided_by_zero

    ! again shamelessly stolen. Jesus forgives!
    subroutine save_array(file_name, array, unit)
        character(*), intent(in) :: file_name
        real(wp), intent(in) :: array(:)
        integer, optional, intent(in) :: unit 
        integer :: index 

        if (present(unit)) then 
            index = unit + 6  ! To avoid 5 and 6?
        else 
            index = 1
        end if

        open(file=file_name, access='stream', unit=index, form='unformatted')
        ! open(file=file_name, unit=index, form='formatted') ! if use formatted write needs a format argument
        ! The "access='stream'" statement is needed if the files need to be loaded in python
        ! This is because Fortran add an integer if FORM='UNFORMATTED'.
        ! https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnaf/index.html
        write(1) array 
        close(1)
    end subroutine save_array
    
end module Basics