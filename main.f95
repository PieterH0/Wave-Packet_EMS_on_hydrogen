! Main program
! Here all simulation parameter are loaded from the PARAMETER.nlm file or
! from the Parameters.f95 file that contains the default values
program main2
    use SpecialFunctions
    use Basics
    use ModulatedPulses
    use FourierTransforms
    use Tmatrix
    use ProbabilityCalculation
    use IntegrationModPulses
    use Parameters
    implicit none

    ! define working precision
    integer, parameter :: wp = dp 

    ! simulation specific variables are defined here
    real(wp) :: period

    ! Folder for saveing data and storing PARAMETER.nml file
    character(len=100) :: folder_character = '' 

    ! Load parameters from the settings file !
    ! The default parameter are hidden in the Parameters.f95 module
    namelist /EMS_SETTINGS/ theta, k0, Nphi, minPhi, maxPhi, Nkpoints, mink, maxk, time, Nthetai_1, Nthetai_2, Nphii, &
        & bx, by, pulse_type, modulated, g2, omega, sigmatrans, sigmalong, sigmatheta, multiple, delay_time, Nint, thetai_1_max
    open(file="" // trim(folder_character) // 'PARAMETER.nml', unit=1)
    read(nml=EMS_SETTINGS, unit=1)
    close(1)

    ! multiple types of pulses are defined in the ModulatedPulses.f95 module
    ! here pointer is set for the pulse type used 
    ! Options are: 'gaussian', 'energyga', 'modulatd', 'modulPRA'
    ! pulse_type must be character of length equal to 8
    if (pulse_type .eq. 'gaussian') then
        pulse_to_use => gaussian_pulse
    else if (pulse_type .eq. 'energyga') then
        pulse_to_use => energy_gaussian_pulse
    else if (pulse_type .eq. 'modulatd') then
        pulse_to_use => modulated_pulse
    else if (pulse_type .eq. 'modulPRA') then
        pulse_to_use => modulated_pulse
    else
        print*,"Error: Unknown Pulse type. Must be either gaussian, energyga, modulatd, modulpra"
    end if

    ! based on the values of g|2|, omega and k0 we calculate dk and tbunch (both set to 0._wp in Parameters as default)
    dk     = omega/k0 
    tbunch = 0._wp    ! k0**2/(omega**2*dble(g2)) ! here dble() converts integer to double precision real type. Outcomment when 2|g|=0! 
    period = 288._wp *2._wp*pi/7._wp

    
    call mainroutine(trim(folder_character), period*0._wp, '', Nkpoints, Nphi, bx, by) 


    !!!! TO COMPILE !!!!
    ! gfortran Basics.f95 FourierTransforms.f95 ModulatedPulses.f95 SpecialFunctions.f95 Parameters.f95 Tmatrix.f95 ProbabilityCalculation.f95 IntegratingModPulses.f95 main2.f95 -march=native -fopenmp -O3 -ftree-vectorize -o runpro

end program main2