module ProbabilityCalculation
    use SpecialFunctions
    use Basics
    use ModulatedPulses
    use FourierTransforms
    use Tmatrix
    use Parameters, only : theta, k0, minPhi, maxPhi, mink, maxk, dk, Nthetai_1, Nthetai_2, Nphii, &
        & modulated, g2, omega, multiple, tbunch, delay_time, Nint, pulse_type ! other parameter set when calling mainroutine
    use IntegrationModPulses, only: set_up_kf_integrals
    implicit none
    private

    integer, parameter :: wp = dp 

    public mainroutine

    
    contains

    subroutine mainroutine(folder, time, titlestring, Nkpoints, Nphi, bx, by)
        ! The the actual simulation based on the setting in the main.f95 and Parameters.nlm
        real(wp), intent(in)         :: time           ! delay time (coefficients in superposition t_c)
        character(len=*), intent(in) :: folder         ! folder to save results
        character(len=*), intent(in) :: titlestring    ! character added to the data file title
        integer, intent(in)   :: Nkpoints, Nphi        ! integration points
        
        integer :: i,j                              ! do loop index
        real(wp), allocatable :: kfpoints(:)        ! kf points
        real(wp), allocatable :: phipoints(:)       ! phi points

        !
        ! For saving evaluations in integrations
        complex(wp) :: amplitude_single_sym, amplitude_single_asym ! if only need single amplitude variable
        real(wp), allocatable :: dprobabilitydk(:)  ! Differential probabilities (depends on kf and phi)
        real(wp) :: probability_single_sym, probability_single_asym 
        real(wp), allocatable :: probability(:)     ! Final probabilities
        !
        real(wp) ::  ka(3), kb(3), kahat(3), kbhat(3), dkf
        real(wp), intent(in) :: bx, by ! impact parameter.
        real :: start, finish   ! taking time

        
        ! for saving
        character(200) :: saving_string_phi, saving_string_prob
        character(20)  :: number_of_k, number_of_p

        ! Some print statements to catch possible errors
        write(*,*)" ####  Starting Calculation  #### "
        write(*,*)"Nk, Np, Nt_i, Np_i, tc=", Nkpoints, Nphi, Nthetai_1, Nphii, time
        write(*,*)"Pulse type = ", pulse_type
        write(*,*)"Modulated, 2|g|, omega, tbunch = ", Modulated, g2, omega, tbunch
        write(*,*)"Impact parameter bx, by = ", bx, by
        write(*,*)"Time delay", delay_time
        write(*,*)"Momentum shift dk = ", dk
        write(*,*)"titlestring : ", titlestring 
        call cpu_time(start)

        ! Allocate for arrays
        allocate(phipoints(Nphi))
        allocate(kfpoints(Nkpoints))
        
        allocate(probability(Nphi))
        allocate(dprobabilitydk(Nkpoints))

        ! Grid over ka and kb
        call linspace(phipoints,  minPhi, maxPhi,   Nphi)

        ! Whether we have multiple seperated kf integral or a single interval
        if (multiple .eqv. .false.) then
            call linspace(kfpoints, mink, maxk, Nkpoints)  ! test for mod.
            !call linspace(kfpoints,    19.16281224699351_wp, 19.176125517950467_wp, Nkpoints) ! works for unmodulated Bojer pulse (see Notebook)
            !call linspace(kfpoints,    19.1_wp, 19.2_wp, Nkpoints)      ! works for Shao pulse
            !call linspace(kfpoints,    19.166_wp, 19.172_wp, Nkpoints)  ! test for mod.
        else
            ! Note that Nkpoints/(2Nint+1) needs to be an integer!
            call set_up_kf_integrals(Nint, Nkpoints, 0.00005_wp, 0.0003_wp, kfpoints)
            call save_array('kfpoints.dat', kfpoints)
        end if

        dkf = kfpoints(2) - kfpoints(1)
        write(*,*)"d_kf = ", dkf
        !!! Loop over grid points !!!
        ! The loop is parallelized over the detection phi_f points in the phipoints vector
        ! As a default the variables are shared
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(kahat, kbhat, ka, kb, amplitude_single_sym, &
        !$OMP amplitude_single_asym, probability_single_sym, probability_single_asym, dprobabilitydk)
        do i=1, Nphi
            ! First loop over "measured" phi_f angles
            ! Evaluate the directions of the scattered particles
            kahat = [sin(theta)*cos(phipoints(i)),      sin(theta)*sin(phipoints(i)),      cos(theta)]
            kbhat = [sin(theta)*cos(pi - phipoints(i)), sin(theta)*sin(pi - phipoints(i)), cos(theta)]
            !
            !
            do j=1, Nkpoints
                ! integrate over measured final momentum length
                ka = kfpoints(j) * kahat
                kb = kfpoints(j) * kbhat
                !
                ! The scattering probability for all phi_f and k_f are calculated with the
                ! "calculalte_T_matrix_element" subroutine with is found in the Tmatrix module
                !
                !
                ! amplitude symmetric wavefunction
                amplitude_single_sym   = 2*pi*calculalte_T_matrix_element(ka, kb, 2, [3, 4], [1, 1], [0, 0], time, &
                                                                        &  bx, by, 1)
                probability_single_sym = (amplitude_single_sym%re)**2 + (amplitude_single_sym%im)**2
                !
                ! amplitude antisymmetric wavefunction
                amplitude_single_asym   = 2*pi*calculalte_T_matrix_element(ka, kb, 2, [3, 4], [1, 1], [0, 0], time, &
                                                                        &  bx, by, -1)
                probability_single_asym = (amplitude_single_asym%re)**2 + (amplitude_single_asym%im)**2
                ! sum norm-square of amplitude and multiply with integration step
                dprobabilitydk(j) = (1._wp/4._wp*probability_single_sym + 3._wp/4._wp*probability_single_asym) * kfpoints(j)**3
                !
            !
            end do
            !
            probability(i) = sum(dprobabilitydk) * dkf
        end do
        !$OMP END PARALLEL DO


        !!! saving in arrays
        ! formating strings for file 
        write(number_of_k, '(I5)') Nkpoints
        write(number_of_p, '(I4)') Nphi
        
        ! concatenating final strings
        saving_string_phi = "C:\Users\piete\Speciale\Data\" // trim(folder) // "\Phipoints.dat"
        !
        saving_string_prob = "C:\Users\piete\Speciale\Data\" // trim(folder) // "\K" &  
                    & // trim(adjustl(number_of_k)) // "P" // trim(adjustl(number_of_p)) // "T" // trim(titlestring) // ".dat"

        ! saving array
        call save_array(trim(saving_string_phi), phipoints)
        call save_array(trim(saving_string_prob), probability)

        call cpu_time(finish)
        print '("Time = ", f10.4," seconds.")', finish-start
        !
    end subroutine mainroutine


end module ProbabilityCalculation