! How to remove ^M in source code on linux %s/^M//g (^M = Ctrl+v+m)
! Note on subroutine random_number
! Returns a single pseudorandom number or an array of pseudorandom numbers from the uniform distribution over the range [0, 1).
! The runtime-library implements the xorshift1024* random number generator(RNG). This generator has a period of 2^{1024} - 1,
! and when using multiple threads up to 2^{512} threads can each generate 2^{512} random numbers before any aliasing occurs.

module parameters
    integer, parameter :: latt_leng = 80
    integer, parameter :: EXCH_CONST = 1 ! exchange constant
    integer, parameter :: n_T_step = 1000
    double precision, parameter :: T_start = 2.0d0
    double precision, parameter :: T_end = 2.5d0
    integer, dimension(-1 : latt_leng, -1 : latt_leng) :: spin_of ! (-1 ... latt_leng) X (-1 ... latt_leng)
    double precision, dimension(-2 : 2) :: accept_prob_chart
end module


program Ising_Model ! computes the Tc of 2D Ising Model, EXCH_CONST = 1
    use parameters
    implicit none ! Require all variables to be explicitly declared
    integer :: T_step
    integer(kind = 8), parameter :: MC_max = 1E11  ! number of iteration of MC simulation per temperature ~n^3 to equilibriate
    double precision, dimension(0 : n_T_step) :: avg_E !(0 ... n_T_increments), Ensemble average of E per site
    double precision, dimension(0 : n_T_step) :: avg_M !(0 ... n_T_increments), Ensemble average of M per site
    double precision, dimension(0 : n_T_step) :: avg_Cv !(0 ... n_T_increments), Ensemble average of Cv per site
    double precision :: T, del_T, Tc_avg_E, Tc_avg_M, Tc_avg_Cv, avg_en, avg_mag, avg_cap ! Unit = kT/J, Tc computed with each of the physical observables
   
    call create_lattice_aligned !Integrate create_lattice routines
    del_T = (T_end - T_start) / n_T_step
    do T_step = 0, n_T_step
        T = T_start + del_T * T_step ! update T, unit = kT/J
        if(EXCH_CONST == -1) then ! Anti-ferromagnetic, EXCH_CONST = -1
            accept_prob_chart = (/ exp(-8d0 / T), exp(-4d0 / T), &
            1d0, exp(4d0 / T), exp(8d0 / T) /) ! Acceptance probability chart
        else ! Ferromagnetic, EXCH_CONST = 1
            accept_prob_chart = (/ exp(8d0 / T), exp(4d0 / T), &
            1d0, exp(-4d0 / T), exp(-8d0 / T) /)
        end if
        call flip_sites_metropolis(T, MC_max, avg_en, avg_mag, avg_cap) ! Iterate MC simulation MC_max times per temperature
        avg_E(T_step) = avg_en ! Ensemble average of E per site
        avg_M(T_step) = avg_mag  ! Ensemble average of M per site
        avg_Cv(T_step) = avg_cap ! (<E**2>-<E>**2)/T**2
    end do
    call compute_Tc_with_E_or_M(avg_E, Tc_avg_E) ! Find maximum slope
    call compute_Tc_with_E_or_M(avg_M, Tc_avg_M)
    call compute_Tc_with_Cv(avg_Cv, Tc_avg_Cv) ! Find maximum value
    call print_data(1, avg_E, Tc_avg_E)
    call print_data(2, avg_M, Tc_avg_M)
    call print_data(3, avg_Cv, Tc_avg_Cv)
    call make_data_file(1, avg_E, Tc_avg_E)
    call make_data_file(2, avg_M, Tc_avg_M)
    call make_data_file(3, avg_Cv, Tc_avg_Cv)
end program


subroutine create_lattice_random
    use parameters
    integer :: x, y
    double precision :: u

    do x = 0, latt_leng - 1
        do y = 0, latt_leng - 1
            call random_number(u) ! Return random real number in [0,1)
            if (u >= 0.5) then
                spin_of(x, y) = 1 ! Assign random spin values
            else
                spin_of(x, y) = -1
            end if
        end do
    end do

    do x = 0, latt_leng - 1 ! Assign spins to top and bottom line of padding sites to make it periodic
        spin_of(x, -1) = spin_of(x, latt_leng - 1)
        spin_of(x, latt_leng) = spin_of(x, 0)
    end do
    do y = 0, latt_leng -1  ! Assign spins to left and right line of padding sites to make it periodic
        spin_of(-1, y) = spin_of(latt_leng - 1, y)
        spin_of(latt_leng, y) = spin_of(0, y)
    end do
end subroutine


subroutine create_lattice_aligned
    use parameters
    integer :: x, y

    do x = 0, latt_leng - 1 ! Assign spins on each site
        do y = 0, latt_leng - 1
            spin_of(x, y) = 1 ! All up spins
        end do
    end do

    do x = 0, latt_leng - 1 ! Assign spins to top and bottom line of padding sites to make it periodic
        spin_of(x, -1) = spin_of(x, latt_leng - 1)
        spin_of(x, latt_leng) = spin_of(x, 0)
    end do
    do y = 0, latt_leng -1  ! Assign spins to left and right line of padding sites to make it periodic
        spin_of(-1, y) = spin_of(latt_leng - 1, y)
        spin_of(latt_leng, y) = spin_of(0, y)
    end do
end subroutine


subroutine flip_sites_metropolis(T, MC_max, avg_en, avg_mag, avg_en2)
    use parameters
    integer :: x, y, sum_spin_nbr, sum_nbr_spin, n_sample 
    integer(kind = 8) :: MC_step
    integer(kind = 8), intent(in) :: MC_max
    double precision, intent(in) :: T
    double precision, intent(out) :: avg_en, avg_en2, avg_mag
    double precision :: accept_prob, tot_en, tot_mag, tot_en2, tot_E, tot_M, tot_E2, sum_tot_E, sum_tot_M, sum_tot_E2
    double precision :: u, v, r

    n_sample = 0 ! Initialize variable to be used in averaging samples
    sum_tot_E = 0d0
    sum_tot_M = 0d0
    sum_tot_E2 = 0d0
    call compute_tot_E_M_E2(tot_en, tot_mag, tot_en2) ! Compute initial values of physical observables
    tot_E = tot_en
    tot_M = tot_mag
    tot_E2 = tot_en2
    do MC_step = 0, MC_max - 1
        call random_number(u) ! Return random real number in [0,1)
        call random_number(v)
        x = FLOOR(u * latt_leng) ! Pick random site, [0, latt_leng - 1]
        y = FLOOR(v * latt_leng)
        sum_nbr_spin = spin_of(x - 1, y) + spin_of(x + 1, y) + spin_of(x, y + 1) + spin_of(x, y - 1)
        accept_prob = accept_prob_chart((spin_of(x, y) * sum_nbr_spin)/2) ! Sigma(SS) = spin_of(x, y) * sum_surrounding = -4, -2 , 0, 2, 4
        call random_number(r) ![0, 1)
        if (r < accept_prob) then
            spin_of(x, y) = - spin_of(x, y) ! If random number < accept_prob => flip site
            
            if (x == 0) then ! Update padding sites
                spin_of(latt_leng, y) = spin_of(x, y) 
            else if (x == latt_leng - 1) then
                spin_of(-1, y) = spin_of(x, y)
            end if 
            if (y == 0) then
                spin_of(x, latt_leng) = spin_of(x, y)
            else if (y == latt_leng - 1) then
                spin_of(x, -1) = spin_of(x, y)
            end if
            
            tot_M = tot_M + 2*spin_of(x, y) ! Update physical observables
            sum_spin_nbr = spin_of(x + 1, y) + spin_of(x - 1, y) + spin_of(x, y + 1) + spin_of(x, y - 1)
            tot_E = tot_E + (- EXCH_CONST) * spin_of(x, y) * sum_spin_nbr
            tot_E2 = tot_E**2
        end if
        
        if (MC_step > MC_step / 2 .and. MOD(MC_step, latt_leng**2) == 0) then
            sum_tot_E = sum_tot_E + tot_E
            sum_tot_M = sum_tot_M + abs(tot_M) ! The value M can flip due to fluctuation
            sum_tot_E2 = sum_tot_E2 + tot_E2
            n_sample = n_sample + 1
        end if
    end do
    avg_en = sum_tot_E / real(n_sample, 8) / latt_leng**2 ! Ensemble average of E per site
    avg_mag = sum_tot_M / real(n_sample, 8) / latt_leng**2 ! Ensemble average of M per site
    avg_en2 = (sum_tot_E2 / real(n_sample, 8) / latt_leng**2 - latt_leng**2 * avg_en**2) / T**2 ! Ensemble average of Cv per site, Cv = (<E**2>-<E>**2)/T**2 
end subroutine


subroutine compute_tot_E_M_E2(tot_E, tot_M, tot_E2)
    use parameters
    integer :: x, y, sum_spin_nbr
    double precision, intent(out) :: tot_E, tot_M, tot_E2
    
    tot_E = 0d0
    tot_M = 0d0
    tot_E2 = 0d0
    do x = 0, latt_leng - 1
        do y = 0, latt_leng - 1
            tot_M = tot_M + spin_of(x, y)
            sum_spin_nbr = spin_of(x + 1, y) + spin_of(x - 1, y) + spin_of(x, y + 1) + spin_of(x, y - 1)
            tot_E = tot_E + (- EXCH_CONST) * spin_of(x, y) * sum_spin_nbr
        end do
    end do
    tot_E = tot_E / 2
    tot_E2 = tot_E**2
end subroutine


subroutine compute_Tc_with_E_or_M(phys_obsv, Tc) ! Compute Tc by calculating T where dp_q / dT is max.
    use parameters
    integer :: T_step
    double precision :: max_slope, slope, del_T, T
    double precision, intent(out) :: Tc
    double precision, dimension(0 : n_T_step), intent(in) :: phys_obsv

    max_slope = 0d0
    slope = 0d0
    del_T = (T_end - T_start) / n_T_step
    do T_step = 0, n_T_step - 1
        T = T_start + del_T * T_step
        slope = (phys_obsv(T_step + 1) - phys_obsv(T_step)) / del_T
        if (slope > max_slope) then 
            max_slope = slope
            Tc = T + 0.5*del_T
        end if
    end do
end subroutine


subroutine compute_Tc_with_Cv(Cv, Tc)
    use parameters
    integer :: T_step
    double precision :: max_value, del_T, T
    double precision, intent(out) :: Tc
    double precision, dimension(0 : n_T_step), intent(in) :: Cv

    max_value = 0d0
    del_T = (T_end - T_start) / n_T_step
    do T_step = 0, n_T_step
        T = T_start + del_T * T_step
        if (Cv(T_step) > max_value) then 
            max_value = Cv(T_step)
            Tc = T
        end if
    end do
end subroutine


subroutine print_data(switch, phys_obsv, Tc)
    use parameters
    integer :: T_step
    integer, intent(in) :: switch
    double precision :: T, del_T
    double precision, intent(in) :: Tc
    double precision, dimension(0 : n_T_step), intent(in) :: phys_obsv

    del_T = (T_end - T_start) / n_T_step
    do T_step = 0, n_T_step
        T = T_start + del_T * T_step
        print *, "kT/J = ", T
        if (switch == 1) then
            print *, "<E> = ", phys_obsv(T_step)
        else if (switch == 2) then 
            print *, "<M> = ", phys_obsv(T_step)
        else if (switch == 3) then 
            print *, "<Cv> = ", phys_obsv(T_step)
        end if
        print *, "--------------"
    end do
    print *, "Tc  = ", Tc, " (kT/J)"
end subroutine


subroutine make_data_file(switch, phys_obsv, Tc)
    use parameters
    integer :: T_step
    integer, intent(in) :: switch
    integer, parameter :: out_unit1 = 10, out_unit2 = 20
    double precision :: del_T, T
    double precision, intent(in) :: Tc
    double precision, dimension(0 : n_T_step), intent(in) :: phys_obsv

    del_T = (T_end - T_start) / n_T_step
    if (switch == 1) then
        open (out_unit1, FILE = 'output_ising_model_avg_E_change.txt', action = "write", STATUS = 'replace')
    else if (switch == 2) then
        open (out_unit1, FILE = 'output_ising_model_avg_M_change.txt', action = "write", STATUS = 'replace')
    else if (switch == 3) then
        open (out_unit1, FILE = 'output_ising_model_avg_Cv_change.txt', action = "write", STATUS = 'replace')
    end if
    
    write (out_unit1,*) "T (kT/J) / avg_" 
    if (switch == 1) then 
        write (out_unit1,*)  "E (no units)"
    else if (switch == 2) then 
        write (out_unit1,*)  "M (no units)"
    else if (switch == 3) then 
        write (out_unit1,*)  "Cv (no units)"
    end if
    do T_step = 0, n_T_step
        T = T_start + del_T * T_step
        write (out_unit1,*) T, " ", phys_obsv(T_step)
    end do
    if (switch == 1) then 
        write (out_unit1,*)  "Tc_computed_with_avg_E = ", Tc,  " (kT/J)"
    else if (switch == 2) then 
        write (out_unit1,*)  "Tc_computed_with_avg_M = ", Tc,  " (kT/J)"
    else if (switch == 3) then 
        write (out_unit1,*)  "Tc_computed_with_avg_Cv = ", Tc,  " (kT/J)"
    end if
    close(out_unit1)
    
    if (switch == 1) then 
        open (out_unit2, FILE = 'output_ising_model_avg_E_change_plot.txt', action = "write", STATUS = 'replace')
    else if (switch == 2) then 
        open (out_unit2, FILE = 'output_ising_model_avg_M_change_plot.txt', action = "write", STATUS = 'replace')
    else if (switch == 3) then 
        open (out_unit2, FILE = 'output_ising_model_avg_Cv_change_plot.txt', action = "write", STATUS = 'replace')
    end if
    do T_step = 0, n_T_step
        T = T_start + del_T * T_step
        write (out_unit2,*) T, " ", phys_obsv(T_step)
    end do
    close(out_unit2)
end subroutine
