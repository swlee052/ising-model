! How to remove ^M in source code on linux %s/^M//g (^M = Ctrl+v+m)
! Note on subroutine random_number
! Returns a single pseudorandom number or an array of pseudorandom numbers from the uniform distribution over the range [0, 1).
! The runtime-library implements the xorshift1024* random number generator(RNG). This generator has a period of 2^{1024} - 1,
! and when using multiple threads up to 2^{512} threads can each generate 2^{512} random numbers before any aliasing occurs.

module parameters
    implicit none
    integer, parameter :: latt_leng = 30
    integer, parameter :: EXCH_CONST = 1 ! exchange constant
    integer, parameter :: n_T_step = 500
    double precision, parameter :: T_start = 0.5d0
    double precision, parameter :: T_end = 3.5d0
    integer, dimension(-1 : latt_leng, -1 : latt_leng) :: spin_of ! (-1 ... latt_leng) X (-1 ... latt_leng)
    double precision, dimension(-2 : 2) :: accept_prob_chart
end module
module Class_Site
    implicit none

    type, public :: Site
    integer :: x, y, spin
    end type
end module
program Ising_Model ! computes the Tc of 2D Ising Model, EXCH_CONST = 1
    use parameters
    implicit none ! Require all variables to be explicitly declared
    integer :: T_step
    integer, parameter :: MC_max = 100* latt_leng**2  ! number of iteration of MC simulation per temperature 
    double precision, dimension(0 : n_T_step) :: avg_E !(0 ... n_T_increments), Ensemble average of E per site
    double precision, dimension(0 : n_T_step) :: avg_M !(0 ... n_T_increments), Ensemble average of M per site
    double precision, dimension(0 : n_T_step) :: avg_Cv !(0 ... n_T_increments), Ensemble average of Cv per site
    double precision :: T, del_T, Tc_avg_E, Tc_avg_M, Tc_avg_Cv, avg_en, avg_mag, avg_cap ! Unit = kT/J, Tc computed with each of the physical observables
    
    call create_lattice(2) ! 1 : random, 2 : aligned
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
        call flip_sites_wolff(T, MC_max, avg_en, avg_mag, avg_cap) ! Iterate MC simulation MC_max times per temperature
        avg_E(T_step) = avg_en ! Ensemble average of E per site
        avg_M(T_step) = avg_mag  ! Ensemble average of M per site
        avg_Cv(T_step) = avg_cap ! (<E**2>-<E>**2)/T**2
    end do
    call compute_Tc_with_E_or_M(avg_E, Tc_avg_E) ! Find maximum slope
    call compute_Tc_with_E_or_M(avg_M, Tc_avg_M)
    call compute_Tc_with_Cv(avg_Cv, Tc_avg_Cv) ! Find maximum value
    call print_data(1, avg_E, Tc_avg_E) ! switch = 1 for <E> data file, 2 for <M>, 3 for <Cv>
    call print_data(2, avg_M, Tc_avg_M)
    call print_data(3, avg_Cv, Tc_avg_Cv)
    call make_data_file(1, avg_E, Tc_avg_E)
    call make_data_file(2, avg_M, Tc_avg_M)
    call make_data_file(3, avg_Cv, Tc_avg_Cv)
end program


subroutine flip_sites_wolff(T, MC_max, avg_en, avg_mag, avg_en2)
    use parameters
    use Class_Site
    implicit none
!    interface 
!        subroutine resize(arr, capacity_new) 
!            use Class_Site
!            integer, intent(in) :: capacity_new
!            integer :: capacity_old
!            type(Site), dimension(:), allocatable, intent(inout) :: arr
!            type(Site), dimension(:), allocatable :: temp
!        end subroutine
        
!        subroutine print_clstr_bag(clstr_bag)
!            use Class_Site
!            implicit none
!            integer :: i
!            type(Site), dimension(:), allocatable, intent(in):: clstr_bag
!        end subroutine

!        subroutine print_bndry_bag(bndry_bag)
!            use Class_Site
!            implicit none
!            integer :: i
!            type(Site), dimension(:), allocatable, intent(in) :: bndry_bag
!        end subroutine
!    end interface
    
    logical, dimension(0 : latt_leng - 1, 0 : latt_leng - 1) :: is_in_clstr_bag
    integer :: i, j, x, y, ran_index, MC_step, n_sample
    integer :: n_site_bndry_bag, n_site_clstr_bag, size_clstr_bag, size_bndry_bag
    integer, intent(in) :: MC_max
    type(Site), dimension(0 : latt_leng**2 - 1) :: clstr_bag
    type(Site), dimension(0 : latt_leng**2 - 1) :: bndry_bag
    type(Site), dimension(0:3) :: nbr_site_bag 
    double precision, intent(in) :: T ! unit = kT/J
    double precision, intent(out) :: avg_en, avg_en2, avg_mag
    double precision :: accept_prob, tot_en, tot_mag, tot_en2, sum_tot_E, sum_tot_M, sum_tot_E2
    double precision :: u, v, w, r

    n_sample = 0 ! Initialize variables and arrays
    n_site_clstr_bag = 0
    n_site_bndry_bag = 0
    size_clstr_bag = 0
    size_bndry_bag = 0
    is_in_clstr_bag (:, :) = .false.
    nbr_site_bag(:) = Site(0, 0, 0)
    clstr_bag(:) = Site(0, 0, 0)
    bndry_bag(:) = Site(0, 0, 0)
    sum_tot_E = 0d0
    sum_tot_M = 0d0
    sum_tot_E2 = 0d0
    accept_prob = 1 - exp(-2*T)
    
    do MC_step = 0, MC_max - 1
!        print *, "                                                        "
!        print *, "++++++++++++++++++MC Iteration Start++++++++++++++++++++"
        call random_number(u) ! Return random real number in [0,1)
        call random_number(v)
        x = FLOOR(u * latt_leng) ! Pick random site, [0, latt_leng - 1]
        y = FLOOR(v * latt_leng)
        !allocate(clstr_bag(0:1))
        !allocate(bndry_bag(0:1))
        !size_clstr_bag = 2
        !size_bndry_bag = 2
        !clstr_bag(:) = Site(0, 0, 0)
        !bndry_bag(:) = Site(0, 0, 0)
        clstr_bag(0) = Site(x, y, spin_of(x, y))
        bndry_bag(0) = Site(x, y, spin_of(x, y))
        is_in_clstr_bag(x, y) = .true.
        n_site_clstr_bag = n_site_clstr_bag + 1
        n_site_bndry_bag = n_site_bndry_bag + 1
!        call print_bndry_bag(bndry_bag, n_site_bndry_bag)
!        call print_clstr_bag(clstr_bag, n_site_clstr_bag)
!        call print_spin_of
!        call print_is_in_clstr_bag(is_in_clstr_bag)
!        call print_nbr_site_bag(nbr_site_bag)
!        print *, "++++++++++++++++++Form Cluster Start++++++++++++++++++++"
        do while (n_site_bndry_bag > 0) 
            call random_number(w)
            ran_index = FLOOR(w * n_site_bndry_bag) ! Pick random site from clstr_bag, [0, n_site_clstr_bag - 1]
            
            nbr_site_bag(0)%x = bndry_bag(ran_index)%x - 1 ! Assign nbr sites to nbr_site_bag
            nbr_site_bag(0)%y = bndry_bag(ran_index)%y
            nbr_site_bag(0)%spin = spin_of(nbr_site_bag(0)%x, nbr_site_bag(0)%y)
            nbr_site_bag(1)%x = bndry_bag(ran_index)%x + 1
            nbr_site_bag(1)%y = bndry_bag(ran_index)%y
            nbr_site_bag(1)%spin = spin_of(nbr_site_bag(1)%x, nbr_site_bag(1)%y)
            nbr_site_bag(2)%x = bndry_bag(ran_index)%x
            nbr_site_bag(2)%y = bndry_bag(ran_index)%y - 1
            nbr_site_bag(2)%spin = spin_of(nbr_site_bag(2)%x, nbr_site_bag(2)%y)
            nbr_site_bag(3)%x = bndry_bag(ran_index)%x 
            nbr_site_bag(3)%y = bndry_bag(ran_index)%y + 1
            nbr_site_bag(3)%spin = spin_of(nbr_site_bag(3)%x, nbr_site_bag(3)%y)
            
            do i = 0, 3 ! Add nbr_sites to cluster with magic probability = 1-exp(-2beta)
                call cnvrt_padng_to_latt(nbr_site_bag(i)%x, nbr_site_bag(i)%y) ! convert padding sites to original sites, two if-else if checks
                if (.not.(is_in_clstr_bag(nbr_site_bag(i)%x, nbr_site_bag(i)%y)) &
                .and. (nbr_site_bag(i)%spin == bndry_bag(ran_index)%spin)) then
                    call random_number(r) ![0, 1)
                    if (r < accept_prob) then
!                        if(n_site_clstr_bag == size_clstr_bag) then ! resize if necessary
!                            call resize(clstr_bag, 2*size_clstr_bag)
!                            size_clstr_bag = 2*size_clstr_bag
!                        else if (n_site_clstr_bag <= size_clstr_bag/4) then
!                            call resize(clstr_bag, size_clstr_bag/2)
!                            size_clstr_bag = size_clstr_bag/2
!                        end if
!                        if(n_site_bndry_bag == size_bndry_bag) then  ! resize if necessary
!                            call resize(bndry_bag, 2*size_bndry_bag)
!                            size_bndry_bag = 2*size_bndry_bag
!                        else if (n_site_bndry_bag <= size_bndry_bag/4) then
!                            call resize(bndry_bag, size_bndry_bag/2)
!                            size_bndry_bag = size_bndry_bag/2
!                        end if
                        
                        clstr_bag(n_site_clstr_bag) = Site(nbr_site_bag(i)%x, nbr_site_bag(i)%y, &
                        spin_of(nbr_site_bag(i)%x, nbr_site_bag(i)%y))
                        bndry_bag(n_site_bndry_bag) = Site(nbr_site_bag(i)%x, nbr_site_bag(i)%y, &
                        spin_of(nbr_site_bag(i)%x, nbr_site_bag(i)%y))
                        is_in_clstr_bag(nbr_site_bag(i)%x, nbr_site_bag(i)%y) = .true.
                        n_site_clstr_bag = n_site_clstr_bag + 1
                        n_site_bndry_bag = n_site_bndry_bag + 1
                    end if
                end if
            end do
            bndry_bag(ran_index) = bndry_bag(n_site_bndry_bag - 1) ! Pop chosen random site from bndry_bag
            bndry_bag(n_site_bndry_bag - 1) = Site(0, 0, 0)
            n_site_bndry_bag = n_site_bndry_bag - 1
        end do
!        print *, "++++++++++++++++++Form Cluster End++++++++++++++++++++"
!        call print_bndry_bag(bndry_bag, n_site_bndry_bag)
!        call print_clstr_bag(clstr_bag, n_site_clstr_bag)
!        call print_spin_of
!        call print_is_in_clstr_bag(is_in_clstr_bag)
!        call print_nbr_site_bag(nbr_site_bag)
!        print *, "++++++++++++++++++Flip Cluster Start++++++++++++++++++++"
        do j = 0, n_site_clstr_bag - 1 ! Flip sites in clstr_bag
            spin_of(clstr_bag(j)%x, clstr_bag(j)%y) = - spin_of(clstr_bag(j)%x, clstr_bag(j)%y)
            if (clstr_bag(j)%x == 0) then ! Update padding sites
                spin_of(latt_leng, clstr_bag(j)%y) = spin_of(clstr_bag(j)%x, clstr_bag(j)%y) 
            else if (clstr_bag(j)%x == latt_leng - 1) then
                spin_of(-1, clstr_bag(j)%y) = spin_of(clstr_bag(j)%x, clstr_bag(j)%y)
            end if 
            if (clstr_bag(j)%y == 0) then
                spin_of(clstr_bag(j)%x, latt_leng) = spin_of(clstr_bag(j)%x, clstr_bag(j)%y)
            else if (clstr_bag(j)%y == latt_leng - 1) then
                spin_of(clstr_bag(j)%x, -1) = spin_of(clstr_bag(j)%x, clstr_bag(j)%y)
            end if
        end do
!        print *, "++++++++++++++++++Flip Cluster End++++++++++++++++++++"
!        call print_bndry_bag(bndry_bag, n_site_bndry_bag)
!        call print_clstr_bag(clstr_bag, n_site_clstr_bag)
!        call print_spin_of
!        call print_is_in_clstr_bag(is_in_clstr_bag)
!        call print_nbr_site_bag(nbr_site_bag)
!        print *, "++++++++++++++++++Deallocate Start++++++++++++++++++++"
        !deallocate(clstr_bag)
        !deallocate(bndry_bag)
        !size_clstr_bag = 0
        !size_bndry_bag = 0
        nbr_site_bag(:) = Site(0, 0, 0)
        is_in_clstr_bag(:, :) = .false.
        n_site_clstr_bag = 0
        n_site_bndry_bag = 0
        do i = 0, n_site_clstr_bag - 1
            clstr_bag(:) = Site(0, 0, 0)
            bndry_bag(:) = Site(0, 0, 0)
        end do
!        print *, "++++++++++++++++++Deallocate End++++++++++++++++++++"
!        call print_bndry_bag(bndry_bag, n_site_bndry_bag)
!        call print_clstr_bag(clstr_bag, n_site_clstr_bag)
!        call print_spin_of
!        call print_is_in_clstr_bag(is_in_clstr_bag)
!        call print_nbr_site_bag(nbr_site_bag)
!        print *, "++++++++++++++++++Add to Ensemble Start++++++++++++++++++++" 
        if (MC_step > MC_step / 2 .and. MOD(MC_step, latt_leng**2) == 0) then ! Create ensemble for sampled microstates
            tot_en = 0d0
            tot_mag = 0d0
            tot_en2 = 0d0
            call compute_tot_E_M_E2(tot_en, tot_mag, tot_en2)
            sum_tot_E = sum_tot_E + tot_en
            sum_tot_M = sum_tot_M + tot_mag ! The value M can flip due to fluctuation
            sum_tot_E2 = sum_tot_E2 + tot_en2
            n_sample = n_sample + 1
!            print *, "sum_tot_E =", sum_tot_E
!            print *, "sum_tot_M =", sum_tot_M
!            print *, "sum_tot_E2 =", sum_tot_E2
!            print *, "tot_en = ", tot_en
!            print *, "tot_mag = ", tot_mag
!            print *, "tot_en2 = ", tot_en2
!            print *, "n_sample = ", n_sample
        end if
!        print *, "++++++++++++++++++Add to Ensemble End++++++++++++++++++++"
!        print *, "++++++++++++++++++MC Iteration End++++++++++++++++++++"
    end do
!    print *, "++++++++++++++++++Ensemble averaging Start++++++++++++++++++++" 
    avg_en = sum_tot_E / real(n_sample, 8) / latt_leng**2 ! Ensemble average of E per site
    avg_mag = sum_tot_M / real(n_sample, 8) / latt_leng**2 ! Ensemble average of M per site
    avg_en2 = (sum_tot_E2 / real(n_sample, 8) / latt_leng**2 - latt_leng**2 * avg_en**2) / T**2 ! Ensemble average of Cv per site, Cv = (<E**2>-<E>**2)/T**2 
end subroutine

subroutine create_lattice(switch)
    use parameters
    implicit none
    integer :: x, y
    integer, intent(in) :: switch
    double precision :: u

    if (switch == 1) then ! random
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
    
    else if (switch == 2) then ! Aligned, i.e., All up spins
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
    end if
end subroutine
subroutine resize(arr, capacity_new) 
    use Class_Site
    implicit none 
    
    integer, intent(in) :: capacity_new
    integer :: i, capacity_old
    type(Site), dimension(:), allocatable, intent(inout) :: arr
    type(Site), dimension(:), allocatable :: temp
    
    capacity_old = size(arr)
    allocate(temp(0:capacity_new-1))
    temp(:) = Site(0, 0, 0)
    do i = 0, min(capacity_old, capacity_new)-1
        temp(i) = arr(i)
    end do
    deallocate(arr)
    allocate(arr(0:capacity_new-1))
    arr(:) = Site(0, 0, 0)
    do i = 0, min(capacity_old, capacity_new)-1
        arr(i) = temp(i)
    end do
    deallocate(temp)
end subroutine
subroutine cnvrt_padng_to_latt(x, y)
    use parameters
    implicit none
    integer, intent(inout) :: x, y

    if (x == -1) then ! Check if the nbr site is a padding site
        x = latt_leng - 1
    else if (x == latt_leng) then
        x = 0
    end if
    if (y == -1) then
        y = latt_leng - 1
    else if (y == latt_leng) then
        y = 0
    end if   
end subroutine
subroutine compute_tot_E_M_E2(tot_E, tot_M, tot_E2)
    use parameters
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
    integer :: T_step
    integer, intent(in) :: switch
    integer, parameter :: out_unit1 = 10, out_unit2 = 20
    double precision :: del_T, T
    double precision, intent(in) :: Tc
    double precision, dimension(0 : n_T_step), intent(in) :: phys_obsv

    del_T = (T_end - T_start) / n_T_step
    if (switch == 1) then
        open (out_unit1, FILE = 'output_ising_model_clstr_avg_E_change.txt', action = "write", STATUS = 'replace')
    else if (switch == 2) then
        open (out_unit1, FILE = 'output_ising_model_clstr_avg_M_change.txt', action = "write", STATUS = 'replace')
    else if (switch == 3) then
        open (out_unit1, FILE = 'output_ising_model_clstr_avg_Cv_change.txt', action = "write", STATUS = 'replace')
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
        open (out_unit2, FILE = 'output_ising_model_clstr_avg_E_change_plot.txt', action = "write", STATUS = 'replace')
    else if (switch == 2) then 
        open (out_unit2, FILE = 'output_ising_model_clstr_avg_M_change_plot.txt', action = "write", STATUS = 'replace')
    else if (switch == 3) then 
        open (out_unit2, FILE = 'output_ising_model_clstr_avg_Cv_change_plot.txt', action = "write", STATUS = 'replace')
    end if
    do T_step = 0, n_T_step
        T = T_start + del_T * T_step
        write (out_unit2,*) T, " ", phys_obsv(T_step)
    end do
    close(out_unit2)
end subroutine


!For debugging
subroutine test_cpu_time
    real :: start, finish
    call cpu_time(start)
    ! put code to test here
    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start
end subroutine
subroutine print_clstr_bag(clstr_bag)
    use Class_Site
    implicit none
    integer :: i
    type(Site), dimension(:), allocatable, intent(in):: clstr_bag

    do i = 0, size(clstr_bag) - 1
        print *, clstr_bag(i)%x, clstr_bag(i)%y, clstr_bag(i)%spin
    end do
    !print *, "n_site_clstr_bag = ", n_site_clstr_bag
    !print *, "size_clstr_bag = ", size(clstr_bag)
    print *, "======clstr_bag======"
end subroutine
subroutine print_bndry_bag(bndry_bag)
    use Class_Site
    implicit none
    integer :: i
    type(Site), dimension(:), allocatable, intent(in) :: bndry_bag

    do i = 0, size(bndry_bag) - 1
    print *, bndry_bag(i)%x, bndry_bag(i)%y, bndry_bag(i)%spin
    end do
    !print *, "n_site_bndry_bag = ", n_site_bndry_bag
    !print *, "size_bndry_bag = ", size_bndry_bag
    print *, "======bndry_bag======"
end subroutine
subroutine print_spin_of
    use Class_Site
    use parameters
    implicit none
    integer :: i
    
    do i = -1, latt_leng
        print *, spin_of(:, i)
    end do
    print *, "======spin_of======"
end subroutine
subroutine print_is_in_clstr_bag(arr)
    use Class_Site
    use parameters
    implicit none
    integer :: i
    logical, dimension(0 : latt_leng - 1, 0 : latt_leng - 1) :: arr
    
    do i = 0, latt_leng - 1
        print *, arr(0:, i)
    end do
    print *, "======is_in_clstr_bag======"
end subroutine
subroutine print_nbr_site_bag(arr)
    use Class_Site
    implicit none
    integer :: i
    type(Site), dimension(0:3), intent(in) :: arr

    do i = 0, 3
    print *, arr(i)%x, arr(i)%y, arr(i)%spin
    end do
    print *, "======nbr_site_bag======"
end subroutine
