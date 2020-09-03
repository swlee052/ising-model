module parameterstemp
    integer, parameter :: latt_leng = 5
    integer, parameter :: n_flips = 1 * latt_leng**2
    integer, parameter :: n_T_increment = 500
    integer, parameter :: EXCH_CONST = 1 ! exchange constant

    integer, dimension(-1 : latt_leng, -1 : latt_leng) :: spin_of ! (-1 ... latt_leng) X (-1 ... latt_leng)
    double precision, dimension(-2 : 2) :: accept_prob_chart
    double precision, dimension(0 : n_T_increment) :: T ! (0 ... n_T_increment), Unit = kT/J
end module


program Ising_Model ! computes the critical temperature of 2D Ising Model, EXCH_CONST = 1
    use parameterstemp
    implicit none       ! Require all variables to be explicitly declared
    integer :: i

    call create_lattice_sorted
    do i = 0, n_T_increment
        T(i) = 5d-1 + 3d0 * real(i, 8) / real(n_T_increment, 8) ! update temperature, unit = kT/abs(EXCH_CONST)
        if(EXCH_CONST == -1) then ! Anti-ferromagnetic, EXCH_CONST = -1
            accept_prob_chart = (/ exp(-8d0 / T(i)), exp(-4d0 / T(i)), &
            1d0, exp(4d0 / T(i)), exp(8d0 / T(i)) /) ! Acceptance probability chart
        else ! Ferromagnetic, EXCH_CONST = 1
            accept_prob_chart = (/ exp(8d0 / T(i)), exp(4d0 / T(i)), &
            1d0, exp(-4d0 / T(i)), exp(-8d0 / T(i)) /)
        end if
        call flip_sites_metropolis ! Iterate n_flips times
    end do
    
    
end program


subroutine create_lattice_sorted
    use parameterstemp
    integer :: i, j, k

    do i = 0, latt_leng - 1, 1 ! Assign spins on each site
        do j = 0, latt_leng - 1, 1
            spin_of(i, j) = 1 ! All up spins
        end do
    end do

    do k = 0, latt_leng - 1, 1 ! Assign spins to padding sites to make it periodic
        spin_of(-1, k) = spin_of(latt_leng - 1, k)
        spin_of(latt_leng, k) = spin_of(0, k)
        spin_of(k, -1) = spin_of(k, latt_leng - 1)
        spin_of(k, latt_leng) = spin_of(k, 0)
    end do
end subroutine


subroutine flip_sites_metropolis
    use parameterstemp
    integer :: x, y, i, sum_spin_nbr, n_measures, row, col
    double precision :: accept_prob
    double precision :: u, v, r

    n_measures = 0
    do i = 0, n_flips - 1
        call random_number(u) ! Return random real number in [0,1)
        call random_number(v)
        x = FLOOR(u * latt_leng) ! Pick random site, [0, latt_leng - 1]
        y = FLOOR(v * latt_leng)
        call random_number(r) ![0, 1)
        sum_spin_nbr = spin_of(x - 1, y) + spin_of(x + 1, y) + spin_of(x, y + 1) + spin_of(x, y - 1)
        accept_prob = accept_prob_chart ((spin_of(x, y) * sum_spin_nbr)/2) ! Sigma(SS) = spin_of(x, y) * sum_nbr = -4, -2 , 0, 2, 4
        if (r < accept_prob) then
            spin_of(x, y) = - spin_of(x, y) ! If random number < accept_prob => flip site
            if (x == 0) then
                spin_of(latt_leng, y) = spin_of(x, y) ! Update padding sites
            else if (y == 0) then
                spin_of(x, latt_leng) = spin_of(x, y)
            else if (x == latt_leng - 1) then
                spin_of(-1, y) = spin_of(x, y)
            else if (y == latt_leng - 1) then
                spin_of(x, -1) = spin_of(x, y)
            end if
        print *, "-------"
        do row = -1, latt_leng
            print '(1x, 7I2)', (spin_of(row,col), col = -1, latt_leng) 
        end do
        print *, "-------"
        end if
    end do
end subroutine
