subroutine compute_Cv(heat_cap) ! compute the heat capacity of the lattice using density of state
    use parameters
    integer :: tot_E
    double precision :: Z, E_avg, M_avg, E2_avg, weight
    double precision, intent(out) :: heat_cap
    
    Z = 0, E_avg = 0, M_avg = 0, E2_avg = 0, weight = 0
    do tot_E = -2*(latt_leng**2), 2*(latt_leng**2), 4
    weight = exp(-tot_E / T) * dens_of_tot_E
    Z = Z + weight
    E_avg = E_avg + weight * tot_E
    E2_avg = E2_avg + weight * tot_E ** 2
    end do
    E_avg = E_avg / Z
    E2_avg = E2_avg/ Z
    heat_cap = (E2_avg - E_avg ** 2) / (real (latt_leng**2, 8) * T**2)
end subroutine


subroutine compute_density_of_state_E(dens_of_tot_E)
   use parameters
   integer :: tot_E
   integer, dimension(-2*(latt_leng**2) : 2*(latt_leng**2) : 4), intent(out) :: dens_of_tot_E
    
   do tot_E = -2*(latt_leng**2), 2*(latt_leng**2), 4 ! compute density of state for each possible E
       dens_of_tot_E(tot_E) = 
end subroutine
