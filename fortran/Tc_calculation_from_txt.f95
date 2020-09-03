module parameters
    integer, parameter :: latt_leng = 30
    integer, parameter :: EXCH_CONST = 1 ! exchange constant
    integer, parameter :: n_T_step = 500
    double precision, parameter :: T_start = 0.5d0
    double precision, parameter :: T_end = 3.5d0
    integer, dimension(-1 : latt_leng, -1 : latt_leng) :: spin_of ! (-1 ... latt_leng) X (-1 ... latt_leng)
end module

program Tc_calculate
    use parameters
    implicit none
    integer i
    double precision, dimension(0 : n_T_step):: M
    double precision, dimension(0 : n_T_step):: Cv
    double precision, dimension(0 : n_T_step):: E
    double precision :: T, Tc_E, Tc_M, Tc_Cv
    
    open (unit=15, file='output_ising_model_avg_E_change_plot_5.txt', status='old',    &
             access='sequential', form='formatted', action='read' )
    open (unit=20, file='output_ising_model_avg_M_change_plot_5.txt', status='old',    &
             access='sequential', form='formatted', action='read' )
    open (unit=25, file='output_ising_model_avg_Cv_change_plot_5.txt', status='old',    &
             access='sequential', form='formatted', action='read' )
    
    do i = 0, n_T_step
        read (15, *) T, E(i)
        read (20, *) T, M(i)
        read (25, *) T, Cv(i)
    end do
    
    call compute_Tc_with_E_or_M(E, Tc_E)
    call compute_Tc_with_E_or_M(M, Tc_M)
    call compute_Tc_with_Cv(Cv, Tc_Cv)
    print *, "Tc_E = ", Tc_E
    print *, "Tc_M = ", Tc_M 
    print *, "Tc_Cv = ", Tc_Cv 
    close(15)
    close(20)
    close(25)
end program


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
