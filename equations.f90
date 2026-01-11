module equations
    use parameters
    use particles
    implicit none

    !> Physical equations for SPH fluid dynamics
    !> Optimized for Structure of Arrays (SoA) and vectorization.

contains

    ! --- EQUATION OF STATE (EoS) ---

    !> Ideal Gas Equation of State: P = rho * R * T / M
    !> Simplified here as P = (gamma-1) * rho * u
    pure subroutine pressure_ideal_gas(sys, gamma_eos, step)
        type(ParticleSystem), intent(inout) :: sys
        real(prec), intent(in)              :: gamma_eos
        integer,    intent(in)              :: step
        integer :: p

        ! Vectorizable loop
        do p = 1, sys%nPart
            sys%pressure(step, p) = (gamma_eos - 1.0_prec) * &
                                    sys%density(step, p) * &
                                    sys%u_therm(step, p)
            
            ! Sound speed for Ideal Gas: c = sqrt(gamma * P / rho)
            sys%c_sound(step, p) = sqrt(gamma_eos * sys%pressure(step, p) / &
                                        sys%density(step, p))
        end do
    end subroutine pressure_ideal_gas


    !> Weakly Compressible Equation of State (Tait's Equation)
    !> P = B * ((rho/rho_0)^gamma - 1)
    pure subroutine pressure_wc(sys, rho_0, c_0, gamma_eos, step)
        type(ParticleSystem), intent(inout) :: sys
        real(prec), intent(in)              :: rho_0, c_0, gamma_eos
        integer,    intent(in)              :: step
        real(prec) :: B
        integer :: p

        ! Stiffening parameter B
        B = rho_0 * c_0**2 / gamma_eos

        do p = 1, sys%nPart
            sys%pressure(step, p) = B * ((sys%density(step, p) / rho_0)**gamma_eos - 1.0_prec)
            
            ! Sound speed for WC: c = c_0 * (rho/rho_0)^((gamma-1)/2)
            sys%c_sound(step, p) = c_0 * (sys%density(step, p) / rho_0)**((gamma_eos - 1.0_prec) / 2.0_prec)
        end do
    end subroutine pressure_wc


    ! --- ARTIFICIAL VISCOSITY ---

    !> Monaghan-type Artificial Viscosity
    !> Computes the term Pi_ab for the momentum equation
    pure function compute_art_visc(sys, a, b, alpha, beta, step) result(pi_ab)
        type(ParticleSystem), intent(in) :: sys
        integer, intent(in)              :: a, b ! Indices of particles
        real(prec), intent(in)           :: alpha, beta
        integer, intent(in)              :: step
        real(prec)                       :: pi_ab

        real(prec) :: v_dot_x, x_dot_x, h_ab, c_ab, rho_ab, mu_ab
        real(prec) :: v_ab(nDim), x_ab(nDim)
        
        pi_ab = 0.0_prec
        
        ! Relative velocity and position
        v_ab = sys%velocity(:, step, a) - sys%velocity(:, step, b)
        x_ab = sys%coords(:, step, a) - sys%coords(:, step, b)
        
        v_dot_x = dot_product(v_ab, x_ab)
        
        ! Only apply viscosity for particles approaching each other
        if (v_dot_x < 0.0_prec) then
            x_dot_x = dot_product(x_ab, x_ab)
            h_ab   = 0.5_prec * (sys%h_part(a) + sys%h_part(b))
            c_ab   = 0.5_prec * (sys%c_sound(step, a) + sys%c_sound(step, b))
            rho_ab = 0.5_prec * (sys%density(step, a) + sys%density(step, b))
            
            ! Viscosity parameter mu_ab
            mu_ab = (h_ab * v_dot_x) / (x_dot_x + 0.01_prec * h_ab**2)
            
            pi_ab = (-alpha * c_ab * mu_ab + beta * mu_ab**2) / rho_ab
        end if
    end function compute_art_visc

end module equations