module forces
    use parameters
    use particles
    use kernels
    use equations
    implicit none

    !> Module to compute the momentum equation and accelerations
    !> Standard SPH formulation for pressure and artificial viscosity

contains

    !> Compute acceleration for all particles
    !> @param sys: The particle system (SoA)
    !> @param alpha: Artificial viscosity coefficient (linear)
    !> @param beta: Artificial viscosity coefficient (quadratic)
    !> @param step: The time slot index (1, 2, or 3)
    subroutine compute_all_forces(sys, alpha, beta, step)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: alpha, beta
        integer,              intent(in)    :: step
        
        integer    :: i, k, d, neighbor_idx
        real(prec) :: r, r2, h_i, W, dWdr, pi_ab, p_rho2_ab
        real(prec) :: dx(nDim), gradW(nDim)
        real(prec) :: g_europa(3)

        ! Gravity on Europa (negative Y-direction)
        g_europa = 0.0_prec
        g_europa(2) = -1.314_prec 

        !$OMP PARALLEL DO PRIVATE(i, k, d, neighbor_idx, dx, r2, r, W, dWdr, &
        !$OMP                     h_i, gradW, pi_ab, p_rho2_ab)
        do i = 1, sys%nPart
            ! 1. Initialize acceleration with gravity (only for mobile particles)
            if (sys%mobile(i)) then
                ! Initialize with gravity vector directly
                sys%acceler(1:nDim, i) = G_VEC(1:nDim)
            else
                sys%acceler(1:nDim, i) = 0.0_prec
                cycle 
            end if
            
            h_i = sys%h_part(i)
            
            ! 2. Loop over neighbors
            do k = 1, sys%neigh_count(i)
                neighbor_idx = sys%neigh_list(sys%neigh_ptr(i) + k - 1)
                
                ! Distance vector dx = x_i - x_j
                r2 = 0.0_prec
                do d = 1, nDim
                    dx(d) = sys%coords(d, step, i) - sys%coords(d, step, neighbor_idx)
                    r2 = r2 + dx(d)**2
                end do
                r = sqrt(r2)
                
                ! Avoid calculation if particles are exactly at the same spot
                if (r < 1.0e-12_prec) cycle
                
                ! 3. Kernel derivative
                call cubic_spline(r, h_i, W, dWdr)
                
                ! 4. Compute Grad W (vector)
                ! gradW = (dW/dr) * (vec_x_ab / r)
                do d = 1, nDim
                    gradW(d) = (dx(d) / r) * dWdr
                end do
                
                ! 5. Pressure term: (P_i/rho_i^2 + P_j/rho_j^2)
                p_rho2_ab = (sys%pressure(step, i) / sys%density(step, i)**2) + &
                            (sys%pressure(step, neighbor_idx) / sys%density(step, neighbor_idx)**2)
                
                ! 6. Artificial Viscosity term (Pi_ab)
                pi_ab = compute_art_visc(sys, i, neighbor_idx, alpha, beta, step)
                
                ! 7. Momentum Equation update
                ! dv/dt = -sum( m_j * (P/rho^2 + Pi_ab) * gradW )
                do d = 1, nDim
                    sys%acceler(d, i) = sys%acceler(d, i) - &
                                        sys%mass(neighbor_idx) * &
                                        (p_rho2_ab + pi_ab) * gradW(d)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine compute_all_forces

end module forces