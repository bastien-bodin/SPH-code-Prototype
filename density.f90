module density
    use parameters
    use particles
    use kernels
    implicit none

contains

    !> Unified density update routine
    subroutine update_density(sys, step)
        type(ParticleSystem), intent(inout) :: sys
        integer,              intent(in)    :: step

        if (DENSITY_METHOD == DENSITY_SUMMATION) then
            call compute_density_summation(sys, step)
        else
            call compute_density_continuity(sys, step)
        end if
    end subroutine update_density

    !> Method 1: Direct Summation (rho_i = sum m_j * W_ij)
    subroutine compute_density_summation(sys, step)
        type(ParticleSystem), intent(inout) :: sys
        integer,              intent(in)    :: step
        integer    :: i, k, neighbor_idx
        real(prec) :: r, r2, h_i, W, dWdr
        real(prec) :: dx(nDim)

        !$OMP PARALLEL DO PRIVATE(i, k, neighbor_idx, dx, r2, r, W, dWdr, h_i)
        do i = 1, sys%nPart
            h_i = sys%h_part(i)
            ! Self contribution
            call cubic_spline(0.0_prec, h_i, W, dWdr)
            sys%density(step, i) = sys%mass(i) * W
            
            do k = 1, sys%neigh_count(i)
                neighbor_idx = sys%neigh_list(sys%neigh_ptr(i) + k - 1)
                r2 = sum((sys%coords(:, step, i) - sys%coords(:, step, neighbor_idx))**2)
                r = sqrt(r2)
                
                call cubic_spline(r, h_i, W, dWdr)
                sys%density(step, i) = sys%density(step, i) + sys%mass(neighbor_idx) * W
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine compute_density_summation

    !> Method 2: Continuity Equation (drho/dt)
    subroutine compute_density_continuity(sys, step)
        type(ParticleSystem), intent(inout) :: sys
        integer,              intent(in)    :: step
        integer    :: i, k, j, neighbor_idx
        real(prec) :: r, r2, h_i, W, dWdr
        real(prec) :: dx(nDim), dv(nDim), gradW(nDim)

        !$OMP PARALLEL DO PRIVATE(i, k, neighbor_idx, dx, dv, r2, r, W, dWdr, h_i, gradW, j)
        do i = 1, sys%nPart
            sys%drhodt(i) = 0.0_prec
            h_i = sys%h_part(i)
            
            do k = 1, sys%neigh_count(i)
                neighbor_idx = sys%neigh_list(sys%neigh_ptr(i) + k - 1)
                dx = sys%coords(:, step, i) - sys%coords(:, step, neighbor_idx)
                dv = sys%velocity(:, step, i) - sys%velocity(:, step, neighbor_idx)
                r = sqrt(sum(dx**2))
                
                call cubic_spline(r, h_i, W, dWdr)
                gradW = (dx / r) * dWdr
                
                sys%drhodt(i) = sys%drhodt(i) + sys%mass(neighbor_idx) * dot_product(dv, gradW)
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine compute_density_continuity

end module density