module integrator
    use parameters
    use particles
    implicit none

contains

    !> Phase 1: Initial update (Predictor or First Kick)
    subroutine integration_phase_1(sys, dt)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: dt
        
        select case (SELECTED_INT)
        case (INT_RK22)
            call predict_rk22(sys, dt)
        case (INT_VERLET)
            call kick_drift_verlet(sys, dt)
        end select
    end subroutine integration_phase_1

    !> Phase 2: Final update (Corrector or Final Kick)
    subroutine integration_phase_2(sys, dt)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: dt
        
        select case (SELECTED_INT)
        case (INT_RK22)
            call correct_rk22(sys, dt)
        case (INT_VERLET)
            call final_kick_verlet(sys, dt)
        case (INT_EULER)
            call step_euler(sys, dt)
        end select
    end subroutine integration_phase_2

    ! --- Simple Euler ---

    subroutine step_euler(sys, dt)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: dt
        integer :: i
        !$OMP PARALLEL DO
        do i = 1, sys%nPart
            if (.not. sys%mobile(i)) cycle
            if (DENSITY_METHOD == DENSITY_CONTINUITY) &
                sys%density(2, i) = sys%density(1, i) + sys%drhodt(i) * dt
            sys%velocity(:, 2, i) = sys%velocity(:, 1, i) + sys%acceler(:, i) * dt
            sys%coords(:, 2, i)   = sys%coords(:, 1, i)   + sys%velocity(:, 1, i) * dt
        end do
        !$OMP END PARALLEL DO
    end subroutine step_euler

    ! --- RK22 (Midpoint Method) ---
    
    ! Predictor: Move from Slot 1 to Slot 2 (t -> t + dt/2)
    subroutine predict_rk22(sys, dt)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: dt
        integer :: i
        real(prec) :: dt2
        dt2 = dt * 0.5_prec

        !$OMP PARALLEL DO
        do i = 1, sys%nPart
            if (.not. sys%mobile(i)) cycle
            if (DENSITY_METHOD == DENSITY_CONTINUITY) &
                sys%density(2, i) = sys%density(1, i) + sys%drhodt(i) * dt2
            
            sys%velocity(:, 2, i) = sys%velocity(:, 1, i) + sys%acceler(:, i) * dt2
            sys%coords(:, 2, i)   = sys%coords(:, 1, i)   + sys%velocity(:, 2, i) * dt2
        end do
        !$OMP END PARALLEL DO
    end subroutine predict_rk22

    ! Corrector: Use derivatives from Slot 2 to move Slot 1 to Slot 3 (t -> t + dt)
    subroutine correct_rk22(sys, dt)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: dt
        integer :: i

        !$OMP PARALLEL DO
        do i = 1, sys%nPart
            if (.not. sys%mobile(i)) cycle
            if (DENSITY_METHOD == DENSITY_CONTINUITY) &
                sys%density(3, i) = sys%density(1, i) + sys%drhodt(i) * dt
            
            sys%velocity(:, 3, i) = sys%velocity(:, 1, i) + sys%acceler(:, i) * dt
            sys%coords(:, 3, i)   = sys%coords(:, 1, i)   + sys%velocity(:, 2, i) * dt
        end do
        !$OMP END PARALLEL DO
    end subroutine correct_rk22

    ! --- Symplectic Velocity Verlet (3-slot implementation) ---

    ! Phase 1: v(t+1/2) = v(t) + a(t)*dt/2 ; x(t+1) = x(t) + v(t+1/2)*dt
    subroutine kick_drift_verlet(sys, dt)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: dt
        integer :: i, d
        real(prec) :: dt_half

        dt_half = 0.5_prec * dt

        !$OMP PARALLEL DO PRIVATE(i, d)
        do i = 1, sys%nPart
            if (.not. sys%mobile(i)) cycle
            
            ! 1. Kick velocity to Slot 2 (t + dt/2)
            do d = 1, nDim
                sys%velocity(d, 2, i) = sys%velocity(d, 1, i) + sys%acceler(d, i) * dt_half
            end do
            
            ! 2. Drift position to Slot 3 (t + dt) using half-step velocity
            do d = 1, nDim
                sys%coords(d, 3, i) = sys%coords(d, 1, i) + sys%velocity(d, 2, i) * dt
            end do
            
            ! 3. Density prediction at t+dt (for Slot 3 forces)
            if (DENSITY_METHOD == DENSITY_CONTINUITY) then
                sys%density(3, i) = sys%density(1, i) + sys%drhodt(i) * dt
            end if
        end do
        !$OMP END PARALLEL DO
    end subroutine kick_drift_verlet

    ! Phase 2: v(t+1) = v(t+1/2) + a(t+1)*dt/2
    subroutine final_kick_verlet(sys, dt)
        type(ParticleSystem), intent(inout) :: sys
        real(prec),           intent(in)    :: dt
        integer :: i, d
        real(prec) :: dt_half

        dt_half = 0.5_prec * dt

        !$OMP PARALLEL DO PRIVATE(i, d)
        do i = 1, sys%nPart
            if (.not. sys%mobile(i)) cycle
            
            ! Update velocity from Slot 2 (t+dt/2) to Slot 3 (t+dt)
            do d = 1, nDim
                sys%velocity(d, 3, i) = sys%velocity(d, 2, i) + sys%acceler(d, i) * dt_half
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine final_kick_verlet

    ! --- Data Management ---

    !> Shift data from final slot to Slot 1 for the next time step
    subroutine shift_data(sys)
        type(ParticleSystem), intent(inout) :: sys
        integer :: last_slot, i 

        last_slot = merge(2, 3, SELECTED_INT == INT_EULER)

        !$OMP PARALLEL DO PRIVATE(i)
        do i = 1, sys%nPart
            if (sys%mobile(i)) then
                sys%coords(:, 1, i)   = sys%coords(:, last_slot, i)
                sys%velocity(:, 1, i) = sys%velocity(:, last_slot, i)
            end if
            
            sys%density(1, i)  = sys%density(last_slot, i)
            sys%pressure(1, i) = sys%pressure(last_slot, i)
        end do
        !$OMP END PARALLEL DO
    end subroutine shift_data

end module integrator