module application
    use parameters
    use particles
    use sort_parts
    use get_neighbours
    use kernels
    use equations
    use density
    use forces
    use integrator
    use geometries
    implicit none

    !> Principal controller for the SPH simulation
    !> Handles initialization, time-stepping, and data management
    type :: SPH_App
        type(ParticleSystem) :: sys
        type(ParticleSort)   :: srt
        
        ! Simulation timing and control
        real(prec) :: dt          = 0.0005_prec
        real(prec) :: maxTime     = 2.0_prec
        real(prec) :: currentTime = 0.0_prec
        integer    :: currentIt   = 0
        integer    :: saveInt     = 50    ! Save every N iterations
        
        ! Physics parameters for Europa cryovolcanism
        real(prec) :: alpha       = 0.1_prec  ! Artificial Viscosity alpha
        real(prec) :: beta        = 0.0_prec  ! Artificial Viscosity beta
        real(prec) :: rho0        = 1000.0_prec
        real(prec) :: c0          = 15.0_prec ! Reference speed of sound
        real(prec) :: gamma_eos   = 7.0_prec  ! EoS stiffening parameter
        real(prec) :: dom_dim     = 2.0_prec  ! Domain dimension (square/cube)
        
    contains
        procedure :: run_simulation
        procedure :: solve_one_step
        procedure :: setup_scene
        procedure :: save_results
    end type SPH_App

contains

    !> 1. Setup the initial scene (Dam Break on Europa)
    subroutine setup_scene(self)
        class(SPH_App) :: self
        real(prec) :: spacing, h0, y_top, y_part, p_init, rho_init, B, g_abs
        integer    :: i, s
        
        ! --- Physical Configuration (Earth Dam-Break) ---
        spacing = 0.04_prec
        h0      = 1.2_prec * spacing
        self%rho0      = 1000.0_prec
        self%c0        = 65.0_prec
        self%gamma_eos = 7.0_prec
        self%dt        = 0.0001_prec
        
        call self%sys%init(40000)
        call self%srt%init(self%dom_dim, h0, KAPPA_CUBIC)
        
        ! Generate geometry
        call setup_dam_break(self%sys, 4.0_prec, 3.0_prec, &
                             1.0_prec, 2.0_prec, spacing, h0)
                             
        ! --- Hydrostatic Initialization ---
        ! Top of the fluid column is at y_start + H_fluid
        y_top = (0.5_prec * spacing) + 2.0_prec
        ! Tait Equation constant B = (rho0 * c0^2) / gamma
        B = (self%rho0 * self%c0**2) / self%gamma_eos
        g_abs = abs(G_VEC(nDim))

        do i = 1, self%sys%nPart
            y_part = self%sys%coords(nDim, 1, i)
            
            ! 1. Compute Hydrostatic Pressure: P = rho * g * depth
            if (y_part < y_top) then
                p_init = self%rho0 * g_abs * (y_top - y_part)
            else
                p_init = 0.0_prec
            end if
            
            ! 2. Compute Consistent Density from EoS: 
            ! rho = rho0 * (P/B + 1)^(1/gamma)
            rho_init = self%rho0 * ((p_init / B) + 1.0_prec)**(1.0_prec / self%gamma_eos)
            
            ! 3. Apply to all time slots to prevent initial oscillations
            do s = 1, nSteps
                self%sys%pressure(s, i) = p_init
                self%sys%density(s, i)  = rho_init
            end do
            
            ! 4. Anchor fixed particles (DBC)
            if (.not. self%sys%mobile(i)) then
                do s = 2, nSteps
                    self%sys%coords(:, s, i) = self%sys%coords(:, 1, i)
                end do
            end if
        end do

        ! Mass calculation for each particle
        self%sys%mass(:) = self%rho0 * (spacing**nDim)
        self%sys%h_part(:) = h0
        
        print *, "Hydrostatic setup complete. Max Pressure: ", &
                 self%rho0 * g_abs * 2.0_prec, " Pa"
    end subroutine setup_scene


    !> 2. Main simulation time loop
    subroutine run_simulation(self)
        class(SPH_App) :: self
        integer, allocatable :: next_array(:)
        
        print *, "Starting SPH Solver | Int: ", SELECTED_INT, &
                 " | Density: ", DENSITY_METHOD
        
        do while (self%currentTime < self%maxTime)
            self%currentIt = self%currentIt + 1
            
            ! Perform a complete physical step
            call self%solve_one_step(next_array)
            
            self%currentTime = self%currentTime + self%dt
            
            ! Data output at specified intervals
            if (mod(self%currentIt, self%saveInt) == 0) then
                call self%save_results()
                print '(A,F8.4,A,I6)', " Time: ", self%currentTime, &
                      " | Iter: ", self%currentIt
            end if
        end do
        
        if (allocated(next_array)) deallocate(next_array)
    end subroutine run_simulation


    !> 3. Orchestrates one complete SPH step (Workflow)
    subroutine solve_one_step(self, next_array)
        class(SPH_App) :: self
        integer, allocatable, intent(inout) :: next_array(:)
        integer :: target_slot

        ! --- STEP 1: Evaluation at t (Slot 1) ---
        call self%srt%build_grid()
        call self%srt%sort(self%sys, self%dom_dim, next_array)
        call compute_all_neighbors(self%sys, self%srt, next_array, KAPPA_CUBIC)
        
        call update_density(self%sys, 1)
        call pressure_wc(self%sys, self%rho0, self%c0, self%gamma_eos, 1)
        call compute_all_forces(self%sys, self%alpha, self%beta, 1)

        ! --- STEP 2: Predictor / First Kick ---
        call integration_phase_1(self%sys, self%dt)

        ! --- STEP 3: Evaluation at intermediate point ---
        if (SELECTED_INT /= INT_EULER) then
            ! For RK22, we evaluate at Slot 2 (t + dt/2)
            ! For Verlet, we evaluate at Slot 3 (t + dt)
            target_slot = merge(2, 3, SELECTED_INT == INT_RK22)
            
            call self%srt%build_grid()
            call self%srt%sort(self%sys, self%dom_dim, next_array)
            call compute_all_neighbors(self%sys, self%srt, next_array, KAPPA_CUBIC)
            
            call update_density(self%sys, target_slot)
            call pressure_wc(self%sys, self%rho0, self%c0, self%gamma_eos, target_slot)
            call compute_all_forces(self%sys, self%alpha, self%beta, target_slot)
        end if

        ! --- STEP 4: Corrector / Final Kick ---
        call integration_phase_2(self%sys, self%dt)
        call shift_data(self%sys)
        
    end subroutine solve_one_step


    !> 4. Export simulation data to CSV format
    subroutine save_results(self)
        class(SPH_App) :: self
        integer :: u, i
        character(len=32) :: filename
        
        write(filename, '(A,I0.4,A)') "output_", self%currentIt, ".csv"
        open(newunit=u, file=filename, status='replace')
        
        ! Header for CSV post-processing (ParaView/Python)
        write(u, '(A)') "x,y,v_mag,rho,p,mobile"
        
        do i = 1, self%sys%nPart
            write(u, '(5(F12.6,A),L2)') &
                self%sys%coords(1, 1, i), ",", &
                self%sys%coords(2, 1, i), ",", &
                sqrt(sum(self%sys%velocity(:, 1, i)**2)), ",", &
                self%sys%density(1, i), ",", &
                self%sys%pressure(1, i), ",", &
                self%sys%mobile(i)
        end do
        close(u)
    end subroutine save_results

end module application