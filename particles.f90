module particles
    use parameters
    implicit none

    !> Structure of Arrays (SoA) for high performance SPH
    type :: ParticleSystem
        integer :: nPart          = 0
        integer :: maxPart        = 0
        
        ! --- Vector Properties (dimension, time_slot, particle_index) ---
        real(prec), allocatable :: coords(:,:,:)   
        real(prec), allocatable :: velocity(:,:,:) 
        real(prec), allocatable :: acceler(:,:)    ! dV/dt
        
        ! --- Scalar Properties (time_slot, particle_index) ---
        real(prec), allocatable :: density(:,:)    
        real(prec), allocatable :: pressure(:,:)   
        real(prec), allocatable :: c_sound(:,:)    
        real(prec), allocatable :: u_therm(:,:)    
        
        ! --- Time Derivatives (particle_index) ---
        ! These are needed for the Continuity and Energy equations
        real(prec), allocatable :: drhodt(:)       ! dRho/dt
        real(prec), allocatable :: dudt(:)         ! du/dt (Internal Energy rate)
        
        ! --- Constant/Slower Properties (particle_index) ---
        real(prec), allocatable :: mass(:)         
        real(prec), allocatable :: h_part(:)       
        logical,    allocatable :: mobile(:)       
        
        ! --- Neighbor Management ---
        integer, allocatable    :: neigh_list(:)   
        integer, allocatable    :: neigh_ptr(:)    
        integer, allocatable    :: neigh_count(:)  
        integer                 :: max_neighbors   
        
    contains
        procedure :: init => init_system
        procedure :: add_particle => add_single_particle
    end type ParticleSystem

contains

    subroutine init_system(self, capacity)
        class(ParticleSystem) :: self
        integer, intent(in)   :: capacity

        self%maxPart = capacity
        
        ! Vectors
        allocate(self%coords(nDim, nSteps, self%maxPart))
        allocate(self%velocity(nDim, nSteps, self%maxPart))
        allocate(self%acceler(nDim, self%maxPart))
        
        ! Scalars
        allocate(self%density(nSteps, self%maxPart))
        allocate(self%pressure(nSteps, self%maxPart))
        allocate(self%c_sound(nSteps, self%maxPart))
        allocate(self%u_therm(nSteps, self%maxPart))
        
        ! Derivatives (Missing parts added here)
        allocate(self%drhodt(self%maxPart))
        allocate(self%dudt(self%maxPart))
        
        ! Constants
        allocate(self%mass(self%maxPart))
        allocate(self%h_part(self%maxPart))
        allocate(self%mobile(self%maxPart))

        ! Neighbors
        if (nDim == 2) then
            self%max_neighbors = self%maxPart * 64
        else
            self%max_neighbors = self%maxPart * 256
        end if
        
        allocate(self%neigh_list(self%max_neighbors))
        allocate(self%neigh_ptr(self%maxPart))
        allocate(self%neigh_count(self%maxPart))
        
        ! Initialize derivatives to zero
        self%drhodt = 0.0_prec
        self%dudt   = 0.0_prec
    end subroutine init_system

    subroutine add_single_particle(self, x_vec, v_vec, m, h, is_mobile)
        class(ParticleSystem) :: self
        real(prec), dimension(nDim), intent(in) :: x_vec, v_vec
        real(prec), intent(in) :: m, h
        logical, intent(in)    :: is_mobile
        
        if (self%nPart >= self%maxPart) then
            print *, "Error: ParticleSystem capacity exceeded. Resize needed."
            return
        end if

        self%nPart = self%nPart + 1
        self%coords(:, 1, self%nPart)   = x_vec
        self%velocity(:, 1, self%nPart) = v_vec
        self%mass(self%nPart)           = m
        self%h_part(self%nPart)         = h
        self%mobile(self%nPart)         = is_mobile
    end subroutine add_single_particle

end module particles