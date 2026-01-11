module sort_parts
    use parameters
    use particles
    implicit none

    !> Standard Linked-Cell sorting for SPH
    !> This module handles the spatial partitioning of particles into a grid
    type :: ParticleSort
        real(prec)    :: h_max     = 0.0_prec ! Maximum smoothing length
        real(prec)    :: cellSize  = 0.0_prec ! Size of a cell (usually kappa * h_max)
        integer       :: nCellsRow = 0        ! Number of cells along one axis
        integer       :: nCellsTot = 0        ! Total number of cells (Row^nDim)
        
        ! --- Linked-Cell Arrays ---
        ! head(i) stores the index of the first particle in cell i
        ! next(p) stores the index of the next particle in the same cell as p
        integer, allocatable :: head(:)
        
        logical :: is_initialized = .false.
        
    contains
        procedure :: init      => init_sort
        procedure :: build_grid => build_spatial_grid
        procedure :: sort      => sort_particles
    end type ParticleSort

contains

    !> Setup the grid dimensions based on the domain size and smoothing length
    subroutine init_sort(self, dom_dim, h_max_val, kappa)
        class(ParticleSort) :: self
        real(prec), intent(in) :: dom_dim   ! Domain size (assumed cubic/square)
        real(prec), intent(in) :: h_max_val ! Max smoothing length
        integer,    intent(in) :: kappa     ! Kernel support factor

        self%h_max    = h_max_val
        self%cellSize = real(kappa, prec) * self%h_max
        
        ! Number of cells per axis (at least 1)
        self%nCellsRow = max(1, int(dom_dim / self%cellSize))
        
        ! Total cells depends on nDim (from parameters)
        if (nDim == 2) then
            self%nCellsTot = self%nCellsRow**2
        else
            self%nCellsTot = self%nCellsRow**3
        end if

        if (allocated(self%head)) deallocate(self%head)
        allocate(self%head(self%nCellsTot))
        
        self%is_initialized = .true.
    end subroutine init_sort

    !> Reset the grid for the current time step
    subroutine build_spatial_grid(self)
        class(ParticleSort) :: self
        if (.not. self%is_initialized) return
        self%head = 0  ! Reset all cells (0 means empty)
    end subroutine build_spatial_grid

    !> Assign particles to cells using the Head-Next method
    !> Complexity: O(N), no dynamic allocation
    subroutine sort_particles(self, sys, dom_dim, next_array)
        class(ParticleSort) :: self
        type(ParticleSystem), intent(in) :: sys
        real(prec),           intent(in) :: dom_dim
        integer, allocatable, intent(out) :: next_array(:)

        integer :: p, i, idx
        ! Always dimension to 3 to avoid compile-time bounds errors in 3D logic
        integer    :: iCell(3) 
        real(prec) :: pos(3)

        if (allocated(next_array)) deallocate(next_array)
        allocate(next_array(sys%nPart))
        next_array = 0
        iCell = 1 ! Initialize default cell indices

        do p = 1, sys%nPart
            ! Copy current position (time slot 1)
            ! We only copy up to nDim, others remain 1
            pos(1:nDim) = sys%coords(1:nDim, 1, p)
            
            ! 1. Calculate grid coordinates (ix, iy, iz)
            do i = 1, nDim
                ! Use dom_dim to clamp positions and avoid out-of-grid indices
                ! This also resolves the "unused dummy argument" warning
                if (pos(i) < 0.0_prec) pos(i) = 0.0_prec
                if (pos(i) >= dom_dim) pos(i) = dom_dim - 0.0001_prec * self%cellSize
                
                iCell(i) = int(pos(i) / self%cellSize) + 1
                iCell(i) = min(max(1, iCell(i)), self%nCellsRow)
            end do
            
            ! 2. Compute linear cell index (Generic 2D/3D)
            if (nDim == 2) then
                idx = iCell(1) + (iCell(2) - 1) * self%nCellsRow
            else ! nDim == 3
                idx = iCell(1) + (iCell(2) - 1) * self%nCellsRow + &
                      (iCell(3) - 1) * (self%nCellsRow**2)
            end if
            
            ! 3. Update Head and Next (Linked-Cell core)
            next_array(p) = self%head(idx)
            self%head(idx) = p
        end do
    end subroutine sort_particles

end module sort_parts