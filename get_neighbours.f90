module get_neighbours
    use parameters
    use particles
    use sort_parts
    implicit none

contains

    !> Find neighbors for all particles using the Linked-Cell structure
    !> Optimized for SoA and memory contiguity
    subroutine compute_all_neighbors(sys, sort, next_array, kappa)
        type(ParticleSystem), intent(inout) :: sys
        type(ParticleSort),   intent(in)    :: sort
        integer,              intent(in)    :: next_array(:)
        integer,              intent(in)    :: kappa ! Kernel support factor
        
        integer :: i, p, neighbor_idx, current_pos
        integer :: ic, jc, kc, cell_idx
        integer :: iCell(3), nCell(3)
        real(prec) :: r2, search_dist2
        real(prec) :: dx(3)

        ! 1. Reset neighbor counts and pointers
        sys%neigh_count = 0
        sys%neigh_ptr   = 0
        current_pos     = 1 ! Position in the global sys%neigh_list

        ! 2. Loop over all particles to find their neighbors
        do p = 1, sys%nPart
            sys%neigh_ptr(p) = current_pos
            search_dist2 = (real(kappa, prec) * sys%h_part(p))**2
            
            ! Find current particle cell indices
            do i = 1, nDim
                iCell(i) = int(sys%coords(i, 1, p) / sort%cellSize) + 1
                iCell(i) = min(max(1, iCell(i)), sort%nCellsRow)
            end do
            if (nDim == 2) iCell(3) = 1 ! Safety for 2D

            ! 3. Check neighboring cells (-1 to +1 in each dimension)
            ! This triple loop is generic: kc will stay at 0 for 2D
            do kc = merge(-1, 0, nDim == 3), merge(1, 0, nDim == 3)
            do jc = -1, 1
            do ic = -1, 1
                
                nCell(1) = iCell(1) + ic
                nCell(2) = iCell(2) + jc
                nCell(3) = iCell(3) + kc
                
                ! Boundary check for the grid
                if (any(nCell(1:nDim) < 1) .or. any(nCell(1:nDim) > sort%nCellsRow)) cycle
                
                ! Compute neighbor cell linear index
                if (nDim == 2) then
                    cell_idx = nCell(1) + (nCell(2) - 1) * sort%nCellsRow
                else
                    cell_idx = nCell(1) + (nCell(2) - 1) * sort%nCellsRow + &
                               (nCell(3) - 1) * (sort%nCellsRow**2)
                end if
                
                ! 4. Iterate through particles in this cell using Head-Next
                neighbor_idx = sort%head(cell_idx)
                
                do while (neighbor_idx > 0)
                    ! Avoid self-interaction
                    if (neighbor_idx /= p) then
                        ! Square distance calculation (faster than sqrt)
                        r2 = 0.0_prec
                        do i = 1, nDim
                            dx(i) = sys%coords(i, 1, p) - sys%coords(i, 1, neighbor_idx)
                            r2 = r2 + dx(i)**2
                        end do
                        
                        ! 5. If within kernel support, add to list
                        if (r2 <= search_dist2) then
                            ! Check for neigh_list overflow
                            if (current_pos <= sys%max_neighbors) then
                                sys%neigh_list(current_pos) = neighbor_idx
                                sys%neigh_count(p) = sys%neigh_count(p) + 1
                                current_pos = current_pos + 1
                            end if
                        end if
                    end if
                    ! Move to next particle in the same cell
                    neighbor_idx = next_array(neighbor_idx)
                end do
                
            end do
            end do
            end do
        end do
        
    end subroutine compute_all_neighbors

end module get_neighbours