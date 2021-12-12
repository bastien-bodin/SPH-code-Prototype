module get_neighbours
    use parameters
    use particles
    use sort_parts
    use equations, only : eval_r
    implicit none

contains
    subroutine getNeighbours(PA, sort,dom_dim,kappa, integStep)

        implicit none

        ! I/O parameters
        type(ParticleArray), target                :: PA
        type(ParticleSort), intent(in)             :: sort
        real(kind=prec), intent(in)                :: dom_dim               ! length of a side of the domain
        integer, intent(in)                        :: kappa                 ! from the kernel

        integer :: integStep

        ! Internal parameters
        integer, dimension(1:3)                            :: iCell         ! Identifies the cell of the particle
        integer                                            :: i,j,k         ! loop counter
        integer, dimension(1:27)                           :: cellsToCheck  ! number of the cells to check for neighbours
        real(kind=prec), dimension(:), pointer, contiguous :: part_xyz      ! position of the particle
        real(kind=prec), dimension(:), pointer, contiguous :: neig_xyz      ! position of a neighbour
        class(Particle), pointer                     :: cur_ptr             ! pointer toward a particle
        type(ParticleList), pointer                        :: cur_nei       ! current list of neighbours
        real(kind=prec)                                    :: r             ! distance


        part_xyz => PA%Part%coords(:,IntegStep)

        !> compute the location of the particle in the cells
        do j = 1,3
            if (part_xyz(j) == dom_dim) then
                iCell(j) = sort%nCells_Row
            else
                iCell(j) = nint( ( part_xyz(j) - mod(part_xyz(j), sort%CellSize) )/sort%CellSize ) + 1
            end if
        end do

        !> compute the number of the neighbouring cells
        do i=-1,1
            do j=-1,1
                do k=-1,1
                    if ( (iCell(1) + i > 0) .and. (iCell(2) + j > 0) .and. (iCell(3) + k > 0) &
                    & .and. &
                    & (iCell(1) + i <= sort%nCells_Row) .and. (iCell(2) + j <= sort%nCells_Row) &
                    & .and. (iCell(3) + k <= sort%nCells_Row)) then
                        cellsToCheck((i+1) * 9 + (j+1) * 3 + (k + 2)) = &
                        & ((iCell(1) + i - 1) * sort%nCells_Row + (iCell(2) + j - 1)) * sort%nCells_Row + (iCell(3) + k)
                    else
                        cellsToCheck((i+1) * 9 + (j+1) * 3 + (k + 2)) = 0
                    end if
                end do
            end do
        end do


        !> stores the neighbours of the particle in Neigh_List
        cur_nei => PA%lst_neigh
        if (cur_nei%nb_elts /= 0) then
            call cur_nei%resetList()
        end if
        do i = 1,27
            if (cellsToCheck(i) > 0) then
                ! storage => sort%LCells(cellsToCheck(i))
                do j = 1, sort%LCells(cellsToCheck(i))%nb_elts
                    cur_ptr => sort%LCells(cellsToCheck(i))%lst_PA(j)%PTR%Part
                    neig_xyz => cur_ptr%coords(:,integStep)
                    r = eval_r(part_xyz,neig_xyz)
                    if ( (r <= kappa*PA%Part%h_part) .and. (r > 1E-12) ) then
                        call cur_nei%addElt(cur_ptr)
                    end if
                end do
            end if
        end do

    end subroutine getNeighbours
end module get_neighbours