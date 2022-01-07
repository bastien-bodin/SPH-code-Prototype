module sort_parts
    use parameters
    use particles
    implicit none

    type :: ParticleSort
        real(kind=prec)                             :: h_max
        real(kind=prec)                             :: CellSize
        integer                                     :: nCells_Row = 1
        integer                                     :: nCells     = 1
        class(ListPA), dimension(:), allocatable    :: LCells
        logical                                     :: init
    contains
        procedure :: get_h_max
        procedure :: setCells
        procedure :: particle_sort
    end type ParticleSort

contains
    subroutine get_h_max(self, lst_parts)
        class(ParticleSort)                         :: self
        type(ListPA)                                :: lst_parts

        integer                                     :: i

        self%h_max = 0.0d0

        do i =1,lst_parts%nb_elts
            if ( lst_parts%lst_PA(i)%PTR%Part%h_part > self%h_max ) then
                self%h_max = lst_parts%lst_PA(i)%PTR%Part%h_part
            end if
        end do
    end subroutine get_h_max

    subroutine setCells(self,dom_dim, kappa,LP)
        class(particleSort)                         :: self
        real(kind=prec), intent(in)                 :: dom_dim
        integer, intent(in)                         :: kappa
        class(ListPA), intent(in)                   :: LP

        call self%get_h_max(LP)

        self%nCells_Row = 0

        do while ( dom_dim / (self%nCells_Row+1) > kappa * self%h_max )
            self%nCells_Row = self%nCells_Row + 1
        end do
        !self%nCells_Row = self%nCells_Row +1


        self%nCells = self%nCells_Row**2 + 1
        self%CellSize = dom_dim / self%nCells_Row

        allocate(self%LCells(1:self%nCells))
    end subroutine setCells

    subroutine particle_sort(self, dom_dim, LP)
        class(ParticleSort) :: self
        real(kind=prec), intent(in) :: dom_dim
        type(ListPA), intent(in) :: LP


        integer :: i, j
        integer, dimension(1:2) :: iCell
        integer :: part_position
        real(kind=prec), dimension(1:2) :: xz
        integer :: numPart

        numPart = LP%nb_elts
        ! do i = 1, self%nCells
        !     if (self%LCells(i)%nb_elts /= 0) then
        !         call self%LCells(i)%resetListPA
        !     end if
        ! end do
        do i = 1, numPart
            xz = LP%lst_PA(i)%PTR%Part%coords(:,1)
            do j = 1,2
                if (xz(j) == dom_dim) then
                    iCell(j) = self%nCells_Row
                else
                    iCell(j) = nint( ( xz(j) - mod(xz(j), self%CellSize) )/self%CellSize ) + 1
                end if
            end do
            part_position = (iCell(1) - 1) * self%nCells_Row + iCell(2) 
            !print*,'nCells', self%nCells,'part pos:', part_position
            call self%LCells(part_position)%addEltPA(LP%lst_PA(i)%PTR)
        end do
    end subroutine particle_sort
end module sort_parts