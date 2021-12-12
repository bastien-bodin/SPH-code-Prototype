module particles
    use parameters
    implicit none

    type :: Particle
        real(kind=prec), dimension(3,3)  :: coords                  ! coordinates of the particle
        real(kind=prec), dimension(3,3)  :: velocity                ! velocity of the particle
        real(kind=prec), dimension(1:3)  :: acceler                 ! acceleration of the particle
        real(kind=prec), dimension(1:3)  :: density            ! density of the particle
        real(kind=prec), dimension(1:3)  :: pressure           ! pressure at the particle
        real(kind=prec), dimension(1:3)  :: c_sound            ! sound speed
        real(kind=prec), dimension(1:3)  :: u_therm                 ! Thermal Energy
        real(kind=prec)                  :: mass = 1.               ! mass of the particle
        real(kind=prec)                  :: h_part = 0.5            ! smoothing length
        real(kind=prec)                  :: max_mu_ab               ! from the artificial viscosity
        real(kind=prec)                  :: molMass                 ! Molar Mass
        logical                          :: mobile = .false.        ! if true, mobile part -> solve eq. of motion
    end type Particle

    type :: PointPart
        type(Particle), pointer :: PTR => null()
    end type PointPart

    type :: ParticleList
        integer                          :: nb_elts = 0
        integer                          :: incr = 35
        integer                          :: max_nb_elts = 0
        type(PointPart), dimension(:), allocatable :: lst_parts ! list of particles
    contains
        procedure :: initList
        procedure :: addElt
        procedure :: resetList
    end type ParticleList

    type :: ParticleArray
        type(Particle) :: Part
        type(ParticleList) :: lst_neigh
    end type ParticleArray

    type :: PointPA
        type(ParticleArray), pointer :: PTR
    end type PointPA

    type :: ListPA
        integer                          :: nb_elts = 0
        integer                          :: incr = 35
        integer                          :: max_nb_elts = 0
        type(PointPA), dimension(:), allocatable :: lst_PA
    contains
        procedure :: initListPA
        procedure :: addEltPA
        procedure :: resetListPA
    end type ListPA

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Particle List
    subroutine initList(self)
        class(ParticleList) :: self

        if ( self%incr < 1 ) then
            self%incr = 1
        end if

        self%max_nb_elts = self%incr                    ! put the number of elements to incr
        allocate(self%lst_parts(1:self%max_nb_elts))    ! allocation
    end subroutine initList

    subroutine addElt(self, Part_Ptr)
        class(ParticleList) :: self
        type(Particle), pointer :: Part_Ptr

        type(PointPart), dimension(:), allocatable :: temp_lst

        if (self%nb_elts == 0) then
            call self%initList()
        else if (self%nb_elts == self%max_nb_elts) then
            allocate(temp_lst(1:self%max_nb_elts + self%incr))
            temp_lst(1:self%max_nb_elts) = self%lst_parts(1:self%max_nb_elts)
            call move_alloc(temp_lst,self%lst_parts)
            self%max_nb_elts = self%max_nb_elts + self%incr
        end if

        self%nb_elts = self%nb_elts + 1
        self%lst_parts(self%nb_elts)%PTR => Part_Ptr
    end subroutine addElt

    subroutine resetList(self)
        class(ParticleList) :: self

        self%nb_elts = 0
        deallocate(self%lst_parts)
    end subroutine resetList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> ParticleArray List
    subroutine initListPA(self)
        class(ListPA) :: self

        if ( self%incr < 1 ) then
            self%incr = 1
        end if

        self%max_nb_elts = self%incr                    ! put the number of elements to incr
        allocate(self%lst_PA(1:self%max_nb_elts))    ! allocation
    end subroutine initListPA

    subroutine addEltPA(self, PA_Ptr)
        class(ListPA) :: self
        type(ParticleArray), pointer :: PA_Ptr

        type(PointPA), dimension(:), allocatable :: temp_lst

        if ((self%nb_elts == 0) .and. (self%max_nb_elts ==0)) then
            call self%initListPA()
        else if (self%nb_elts == self%max_nb_elts) then
            allocate(temp_lst(1:self%max_nb_elts + self%incr))
            temp_lst(1:self%max_nb_elts) = self%lst_PA(1:self%max_nb_elts)
            call move_alloc(temp_lst,self%lst_PA)
            self%max_nb_elts = self%max_nb_elts + self%incr
        end if

        self%nb_elts = self%nb_elts + 1
        self%lst_PA(self%nb_elts)%PTR => PA_Ptr
    end subroutine addEltPA

    subroutine resetListPA(self)
        class(ListPA) :: self

        self%nb_elts = 0
        deallocate(self%lst_PA)
    end subroutine resetListPA
end module particles