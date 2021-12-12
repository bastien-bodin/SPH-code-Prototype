module application
    use parameters
    use particles
    use sort_parts
    use get_neighbours
    use kernels
    use equations
    implicit none

    type :: SPH_App
        type(ParticleSort)  :: srt ! sorting machine
        type(listPA) :: part ! array of pointers toward particles
        integer :: numFP ! number of fixed particles
        integer :: numMP ! number of mobile particles
        integer :: numPart ! total number of particles
        real(kind=prec) :: alpha ! weight. fact. in the art. viscosity formulation
        real(kind=prec) :: beta  ! weight. fact. in the art. viscosity formulation
        integer :: eqnstate ! -> 1 : ideal gas, -> 2 : WC fluid
        real(kind=prec) :: state_gamma ! power in EoS for WC fluids (often taken as 7)
        real(kind=prec) :: molMass ! molar mass for ideal gas
        logical :: kernelCorrection ! -> false : no correction, -> true : correction
        real(kind=prec) :: maxTime  ! simulation time (in seconds)
        real(kind=prec) :: saveInt  ! saving interval (in seconds)
        real(kind=prec) :: timeStep ! time step (adaptive, in seconds)
        real(kind=prec) :: currentTime ! current time (in seconds)
        integer :: neededStep ! steps needed for computation (Euler : 1, Verlet/RK : 2/4)
        real(kind=prec) :: h_0 ! initial smoothing length
        real(kind=prec) :: rho_0 ! density of the fluid at free surface
        real(kind=prec) :: c_0 ! sound speed in normal conditions
        real(kind=prec) :: dom_dim ! dimensions of the domains
        integer :: ite
        integer :: dim ! dimension of the problem (1D,2D or 3D)
    contains
        procedure, public :: initialisation
        procedure, public :: solver
        procedure, private :: timeStepUpdate
        procedure, private :: slUpdate
        procedure, private :: RK22
        procedure, private :: SyplVerlet
    end type SPH_App
    
contains
    subroutine initialisation(self)
        !! read param/fp/mp files and initialises the working environment
        class(SPH_App), target :: self
        character(len=250) :: param_path ! path of the parameters file
        character(len=250) :: fp_path    ! path of the fixed particles file
        character(len=250) :: mp_path    ! path of the mobile particles file
        real(kind=prec) :: x,y,z         ! coords of a particle
        real(kind=prec) :: v_x,v_y,v_z   ! velocity of a particle
        real(kind=prec) :: rho, m        ! density/mass of a particle
        real(kind=prec) :: u             ! thermal energy of a particle
        class(ParticleArray), pointer :: p_part ! pointer toward a particle

        integer :: i ! loop counter

        self%timeStep = 1E-15 ! initialise time step
        self%currentTime = 0.0d0 ! initialise current time

        !> reading paths for init files
        open(unit=1,file='paths.txt')
        read(1,*) param_path
        read(1,*) fp_path
        read(1,*) mp_path
        close(unit=1)

        !> reading and storing parameters for the simulation
        open(unit=2, file=trim(param_path))
        read(2,*) self%dim
        read(2,*) self%numFP
        read(2,*) self%numMP
        read(2,*) self%h_0
        read(2,*) self%c_0
        read(2,*) self%rho_0
        read(2,*) self%dom_dim
        read(2,*) self%alpha
        read(2,*) self%beta
        read(2,*) self%eqnstate
        read(2,*) self%state_gamma
        read(2,*) self%molMass
        read(2,*) self%kernelCorrection
        read(2,*) self%maxTime
        read(2,*) self%saveInt
        close(unit=2)

        !> allocate space for the list of particles
        self%numPart = self%numFP + self%numMP
        allocate(self%part%lst_PA(1:self%numPart))
        self%part%nb_elts = self%numPart

        !> read file for fixed particles
        open(unit=3,file=trim(fp_path))
        do i = 1, self%numFP
            read(3,*) x,y,z,v_x,v_y,v_z,rho,m,u
            allocate(self%part%lst_PA(i)%PTR)
            p_part => self%part%lst_PA(i)%PTR
            p_part%Part%coords = 0.0d0
            p_part%Part%coords(1,1) = x
            p_part%Part%coords(2,1) = y
            p_part%Part%coords(3,1) = z
            p_part%Part%velocity = 0.0d0
            p_part%Part%velocity(1,1) = v_x
            p_part%Part%velocity(2,1) = v_y
            p_part%Part%velocity(3,1) = v_z
            p_part%Part%density = 0.0d0
            p_part%Part%density(1) = rho
            p_part%Part%mass = m
            p_part%Part%h_part = self%h_0
            p_part%Part%u_therm = 0.0d0
            p_part%Part%u_therm(1) = u
            p_part%Part%pressure = 0.0d0
            p_part%Part%c_sound = 0.0d0
            if (self%eqnstate == 1) then
                call P_IG(p_part%Part,self%rho_0,self%molMass, 1)
                call CS_IG(p_part%Part, 1)
            else
                call P_WC(p_part%Part,self%c_0,self%rho_0, self%state_gamma, 1)
                call CS_WC(p_part%Part,self%rho_0, self%c_0, self%state_gamma, 1)
            end if
            p_part%Part%mobile = .false.
            !call p_part%lst_neigh%initList()
        end do
        close(unit=3)

        !> read file for mobile particles
        open(unit=4,file=trim(mp_path))
        do i = 1, self%numMP
            read(4,*) x,y,z,v_x,v_y,v_z,rho,m,u
            allocate(self%part%lst_PA(self%numFP+i)%PTR)
            p_part => self%part%lst_PA(self%numFP+i)%PTR
            p_part%Part%coords = 0.0d0
            p_part%Part%coords(1,1) = x
            p_part%Part%coords(2,1) = y
            p_part%Part%coords(3,1) = z
            p_part%Part%velocity = 0.0d0
            p_part%Part%velocity(1,1) = v_x
            p_part%Part%velocity(2,1) = v_y
            p_part%Part%velocity(3,1) = v_z
            p_part%Part%density = 0.0d0
            p_part%Part%density(1) = rho
            p_part%Part%mass = m
            p_part%Part%h_part = self%h_0
            p_part%Part%u_therm = 0.0d0
            p_part%Part%u_therm(1) = u
            p_part%Part%pressure = 0.0d0
            p_part%Part%c_sound = 0.0d0
            if (self%eqnstate == 1) then
                call P_IG(p_part%Part,self%rho_0,self%molMass, 1)
                call CS_IG(p_part%Part, 1)
            else
                call P_WC(p_part%Part,self%c_0,self%rho_0, self%state_gamma, 1)
                call CS_WC(p_part%Part,self%rho_0, self%c_0, self%state_gamma, 1)
            end if
            p_part%Part%mobile = .true.
            !call p_part%lst_neigh%initList()
        end do
        close(unit=4)
        
        !call self%srt%setCells(self%dom_dim,2,self%part)
        self%ite = 0

        open(unit=24,file='log.txt',form='formatted',access='stream',status='replace')
        write(24,*) 'Initialisation finished'
        print*, 'Initialisation finished'
    end subroutine initialisation

    subroutine solver(self)
        !! solves the problem, given the initial setup
        class(SPH_App) :: self

        integer :: i,j
        integer :: ite      ! iteration counter
        logical :: to_save  ! saving flag; if true -> save

        ite = 0

        open(unit=22, file='result.out', form='unformatted',status='replace')! ,access='stream'
        open(unit=23, file='time.out', form='unformatted',status='replace') ! , access='stream'

        !> time increment + saving status
        do while (self%currentTime <= self%maxTime)
            if ((floor(self%currentTime/self%saveInt) /= floor((self%currentTime+self%timeStep)/self%saveInt)) .or. ite == 0) then
                to_save = .true.
            else
                to_save = .false.
            end if
            self%currentTime = self%currentTime + self%timeStep
            !> update X,V,rho,u,cs,p
            self%neededStep = 2
            do j=1,self%neededStep
                call self%srt%setCells(self%dom_dim,1,self%part)
                call self%srt%particle_sort(self%dom_dim,self%part)
                do i = 1, self%numPart
                    ! solving method (Euler, Verlet, RK)
                    !call RK22(self,self%part%lst_PA(i)%PTR,j)
                    call SyplVerlet(self,self%part%lst_PA(i)%PTR,j)
                end do
                do i=1, self%srt%nCells
                    if (self%srt%LCells(i)%nb_elts /= 0) then
                        call self%srt%LCells(i)%resetListPA()
                    end if
                end do
                deallocate(self%srt%LCells)
            end do

            !> update the current time variables (current time = next time)
            do i=1,self%numPart
                self%part%lst_PA(i)%PTR%Part%density(1) = self%part%lst_PA(i)%PTR%Part%density(self%neededStep+1)
                self%part%lst_PA(i)%PTR%Part%pressure(1) = self%part%lst_PA(i)%PTR%Part%pressure(self%neededStep+1)
                self%part%lst_PA(i)%PTR%Part%c_sound(1) = self%part%lst_PA(i)%PTR%Part%c_sound(self%neededStep+1)
                self%part%lst_PA(i)%PTR%Part%u_therm(1) = self%part%lst_PA(i)%PTR%Part%u_therm(self%neededStep+1)
                self%part%lst_PA(i)%PTR%Part%coords(:,1) = self%part%lst_PA(i)%PTR%Part%coords(:,self%neededStep+1)
                self%part%lst_PA(i)%PTR%Part%velocity(:,1) = self%part%lst_PA(i)%PTR%Part%velocity(:,self%neededStep+1)
            end do



            !> saving data at given intervals
            if (to_save) then
                !> saving parts
                do i = 1,self%numMP
                    write(unit=22) real(self%part%lst_PA(self%numFP+i)%PTR%Part%coords(1,1),4), &
                    & real(self%part%lst_PA(self%numFP+i)%PTR%Part%coords(2,1),4), &
                    & real(self%part%lst_PA(self%numFP+i)%PTR%Part%coords(3,1),4), &
                    & real(self%part%lst_PA(self%numFP+i)%PTR%Part%velocity(1,1),4), &
                    & real(self%part%lst_PA(self%numFP+i)%PTR%Part%velocity(2,1),4), &
                    &  real(self%part%lst_PA(self%numFP+i)%PTR%Part%velocity(3,1),4), &
                    & real(self%part%lst_PA(self%numFP+i)%PTR%Part%pressure(1),4), &
                    & real(self%part%lst_PA(self%numFP+i)%PTR%Part%density(1),4), &
                    & real(self%part%lst_PA(self%numFP+i)%PTR%Part%u_therm(1),4)
                end do
                !> saving time
                write(unit=23) real(self%currentTime,4)
                print*, 'Iteration nb:   ',self%ite
                print*, 'Time (s) =      ',self%currentTime
                print*, 'Time step (s) = ',self%timeStep
                write(24,*) 'Iteration nb:   ',self%ite
                write(24,*) 'Time (s) =      ',self%currentTime
                write(24,*) 'Time step (s) = ',self%timeStep
            end if
            call self%timeStepUpdate
            !call self%slUpdate
            self%ite = self%ite +1
        end do

        close(22)
        close(23)
        close(24)

        call self%part%resetListPA()
    end subroutine solver

    subroutine timeStepUpdate(self)
        !> compute the adaptive time step (c.f. dualSPHysics)
        class(SPH_App) :: self

        real(kind=prec) :: dt_f, dt_ftemp
        real(kind=prec) :: dt_cv, dt_cvtemp
        integer :: i
        class(ParticleArray), pointer :: p_part
        real(kind=prec), dimension(1:3) :: a !  acceleration
        real(kind=prec) :: f_m ! F/M -> acceleration norm
        real(kind=prec) :: h
        real(kind=prec) :: c_s, max_mu

        a = self%part%lst_PA(self%numFP+1)%PTR%Part%acceler(:)
        h = self%part%lst_PA(self%numFP+1)%PTR%Part%h_part
        f_m = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
        c_s = self%part%lst_PA(self%numFP+1)%PTR%Part%c_sound(1)
        max_mu = self%part%lst_PA(self%numFP+1)%PTR%Part%max_mu_ab

        !> compute the time_step
        dt_f = sqrt(h/f_m)
        dt_cv = h/(c_s+max_mu)
        do i = self%numFP+2, self%numPart
            p_part => self%part%lst_PA(i)%PTR
            h = p_part%Part%h_part
            a = p_part%Part%acceler(:)
            f_m = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
            dt_ftemp = sqrt(h/f_m)
            if (dt_ftemp < dt_f) then
                dt_f = dt_ftemp
            end if
            c_s = p_part%Part%c_sound(1)
            max_mu = p_part%Part%max_mu_ab
            dt_cvtemp = h/(c_s+max_mu)
            if (dt_cvtemp< dt_cv) then
                dt_cv = dt_cvtemp
            end if
        end do
        self%timeStep = 0.3 * min(dt_f,dt_cv)
    end subroutine timeStepUpdate

    subroutine slUpdate(self)
        !> computes the new smoothing length at each time step (same h for every particle)
        class(SPH_App) :: self

        real(kind=prec) :: rho_mean
        real(kind=prec) :: h_new
        integer :: i

        rho_mean = 0.0d0

        do i=1,self%numPart
            rho_mean = rho_mean + self%part%lst_PA(i)%PTR%Part%density(1)
        end do
        rho_mean = rho_mean/self%numPart

        h_new = self%h_0*(self%rho_0/rho_mean)**(1.0d0/self%dim)
        if (h_new > 0.5d0 * self%srt%CellSize) then
            h_new = 0.5d0 * self%srt%CellSize
            print*, 'Warning: the new h has been limited'
            write(24,*) 'Warning: the new h has been limited'
        end if
        do i=1,self%numPart
            self%part%lst_PA(i)%PTR%Part%h_part = h_new
        end do
    end subroutine slUpdate

    subroutine RK22(self, PartArray, integStep)
        !> Updates time dependent variables with the RK22 time integration
        class(SPH_App) :: self
        class(ParticleArray) :: PartArray
        integer :: integStep

        type(Quintic) :: KernL
        real(kind=prec) :: r_ab

        real(kind=prec), dimension(:,:), allocatable :: vec_gradW

        real(kind=prec) :: D_rho
        real(kind=prec) :: D_u
        real(kind=prec), dimension(1:3) :: D_v
        real(kind=prec), dimension(1:3) :: D_x
        real(kind=prec), dimension(1:3) :: v_ab
        real(kind=prec), dimension(1:3) :: F_m
        real(kind=prec) :: pi_ab
        integer :: i
        class(Particle), pointer :: p_neigh

        D_rho = 0.0d0
        D_u = 0.0d0
        D_v = (/0.0d0,0.0d0,0.0d0/)
        D_x = (/0.0d0,0.0d0,0.0d0/)
        F_m = (/real(0,prec),real(0,prec),real(-9.81,prec)/)

        call KernL%init_quint(PartArray%Part%h_part,self%dim)

        call getNeighbours(PartArray,self%srt,self%dom_dim,KernL%kappa, integStep)

        allocate(vec_gradW(1:3,1:PartArray%lst_neigh%nb_elts))
        do i = 1, PartArray%lst_neigh%nb_elts
            p_neigh => PartArray%lst_neigh%lst_parts(i)%PTR
            r_ab = eval_r(PartArray%Part%coords(:,integStep),p_neigh%coords(:,integStep))
            KernL%r_ab = r_ab
            KernL%q = KernL%r_ab/KernL%h
            call KernL%kernel_quint()
            call KernL%dwdq_quint()
            vec_gradW(:,i) = (PartArray%Part%coords(:,integStep) - p_neigh%coords(:,integStep))/r_ab * KernL%dw
        end do
        if (integStep == 1) then
            PartArray%Part%max_mu_ab = 0.0d0
        end if

        do i = 1, PartArray%lst_neigh%nb_elts
            p_neigh => PartArray%lst_neigh%lst_parts(i)%PTR
            v_ab(:) = PartArray%Part%velocity(:,integStep) - p_neigh%velocity(:,integStep)
            pi_ab = ArtVisc(PartArray%Part,p_neigh, self%alpha, self%beta, integStep)
            D_rho = D_rho + p_neigh%mass * dot_product(v_ab,vec_gradW(:,i))
            if (PartArray%Part%mobile) then
                D_v = D_v + PartArray%Part%mass * &
                    & (p_neigh%pressure(integStep)/(p_neigh%density(integStep)**2.0d0) + &
                    & PartArray%Part%pressure/(PartArray%Part%density(integStep)**2.0d0) + &
                    & pi_ab) * vec_gradW(:,i)
            end if
            D_u = D_u + 0.5d0 * PartArray%Part%mass * &
                & (p_neigh%pressure(integStep)/(p_neigh%density(integStep)**2.0d0) + &
                & PartArray%Part%pressure(integStep)/(PartArray%Part%density(integStep)**2.0d0) + &
                & pi_ab) * dot_product(v_ab,vec_gradW(:,i))
        end do
        if (PartArray%Part%mobile) then
            D_v = -D_v + F_m
        end if
        if (integStep == 1) then
            PartArray%Part%density(2) = PartArray%Part%density(1) + D_rho * self%timeStep
            PartArray%Part%density(3) = PartArray%Part%density(1) + D_rho * self%timeStep/2.0d0
            PartArray%Part%u_therm(2) = PartArray%Part%u_therm(1) + D_u * self%timeStep
            PartArray%Part%u_therm(3) = PartArray%Part%u_therm(1) + D_u * self%timeStep/2.0d0
            if (PartArray%Part%mobile) then
                PartArray%Part%acceler(:) = D_v
                PartArray%Part%velocity(:,2) = PartArray%Part%velocity(:,1) + D_v * self%timeStep
                PartArray%Part%velocity(:,3) = PartArray%Part%velocity(:,1) + D_v * self%timeStep/2.0d0
                D_x = PartArray%Part%velocity(:,1)
                PartArray%Part%coords(:,2) = PartArray%Part%coords(:,1) + D_x * self%timeStep
                PartArray%Part%coords(:,3) = PartArray%Part%coords(:,1) + D_x * self%timeStep/2.0d0
            else
                PartArray%Part%velocity(:,2) = PartArray%Part%velocity(:,1)
                PartArray%Part%coords(:,2) = PartArray%Part%coords(:,1)
            end if
        else
            PartArray%Part%density(3) = PartArray%Part%density(3) + D_rho * self%timeStep/2.0d0
            PartArray%Part%u_therm(3) = PartArray%Part%u_therm(3) + D_u * self%timeStep/2.0d0
            if (PartArray%Part%mobile) then
                PartArray%Part%acceler(:) = (D_v + PartArray%Part%acceler(:))/2.0d0
                PartArray%Part%velocity(:,3) = PartArray%Part%velocity(:,3) + D_v * self%timeStep/2.0d0
                D_x = PartArray%Part%velocity(:,2)
                PartArray%Part%coords(:,3) = PartArray%Part%coords(:,3) + D_x * self%timeStep/2.0d0
            else
                PartArray%Part%velocity(:,3) = PartArray%Part%velocity(:,2)
                PartArray%Part%coords(:,3) = PartArray%Part%coords(:,2)
            end if
            
        end if

        call P_WC(PartArray%Part,self%c_0,self%rho_0,self%state_gamma, integStep)
        call CS_WC(PartArray%Part,self%rho_0,self%c_0,self%state_gamma, integStep)

        if (PartArray%lst_neigh%nb_elts /= 0) then
            call PartArray%lst_neigh%resetList()
        end if
        deallocate(vec_gradW)

    end subroutine RK22

    subroutine SyplVerlet(self, PartArray, integStep)
        !> Updates time dependent variables with the Symplectic Verlet algorithm
        class(SPH_App) :: self
        class(ParticleArray) :: PartArray
        integer :: integStep

        type(Quintic) :: KernL
        real(kind=prec) :: r_ab

        real(kind=prec), dimension(:,:), allocatable :: vec_gradW

        real(kind=prec) :: D_rho
        real(kind=prec) :: D_u
        real(kind=prec), dimension(1:3) :: D_v
        real(kind=prec), dimension(1:3) :: D_x
        real(kind=prec), dimension(1:3) :: v_ab
        real(kind=prec), dimension(1:3) :: F_m
        real(kind=prec) :: pi_ab
        integer :: i
        class(Particle), pointer :: p_neigh

        D_rho = 0.0d0
        D_u = 0.0d0
        D_v = (/0.0d0,0.0d0,0.0d0/)
        D_x = (/0.0d0,0.0d0,0.0d0/)
        F_m = (/real(0,prec),real(0,prec),real(-9.81,prec)/)

        call KernL%init_quint(PartArray%Part%h_part,self%dim)

        call getNeighbours(PartArray,self%srt,self%dom_dim,KernL%kappa, integStep)

        allocate(vec_gradW(1:3,1:PartArray%lst_neigh%nb_elts))
        do i = 1, PartArray%lst_neigh%nb_elts
            p_neigh => PartArray%lst_neigh%lst_parts(i)%PTR
            r_ab = eval_r(PartArray%Part%coords(:,integStep),p_neigh%coords(:,integStep))
            KernL%r_ab = r_ab
            KernL%q = KernL%r_ab/KernL%h
            call KernL%kernel_quint()
            call KernL%dwdq_quint()
            vec_gradW(:,i) = (PartArray%Part%coords(:,integStep) - p_neigh%coords(:,integStep))/r_ab * KernL%dw
        end do
        
        if (integStep == 1) then
            PartArray%Part%max_mu_ab = 0.0d0
        end if

        do i = 1, PartArray%lst_neigh%nb_elts
            p_neigh => PartArray%lst_neigh%lst_parts(i)%PTR
            v_ab(:) = PartArray%Part%velocity(:,integStep) - p_neigh%velocity(:,integStep)
            pi_ab = ArtVisc(PartArray%Part,p_neigh, self%alpha, self%beta, integStep)
            D_rho = D_rho + p_neigh%mass * dot_product(v_ab,vec_gradW(:,i))
            if (PartArray%Part%mobile) then
                D_v = D_v + PartArray%Part%mass * &
                    & (p_neigh%pressure(integStep)/(p_neigh%density(integStep)**2.0d0) + &
                    & PartArray%Part%pressure/(PartArray%Part%density(integStep)**2.0d0) + &
                    & pi_ab) * vec_gradW(:,i)
            end if
            D_u = D_u + 0.5d0 * PartArray%Part%mass * &
                & (p_neigh%pressure(integStep)/(p_neigh%density(integStep)**2.0d0) + &
                & PartArray%Part%pressure(integStep)/(PartArray%Part%density(integStep)**2.0d0) + &
                & pi_ab) * dot_product(v_ab,vec_gradW(:,i))
        end do
        if (PartArray%Part%mobile) then
            D_v = -D_v + F_m
        end if

        if (integStep == 1) then
            !> half time step
            PartArray%Part%density(2)  = PartArray%Part%density(1) + self%timeStep/2.0d0 * D_rho
            PartArray%Part%u_therm(2)  = PartArray%Part%u_therm(1) + self%timeStep/2.0d0 * D_u
            if (PartArray%Part%mobile) then
                PartArray%Part%acceler(:) = D_v
                PartArray%Part%velocity(:,2) = PartArray%Part%velocity(:,1) + self%timeStep/2.0d0 * D_v
                PartArray%Part%coords(:,2) = PartArray%Part%coords(:,1) + self%timeStep/2.0d0 * PartArray%Part%velocity(:,1)
            else
                PartArray%Part%velocity(:,2) = PartArray%Part%velocity(:,1)
                PartArray%Part%coords(:,2) = PartArray%Part%coords(:,1)
            end if
        else
            !> full time step
            PartArray%Part%density(3) = PartArray%Part%density(1)*(2.0d0 + (D_rho/PartArray%Part%density(2))*self%timeStep )/ &
                            & (2.0d0 - (D_rho/PartArray%Part%density(2))*self%timeStep )

            PartArray%Part%u_therm(3) = PartArray%Part%u_therm(2) + self%timeStep * D_u/2.0d0

            if (PartArray%Part%mobile) then
                PartArray%Part%acceler(:) = D_v
                PartArray%Part%velocity(:,3) = PartArray%Part%velocity(:,1) + self%timeStep * D_v
                PartArray%Part%coords(:,3) = PartArray%Part%coords(:,1) &
                    & + self%timeStep * (PartArray%Part%velocity(:,3) + PartArray%Part%velocity(:,1))/2.0d0
            else
                PartArray%Part%velocity(:,3) = PartArray%Part%velocity(:,1)
                PartArray%Part%coords(:,3) = PartArray%Part%coords(:,1)
            end if
        end if

        call P_WC(PartArray%Part,self%c_0,self%rho_0,self%state_gamma, integStep)
        call CS_WC(PartArray%Part,self%rho_0,self%c_0,self%state_gamma, integStep)

        if (PartArray%lst_neigh%nb_elts /= 0) then
            call PartArray%lst_neigh%resetList()
        end if
        deallocate(vec_gradW)
    end subroutine SyplVerlet

    !> Applies a correction to the kernel gradient
    subroutine kernelCor(PA, gradW, gradW_mod, integStep)
        implicit none
        class(ParticleArray) :: PA

        real(kind=prec), dimension(:,:) :: gradW
        real(kind=prec), dimension(:,:) :: gradW_mod

        real(kind=prec), dimension(3,3) :: M
        real(kind=prec), dimension(3,3) :: L
        integer :: i,j,k
        real(kind=prec) :: detM
        real(kind=prec) :: MDivRho
        integer :: integStep
        class(Particle), pointer :: p_neigh

        M(1,1) = 0.0d0
        M(2,2) = 0.0d0
        M(3,3) = 0.0d0
        M(1,2) = 0.0d0
        M(1,3) = 0.0d0
        M(2,3) = 0.0d0

        do i = 1, PA%lst_neigh%nb_elts
            p_neigh => PA%lst_neigh%lst_parts(i)%PTR
            MDivRho = p_neigh%mass / p_neigh%density(integStep)
            M(1,1) = M(1,1) + MDivRho * (p_neigh%coords(1,integStep) - PA%Part%coords(1,integStep)) * gradW(1,i)
            M(2,2) = M(2,2) + MDivRho * (p_neigh%coords(2,integStep) - PA%Part%coords(2,integStep)) * gradW(2,i)
            M(3,3) = M(3,3) + MDivRho * (p_neigh%coords(3,integStep) - PA%Part%coords(3,integStep)) * gradW(3,i)
            M(1,2) = M(1,1) + MDivRho * (p_neigh%coords(1,integStep) - PA%Part%coords(1,integStep)) * gradW(2,i)
            M(1,3) = M(1,1) + MDivRho * (p_neigh%coords(1,integStep) - PA%Part%coords(1,integStep)) * gradW(3,i)
            M(2,3) = M(1,1) + MDivRho * (p_neigh%coords(2,integStep) - PA%Part%coords(2,integStep)) * gradW(3,i)
        end do

            !> M is symmetric
        M(2,1) = M(1,2)
        M(3,1) = M(1,3)
        M(3,2) = M(2,3)

        !> computes the determinant of M
        detM = M(1,1) * (M(2,2)*M(3,3)-M(3,2)*M(2,3)) &
            &- M(1,2) * (M(2,1)*M(3,3)-M(3,1)*M(2,3)) &
            &+ M(1,3) * (M(2,1)*M(3,2)-M(3,1)*M(2,3))
        
        !> compute M^-1 = L
        L(1,1) = M(2,2)*M(3,3)-M(3,2)*M(2,3)
        L(2,2) = M(1,1)*M(3,3)-M(3,1)*M(1,3)
        L(3,3) = M(1,1)*M(2,2)-M(2,1)*M(1,2)
        L(1,2) = M(3,1)*M(2,3)-M(2,1)*M(3,3)
        L(1,3) = M(2,1)*M(3,2)-M(3,1)*M(2,2)
        L(2,3) = M(3,1)*M(1,2)-M(1,1)*M(3,2)

            !> L is symmetric, because M is symmetric
        L(2,1) = L(1,2)
        L(3,1) = L(1,3)
        L(3,2) = L(2,3)

        L = (1.0d0/detM) * L

        do i=1, PA%lst_neigh%nb_elts
            gradW_mod(1,i) = 0.0d0
            gradW_mod(2,i) = 0.0d0
            gradW_mod(2,i) = 0.0d0
            do j=1,3
                do k=1,3
                    gradW_mod(j,i) = gradW_mod(j,i) + L(j,k) * gradW(k,i)
                end do
            end do
        end do
    
    end subroutine kernelCor
end module application