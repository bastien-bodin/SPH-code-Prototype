module kernels
    use parameters
    implicit none

    type :: kernel
        integer         :: dim    ! dimension of the problem
        integer         :: kappa  ! set the radius scale
        real(kind=prec) :: r_ab
        real(kind=prec) :: h      ! smoothing length
        real(kind=prec) :: q      ! q = r/h
        real(kind=prec) :: W      ! kernel
        real(kind=prec) :: dw     ! derivative of the kernel
    end type kernel

    type, extends(kernel) :: Gaussian
    contains
        procedure :: init_g
        procedure :: kernel_g
        procedure :: dwdq_g
    end type gaussian

    type, extends(kernel) :: Quartic
    contains
        procedure :: init_quart
        procedure :: kernel_quart
        procedure :: dwdq_quart
    end type Quartic

    type, extends(kernel) :: Cubic
    contains
        procedure :: init_cub
        procedure :: kernel_cub
        procedure :: dwdq_cub
    end type Cubic

    type, extends(kernel) :: Quadratic
    contains
        procedure :: init_quadra
        procedure :: kernel_quadra
        procedure :: dwdq_quadra
    end type Quadratic

    type, extends(kernel) :: Quintic
    contains
        procedure :: init_quint
        procedure :: kernel_quint
        procedure :: dwdq_quint
    end type Quintic

    type, extends(kernel) :: QuinticSpline
    contains
        procedure :: init_QS
        procedure :: kernel_QS
        procedure :: dwdq_QS
    end type QuinticSpline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Procedures for Gaussian kernel
    subroutine init_g(self, h, dim)
        class(Gaussian), intent(inout) :: self

        real(kind=prec), intent(in)    :: h
        integer        , intent(in)    :: dim

        self%h    = h
        self%dim  = dim

        self%kappa = 1

    end subroutine init_g

    subroutine kernel_g(self)
        class(Gaussian):: self
        real(kind=prec) :: a_d
        a_d = 1.0d0/(PI**(self%dim/2) * self%h**self%dim)

        !> compute W(r,h)
        self%W = a_d * EXP(-self%q**2)
    end subroutine kernel_g

    subroutine dwdq_g(self)
        class(Gaussian) :: self
        real(kind=prec) :: a_d

        a_d = 1.0d0/(PI**(self%dim/2) * self%h**self%dim)

        !> compute dW/dr
        self%dW = -2.0d0 * a_d/self%h * self%q * EXP(-self%q**2)
    end subroutine dwdq_g


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Procedures for Quartic (Bell-shaped) kernel
    subroutine init_quart(self, h, dim)
        class(Quartic), intent(inout) :: self

        real(kind=prec), intent(in)    :: h
        integer        , intent(in)    :: dim

        self%h    = h
        self%dim  = dim

        self%kappa = 2

    end subroutine init_quart

    subroutine kernel_quart(self)
        class(Quartic) :: self
        real(kind=prec) :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 5.0d0/(4.0d0 * self%h)
        else if (self%dim == 2) then
            a_d = 5.0d0/(PI * self%h**2)
        else
            a_d = 105.0d0/(16.0d0 * PI * self%h)
        end if

        !> compute W(r,h)
        if (self%q <= 1) then
            self%W = a_d * (1 + 3 * self%q) * (1 - self%q)**3
        else
            self%W = 0
        end if
    end subroutine kernel_quart

    subroutine dwdq_quart(self)
        class(Quartic)  :: self
        real(kind=prec) :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 5.0d0/(4.0d0 * self%h)
        else if (self%dim == 2) then
            a_d = 5.0d0/(PI * self%h**2)
        else
            a_d = 105.0d0/(16.0d0 * PI * self%h**3)
        end if

        !> compute dW/dr
        if (self%q <= 1) then
            self%dW = a_d/self%h * 3 * ((1 - self%q)**3 - (1 + 3 * self%q) * (1 - self%q)**2)
        else
            self%dW = 0
        end if
    end subroutine dwdq_quart


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Procedures for Cubic spline kernel
    subroutine init_cub(self, h, dim)
        class(Cubic), intent(inout) :: self
        real(kind=prec), intent(in)    :: h
        integer        , intent(in)    :: dim

        self%h    = h
        self%dim  = dim

        self%kappa = 2

    end subroutine init_cub

    subroutine kernel_cub(self)
        class(Cubic)    :: self
        real(kind=prec) :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 1.0d0/self%h
        else if (self%dim == 2) then
            a_d = 15.0d0/(7 * PI * self%h**2)
        else
            a_d = 3.0d0/(2.0d0 * PI * self%h**3)
        end if

        !> compute W(r,h)
        if (self%q < 1) then
            self%dW = a_d * (1.5d0 - self%q**2 + 0.5d0 * self%q**3)
        else if (self%q <2) then
            self%dW = a_d * 1.0d0/6.0d0 * (2 - self%q)**3
        else
            self%dW = 0
        end if
    end subroutine kernel_cub

    subroutine dwdq_cub(self)
        class(Cubic)    :: self
        real(kind=prec) :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 1.0d0/self%h
        else if (self%dim == 2) then
            a_d = 15.0d0/(7 * PI * self%h**2)
        else
            a_d = 3.0d0/(2.0d0 * PI * self%h**3)
        end if

        !> compute dW/dr
        if (self%q < 1) then
            self%dW = a_d/self%h * (1.5d0 * self%q**2 - 2 * self%q)
        else if (self%q <2) then
            self%dW = - a_d/self%h * 0.5d0 * (2 - self%q)**2
        else
            self%dW = 0
        end if
    end subroutine dwdq_cub


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Procedures for Quadratic kernel
    subroutine init_quadra(self, h, dim)
        class(Quadratic), intent(inout) :: self
        real(kind=prec), intent(in)    :: h
        integer        , intent(in)    :: dim

        self%h    = h
        self%dim  = dim

        self%kappa = 2

    end subroutine init_quadra

    subroutine kernel_quadra(self)
        class(Quadratic) :: self
        real(kind=prec)  :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 1.0d0/self%h
        else if (self%dim == 2) then
            a_d = 2.0d0/(PI * self%h**2)
        else
            a_d = 5.0d0/(4.0d0 * PI * self%h**3)
        end if

        !> compute W(r,h)
        if (self%q <=2) then
            self%W = a_d * (3.0d0/16.0d0 * self%q**2 - 0.75d0 * self%q + 0.75d0)
        else
            self%W = 0
        end if
    end subroutine kernel_quadra

    subroutine dwdq_quadra(self)
        class(Quadratic) :: self
        real(kind=prec)  :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 1.0d0/self%h
        else if (self%dim == 2) then
            a_d = 2.0d0/(PI * self%h**2)
        else
            a_d = 5.0d0/(4.0d0 * PI * self%h**3)
        end if

        !> compute dW/dr
        if (self%q <=2) then
            self%dW = a_d/self%h * (3.0d0/8.0d0 * self%q - 0.75d0)
        else
            self%dW = 0
        end if
    end subroutine dwdq_quadra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Procedures for Quintic kernel
    subroutine init_quint(self, h, dim)
        class(Quintic), intent(inout) :: self

        real(kind=prec), intent(in)    :: h
        integer        , intent(in)    :: dim

        self%h    = h
        self%dim  = dim

        self%kappa = 2

    end subroutine init_quint

    subroutine kernel_quint(self)
        class(Quintic)  :: self
        real(kind=prec) :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 0.75d0/self%h
        else if (self%dim == 2) then
            a_d = 7.0d0/(4.0d0 * PI * self%h**2)
        else
            a_d = 21.0d0/(16.0d0 * PI * self%h**3)
        end if

        !> compute W(r,h)
        if (self%q <=2) then
            self%W = a_d * ((1 - 0.5d0 * self%q)**4 * (2 * self%q + 1))
        else
            self%W = 0
        end if
    end subroutine kernel_quint

    subroutine dwdq_quint(self)
        class(Quintic)  :: self
        real(kind=prec) :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 1.0d0/self%h
        else if (self%dim == 2) then
            a_d = 2.0d0/(PI * self%h**2)
        else
            a_d = 5.0d0/(4.0d0 * PI * self%h**3)
        end if

        !> compute dW/dr
        if (self%q <=2) then
            self%dW = -5.0d0 * a_d/self%h * self%q * (1 - 0.5d0 * self%q)**3
        else
            self%dW = 0
        end if
    end subroutine dwdq_quint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Procedures for Quintic Spline kernel
    subroutine init_QS(self, h, dim)
        class(QuinticSpline), intent(inout) :: self

        real(kind=prec), intent(in)    :: h
        integer        , intent(in)    :: dim

        self%h    = h
        self%dim  = dim

        self%kappa = 3

    end subroutine init_QS

    subroutine kernel_QS(self)
        class(QuinticSpline) :: self
        real(kind=prec)      :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 1.0d0/(120.0d0 * self%h)
        else if (self%dim == 2) then
            a_d = 7.0d0/(478.0d0 * PI * self%h**2)
        else
            a_d = 3.0d0/(359.0d0 * PI * self%h**3)
        end if

        !> compute W(r,h)
        if (self%q < 1) then
            self%W = a_d * ((3.0d0 - self%q)**5 - 6*(2.0d0 - self%q)**5 + 15*(1.0d0 - self%q)**5)
        else if (self%q < 2) then
            self%W = a_d * ((3.0d0 - self%q)**5 - 6.0d0*(2.0d0 - self%q)**5)
        else if (self%q < 3) then
            self%W = a_d * (3.0d0 - self%q)**5
        else
            self%W = 0
        end if
    end subroutine kernel_QS

    subroutine dwdq_QS(self)
        class(QuinticSpline) :: self
        real(kind=prec)      :: a_d

        self%q = self%r_ab/self%h

        if (self%dim == 1) then
            a_d = 1.0d0/(120.0d0 * self%h)
        else if (self%dim == 2) then
            a_d = 7.0d0/(478.0d0 * PI * self%h**2)
        else
            a_d = 3.0d0/(359.0d0 * PI * self%h**3)
        end if

        !> compute dW/dr
        if (self%q < 1) then
            self%W = a_d/self%h * (-5.0d0*(3 - self%q)**4 + 30.0d0*(2 - self%q)**4 + 75.0d0*(1.0d0 - self%q)**5)
        else if (self%q < 2) then
            self%W = a_d * (-5.0d0*(3 - self%q)**4 + 30.0d0*(2 - self%q)**4)
        else if (self%q < 3) then
            self%W = a_d * (-5.0d0*(3 - self%q)**4)
        else
            self%W = 0
        end if
    end subroutine dwdq_QS

end module kernels