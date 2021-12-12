module equations
    use parameters
    use particles
    implicit none

    ! here :
    ! Miscellaneous
    !   > eval_r -> Returns the distance between two points
    ! Sound speed
    !   > CS_IG  -> Compute the sound speed for ideal gas
    !   > CS_WC  -> Compute the sound speed for weakly compressible fluids
    ! Pressure
    !   > P_IG   -> Compute the pressure for ideal gas
    !   > P_IG   -> Compute the pressure for weakly compressible fluids

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Miscellaneous
    function eval_r(xyz1,xyz2) result(r_ab)
    !> compute the distance between two points
        real(kind=prec), dimension(1:3) :: xyz1
        real(kind=prec), dimension(1:3) :: xyz2
        real(kind=prec)                 :: r_ab

        r_ab = sqrt(sum( (xyz1(:) - xyz2(:)) * (xyz1(:) - xyz2(:)) ))
    end function eval_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sound speed
    subroutine CS_IG(Part, integStep)
    !> Compute the sound speed for ideal gas
        class(Particle) :: Part
        integer :: IntegStep

        Part%c_sound(IntegStep+1) = Part%c_sound(IntegStep)
    end subroutine CS_IG

    subroutine CS_WC(Part,rho_0,c_0,state_gamma, integStep)
    !> Compute the sound speed for weakly compressible fluids
        class(Particle) :: Part
        real(kind=prec), intent(in)       :: rho_0
        real(kind=prec), intent(in)       :: c_0
        real(kind=prec), intent(in)       :: state_gamma
        integer :: IntegStep

        Part%c_sound(IntegStep+1) = c_0 * SQRT( (Part%density(IntegStep+1)/rho_0)**(state_gamma-1) )
    end subroutine CS_WC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pressure
    subroutine P_IG(Part,rho_0,molMass, integStep)
    !> Compute the pressure for ideal gas
        class(Particle) :: Part
        real(kind=prec), intent(in) :: rho_0
        real(kind=prec), intent(in) :: molMass
        integer :: IntegStep

        Part%pressure(IntegStep) = (part%density(IntegStep) / rho_0 - 1) * idealGasCst * 293.15d0 / molMass
    end subroutine P_IG

    subroutine P_WC(Part,c_0,rho_0,state_gamma, integStep)
    !> Compute the pressure for weakly compressible fluid
        class(Particle) :: Part
        real(kind=prec), intent(in)       :: rho_0
        real(kind=prec), intent(in)       :: c_0
        real(kind=prec), intent(in)       :: state_gamma
        integer :: IntegStep

        real(kind=prec) :: B

        B = c_0**2 * rho_0/state_gamma

        Part%pressure(IntegStep+1) = B * ((Part%density(IntegStep+1)/rho_0)**state_gamma - 1)
    end subroutine P_WC

    function ArtVisc(Part, Neigh, alpha, beta, IntegStep) result(pi_ab)
        class(Particle) :: Part
        class(Particle) :: Neigh
        real(kind=prec) :: alpha
        real(kind=prec) :: beta
        integer :: IntegStep
        real(kind=prec) :: pi_ab

        real(kind=prec) :: mu_ab = 0.0d0
        real(kind=prec), dimension(1:3) :: v_ab
        real(kind=prec), dimension(1:3) :: x_ab
        real(kind=prec) :: v_dot_x
        real(kind=prec) :: x_dot_x
        real(kind=prec) :: eta2
        real(kind=prec) :: c_ab
        real(kind=prec) :: rho_ab

        v_ab = Part%velocity(:,IntegStep) - Neigh%velocity(:,IntegStep)
        x_ab = Part%coords(:,IntegStep) - Neigh%coords(:,IntegStep)
        v_dot_x = dot_product(v_ab,x_ab)
        x_dot_x = dot_product(x_ab,x_ab)
        eta2 = 0.01d0 * Part%h_part**2.0d0

        if (v_dot_x < 0) then
            mu_ab = Part%h_part * v_dot_x/(x_dot_x + eta2)
            c_ab = 0.5d0 * (Part%c_sound(IntegStep) + Neigh%c_sound(IntegStep))
            rho_ab = 0.5d0 * (Part%density(IntegStep) + Neigh%density(IntegStep))
            pi_ab = (-alpha * c_ab * mu_ab + beta * mu_ab**2.0d0)/rho_ab
        else
            pi_ab = 0.0d0
        end if

        if ((IntegStep == 1) .and. (mu_ab > Part%max_mu_ab)) then
            Part%max_mu_ab = mu_ab
        end if

    end function ArtVisc

end module equations