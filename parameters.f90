module parameters
    implicit none

    integer, parameter :: prec = selected_real_kind(15, 307)
    integer, parameter :: nDim = 2 

    ! --- Earth Gravity Vector ---
    ! G_ABS: Gravity magnitude
    ! G_VEC: Array constructor that puts 0.0 in the first (nDim-1) components
    !        and -9.81 in the last component (Y in 2D, Z in 3D).
    real(prec), parameter :: G_ABS = -9.81_prec
    ! Reference 3D vector [0, 0, G_ABS]
    real(prec), parameter :: G_REF(3) = [0.0_prec, 0.0_prec, G_ABS]
    
    ! We "slice" the end of the reference vector to match nDim.
    ! If nDim=2, we take the last 2 elements: [0, G_ABS]
    ! If nDim=3, we take the last 3 elements: [0, 0, G_ABS]
    real(prec), parameter :: G_VEC(nDim) = G_REF(3-nDim+1 : 3)

    ! --- Integration Methods ---
    integer, parameter :: INT_EULER    = 1
    integer, parameter :: INT_VERLET   = 2
    integer, parameter :: INT_RK22     = 3
    
    ! Chosen method for the simulation
    integer, parameter :: SELECTED_INT = INT_RK22

    ! Euler needs 2 slots (t, t+dt)
    ! RK22 and Verlet need 3 slots (t, t+dt/2, t+dt)
    integer, parameter :: nSteps = merge(2, 3, SELECTED_INT == INT_EULER)

    ! --- Density Methods ---
    integer, parameter :: DENSITY_SUMMATION  = 1
    integer, parameter :: DENSITY_CONTINUITY = 2
    integer, parameter :: DENSITY_METHOD     = DENSITY_CONTINUITY

    real(prec), parameter :: PI = 3.141592653589793238_prec
end module parameters