program main
    use application
    implicit none

    !> Main entry point for the SPH Cryolava Simulator
    type(SPH_App) :: MyApp

    print *, "=========================================="
    print *, "   EUROPA CRYOLAVA SPH PROTOTYPE v2.0    "
    print *, "=========================================="

    ! 1. Initialize the simulation environment
    ! This allocates memory and generates the Cabrera-Crespo geometry
    call MyApp%setup_scene()

    ! 2. Launch the solver
    ! This runs the time-stepping loop (RK22 or Verlet)
    call MyApp%run_simulation()

    print *, "=========================================="
    print *, "      SIMULATION COMPLETED SUCCESSFULLY   "
    print *, "=========================================="

end program main