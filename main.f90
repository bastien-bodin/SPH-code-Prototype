program main
    use application
    implicit none

    type(SPH_App) :: MyApp

    call MyApp%initialisation()
    call MyApp%solver()
    
end program main