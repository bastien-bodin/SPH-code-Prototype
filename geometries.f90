module geometries
    use parameters
    use particles
    implicit none

contains

    subroutine setup_dam_break(sys, L_tank, H_tank, L_fluid, H_fluid, &
                               spacing, h0)
        type(ParticleSystem), intent(inout) :: sys
        real(prec), intent(in) :: L_tank, H_tank, L_fluid, H_fluid, &
                                  spacing, h0
        
        print *, "--- Geometry: Cabrera-Crespo (2007) DBC ---"
        print *, "Gap: spacing/2 in X and Y between layers"
        
        ! 1. Create Tank (Layers at spacing/2)
        call add_cabrera_tank(sys, L_tank, H_tank, spacing, h0)
        
        ! 2. Create Fluid Block (Origin centered in gap)
        ! Fluid begins at x = +spacing/2, y = +spacing/2
        call add_fluid_block(sys, 0.5_prec*spacing, 0.5_prec*spacing, &
                             L_fluid, H_fluid, spacing, h0)
                             
        print *, "Total particles: ", sys%nPart
    end subroutine setup_dam_break


    subroutine add_cabrera_tank(sys, L_tank, H_tank, spacing, h0)
        type(ParticleSystem), intent(inout) :: sys
        real(prec), intent(in) :: L_tank, H_tank, spacing, h0
        
        integer    :: i, j, layer
        real(prec) :: x, y, dx_half, dy_half, wall_inner
        real(prec) :: v0(nDim)
        
        v0 = 0.0_prec
        dx_half = 0.5_prec * spacing
        dy_half = 0.5_prec * spacing
        ! wall_inner defines the contact face at -spacing/2
        wall_inner = -0.5_prec * spacing 

        ! --- BOTTOM WALL (2 Layers, dy = spacing/2) ---
        do layer = 0, 1
            y = wall_inner - (layer * dy_half)
            do i = -1, int(L_tank/spacing) + 1
                ! Staggering shift in X for the second layer
                x = (i * spacing) + wall_inner + (layer * dx_half)
                call sys%add_particle([x, y], v0, 1.0_prec, h0, .false.)
            end do
        end do

        ! --- SIDE WALLS (2 Layers, dx = spacing/2) ---
        do layer = 0, 1
            do j = 0, int(H_tank/spacing)
                ! Staggering shift in Y for the second layer
                y = (j * spacing) + wall_inner + (layer * dy_half)
                
                ! Left Wall: offsets to the left (negative x)
                x = wall_inner - (layer * dx_half)
                call sys%add_particle([x, y], v0, 1.0_prec, h0, .false.)
                
                ! Right Wall: offsets to the right (positive x)
                x = (L_tank + wall_inner) + (layer * dx_half)
                call sys%add_particle([x, y], v0, 1.0_prec, h0, .false.)
            end do
        end do
    end subroutine add_cabrera_tank


    subroutine add_fluid_block(sys, x0, y0, width, height, spacing, h0)
        type(ParticleSystem), intent(inout) :: sys
        real(prec), intent(in) :: x0, y0, width, height, spacing, h0
        integer :: i, j
        real(prec) :: x, y
        real(prec) :: v0(nDim)
        
        v0 = 0.0_prec
        do i = 0, int(width/spacing) - 1
            do j = 0, int(height/spacing) - 1
                x = x0 + i * spacing
                y = y0 + j * spacing
                call sys%add_particle([x, y], v0, 1.0_prec, h0, .true.)
            end do
        end do
    end subroutine add_fluid_block

end module geometries