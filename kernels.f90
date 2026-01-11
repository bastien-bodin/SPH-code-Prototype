module kernels
    use parameters
    implicit none

    !> Support radius factors for different kernels
    !> The kernel vanishes for r > kappa * h
    integer, parameter :: KAPPA_CUBIC  = 2
    integer, parameter :: KAPPA_QUINTIC = 3

contains

    ! --- CUBIC SPLINE KERNEL ---

    !> Compute Cubic Spline W and dW/dr
    !> Standard B-Spline kernel for SPH simulations.
    pure subroutine cubic_spline(r, h, W, dWdr)
        real(prec), intent(in)  :: r, h
        real(prec), intent(out) :: W, dWdr
        
        real(prec) :: q, alpha_d, invh
        
        q = r / h
        invh = 1.0_prec / h

        ! Normalization constant alpha_d depends on nDim
        if (nDim == 2) then
            alpha_d = 10.0_prec / (7.0_prec * PI * h**2)
        else ! nDim == 3
            alpha_d = 1.0_prec / (PI * h**3)
        end if

        if (q < 1.0_prec) then
            W = alpha_d * (1.0_prec - 1.5_prec * q**2 + 0.75_prec * q**3)
            dWdr = alpha_d * invh * (-3.0_prec * q + 2.25_prec * q**2)
        else if (q < 2.0_prec) then
            W = alpha_d * 0.25_prec * (2.0_prec - q)**3
            dWdr = -alpha_d * invh * 0.75_prec * (2.0_prec - q)**2
        else
            W = 0.0_prec
            dWdr = 0.0_prec
        end if
    end subroutine cubic_spline


    ! --- QUINTIC SPLINE KERNEL ---

    !> Compute Quintic Spline W and dW/dr
    !> More stable than cubic spline for high-pressure gradients.
    pure subroutine quintic_spline(r, h, W, dWdr)
        real(prec), intent(in)  :: r, h
        real(prec), intent(out) :: W, dWdr
        
        real(prec) :: q, alpha_d, invh
        
        q = r / h
        invh = 1.0_prec / h

        if (nDim == 2) then
            alpha_d = 7.0_prec / (478.0_prec * PI * h**2)
        else ! nDim == 3
            alpha_d = 3.0_prec / (359.0_prec * PI * h**3)
        end if

        if (q < 1.0_prec) then
            W = alpha_d * ((3.0_prec-q)**5 - 6.0_prec*(2.0_prec-q)**5 + &
                15.0_prec*(1.0_prec-q)**5)
            dWdr = alpha_d * invh * (-5.0_prec*(3.0_prec-q)**4 + &
                   30.0_prec*(2.0_prec-q)**4 - 75.0_prec*(1.0_prec-q)**4)
        else if (q < 2.0_prec) then
            W = alpha_d * ((3.0_prec-q)**5 - 6.0_prec*(2.0_prec-q)**5)
            dWdr = alpha_d * invh * (-5.0_prec*(3.0_prec-q)**4 + &
                   30.0_prec*(2.0_prec-q)**4)
        else if (q < 3.0_prec) then
            W = alpha_d * (3.0_prec-q)**5
            dWdr = alpha_d * invh * (-5.0_prec*(3.0_prec-q)**4)
        else
            W = 0.0_prec
            dWdr = 0.0_prec
        end if
    end subroutine quintic_spline

end module kernels