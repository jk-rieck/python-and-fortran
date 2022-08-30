subroutine solve(D, phi, delta, nx, ny)
! =====================================================
! Solve poisson equation with Gauss-Seidel
! =====================================================
    integer, intent(in)                            :: nx, ny
    real, intent(in)                               :: delta
    real(kind=8), intent(in), dimension(ny, nx)    :: D(ny, nx)
    real(kind=8), intent(inout), dimension(ny, nx) :: phi(ny, nx)
    integer                                        :: i, j
    real(kind=8)                                   :: tmp
 
    do i = 2, nx-1
        do j = 2, ny-1
            phi(j, i) = 0.25 * (phi(j, i+1) + phi(j, i-1) &
                            & + phi(j+1, i) + phi(j-1, i) &
                            & - ((delta * delta) * D(j, i)))
        end do
    end do
 
end subroutine solve