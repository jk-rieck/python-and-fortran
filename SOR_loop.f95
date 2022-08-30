subroutine solve(D, phi, om, delta, nx, ny)
! =====================================================
! Solve poisson equation with SOR
! =====================================================
    integer, intent(in)                            :: nx, ny
    real, intent(in)                               :: delta, om
    real(kind=8), intent(in), dimension(ny, nx)    :: D(ny, nx)
    real(kind=8), intent(inout), dimension(ny, nx) :: phi(ny, nx)
    integer                                        :: i, j
    real(kind=8)                                   :: tmp
 
    do i = 2, nx-1
        do j = 2, ny-1
            tmp = (phi(j, i+1) + phi(j, i-1) &
               & + phi(j+1, i) + phi(j-1, i) &
               & - ((delta * delta) * D(j, i)))
            phi(j, i) = (((1 - om) * phi(j, i)) &
                            & + om * (tmp / 4))
        end do
    end do
 
end subroutine solve