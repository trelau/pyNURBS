module map_surface
use math
use evaluate
implicit none

private
public :: surface_distance_map

contains

subroutine surface_distance_map(n, p, uk, m, q, vk, cpw, nu, nv, ugrid, vgrid, &
                                udist, vdist, uparam, vparam, xy)
    !> Build a distance map of the surface.
    
    !f2py intent(in) n, p, uk, q, m, vk, cpw, nu, nv, ugrid, vgrid
    !f2py intent(out) udist, vdist, uparam, vparam, xy
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    !f2py depend(nu) ugrid
    !f2py depend(nv) vgrid
    !f2py depend(nu, nv) udist, vdist, xy
    
    ! Input
    integer, intent(in) :: n, p, m, q, nu, nv
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    double precision, intent(in) :: ugrid(nu), vgrid(nv)
    
    ! Ouptut
    double precision, intent(out) :: udist(nu, nv), vdist(nu, nv)
    double precision, intent(out) :: uparam(nu * nv), vparam(nu * nv)
    double precision, intent(out) :: xy(nu * nv, 2)

    ! Working
    integer :: i, j, k
    double precision :: d, dv, du, du_0, du_1, dv_0, dv_1
    double precision :: p0(3), p1(3), dp(3)
    
    udist(:, :) = 0.0d0
    vdist(:, :) = 0.0d0
    uparam(:) = 0.0d0
    vparam(:) = 0.0d0
    xy(:, :) = 0.0d0
    
    ! Map u-direction.
    do j = 1, nv
        dv = vgrid(j)
        d = 0.0d0
        do i = 2, nu
            du_0 = ugrid(i - 1)
            du_1 = ugrid(i)
            call surface_point(n, p, uk, m, q, vk, cpw, du_0, dv, p0)
            call surface_point(n, p, uk, m, q, vk, cpw, du_1, dv, p1)
            dp = p1 - p0
            d = d + norm(dp)
            udist(i, j) = d
        end do
    end do
    
    ! Map v-direction.
    do i = 1, nu
        du = ugrid(i)
        d = 0.0d0
        do j = 2, nv
            dv_0 = vgrid(j - 1)
            dv_1 = vgrid(j)
            call surface_point(n, p, uk, m, q, vk, cpw, du, dv_0, p0)
            call surface_point(n, p, uk, m, q, vk, cpw, du, dv_1, p1)
            dp = p1 - p0
            d = d + norm(dp)
            vdist(i, j) = d
        end do
    end do
    
    ! Build data for inversion.
    k = 1
    do j = 1, nv
        do i = 1, nu
            xy(k, :) = (/ udist(i, j), vdist(i, j) /)
            uparam(k) = ugrid(i)
            vparam(k) = vgrid(j)
            k = k + 1
        end do
    end do

end subroutine surface_distance_map

end module map_surface
