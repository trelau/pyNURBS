module calculate
use math
use geom_utils
use divide
implicit none

private
public :: triangle_area
public :: cp_lengths, cp_net_areas
public :: arc_length_bezier, surface_area_bezier

contains

function triangle_area(m, tri) result(area)
    !> Calculate the area of a triangle.
    !> m - Number columns in triangle array.
    !> tri - Triangle array.
    !> area - Area of triangle.
    
    ! Input
    integer, intent(in) :: m
    double precision, intent(in) :: tri(3, m)
    
    ! Output
    double precision :: area
    
    ! Working
    double precision :: mag
    double precision :: v1(m), v2(m), vc(m)
    
    v1 = tri(2, :) - tri(1, :)
    v2 = tri(3, :) - tri(1, :)
    
    vc = cross(v1, v2)
    mag = norm(vc)
    area = 0.50d0 * mag
    
end function triangle_area

subroutine cp_lengths(n, cp, lp, lc)
    !> Calculate the linear length of the control polygon and its chord.
    !> n - Number of control points - 1.
    !> cp - Dehomogenized control points.
    !> lp - Length of control polygon.
    !> lc - Chord length.
    
    !f2py intent(in) n, cp
    !f2py intent(out) lp, lc
    !f2py depend(n) cp
    
    ! Input
    integer, intent(in) :: n
    
    ! Output
    double precision, intent(out) :: lp, lc
    double precision, intent(in) :: cp(0:n, 0:2)    
    
    ! Working
    integer :: i
    double precision :: vec(0:2)
    
    vec = cp(n, :) - cp(0, :)
    lc = norm(vec)
    lp = 0.0d0
    do i = 0, n - 1
        vec = cp(i + 1, :) - cp(i, :)
        lp = lp + norm(vec)
    end do
    
end subroutine cp_lengths

subroutine cp_net_areas(n, m, cp, anet, ap)
    !> Calculate the surface areas of the control net.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cp - Dehomogenized control points.
    !> anet - Area of control net.
    !> ap - Area of the plane formed by the four corners of the
    !> control net.
    
    !f2py intent(in) n, m, cp
    !f2py intent(out) anet, ap
    !f2py depend(n, m) cp
    
    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: cp(0:n, 0:m, 0:2)
    
    ! Output
    double precision, intent(out) :: anet, ap
    
    ! Working
    integer :: i, ntri
    integer :: triangles(0:2 * n * m - 1, 0:2)
    double precision :: v1(0:2), v2(0:2), v3(0:2)
    double precision :: verts(0:(n + 1) * (m + 1) - 1, 0:2)
    double precision :: c1(0:2), c2(0:2), c3(0:2)
                        
    call tessellate_cp_net(n, m, cp, verts, triangles)
    ntri = size(triangles, 1)
    anet = 0.0d0
    do i = 0, ntri - 1
        v1 = verts(triangles(i, 0), :)
        v2 = verts(triangles(i, 1), :)
        v3 = verts(triangles(i, 2), :)
        c1 = v2 - v1
        c2 = v3 - v1
        c3 = cross(c1, c2)
        anet = anet + 0.50d0 * norm(c3)
    end do
    c1 = cp(n, 0, :) - cp(0, 0, :)
    c2 = cp(n, m, :) - cp(0, 0, :)
    c3 = cross(c1, c2)
    ap = 0.50d0 * norm(c3)
    c1 = cp(n, m, :) - cp(0, 0, :)
    c2 = cp(0, m, :) - cp(0, 0, :)
    c3 = cross(c1, c2)
    ap = ap + 0.50d0 * norm(c3)
end subroutine cp_net_areas

recursive subroutine arc_length_bezier(n, cpw, u0, u1, tol, l)
    !> Estimate the arc length of a Bezier curve.
    !> n - Number of control points - 1.
    !> cpw - Control points.
    !> u0 - Starting parameter.
    !> u1 - Ending parameter.
    !> tol - Tolerance for comparing control polygon length and chord
    !> length.
    !> l - Approximate arc length.
    
    !f2py intent(in) n, cpw, tol, u0, u1
    !f2py intent(in,out) l
    !f2py depend(n) cpw
    
    ! Input
    integer, intent(in) :: n
    double precision, intent(in) :: u0, u1, tol
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Input/Output
    double precision, intent(inout) :: l(1)
    
    ! Working
    integer :: i
    double precision :: lp, lc, a1, b1, a2, b2
    double precision :: cp(0:n, 0:2)
    double precision :: qw1(0:n, 0:3), qw2(0:n, 0:3)
    double precision :: temp2(0:n, 0:3), w(0:n)
                        
    ! Extract curve if needed
    if ((u0 .gt. 0.0d0) .or. (u1 .lt. 1.0d0)) then
        call extract_bezier_curve(cpw, n, u0, u1, 0.0d0, 1.0d0, temp2, a1, b1)
    else
        temp2 = cpw
    end if

    call dehomogenize_array1d(n, temp2, cp, w)
    call cp_lengths(n, cp, lp, lc)
    
    if (dabs(lp - lc) .le. tol) then
        l(1) = l(1) + (2.0d0 * lc + (dble(n) - 1.0d0) * lp) / (dble(n) + 1.0d0)
    else
        call split_bezier_curve(temp2, n, 0.50d0, 0.0d0, 1.0d0, qw1, a1, b1, &
                                qw2, a2, b2)
        call arc_length_bezier(n, qw1, 0.0d0, 1.0d0, tol, l)
        call arc_length_bezier(n, qw2, 0.0d0, 1.0d0, tol, l)
    end if
    
end subroutine arc_length_bezier

recursive subroutine surface_area_bezier(n, m, cpw, tol, a)
    !> Estimate the surface area of a Bezier surface.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cpw - Control points.
    !> tol - Tolerance for comparing area of control net to area of
    !> base area of formed by the four corners of the control net.
    !> a - Approximate area of surface.
    
    !f2py intent(in) n, m, cpw, tol
    !f2py intent(in,out) a
    !f2py depend(n, m) cpw
    
    ! Input 
    integer, intent(in) :: n, m
    double precision, intent(in) :: tol
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Input/Output
    double precision, intent(inout) :: a(1)
    
    integer :: i
    double precision :: ap, anet
    double precision :: cp(0:n, 0:m, 0:2)
    double precision :: qw1(0:n, 0:m, 0:3), qw2(0:n, 0:m, 0:3)
    double precision :: qw3(0:n, 0:m, 0:3), qw4(0:n, 0:m, 0:3)
    double precision :: qw5(0:n, 0:m, 0:3), qw6(0:n, 0:m, 0:3)
    double precision :: w(0:n, 0:m)
    double precision :: ab1(0:1, 0:1), ab2(0:1, 0:1)

    call dehomogenize_array2d(n, m, cpw, cp, w)
    call cp_net_areas(n, m, cp, anet, ap)
    
    if (abs(anet - ap) .le. tol) then
        a(1) = a(1) + (4.0d0 * ap + dble(n + m + n * m - 3) * anet) / &
               dble((n + 1) * (m + 1))
    else
        call split_bezier_surface(cpw, n, m, 0.50d0, -1.0d0, &
                                  0.0d0, 1.0d0, 0.0d0, 1.0d0, &
                                  qw1, qw2, ab1, ab2)
        call split_bezier_surface(qw1, n, m, -1.0d0, 0.50d0, &
                                  0.0d0, 1.0d0, 0.0d0, 1.0d0, &
                                  qw3, qw4, ab1, ab2)
        call split_bezier_surface(qw2, n, m, -1.0d0, 0.50d0, &
                                  0.0d0, 1.0d0, 0.0d0, 1.0d0, &
                                  qw5, qw6, ab1, ab2)
        call surface_area_bezier(n, m, qw3, tol, a)
        call surface_area_bezier(n, m, qw4, tol, a)
        call surface_area_bezier(n, m, qw5, tol, a)
        call surface_area_bezier(n, m, qw6, tol, a)
    end if
    
end subroutine surface_area_bezier

end module calculate