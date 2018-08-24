module geom_utils
use config, only: gtol
use math
implicit none

private
public :: dehomogenize_array1d, dehomogenize_array2d, tessellate_cp_net
public :: local_to_global_param, global_to_local_param
public :: is_curve_flat, is_surface_flat
public :: equivalence_array1d, barycentric_params, are_points_equal
public :: is_surface_closed
public :: invert_point_in_triangle
public :: angle_between_vecs

contains

subroutine dehomogenize_array1d(n, cpw, cp, w)
    !> Convert the 1D homogeneous control points cpw to
    !> dehomogenized coordinates (x/w, y/w, z/w).
    !> n - Number of control points - 1.
    !> cpw - Homogeneous control points.
    !> cp - Dehomogenized control points.
    !> w - Weights.
    
    !f2py intent(in) n, cpw
    !f2py intent(out) cp, w
    !f2py depend(n) cpw, cp, w
    
    ! Input
    integer, intent(in) :: n
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: cp(0:n, 0:2), w(0:n)
    
    ! Working
    integer :: i
    
    do i = 0, n
        w(i) = cpw(i, 3)
        cp(i, :) = cpw(i, :2) / w(i)
    end do
    
end subroutine dehomogenize_array1d

subroutine dehomogenize_array2d(n, m, cpw, cp, w)
    !> Convert the 2D homogeneous control points cpw to
    !> dehomogenized coordinates (x/w, y/w, z/w).
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cpw - Homogeneous control points.
    !> cp - Dehomogenized control points.
    !> w - Weights.
    
    !f2py intent(in) n, m, cpw
    !f2py intent(out) cp, w
    !f2py depend(n, m) cpw, cp, w
    
    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: cp(0:n, 0:m, 0:2)
    double precision, intent(out) :: w(0:n, 0:m)
    
    ! Working
    integer :: i, j
    
    do j = 0, m
        do i = 0, n
            w(i, j) = cpw(i, j, 3)
            cp(i, j, :) = cpw(i, j, :2) / w(i, j)
        end do
    end do
    
end subroutine dehomogenize_array2d

subroutine tessellate_cp_net(n, m, cp, verts, triangles)
    !> Tessellate the control net of the surface as triangles.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cp - Dehomogenized control points.
    !> verts - Vertices of triangles.
    !> triangles - Connectivity of triangles where each row is a
    !> triangle and each column is a vertex id.
    
    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: cp(0:n, 0:m, 0:2)
    
    ! Output
    integer, intent(out) :: triangles(0:2 * n * m - 1, 0:2)
    double precision, intent(out) :: verts(0:(n + 1) * (m + 1) - 1, 0:2)
                                     
    ! Working
    integer :: i, j, t, v, v1, v2, v3, v4
    integer :: vert_id(0:n, 0:m)
    
    v = 0
    do j = 0, m
        do i = 0, n
            verts(v, :) = cp(i, j, :)
            vert_id(i, j) = v
            v = v + 1
        end do
    end do
    t = 0
    do j = 0, m - 1
        do i = 0, n - 1
            v1 = vert_id(i, j)
            v2 = vert_id(i + 1, j)
            v3 = vert_id(i + 1, j + 1)
            v4 = vert_id(i, j + 1)
            triangles(t, :) = (/ v1, v2, v3 /)
            triangles(t + 1, :) = (/ v1, v3, v4 /)
            t = t + 2
        end do
    end do
end subroutine tessellate_cp_net

function local_to_global_param(a, b, ul) result(ug)
    !> Convert the local parameter to the global domain
    !> (a <= u <= b).
    !> a - Lower global parameter.
    !> b - Upper global parameter.
    !> ul - Local parameter.
    !> ug - Global parameter.
    
    ! Input
    double precision, intent(in) :: a, b, ul
    
    ! Output
    double precision :: ug
    
    ! Working
    double precision :: ui
    
    ui = ul
    if (ui .lt. 0.0d0) then
        ui = 0.0d0
    else if (ui .gt. 1.0d0) then
        ui = 1.0d0
    end if
    ug = a + ui * (b - a)
    
end function local_to_global_param

function global_to_local_param(a, b, ug) result(ul)
    !> Convert the global parameter to the local domain
    !> (0 <= u <= 1).
    !> a - Lower global parameter.
    !> b - Upper global parameter.
    !> ug - Global parameter.
    !> ul - Local parameter.
    
    ! Input
    double precision, intent(in) :: a, b, ug
    
    ! Output
    double precision :: ul
    
    ul = (ul - a) / (b - a)
    if (ul .lt. 0.0d0) then
        ul = 0.0d0
    else if (ul .gt. 1.0d0) then
        ul = 1.0d0
    end if
    
end function global_to_local_param

subroutine is_curve_flat(n, cp, tol, is_flat)
    !> Check flatness of curve.
    !> n - Number of control points - 1.
    !> cp - Dehomogenized control points.
    !> tol - Tolerance.
    !> is_flat - True if curve is flat, False if not.
    
    !f2py intent(in) n, cp, tol
    !f2py intent(out) is_flat
    !f2py depend(n) cp
    
    ! Input    
    integer, intent(in) :: n
    double precision, intent(in) ::  tol
    double precision, intent(in) :: cp(0:n, 0:2)
    
    ! Output
    logical, intent(out) :: is_flat
    
    ! Working
    integer :: i
    double precision :: d, dline
    double precision :: v0(3), v1(3), v2(3), v3(3), vec(3)
    
    is_flat = .true.
    vec = cp(n, :) - cp(0, :)
    dline = norm(vec)
    
    if (dline .le. 0.0d0) then
        ! Curve must be closed. Use the distnace between points.
        do i = 1, n - 1
            v0 = cp(i, :) - cp(0, :)
            d = norm(v0)
            if (d .gt. tol) then
                is_flat = .false.
                return
            end if
        end do
        is_flat = .true.
        return
    end if
    
    do i = 1, n - 1
        v0 = cp(i, :) - cp(0, :)
        v1 = cp(i, :) - cp(n, :)
        v2 = cross(v0, v1)
        d = norm(v2) / dline
        if (d .gt. tol) then
            is_flat = .false.
            return
        end if
    end do
end subroutine is_curve_flat

subroutine is_surface_flat(n, m, cp, tol, is_flat)
    !> Check flatness of surface.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cp - Dehomogenized control points.
    !> tol - Tolerance.
    !> is_flat - True if surface is flat, False if not.

    !f2py intent(in) n, m, cp, tol
    !f2py intent(out) is_flat
    !f2py depend(n, m) cp
    
    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: tol
    double precision, intent(in) :: cp(0:n, 0:m, 0:2)
    
    ! Output
    logical, intent(out) :: is_flat

    ! Working
    integer :: i, j
    double precision :: row(0:m, 0:2), col(0:n, 0:2)
    
    is_flat = .true.
    
    ! Columns
    do j = 0, m
        col = cp(:, j, :)
        call is_curve_flat(n, col, tol, is_flat)
        if (is_flat .eqv. .false.) then
            return
        end if
    end do
    
    ! Rows
    do i = 0, n
        row = cp(i, :, :)
        call is_curve_flat(m, row, tol, is_flat)
        if (is_flat .eqv. .false.) then
            return
        end if
    end do
    
end subroutine is_surface_flat

subroutine equivalence_array1d(d0, d1, array, tol, npts, eq_array)
    !> Equivalence points of an array.
    !> d0 - Dimension 1 of array (number of points).
    !> d1 - Dimension 2 of array (3 for dehomogenized array,
    !> 4 for homogeneous array).
    !> array - Input array.
    !> tol - Tolerance.
    !> npts - Number of unique points.
    !> eq_array - Array of unique points.
    
    !f2py intent(in) d0, d1, array, tol
    !f2py intent(out) npts, eq_array
    !f2py depend(d0, d1) array, eq_array
    
    ! Input
    integer, intent(in) :: d0, d1
    double precision, intent(in) ::  tol
    double precision, intent(in) :: array(d0, d1)
    
    ! Output
    integer, intent(out) :: npts
    double precision, intent(out) :: eq_array(d0, d1)
    
    ! Working
    logical :: unique
    integer :: i, j, k
    double precision :: dist, v(d1)
    
    eq_array(:, :) = 0.0d0
    eq_array(1, :) = array(1, :)
    npts = 1
    do i = 2, d0
        unique = .true.
        do j = 1, npts
            v = eq_array(j, :) - array(i, :)
            dist = 0.0d0
            do k = 1, d1
                dist = dist + v(k) * v(k)
            end do
            dist = dsqrt(dist)
            if (dist .le. tol) then
                unique = .false.
                exit
            end if
        end do
        if (unique .eqv. .true.) then
            npts = npts + 1
            eq_array(npts, :) = array(i, :)
        end if
    end do
end subroutine equivalence_array1d

subroutine barycentric_params(p, t, u, v)
    !> Determine the barycentric coordinates of a point with respect
    !> to the triangle defined by three vertices.
    !> p - Point.
    !> t - Vertices of triangle.
    !> u - Barycentric coordinate.
    !> v - Barycentric coordinate.
    
    ! Input
    double precision, intent(in) :: p(0:2), t(0:2, 0:2)
    
    ! Output
    double precision, intent(out) :: u, v
    
    ! Working
    double precision :: d01, d12, d02, dp1, dp2, denom
    double precision :: v01(0:2), v02(0:2), vp(0:2)                        
                        
    v01 = t(1, :) - t(0, :)
    v02 = t(2, :) - t(0, :)
    vp = p - t(0, :)
    d01 = dot_product(v01, v01)
    d12 = dot_product(v01, v02)
    d02 = dot_product(v02, v02)
    dp1 = dot_product(vp, v01)
    dp2 = dot_product(vp, v02)
    denom = d01 * d02 - d12 * d12
    u = (d02 * dp1 - d12 * dp2) / denom
    v = (d01 * dp2 - d12 * dp1) / denom
    
    if (u .lt. 0.0d0) then
        u = 0.0d0
    elseif (u .gt. 1.0d0) then
        u = 1.0d0
    end if
    
    if (v .lt. 0.0d0) then
        v = 0.0d0
    elseif (v .gt. 1.0d0) then
        v = 1.0d0
    end if
    
end subroutine barycentric_params

function are_points_equal(p1, p2, tol) result(are_equal)
    !> Check to see if two points are coincident within a tolerance.
    !> p1 - A point.
    !> p2 - Other point.
    !> tol - Tolerance.
    !> are_equal - True if points are equal, False if not.
    
    ! Input
    double precision, intent(in) :: tol
    double precision, intent(in) :: p1(:), p2(:)
    
    ! Output
    logical :: are_equal
    
    ! Working
    double precision :: d
    
    are_equal = .false.
    d = norm(p1 - p2)
    if (d .le. tol) then
        are_equal = .true.
    end if
    
end function are_points_equal

subroutine is_surface_closed(n, m, cpw, uclosed, vclosed)
    !> Check to see if a surface is closed.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cpw - Control points.
    !> uclosed - Is surface closed in u-direction.
    !> vclosed - Is surface closed in v-direction.

    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    logical, intent(out) :: uclosed, vclosed
    
    ! Working
    logical :: are_equal
    integer :: i
    double precision :: p1(0:3), p2(0:3)
    
    uclosed = .true.
    vclosed = .true.
    
    ! Check for u-direction by comparing control points at first and last rows.
    do i = 0, m
        p1 = cpw(0, i, :)
        p2 = cpw(n, i, :)
        are_equal = are_points_equal(p1, p2, gtol)
        if (are_equal .eqv. .false.) then
            uclosed = .false.
            exit
        end if
    end do
    
    ! Check for v-direction by comparing control points at first and last columns.
    do i = 0, n
        p1 = cpw(i, 0, :)
        p2 = cpw(i, m, :)
        are_equal = are_points_equal(p1, p2, gtol)
        if (are_equal .eqv. .false.) then
            vclosed = .false.
            exit
        end if
    end do
    
end subroutine is_surface_closed

subroutine invert_point_in_triangle(p, t3d, t2d, u, v)
    !> Find parameters for a point inside a triangle.
    !> p - Point.
    !> t3d - Vertices of triangle in 3-D.
    !> t2d - Vertices of triangle in 2-D.
    !> u - Point parameter in u-direction.
    !> v - Point parameter in v-direction.
    
    ! Input
    double precision, intent(in) :: p(0:2), t3d(0:2, 0:2), t2d(0:2, 0:1)
    
    ! Output
    double precision, intent(out) :: u, v
    
    ! Working
    double precision :: d01, d12, d02, dp1, dp2, denom, ub, vb
    double precision :: v01(0:2), v02(0:2), vp(0:2)                        
                        
    ! Find barycentric parameters.
    v01 = t3d(1, :) - t3d(0, :)
    v02 = t3d(2, :) - t3d(0, :)
    vp = p - t3d(0, :)
    d01 = dot_product(v01, v01)
    d12 = dot_product(v01, v02)
    d02 = dot_product(v02, v02)
    dp1 = dot_product(vp, v01)
    dp2 = dot_product(vp, v02)
    denom = d01 * d02 - d12 * d12
    ub = (d02 * dp1 - d12 * dp2) / denom
    vb = (d01 * dp2 - d12 * dp1) / denom
    
    ! Find 2-D parameters.
    u = (1.0d0 - ub - vb) * t2d(0, 0) + ub * t2d(1, 0) + vb * t2d(2, 0)
    v = (1.0d0 - ub - vb) * t2d(0, 1) + ub * t2d(1, 1) + vb * t2d(2, 1)

end subroutine invert_point_in_triangle

function angle_between_vecs(v1, v2) result(angle)
    !> Calculate the angle between two vectors.
    !> v1 - Vector 1.
    !> v2 - Vector 2.
    
    ! Input
    double precision, intent(in) :: v1(:), v2(:)
    
    ! Output
    double precision :: angle
    
    ! Working
    double precision :: n1, n2, x, arad, pi
    
    n1 = norm(v1)
    n2 = norm(v2)
    x = dot_product(v1, v2) / (n1 * n2)
    if (x .gt. 1.0d0) then
        x = 1.0d0
    elseif (x .lt. -1.0d0) then
        x = -1.0d0
    end if
    arad = dacos(x)
    pi = 4.0d0 * datan(1.0d0)
    angle = arad * 180.0d0 / pi

end function angle_between_vecs

end module geom_utils