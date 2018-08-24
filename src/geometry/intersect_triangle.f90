module intersect_triangle
use math
implicit none

private
public :: intersect_triangle_ray, intersect_triangle_plane, intersect_triangles

contains

subroutine intersect_triangle_ray(v0, v1, v2, origin, ray, eps, &
                                  test_cull, t, u, v, success)
    !> Compute the intersection of a triangle with a ray.
    !> v0 - Vertex 0 of triangle.
    !> v1 - Vertex 1 of triangle.
    !> v2 - Vertex 2 of triangle.
    !> origin - Ray origin.
    !> ray - Ray direction (unit vector).
    !> eps - Nondimensional tolerance used when a ray intersects near the
    !> edges of a triangle.
    !> test_cull - Option to cull back facing triangles.
    !> t - Intersection parameter on ray.
    !> u - Barycentric parameter on triangle.
    !> v - Barycentric parameter on triangle.
    !> success - Logical to specify if an intersection was found at all.
    
    !f2py intent(in) v0, v1, v2, origin, ray, eps, test_cull
    !f2py intent(out) t, u, v, success
    
    ! Input
    logical, intent(in) :: test_cull
    double precision, intent(in) :: eps
    double precision, intent(in) :: v0(0:2), v1(0:2), v2(0:2)
    double precision, intent(in) :: origin(0:2), ray(0:2)
    
    ! Output
    logical, intent(out) :: success
    double precision, intent(out) :: t, u, v
    
    ! Working
    double precision :: det, inv_det
    double precision :: e1(0:2), e2(0:2), pvec(0:2), tvec(0:2), qvec(0:2)
                        
    t = 0.0d0
    u = 0.0d0
    v = 0.0d0
    success = .true.
    
    e1 = v1 - v0
    e2 = v2 - v0
    pvec = cross(ray, e2)
    det = dot_product(e1, pvec)
    if (test_cull) then
        if (det .le. eps) then
            success = .false.
            return
        end if
        tvec = origin - v0
        u = dot_product(tvec, pvec)
        if ((u .le. -eps) .or. (u .ge. det + eps)) then
            success = .false.
            return
        end if
        qvec = cross(tvec, e1)
        v = dot_product(ray, qvec)
        if ((v .le. -eps) .or. (u + v .ge. det + eps)) then
            success = .false.
            return
        end if
        t = dot_product(e2, qvec)
        inv_det = 1.0d0 / det
        t = t * inv_det
        u = u * inv_det
        v = v * inv_det
    else
        if (abs(det) .lt. eps) then
            success = .false.
            return
        end if
        inv_det = 1.0d0 / det
        tvec = origin - v0
        u = dot_product(tvec, pvec) * inv_det
        if ((u .le. -eps) .or. (u .ge. 1.0d0 + eps)) then
            success = .false.
            return
        end if
        qvec = cross(tvec, e1)
        v = dot_product(ray, qvec) * inv_det
        if ((v .le. -eps) .or. (u + v .ge. 1.0d0 + eps)) then
            success = .false.
            return
        end if
        t = dot_product(e2, qvec) * inv_det
    end if
    
end subroutine intersect_triangle_ray

subroutine intersect_triangles(t1, t2, tol, eps, npts, pnts)
    !> Find the intersection segment between two triangles.
    !> t1 - Vertices of triangle 1.
    !> t2 - Vertices of triangle 2.
    !> tol - Geometric tolerance.
    !> eps - Nondimensional tolerance.
    !> npts - Number of intersection points. This should be two
    !> if a valid intersection segment was found, zero otherwise.
    !> pnts - Intersection points.
    
    !f2py intent(in) t1, t2, tol, eps
    !f2py intent(out) npts, pnts
    
    ! Input
    double precision, intent(in) :: tol, eps
    double precision, intent(in) :: t1(0:2, 0:2), t2(0:2, 0:2)
    
    ! Ouput
    integer, intent(out) :: npts
    double precision, intent(out) :: pnts(0:1, 0:2)
    
    ! Working
    logical :: test
    integer :: i
    double precision :: d0, di, t, u, v
    double precision :: vn1(0:2), vn2(0:2), verts(0:2, 0:2)
    double precision :: edges(0:2, 0:2), e1(0:2), e2(0:2)
    double precision :: v0(0:2), v1(0:2), v2(0:2)
    double precision ::  p0(0:2), ray(0:2), p12(0:2)
                        
    npts = 0
    pnts(:, :) = 0.0d0
    
    e1 = t1(1, :) - t1(0, :)
    e2 = t1(2, :) - t1(0, :)
    vn1 = cross(e1, e2)
    e1 = t2(1, :) - t2(0, :)
    e2 = t2(2, :) - t2(0, :)
    vn2 = cross(e1, e2)
    
    test = .false.
    d0 = dot_product(vn1, t2(0, :) - t1(0, :))
    do i = 1, 2
        di = dot_product(vn1, t2(i, :) - t1(0, :))
        if (d0 * di .le. 0.0d0) then
            test = .true.
            exit
        end if
    end do
    if (test .eqv. .false.) then
        return
    end if
    
    test = .false.
    d0 = dot_product(vn2, t1(0, :) - t2(0, :))
    do i = 1, 2
        di = dot_product(vn2, t1(i, :) - t2(0, :))
        if (d0 * di .le. 0.0d0) then
            test = .true.
            exit
        end if
    end do
    if (test .eqv. .false.) then
        return
    end if
    
    ! Intersect edges of triangle 1 with triangle 2.
    verts(:, :) = 0.0d0
    edges(:, :) = 0.0d0
    edges(0, :) = t1(1, :) - t1(0, :)
    verts(0, :) = t1(0, :)
    edges(1, :) = t1(2, :) - t1(0, :)
    verts(1, :) = t1(0, :)
    edges(2, :) = t1(2, :) - t1(1, :)
    verts(2, :) = t1(1, :)
    v0 = t2(0, :)
    v1 = t2(1, :)
    v2 = t2(2, :)
    do i = 0, 2
        p0 = verts(i, :)
        ray = edges(i, :)
        call intersect_triangle_ray(v0, v1, v2, p0, ray, eps, .false., &
                                    t, u, v, test)
        if ((test .eqv. .true.) .and. &
            (t .ge. -eps) .and. (t .le. 1.0d0 + eps)) then
            pnts(npts, :) = p0 + t * ray
            npts = npts + 1
            if (npts .eq. 2) then
                p12 = pnts(1, :) - pnts(0, :)
                di = norm(p12)
                if ((di .le. tol) .and. (i .lt. 2)) then
                    npts = npts - 1
                elseif (di .gt. tol) then
                    return
                elseif ((di .le. tol) .and. (i .eq. 2)) then
                    npts = npts - 1
                end if
            end if
        end if
    end do
    
    ! Intersect edges of triangle 2 with triangle 1.
    edges(0, :) = t2(1, :) - t2(0, :)
    verts(0, :) = t2(0, :)
    edges(1, :) = t2(2, :) - t2(0, :)
    verts(1, :) = t2(0, :)
    edges(2, :) = t2(2, :) - t2(1, :)
    verts(2, :) = t2(1, :)
    v0 = t1(0, :)
    v1 = t1(1, :)
    v2 = t1(2, :)
    do i = 0, 2
        p0 = verts(i, :)
        ray = edges(i, :)
        call intersect_triangle_ray(v0, v1, v2, p0, ray, eps, &
                                    .false., t, u, v, test)
        if ((test .eqv. .true.) .and. &
            (t .ge. -eps) .and. (t .le. 1.0d0 + eps)) then
            pnts(npts, :) = p0 + t * ray
            npts = npts + 1
            if (npts .eq. 2) then
                p12 = pnts(1, :) - pnts(0, :)
                di = norm(p12)
                if ((di .le. tol) .and. (i .lt. 2)) then
                    npts = npts - 1
                elseif (di .gt. tol) then
                    return
                elseif ((di .le. tol) .and. (i .eq. 2)) then
                    npts = npts - 1
                end if
            end if
        end if
    end do
    ! No intersection found.
    npts = 0
    pnts(:, :) = 0.0d0
end subroutine intersect_triangles

subroutine intersect_triangle_plane(t, p0, pnorm, gtol, itol, ni, pi)
    !> Find the intersection segment of a triangle and a plane.
    !> t - Vertices of triangle.
    !> p0 - Origin of plane.
    !> pnorm - Unit normal vector of plane.
    !> gtol - Geometric tolerance.
    !> itol - Intersection tolerance for checking minimum length of
    !> intersection segment.
    !> ni - Number of intersection points. This should be two if a valid
    !> intersection segment was found, zero otherwise.
    !> pi - Intersection points.
    
    !f2py intent(in) t, p0, pnorm, gtol, itol
    !f2py intent(out) ni, pi
    
    !Input
    double precision, intent(in) :: gtol, itol
    double precision, intent(in) :: t(0:2, 0:2), p0(0:2), pnorm(0:2)
                                    
    ! Output
    integer, intent(out) :: ni
    double precision, intent(out) :: pi(0:1, 0:2)
    
    ! Working
    logical :: test
    integer :: i
    double precision :: d0, d1, di, denom, ui, mag
    double precision :: verts(0:2, 0:2), edges(0:2, 0:2), p(0:2), v(0:2)
    double precision :: v1(0:2), v2(0:2), tri_norm(0:2)
    
    ni = 0
    pi(:, :) = 0.0d0
    
    ! Check for parallel triangle and plane. Return no intersection if
    ! they're parallel.
    v1 = t(1, :) - t(0, :)
    v2 = t(2, :) - t(0, :)
    tri_norm = cross(v1, v2)
    mag = norm(tri_norm)
    tri_norm = tri_norm / mag
    if (1.0d0 - dabs(dot_product(pnorm, tri_norm)) .le. 1.0d-12) then
        return
    end if
    
    ! Test for possible intersection.
    test = .false.
    d0 = dot_product(pnorm, t(0, :) - p0)
    do i = 1, 2
        di = dot_product(pnorm, t(i, :) - p0)
        if (d0 * di .le. 0.0d0) then
            test = .true.
            exit
        end if
    end do
    if (test .eqv. .false.) then
        return
    end if
    
    ! Build array of vertices and edges.
    edges(0, :) = t(1, :) - t(0, :)
    verts(0, :) = t(0, :)
    edges(1, :) = t(2, :) - t(0, :)
    verts(1, :) = t(0, :)
    edges(2, :) = t(2, :) - t(1, :)
    verts(2, :) = t(1, :)
    
    ! Intersect edges with plane.
    do i = 0, 2
        d0 = dot_product(pnorm, verts(i, :) - p0)
        d1 = dot_product(pnorm, verts(i, :) + edges(i, :) - p0)
        if ((dabs(d0) .le. gtol) .and. (dabs(d1) .le. gtol)) then
            pi(0, :) = verts(i, :)
            pi(1, :) = verts(i, :) + edges(i, :)
            ni = 2
            return
        end if
        if (d0 * d1 .le. 0.0d0) then
            ! Find intersection
            denom = dot_product(pnorm, edges(i, :))
            ui = dot_product(pnorm, p0 - verts(i, :)) / denom
            if (ui .lt. 0.0d0) then
                ui = 0.0d0
            end if
            if (ui .gt. 1.0d0) then
                ui = 1.0d0
            end if
            p = verts(i, :) + ui * edges(i, :)
            pi(ni, :) = p
            ni = ni + 1
            if (ni .eq. 2) then
                v = pi(1, :) - pi(0, :)
                di = norm(v)
                if ((di .le. itol) .and. (i .lt. 2)) then
                    ni = ni - 1
                elseif ((di .le. itol) .and. (i .eq. 2)) then
                    ni = 0
                    pi(:, :) = 0.0d0
                    return
                else
                    return
                end if
            end if
        end if
    end do
    ni = 0
end subroutine intersect_triangle_plane

end module intersect_triangle