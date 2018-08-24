module intersect_curve
use config, only: ftol
use math
use generic_list
use geom_utils
use evaluate
use modify
use divide
use bounding_box
use intersect_bbox
use geom_results
use nelder_mead, only: simplex
implicit none

private :: ci_data, ci_ptr

! Curve intersection data.
type ci_data
    double precision :: u0
    double precision :: u1
    double precision :: t, u, v
end type ci_data

type ci_ptr
    type(ci_data), pointer :: ci
end type ci_ptr

contains

subroutine intersect_curve_curve(n1, p1, uk1, cpw1, n2, p2, uk2, cpw2, itol, flag)
    !> Intersect two curves.
    !> n1 - Numer of control points - 1 for curve 1.
    !> p1 - Degree for curve 1.
    !> uk1 - Knot vector for curve 1.
    !> cpw1 - Control points for curve 1.
    !> n2 - Numer of control points - 1 for curve 2.
    !> p2 - Degree for curve 2.
    !> uk2 - Knot vector for curve 2.
    !> cpw2 - Control points for curve 2.
    !> itol - Intersection tolerance.
    !> flag - Status flag.
    
    !f2py intent(in) n1, p1, uk1, cpw1, n2, p2, uk2, cpw2, itol
    !f2py intent(out) flag
    !f2py depend(n1, p1) uk1
    !f2py depend(n2, p2) uk2
    !f2py depend(n1) cpw1
    !f2py depend(n2) cpw2
    
    ! Input
    integer, intent(in) :: n1, p1, n2, p2
    double precision, intent(in) :: uk1(0:n1 + p1 + 1)
    double precision, intent(in) :: uk2(0:n2 + p2 + 1)
    double precision, intent(in) :: cpw1(0:n1, 0:3)
    double precision, intent(in) :: cpw2(0:n2, 0:3)
    double precision, intent(in) :: itol

    ! Output
    integer, intent(out) :: flag
    
    ! Working
    logical :: is_cand
    integer :: i, j
    integer :: nsol
    
    integer :: nbez1, nbez2
    double precision :: ucheck, di, a1, b1, a2, b2
    double precision :: qw1(0:n1 + p1 + 1, 0:p1, 0:3)
    double precision :: qw2(0:n2 + p2 + 1, 0:p2, 0:3)
    double precision :: ab1(0:n1 + p1 + 1, 0:1)
    double precision :: ab2(0:n2 + p2 + 1, 0:1)
    double precision :: qwi1(0:p1, 0:3), abi1(0:1)
    double precision :: qwi2(0:p2, 0:3), abi2(0:1)
    
    logical :: unique, test1, test2
    double precision :: p0(0:2), pi1(0:2), pi2(0:2), vi(0:2), pi12(0:2)
    double precision :: p10(0:2), p11(0:2), p20(0:2), p21(0:2)
    
    ! Nelder-Mead
    logical :: success
    double precision :: tol
    double precision :: x(2)
    
    ! Generic list for storing results.
    type(list_node_t), pointer :: cci_list => null()
    type(ci_ptr) :: cci_ptr
    
    ! Initialize
    nsol = 0
    call list_init(cci_list)
    
    ! Step 1: Check for possible intersection.
    is_cand = candidate_intersect(n1, cpw1, n2, cpw2)
    if (is_cand .eqv. .false.) then
        flag = 0
        call list_free(cci_list)
        call set_empty_results()
        return
    end if
    
    ! Step 2: Decompose curves into Bezier segments.
    call decompose_curve(n1, p1, uk1, cpw1, nbez1, qw1, ab1)
    call decompose_curve(n2, p2, uk2, cpw2, nbez2, qw2, ab2)
    
    ! Step 3: Perform recursive subdivision on candidate curves.
    do i = 0, nbez1 - 1
        qwi1 = qw1(i, :, :)
        abi1 = ab1(i, :)
        do j = 0, nbez2 - 1
            qwi2 = qw2(j, :, :)
            abi2 = ab2(j, :)
            is_cand = candidate_intersect(p1, qwi1, p2, qwi2)
            if (is_cand) then
                call subdivide(p1, qwi1, abi1, p2, qwi2, abi2)
            end if
        end do
    end do
    
    ! Return if no intersections are found.
    if (nsol .eq. 0) then
        flag = 0
        call list_free(cci_list)
        call set_empty_results()
        return
    end if
    
    ! Allocate result arrays.
    call reset_results()
    allocate(points_(nsol + 4, 3))
    allocate(params2d_(nsol + 4, 2))
    points_(:, :) = 0.0d0
    params2d_(:, :) = 0.0d0
    
    ! Manually check curve endpoint intersection to improve robustness.
    call curve_point(n1, p1, uk1, cpw1, uk1(p1), p10)
    call curve_point(n1, p1, uk1, cpw1, uk1(n1 + 1), p11)
    call curve_point(n2, p2, uk2, cpw2, uk2(p2), p20)
    call curve_point(n2, p2, uk2, cpw2, uk2(n2 + 1), p21)
    test1 = are_points_equal(p10, p20, itol)
    test2 = are_points_equal(p10, p21, itol)
    if (test1) then
        pi12 = 0.5d0 * (p10 + p20)
        call add_pnt(uk1(p1), uk2(p2), pi12)
    elseif (test2) then
        pi12 = 0.5d0 * (p10 + p21)
        call add_pnt(uk1(p1), uk2(n2 + 1), pi12)
    end if    
    
    test1 = are_points_equal(p11, p20, itol)
    test2 = are_points_equal(p11, p21, itol)
    if (test1) then
        pi12 = 0.5d0 * (p11 + p20)
        call add_pnt(uk1(n1 + 1), uk2(p2), pi12)
    elseif (test2) then
        pi12 = 0.5d0 * (p11 + p21)
        call add_pnt(uk1(n1 + 1), uk2(n2 + 1), pi12)
    end if

    ! Step 4: Refine intersection points.
    a1 = uk1(0)
    b1 = uk1(n1 + p1 + 1)
    a2 = uk2(0)
    b2 = uk2(n2 + p2 + 1)
    tol = itol / 100.0d0
    do i = 1, nsol
        ! Get next
        cci_list = list_next(cci_list)
        cci_ptr = transfer(list_get(cci_list), cci_ptr)
        x(1) = cci_ptr%ci%u0
        x(2) = cci_ptr%ci%u1
        ! Call Nelder-Mead.
        call simplex(obj, 2, x, tol, success)
        ! Check parameter domains.
        ucheck = x(1)
        call check_param(n1, p1, uk1, ucheck)
        x(1) = ucheck
        ucheck = x(2)
        call check_param(n2, p2, uk2, ucheck)
        x(2) = ucheck
        ! Get point and add to results only if unique.
        call curve_point(n1, p1, uk1, cpw1, x(1), pi1)
        call curve_point(n2, p2, uk2, cpw2, x(2), pi2)
        vi = pi1 - pi2
        di = norm(vi)
        if (di .le. itol) then
            pi12 = 0.50d0 * (pi1 + pi2)
            call add_pnt(x(1), x(2), pi12)
        end if
    end do
    
    ! Return results.
    call list_free(cci_list)
    flag = 1
    
    ! Internal subroutines.
    contains
    
    subroutine add_pnt(u1, u2, pi)
        ! Add point to intersection results.
        
        ! Input
        double precision, intent(in) :: u1, u2
        double precision, intent(in) :: pi(0:2)
        
        ! Working
        logical :: unique, test
        integer :: i
        double precision :: p0(0:2)
        
        unique = .true.
        do i = 1, npts_
            p0 = points_(i, :)
            test = are_points_equal(pi, p0, itol)
            if (test) then
                unique = .false.
                exit
            end if
        end do
        if (unique) then
            npts_ = npts_ + 1
            points_(npts_, :) = pi
            params2d_(npts_, :) = (/ u1, u2 /)
        end if
        
    end subroutine add_pnt
    
    subroutine obj(nx, x, fx)
        integer, intent(in) :: nx
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: fx
        
        ! Working
        double precision :: factor
        double precision :: r(0:2)
        double precision :: pi1(0:2), pi2(0:2)
        
        ! Penalize objective function if variable is outside domain.
        factor = 1.0d0
        if ((x(1) .lt. a1) .or. (x(1) .gt. b1))  then
            factor = 1000.0d0
        end if
        if ((x(2) .lt. a2) .or. (x(2) .gt. b2))  then
            factor = 1000.0d0
        end if
        
        ! Evaluate distance between points.
        call curve_point(n1, p1, uk1, cpw1, x(1), pi1)
        call curve_point(n2, p2, uk2, cpw2, x(2), pi2)
        r = pi1 - pi2
        fx = norm(r) * factor
        
    end subroutine obj
    
    recursive subroutine subdivide(n1, cpw1, ab1, n2, cpw2, ab2)
        ! Recursive subdivision.
        integer, intent(in) :: n1, n2
        double precision, intent(in) :: cpw1(0:n1, 0:3), ab1(0:1), &
                                        cpw2(0:n2, 0:3), ab2(0:1)
                                        
        logical :: flat1, flat2, candidate
        double precision :: cp1(0:n1, 0:2), w1(0:n1), &
                            cp2(0:n2, 0:2), w2(0:n2), &
                            qw11(0:n1, 0:3), qw12(0:n1, 0:3), &
                            qw21(0:n2, 0:3), qw22(0:n2, 0:3), &
                            ab11(0:1), ab12(0:1), ab21(0:1), ab22(0:1)
        
        call dehomogenize_array1d(n1, cpw1, cp1, w1)
        call dehomogenize_array1d(n2, cpw2, cp2, w2)
        
        call is_curve_flat(n1, cp1, ftol, flat1)
        call is_curve_flat(n2, cp2, ftol, flat2)
        
        if (flat1 .and. flat2) then
            call intersect(n1, cp1, ab1, n2, cp2, ab2)
        else
            call split_bezier_curve(cpw1, n1, 0.50d0, ab1(0), ab1(1), &
                                    qw11, ab11(0), ab11(1), &
                                    qw12, ab12(0), ab12(1))
            call split_bezier_curve(cpw2, n2, 0.50d0, ab2(0), ab2(1), &
                                    qw21, ab21(0), ab21(1), &
                                    qw22, ab22(0), ab22(1))
            
            candidate = candidate_intersect(n1, qw11, n2, qw21)
            if (candidate) then
                call subdivide(n1, qw11, ab11, n2, qw21, ab21)
            end if
            candidate = candidate_intersect(n1, qw11, n2, qw22)
            if (candidate) then
                call subdivide(n1, qw11, ab11, n2, qw22, ab22)
            end if
            candidate = candidate_intersect(n1, qw12, n2, qw21)
            if (candidate) then
                call subdivide(n1, qw12, ab12, n2, qw21, ab21)
            end if
            candidate = candidate_intersect(n1, qw12, n2, qw22)
            if (candidate) then
                call subdivide(n1, qw12, ab12, n2, qw22, ab22)
            end if
        end if
        
    end subroutine subdivide
    
    subroutine intersect(n1, cp1, ab1, n2, cp2, ab2)
        ! Intersect lines.
        integer, intent(in) :: n1, n2
        double precision, intent(in) :: cp1(0:n1, 0:2), ab1(0:1), &
                                        cp2(0:n2, 0:2), ab2(0:1)
                                        
        double precision :: p1(0:2), p2(0:2), p3(0:2), p4(0:2), &
                            d1343, d4321, d1321, d4343, d2121, &
                            numer, denom, u12, u34, ui1, ui2, &
                            p12(0:2), p34(0:2)
                            
        p1 = cp1(0, :)
        p2 = cp1(n1, :)
        p3 = cp2(0, :)
        p4 = cp2(n2, :)
        d4321 = dot_product(p4 - p3, p2 - p1)
        d1321 = dot_product(p1 - p3, p2 - p1)
        d4343 = dot_product(p4 - p3, p4 - p3)
        d2121 = dot_product(p2 - p1, p2 - p1)
        denom = d2121 * d4343 - d4321 * d4321
        if (dabs(denom) .le. 1.0d-12) then
            ui1 = local_to_global_param(ab1(0), ab1(1), 0.50d0)
            ui2 = local_to_global_param(ab2(0), ab2(1), 0.50d0)
            ! Store result.
            allocate(cci_ptr%ci)
            cci_ptr%ci%u0 = ui1
            cci_ptr%ci%u1 = ui2
            call list_insert(cci_list, transfer(cci_ptr, list_data))
            nsol = nsol + 1
        else
            d1343 = dot_product(p1 - p3, p4 - p3)
            numer = d1343 * d4321 - d1321 * d4343
            u12 = numer / denom
            u34 = (d1343 + u12 * d4321) / d4343
            ui1 = local_to_global_param(ab1(0), ab1(1), u12)
            ui2 = local_to_global_param(ab2(0), ab2(1), u34)
            ! Store result.
            allocate(cci_ptr%ci)
            cci_ptr%ci%u0 = ui1
            cci_ptr%ci%u1 = ui2
            call list_insert(cci_list, transfer(cci_ptr, list_data))
            nsol = nsol + 1
        end if
        
    end subroutine intersect
    
    function candidate_intersect(n1, cpw1, n2, cpw2) result(candidate)
        ! Check if the two curves intersect.
        integer, intent(in) :: n1, n2
        double precision, intent(in) :: cpw1(0:n1, 0:3), cpw2(0:n2, 0:3)
        logical :: candidate
        
        double precision :: bbox1(0:2, 0:1), bbox2(0:2, 0:1)
        
        bbox1 = curve_bbox(n1, cpw1)
        bbox2 =  curve_bbox(n2, cpw2)
        candidate = bboxes_intersect(3, bbox1, bbox2, itol)
        
    end function candidate_intersect

end subroutine intersect_curve_curve

subroutine intersect_curve_plane(n, p, uk, cpw, origin, normal, itol, flag)
    !> Intersect a curve and a plane.
    !> n - Numer of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> origin - Origin point of plane.
    !> normal - Unit normal vector of plane.
    !> itol - Intersection tolerance.
    !> flag - Status flag.
    
    !f2py intent(in) n, p, uk, cpw, origin, normal, itol
    !f2py intent(out) flag
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    double precision, intent(in) :: origin(0:2), normal(0:2)
    double precision, intent(in) :: itol
    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    logical :: is_cand, unique
    integer :: i, j
    integer :: nsol, nbez
    double precision :: di, ucheck, a, b
    double precision :: pnorm(0:2)
    double precision :: qw(0:n + p + 1, 0:p, 0:3), ab(0:n + p + 1, 0:1)
    double precision :: qwi(0:p, 0:3), abi(0:1)
    double precision :: p0(0:2), p1(0:2), vi(0:2)
    
    ! Nelder-Mead
    logical :: success
    double precision :: tol
    double precision :: x(1)
    
    ! Generic list for storing results.
    type(list_node_t), pointer :: cpi_list => null()
    type(ci_ptr) :: cpi_ptr
    
    ! Initialize
    nsol = 0
    call list_init(cpi_list)
    
    ! Make sure a unit normal vector is used.
    di = norm(normal)
    pnorm = normal / di
    
    ! Step 1: Check bounding box and plane for possible intersection.
    is_cand = candidate_intersect(n, cpw)
    if (is_cand .eqv. .false.) then
        flag = 0
        call list_free(cpi_list)
        call set_empty_results()
        return
    end if
    
    ! Step 2: Decompose curve into Bezier segments.
    call decompose_curve(n, p, uk, cpw, nbez, qw, ab)
    
    ! Step 3: Perform recursive subdivision on candidate curves.
    do i = 0, nbez - 1
        qwi = qw(i, :, :)
        abi = ab(i, :)
        is_cand = candidate_intersect(p, qwi)
        if (is_cand) then
            call subdivide(p, qwi, abi)
        end if
    end do
    
    ! Return if no intersections are found.
    if (nsol .eq. 0) then
        flag = 0
        call list_free(cpi_list)
        call set_empty_results()
        return
    end if
    
    ! Allocate result arrays.
    call reset_results()
    allocate(points_(nsol, 3))
    allocate(params1d_(nsol))
    points_(:, :) = 0.0d0
    params1d_(:) = 0.0d0
    ! Step 4: Refine intersection points.
    a = uk(0)
    b = uk(n + p + 1)
    tol = itol / 100.0d0
    do i = 1, nsol
        ! Get next
        cpi_list = list_next(cpi_list)
        cpi_ptr = transfer(list_get(cpi_list), cpi_ptr)
        x(1) = cpi_ptr%ci%u0
        ! Call Nelder-Mead.
        call simplex(obj, 1, x, tol, success)
        ! Check parameter domains.
        ucheck = x(1)
        call check_param(n, p, uk, ucheck)
        x(1) = ucheck
        ! Get point and add to results only if unique.
        call curve_point(n, p, uk, cpw, x(1), p0)
        ! Project to plane and use midpoint.
        p1 = p0 - pnorm * dot_product(pnorm, p0 - origin)
        p0 = 0.50d0 * (p0 + p1)
        unique = .true.
        do j = 1, npts_
            p1 = points_(j, :)
            vi = p1 - p0
            di = norm(vi)
            if (di .le. itol) then
                unique = .false.
                exit
            end if
        end do
        if (unique) then
            npts_ = npts_ + 1
            points_(npts_, :) = p0
            params1d_(npts_) = x(1)
        end if
    end do
    
    ! Return results.
    call list_free(cpi_list)
    flag = 1
    
    ! Internal subroutines.
    contains
    
    subroutine obj(nx, x, fx)
        ! Objective function for point refinement.
    
        ! Input/Output
        integer, intent(in) :: nx
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: fx
        
        ! Working
        double precision :: factor
        double precision :: pi(0:2), pj(0:2), r(0:2)
        
        ! Penalize objective function if variable is outside domain.
        factor = 1.0d0
        if ((x(1) .lt. a) .or. (x(1) .gt. b))  then
            factor = 1000.0d0
        end if
        
        ! Evaluate curve point and project to plane.
        call curve_point(n, p, uk, cpw, x(1), pi)
        pj = pi - pnorm * dot_product(pnorm, pi - origin)
        r = pi - pj
        fx = norm(r) * factor
        
    end subroutine obj
    
    recursive subroutine subdivide(n, cpw, ab)
        ! Recursive subdivision.
        integer, intent(in) :: n
        double precision, intent(in) :: cpw(0:n, 0:3), ab(0:1)
                                        
        logical :: flat, candidate
        double precision :: cp(0:n, 0:2), w(0:n), &
                            qw1(0:n, 0:3), qw2(0:n, 0:3), &
                            ab1(0:1), ab2(0:1)
        
        call dehomogenize_array1d(n, cpw, cp, w)
        
        call is_curve_flat(n, cp, ftol, flat)
        
        if (flat) then
            call intersect(n, cp, ab)
        else
            call split_bezier_curve(cpw, n, 0.50d0, ab(0), ab(1), &
                                    qw1, ab1(0), ab1(1), &
                                    qw2, ab2(0), ab2(1))

            candidate = candidate_intersect(n, qw1)
            if (candidate) then
                call subdivide(n, qw1, ab1)
            end if
            candidate = candidate_intersect(n, qw2)
            if (candidate) then
                call subdivide(n, qw2, ab2)
            end if
        end if
        
    end subroutine subdivide
    
    subroutine intersect(n, cp, ab)
        ! Intersect line and plane.
        integer, intent(in) :: n
        double precision, intent(in) :: cp(0:n, 0:2), ab(0:1)
                                        
        double precision :: p1(0:2), p2(0:2), ui, ug, denom
                            
        p1 = cp(0, :)
        p2 = cp(n, :)
        denom = dot_product(pnorm, p2 - p1)
        if (dabs(denom) .le. 1.0d-12) then
            ug = local_to_global_param(ab(0), ab(1), 0.50d0)
            ! Store result.
            allocate(cpi_ptr%ci)
            cpi_ptr%ci%u0 = ug
            call list_insert(cpi_list, transfer(cpi_ptr, list_data))
            nsol = nsol + 1
        else
            ui = dot_product(pnorm, origin - p1) / denom
            ug = local_to_global_param(ab(0), ab(1), ui)
            ! Store result.
            allocate(cpi_ptr%ci)
            cpi_ptr%ci%u0 = ug
            call list_insert(cpi_list, transfer(cpi_ptr, list_data))
            nsol = nsol + 1
        end if
    end subroutine intersect
    
    function candidate_intersect(n, cpw) result(candidate)
    ! Check for possible curve-plane intersection.
        integer, intent(in) :: n
        double precision, intent(in) :: cpw(0:n, 0:3)
        
        logical :: candidate
        
        double precision :: bbox(0:2, 0:1)
        
        bbox = curve_bbox(n, cpw)
        candidate = bbox_intersects_plane(bbox, origin, pnorm, itol)
        
    end function candidate_intersect

end subroutine intersect_curve_plane

subroutine intersect_curve_surface(cn, cp, cuk, cpw, sn, sp, suk, sm, sq, &
                                   svk, spw, itol, flag)
    !> Intersect a curve and a surface.
    !> cn - Number of curve control points - 1.
    !> cp - Degree of curve.
    !> cuk - Knot vector of curve.
    !> cpw - Curve control points.
    !> sn - Number of surface control points - 1 in u-direction.
    !> sp - Degree of surface in u-direction.
    !> suk - Knot vector of surface in u-direction.
    !> sm - Number of surface control points - 1 in v-direction.
    !> sq - Degree of surface in v-direction.
    !> svk - Knot vector of surface in v-direction.
    !> spw - Surface control points.
    !> itol - Intersection tolerance.
    !> flag - Status flag.
    
    !f2py intent(in) cn, cp, cuk, cpw, sn, sp, suk, sm, sq, svk, spw, itol
    !f2py intent(out) flag
    !f2py depend(cn, cp) cuk
    !f2py depend(cn) cpw
    !f2py depend(sn, sp) suk
    !f2py depend(sm, sq) svk
    !f2py depend(sn, sm) spw
    
    ! Input
    integer, intent(in) :: cn, cp, sn, sp, sm, sq
    double precision, intent(in) :: cuk(0: cn + cp + 1)
    double precision, intent(in) :: cpw(0:cn, 0:3)
    double precision, intent(in) :: suk(0:sn + sp + 1)
    double precision, intent(in) :: svk(0:sm + sq + 1)
    double precision, intent(in) :: spw(0:sn, 0:sm, 0:3)
    double precision, intent(in) :: itol
    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    integer :: nsol
    
    logical :: is_cand
    integer :: i, j, k, cbez, sbez_u
    double precision :: cqw(0:cn + cp + 1, 0:cp, 0:3)
    double precision :: cab(0:cn + cp + 1, 0:1)
    double precision :: sqw_u(0:sn + sp + 1, 0:sp, 0:sm, 0:3)
    double precision :: sab_u(0:sn + sp + 1, 0:1)
    
    double precision :: suk_ui(0:2 * sp - 1)
    double precision :: sqw_ui(0:sp, 0:sm, 0:3)
    
    integer :: sbez_v
    double precision :: sqw_v(0:sm + sq + 1, 0:sp, 0:sq, 0:3)
    double precision :: sab_v(0:sm + sq + 1, 0:1)
    
    double precision :: sqwi(0:sp, 0:sq, 0:3)
    double precision :: sabi(0:1, 0:1)
    
    double precision :: cqwi(0:cp, 0:3)
    double precision :: cabi(0:1)
    
    logical :: unique
    double precision :: di, ucheck, vcheck, a, b, au, bu, av, bv
    double precision :: pc(0:2), ps(0:2), vi(0:2), p0(0:2), pi(0:2)
    
    
    ! Nelder-Mead
    logical :: success
    double precision :: tol
    double precision :: x(3)
    
    ! Generic list for storing results.
    type(list_node_t), pointer :: csi_list => null()
    type(ci_ptr) :: csi_ptr
    
    ! Initialize
    nsol = 0
    call list_init(csi_list)
    
    ! Step 1: Check for possible intersection.
    is_cand = candidate_intersect(cn, cpw, sn, sm, spw)
    if (is_cand .eqv. .false.) then
        flag = 0
        call list_free(csi_list)
        call set_empty_results()
        return
    end if
    
    ! Step 2: Decompose curve and surface into Bezier segments.
    call decompose_curve(cn, cp, cuk, cpw, cbez, cqw, cab)
    
    ! Decompose surface into u-direction strips.
    call decompose_surface(sn, sp, suk, sm, sq, svk, spw, 'u', sn + sp + 1, &
                           sp, sm, sbez_u, sqw_u, sab_u)
    
    ! Step 3: Perform recursive subdivision to find potential intersection points.
    do i = 0, sbez_u - 1
        ! Decompose each strip in v-direction.
        sqw_ui = sqw_u(i, :, :, :)
        suk_ui(0:sp) = sab_u(i, 0)
        suk_ui(sp + 1:) = sab_u(i, 1)
        call decompose_surface(sp, sp, suk_ui, sm, sq, svk, sqw_ui, 'v', &
                               sm + sq + 1, sp, sq, sbez_v, sqw_v, sab_v)
        ! Check each Bezier surface patch with Bezier curve segment.
        do j = 0, sbez_v - 1
            sqwi = sqw_v(j, :, :, :)
            sabi(0, :) = sab_u(i, :)
            sabi(1, :) = sab_v(j, :)
            do k = 0, cbez - 1
                cqwi = cqw(k, :, :)
                cabi = cab(k, :)
                is_cand = candidate_intersect(cp, cqwi, sp, sq, sqwi)
                if (is_cand) then
                    call subdivide(cp, cqwi, cabi, sp, sq, sqwi, sabi)
                end if
            end do
        end do
    end do
    
    ! Return if not candidate points are found.
    if (nsol .eq. 0) then
        flag = 0
        call list_free(csi_list)
        call set_empty_results()
        return
    end if
    
    ! Allocate result arrays.
    call reset_results()
    allocate(points_(nsol,  3))
    allocate(params2d_(nsol, 3))
    points_(:, :) = 0.0d0
    params2d_(:, :) = 0.0d0
    
    ! Step 4: Refine points.
    a = cuk(0)
    b = cuk(cn + cp + 1)
    au = suk(0)
    bu = suk(sn + sp + 1)
    av = svk(0)
    bv = svk(sm + sq + 1)
    tol = itol / 100.0d0
    do i = 1, nsol
        ! Get next
        csi_list = list_next(csi_list)
        csi_ptr = transfer(list_get(csi_list), csi_ptr)
        x(1) = csi_ptr%ci%t
        x(2) = csi_ptr%ci%u
        x(3) = csi_ptr%ci%v
        ! Call Nelder-Mead
        call simplex(obj, 3, x, tol, success)
        ! Check parameter domains.
        ucheck = x(1)
        call check_param(cn, cp, cuk, ucheck)
        x(1) = ucheck
        ucheck = x(2)
        call check_param(sn, sp, suk, ucheck)
        x(2) = ucheck
        vcheck = x(3)
        call check_param(sm, sq, svk, vcheck)
        x(3) = vcheck
        ! Add to results if point is unique.
        call curve_point(cn, cp, cuk, cpw, x(1), pc)
        call surface_point(sn, sp, suk, sm, sq, svk, spw, x(2), x(3), ps)
        vi = pc - ps
        di = norm(vi)
        if (di .le. itol) then
            pi = 0.50d0 * (pc + ps)
            unique = .true.
            do j = 1, npts_
                p0 = points_(j, :)
                vi = pi - p0
                di = norm(vi)
                if (di .le. itol) then
                    unique = .false.
                    exit
                end if
            end do
            if (unique) then
                npts_ = npts_ + 1
                points_(npts_, :) = pi
                params2d_(npts_, :) = x
            end if
        end if
    end do
    
    ! Normal completion.
    call list_free(csi_list)
    flag = 1
    
    ! Internal subroutines.
    contains
    
    subroutine obj(nx, x, fx)
        ! Subroutine for point refinment.
        ! Input/Output
        integer, intent(in) :: nx
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: fx
        
        ! Working
        double precision :: factor
        double precision :: r(0:2)
        double precision :: pc(0:2)
        double precision :: ps(0:2)
        
        ! Penalize objective function if variable is outside domain.
        factor = 1.0d0
        if ((x(1) .lt. a) .or. (x(1) .gt. b))  then
            factor = 1000.0d0
        end if
        if ((x(2) .lt. au) .or. (x(2) .gt. bu))  then
            factor = 1000.0d0
        end if
        if ((x(3) .lt. av) .or. (x(3) .gt. bv))  then
            factor = 1000.0d0
        end if
        
        ! Evaluate curve and surface point.
        call curve_point(cn, cp, cuk, cpw, x(1), pc)
        call surface_point(sn, sp, suk, sm, sq, svk, spw, x(2), x(3), ps)
        r = pc - ps
        fx = norm(r) * factor
        
    end subroutine obj
    
    recursive subroutine subdivide(cn, cpw, cab, sn, sm, spw, sab)
        ! Recursive subdivision.
        ! Input
        integer, intent(in) :: cn, sn, sm
        double precision, intent(in) :: cpw(0:cn, 0:3)
        double precision, intent(in) :: cab(0:1)
        double precision, intent(in) :: spw(0:sn, 0:sm, 0:3)
        double precision, intent(in) :: sab(0:1, 0:1)
                                        
        ! Working
        logical :: flat1, flat2, candidate
        double precision :: cpp(0:cn, 0:2), wc(0:cn)
        double precision :: spp(0:sn, 0:sm, 0:2), ws(0:sn, 0:sm)
        double precision :: cqw1(0:cn, 0:3), cqw2(0:cn, 0:3)
        double precision :: cab1(0:1), cab2(0:1)
        double precision :: sqw1(0:sn, 0:sm, 0:3), sab1(0:1, 0:1)
        double precision :: sqw2(0:sn, 0:sm, 0:3), sab2(0:1, 0:1)
        double precision :: sqw3(0:sn, 0:sm, 0:3), sab3(0:1, 0:1)
        double precision :: sqw4(0:sn, 0:sm, 0:3), sab4(0:1, 0:1)
        double precision :: sqw5(0:sn, 0:sm, 0:3), sab5(0:1, 0:1)
        double precision :: sqw6(0:sn, 0:sm, 0:3), sab6(0:1, 0:1)
        
        call dehomogenize_array1d(cn, cpw, cpp, wc)
        call dehomogenize_array2d(sn, sm, spw, spp, ws)
        
        call is_curve_flat(cn, cpp, ftol, flat1)
        call is_surface_flat(sn, sm, spp, ftol, flat2)
        
        if (flat1 .and. flat2) then
            call intersect(cn, cpp, cab, sn, sm, spp, sab)
        else
            call split_bezier_curve(cpw, cn, 0.5d0, cab(0), cab(1), &
                                    cqw1, cab1(0), cab1(1), &
                                    cqw2, cab2(0), cab2(1))         
            call split_bezier_surface(spw, sn, sm, 0.5d0, -1.0d0, &
                                      sab(0, 0), sab(0, 1), &
                                      sab(1, 0), sab(1, 1), &
                                      sqw1, sqw2, sab1, sab2)
            call split_bezier_surface(sqw1, sn, sm, -1.0d0, 0.5d0, &
                                      sab1(0, 0), sab1(0, 1), &
                                      sab1(1, 0), sab1(1, 1), &
                                      sqw3, sqw4, sab3, sab4)
            call split_bezier_surface(sqw2, sn, sm, -1.0d0, 0.5d0, &
                                      sab2(0, 0), sab2(0, 1), &
                                      sab2(1, 0), sab2(1, 1), &
                                      sqw5, sqw6, sab5, sab6)
            ! Curve 1 and all surfaces.
            candidate = candidate_intersect(cn, cqw1, sn, sm, sqw3)
            if (candidate) then
                call subdivide(cn, cqw1, cab1, sn, sm, sqw3, sab3)
            end if
            candidate = candidate_intersect(cn, cqw1, sn, sm, sqw4)
            if (candidate) then
                call subdivide(cn, cqw1, cab1, sn, sm, sqw4, sab4)
            end if
            candidate = candidate_intersect(cn, cqw1, sn, sm, sqw5)
            if (candidate) then
                call subdivide(cn, cqw1, cab1, sn, sm, sqw5, sab5)
            end if
            candidate = candidate_intersect(cn, cqw1, sn, sm, sqw6)
            if (candidate) then
                call subdivide(cn, cqw1, cab1, sn, sm, sqw6, sab6)
            end if
            ! Curve 2 and all surfaces.
            candidate = candidate_intersect(cn, cqw2, sn, sm, sqw3)
            if (candidate) then
                call subdivide(cn, cqw2, cab2, sn, sm, sqw3, sab3)
            end if
            candidate = candidate_intersect(cn, cqw2, sn, sm, sqw4)
            if (candidate) then
                call subdivide(cn, cqw2, cab2, sn, sm, sqw4, sab4)
            end if
            candidate = candidate_intersect(cn, cqw2, sn, sm, sqw5)
            if (candidate) then
                call subdivide(cn, cqw2, cab2, sn, sm, sqw5, sab5)
            end if
            candidate = candidate_intersect(cn, cqw2, sn, sm, sqw6)
            if (candidate) then
                call subdivide(cn, cqw2, cab2, sn, sm, sqw6, sab6)
            end if
        end if
        
    end subroutine subdivide
    
    subroutine intersect(cn, cpp, cab, sn, sm, spp, sab)
        ! Intersect line-plane.
        integer, intent(in) :: cn, sn, sm
        double precision, intent(in) :: cpp(0:cn, 0:2)
        double precision, intent(in) :: cab(0:1)
        double precision, intent(in) :: spp(0:sn, 0:sm, 0:2)
        double precision, intent(in) :: sab(0:1, 0:1)
        
        double precision :: ti, t, u, v, mag, denom
        double precision :: p0(0:2), vu(0:2), vv(0:2), p1(0:2)
        double precision :: p2(0:2), pnorm(0:2), pi(0:2)
        double precision :: t3d(0:2, 0:2), t2d(0:2, 0:1)
                            
        p0 = spp(0, 0, :)
        vu = spp(sn, 0, :) - p0
        vv = spp(0, sm, :) - p0
        p1 = cpp(0, :)
        p2 = cpp(cn, :)
        pnorm = cross(vu, vv)
        mag = norm(pnorm)
        pnorm = pnorm / mag
        denom = dot_product(pnorm, p2 - p1)
        if (dabs(denom) .le. 1.0d-12) then
            t = local_to_global_param(cab(0), cab(1), 0.50d0)
            u = local_to_global_param(sab(0, 0), sab(0, 1), 0.50d0)
            v = local_to_global_param(sab(1, 0), sab(1, 1), 0.50d0)
            ! Store result.
            allocate(csi_ptr%ci)
            csi_ptr%ci%t = t
            csi_ptr%ci%u = u
            csi_ptr%ci%v = v
            call list_insert(csi_list, transfer(csi_ptr, list_data))
            nsol = nsol + 1
        else
            ti = dot_product(pnorm, p0 - p1) / denom
            t = local_to_global_param(cab(0), cab(1), ti)
            pi = p1 + ti * (p2 - p1)
            t3d(0, :) = spp(0, 0, :)
            t3d(1, :) = spp(sn, 0, :)
            t3d(2, :) = spp(0, sm, :)
            t2d(0, :) = (/ sab(0, 0), sab(1, 0) /)
            t2d(1, :) = (/ sab(0, 1), sab(1, 0) /)
            t2d(2, :) = (/ sab(0, 0), sab(1, 1) /)
            call invert_point_in_triangle(pi, t3d, t2d, u, v)
            ! Store result.
            allocate(csi_ptr%ci)
            csi_ptr%ci%t = t
            csi_ptr%ci%u = u
            csi_ptr%ci%v = v
            call list_insert(csi_list, transfer(csi_ptr, list_data))
            nsol = nsol + 1
        end if
        
    end subroutine intersect
    
    function candidate_intersect(cn, cpw, sn, sm, spw) result(candidate)
        ! Check for possible curve-surface intersection.
        integer, intent(in) :: cn, sn, sm
        double precision, intent(in) :: cpw(0:cn, 0:3)
        double precision, intent(in) :: spw(0:sn, 0:sm, 0:3)
        
        logical :: candidate
        
        double precision :: bbox1(0:2, 0:1), bbox2(0:2, 0:1)
        
        bbox1 = curve_bbox(cn, cpw)
        bbox2 = surface_bbox(sn, sm, spw)
        candidate = bboxes_intersect(3, bbox1, bbox2, itol)
        
    end function candidate_intersect

end subroutine intersect_curve_surface

end module intersect_curve