module project
use config, only: ftol, gtol, atol, ptol
use math
use generic_list
use geom_utils
use evaluate
use modify
use divide
use geom_results
implicit none

private
public :: project_point_to_curve, project_point_to_surface

! Point projection data.
type proj_data
    double precision :: u0
    double precision :: v0
end type proj_data

type proj_ptr
    type(proj_data), pointer :: proj
end type proj_ptr

contains

subroutine project_point_to_curve(pnt, n, p, uk, cpw, ends, all_pnts, flag)
    !> Project a point to a curve.
    !> pnt - Point to project.
    !> n - Numer of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> ends - Option to project to nearest end point of curve if no
    !> all_pnts - Option to return only the nearest point during the
    !> subdivision process (False), or all points (True).
    !> orthogonal projections are found.
    !> flag - Status flag.
    
    !f2py intent(in) pnt, n, p, uk, cpw, ends, all_pnts
    !f2py intent(out) flag
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    
    ! Input
    logical, intent(in) :: ends, all_pnts
    integer, intent(in) :: n, p
    double precision, intent(in) :: pnt(0:2)
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    integer :: i, j
    integer :: nsol
    double precision :: dmin
    
    integer :: nbez
    double precision :: qw(0:n + p + 1, 0:p, 0:3), ab(0:n + p + 1, 0:1)
    
    logical :: is_cand
    double precision :: qwi(0:p, 0:3), abi(0:1)
    
    logical :: unique
    double precision :: di, d0, d1
    double precision :: p0(0:2), p1(0:2), vi(0:2)
    
    ! MINPACK
    integer :: info
    double precision :: ucheck
    double precision :: x(1), fvec(1), fjac(1, 1), wa(7)
    
    ! Generic list for storing results.
    type(list_node_t), pointer :: pc_list => null()
    type(proj_ptr) :: pc_ptr
    
    ! Initialize
    flag = 0
    nsol = 0
    dmin = 0.0d0
    call list_init(pc_list)
    
    ! Step 1: Decompose curve into Bezier segments.
    call decompose_curve(n, p, uk, cpw, nbez, qw, ab)
    
    ! Step 2: Perform recursive subdivision on candidate curves.
    do i = 0, nbez - 1
        qwi = qw(i, :, :)
        abi = ab(i, :)
        is_cand = candidate_curve(p, qwi)
        if (is_cand) then
            call subdivide(p, qwi, abi)
        end if
    end do
    
    ! Allocate result arrays.
    call reset_results()
    allocate(points_(nsol, 3))
    allocate(params1d_(nsol))
    allocate(dproj_(nsol))
    points_(:, :) = 0.0d0
    params1d_(:) = 0.0d0
    dproj_(:) = 0.0d0
    
    ! Step 3: Refine projection points.
    do i = 1, nsol
        ! Get next
        pc_list = list_next(pc_list)
        pc_ptr = transfer(list_get(pc_list), pc_ptr)
        x(1) = pc_ptr%proj%u0
        ! Call MINPACK hybrj driver.
        call hybrj1(obj, 1, x, fvec, fjac, 1, 1.490120d-8, info, wa, 7)
        ! if (info .ne. 1) then
        !     print *, "Unsuccessful solution in Fortran curve projection."
        ! end if
        ! Check parameter.
        ucheck = x(1)
        call check_param(n, p, uk, ucheck)
        x(1) = ucheck
        if (info .eq. 1) then
            ! Get point and add to results only if unique.
            call curve_point(n, p, uk, cpw, x(1), p0)
            unique = .true.
            do j = 1, npts_
                p1 = points_(j, :)
                vi = p1 - p0
                di = norm(vi)
                if ((di .le. gtol) .and. (dabs(x(1) - params1d_(j)) .le. ptol)) then
                    unique = .false.
                    exit
                end if
            end do
            if (unique) then
                npts_ = npts_ + 1
                points_(npts_, :) = p0
                params1d_(npts_) = x(1)
                vi = pnt - p0
                di = norm(vi)
                dproj_(npts_) = di
            end if
        end if
    end do
    
    if (npts_ .gt. 0) then        
        ! Normal completion.
        flag = 1
        return
    end if
        
    if (ends .eqv. .false.) then
        ! No points and no endpoints.
        flag = 0
        call set_empty_results()
        return
    end if
    
    ! Select nearest endpoint if no solutions are found.
    call reset_results()
    allocate(points_(1, 3))
    allocate(params1d_(1))
    allocate(dproj_(1))
    
    points_(:, :) = 0.0d0
    params1d_(:) = 0.0d0
    dproj_(:) = 0.0d0
    
    call curve_point(n, p, uk, cpw, uk(p), p0)
    call curve_point(n, p, uk, cpw, uk(n + 1), p1)
    vi = pnt - p0
    d0 = norm(vi)
    vi = pnt - p1
    d1 = norm(vi)
    if (d0 .le. d1) then
        npts_ = 1
        points_(1, :) = p0
        params1d_(1) = uk(p)
        dproj_(1) = d0
    else
        npts_ = 1
        points_(1, :) = p1
        params1d_(1) = uk(n + 1)
        vi = pnt - p0
        di = norm(vi)
        dproj_(1) = d1
    end if
    flag = 1

    ! Internal subroutines for recursive subdivision.
    contains

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
            call project(n, cp, ab)
        else
            call split_bezier_curve(cpw, n, 0.50d0, ab(0), ab(1), &
                                    qw1, ab1(0), ab1(1), &
                                    qw2, ab2(0), ab2(1))
            
            candidate = candidate_curve(n, qw1)
            if (candidate) then
                call subdivide(n, qw1, ab1)
            end if

            candidate = candidate_curve(n, qw2)
            if (candidate) then
                call subdivide(n, qw2, ab2)
            end if
        end if
        
    end subroutine subdivide
    
    subroutine project(n, cp, ab)
        ! Project point to line.
        integer, intent(in) :: n
        double precision, intent(in) :: cp(0:n, 0:2), ab(0:1)
        
        double precision :: uil, uig, vdot, di
        double precision :: v0(0:2), v1(0:2), vi(0:2)
                            
        v0 = cp(n, :) - cp(0, :)
        v1 = pnt - cp(0, :)
        vdot = dot_product(v0, v0)
        if (vdot .le. 0.0d0) then
            return
        end if
        uil = dot_product(v1, v0) / vdot
        uig = local_to_global_param(ab(0), ab(1), uil)
        
        if (all_pnts) then
            ! Store result.
            allocate(pc_ptr%proj)
            pc_ptr%proj%u0 = uig
            call list_insert(pc_list, transfer(pc_ptr, list_data))
            nsol = nsol + 1
        else
            vi = pnt - (cp(0, :) + uil * v0)
            di = norm(vi)
            if (nsol .eq. 0) then
                ! Store result.
                allocate(pc_ptr%proj)
                pc_ptr%proj%u0 = uig
                call list_insert(pc_list, transfer(pc_ptr, list_data))
                nsol = 1
                dmin = di
            else
                if (di .lt. dmin) then
                    ! Retrieve second node and store result.
                    pc_ptr = transfer(list_get(list_next(pc_list)), pc_ptr)
                    pc_ptr%proj%u0 = uig
                    nsol = 1
                    dmin = di
                end if
            end if
        end if
        
    end subroutine project
    
    function candidate_curve(n, cpw) result(candidate)
        ! Test for potential projection to curve.
        integer, intent(in) :: n
        double precision, intent(in) :: cpw(0:n, 0:3)
        
        logical :: candidate
        
        double precision :: cp(0:n, 0:2), w(0:n)
        
        call dehomogenize_array1d(n, cpw, cp, w)
        candidate = point_nearest_bezier_curve(n, cp, pnt)
        
    end function candidate_curve
    
    subroutine obj(ni, x, fvec, fjac, ldfjac, iflag)
        ! Objective function for point refinement.
    
        ! Input/Output
        integer, intent(in) :: ni, ldfjac, iflag
        double precision, intent(in) :: x(1)
        double precision, intent(out) :: fvec(1), fjac(1, 1)
        
        ! Working
        double precision :: cders(0:2, 0:2)
        double precision :: pi(0:2), cu(0:2), cuu(0:2)
        
        ! Evaluate curve point and derivatives.
        call rat_curve_derivs(n, p, uk, cpw, x(1), 2, cders)
        pi = cders(0, :)
        cu = cders(1, :)
        cuu = cders(2, :)
        
        ! Calculate function at x and return in fvec.
        if (iflag .eq. 1) then
            fvec(1) = dot_product(cu, pi - pnt)
            return
        end if
        
        ! Calculate the jacobian at x and return in fjac.
        if (iflag .eq. 2) then
            fjac(1, 1) = dot_product(cuu, pi - pnt) + dot_product(cu, cu)
            return
        end if
    end subroutine obj

end subroutine project_point_to_curve

subroutine project_point_to_surface(pnt, n, p, uk, m, q, vk, cpw, ends, all_pnts, flag)
    !> Project a point to a surface.
    !> pnt - Point to project.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> ends - Option to project to nearest end curve if no
    !> orthogonal projections are found.
    !> all_pnts - Option to return only the nearest point during the
    !> subdivision process (False), or all points (True).
    !> flag - Status flag.
    
    !f2py intent(in) pnt, n, p, uk, m, q, vk, cpw, ends
    !f2py intent(in) all_pnts
    !f2py intent(out) flag
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    
    ! Input
    logical, intent(in) :: ends, all_pnts
    integer, intent(in) :: n, p, m, q
    double precision, intent(in) :: pnt(0:2)
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
                                    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    integer :: i, j
    integer :: nsol, nsub
    double precision :: dmin
    
    integer :: nbez_u
    double precision :: qw_u(0:n + p + 1, 0:p, 0:m, 0:3)
    double precision :: ab_u(0:n + p + 1, 0:1)
    
    integer :: nbez_v
    double precision :: qw_v(0:m + q + 1, 0:p, 0:q, 0:3)
    double precision :: ab_v(0:m + q + 1, 0:1)
    
    double precision :: uk_ui(0:2 * p - 1)
    double precision :: qw_ui(0:p, 0:m, 0:3)
    
    logical :: is_cand
    double precision :: qwi(0:p, 0:q, 0:3)
    double precision :: abi(0:1, 0:1)
    
    logical :: unique
    double precision :: di, dp
    double precision :: p0(0:2), p1(0:2), vi(0:2), duv(0:1)
    
    ! MINPACK
    integer :: info
    double precision :: ucheck, vcheck
    double precision :: x(2), fvec(2), fjac(2, 2), wa(15)
    
    ! Projection to boundary curves if needed.
    integer :: cpts
    double precision :: u, v
    double precision :: cpoints(0:2)
    double precision :: cparams(0:1)
    double precision :: cdist
    double precision, allocatable :: qw(:, :)
    
    ! Generic list for storing results.
    type(list_node_t), pointer :: ps_list => null()
    type(proj_ptr) :: ps_ptr
    
    ! Initialize
    nsol = 0
    dmin = 0.0d0
    call list_init(ps_list)
    
    ! Step 1: Decompose surface into u-direction Bezier strips.
    ! u-direction.
    call decompose_surface(n, p, uk, m, q, vk, cpw, 'u', n + p + 1, p, m, &
                           nbez_u, qw_u, ab_u)
    
    ! Step 2: Perform recursive subdivision on candidate Bezier patches.
    nsub = 0
    do i = 0, nbez_u - 1
        ! Decompose each strip in v-direction.
        qw_ui = qw_u(i, :, :, :)
        uk_ui(0:p) = ab_u(i, 0)
        uk_ui(p + 1:) = ab_u(i, 1)
        call decompose_surface(p, p, uk_ui, m, q, vk, qw_ui, 'v', m + q + 1, &
                               p, q, nbez_v, qw_v, ab_v)
                               
        ! Check each patch and subdivide if candidate.
        do j = 0, nbez_v - 1
            qwi = qw_v(j, :, :, :)
            abi(0, :) = ab_u(i, :)
            abi(1, :) = ab_v(j, :)
            is_cand = candidate_surface(p, q, qwi)
            if (is_cand) then
                call subdivide(p, q, qwi, abi)
            end if
        end do
    end do
    
    ! Return if no projection found and end curve option if false.
    if ((nsol .eq. 0) .and. (ends .eqv. .false.)) then
        flag = 0
        call set_empty_results()
        return
    end if
    
    ! Allocate result arrays.
    call reset_results()
    allocate(points_(nsol, 3))
    allocate(params2d_(nsol, 2))
    allocate(dproj_(nsol))
    points_(:, :) = 0.0d0
    params2d_(:, :) = 0.0d0
    dproj_(:) = 0.0d0
    
    ! Step 3: Refine projection points.
    do i = 1, nsol
        ! Get next
        ps_list = list_next(ps_list)
        ps_ptr = transfer(list_get(ps_list), ps_ptr)
        x(1) = ps_ptr%proj%u0
        x(2) = ps_ptr%proj%v0
        ! Call MINPACK hybrj driver.
        call hybrj1(obj, 2, x, fvec, fjac, 2, 1.490120d-8, info, wa, 15)
        ! if (info .ne. 1) then
        !     print *, "Unsuccessful solution in Fortran surface projection.", info
        ! end if
        ! Check parameter.
        ucheck = x(1)
        call check_param(n, p, uk, ucheck)
        x(1) = ucheck
        vcheck = x(2)
        call check_param(m, q, vk, vcheck)
        x(2) = vcheck
        if (info .eq. 1) then
            ! Get point and add to results only if unique.
            call surface_point(n, p, uk, m, q, vk, cpw, x(1), x(2), p0)
            unique = .true.
            do j = 1, npts_
                p1 = points_(j, :)
                vi = p1 - p0
                di = norm(vi)
                duv = x - params2d_(j, :)
                dp = norm(duv)
                if ((di .le. gtol) .and. (dp .le. ptol)) then
                    unique = .false.
                    exit
                end if
            end do
            if (unique) then
                npts_ = npts_ + 1
                points_(npts_, :) = p0
                params2d_(npts_, :) = x(:)
                vi = pnt - p0
                di = norm(vi)
                dproj_(npts_) = di
            end if
        end if
    end do
    
    ! Normal completion.
    if (npts_ > 0) then
        flag = 1
        return
    end if
    
    ! Step 4: Project to boundary curves if no projection if found and end
    ! curve option is true.
    cpts = 1
    cpoints(:) = 0.0d0
    cparams(:) = 0.0d0
    cdist = 0.0d0
    
    ! In u-direction at v=0.
    allocate(qw(0:n, 0:3))
    qw = cpw(:, 0, :)
    call project_point_to_curve(pnt, n, p, uk, qw, .true., .false., flag)
     
    v = vk(q)
    cpoints = points_(1, :)
    cparams = (/ params1d_(1), v /)
    cdist = dproj_(1)
    do i = 1, npts_
        if (dproj_(i) .lt. cdist) then
            cpoints = points_(i, :)
            cparams = (/ params1d_(i), v /)
            cdist = dproj_(i)
        end if
    end do
    
    ! In u-direction at v=1.
    qw = cpw(:, m, :)
    call project_point_to_curve(pnt, n, p, uk, qw, .true., .false., flag)
    v = vk(m + 1)
    do i = 1, npts_
        if (dproj_(i) .lt. cdist) then
            cpoints = points_(i, :)
            cparams = (/ params1d_(i), v /)
            cdist = dproj_(i)
        end if
    end do
    
    ! In v-direction at u=0.
    deallocate(qw)
    allocate(qw(0:m, 0:3))
    qw = cpw(0, :, :)
    call project_point_to_curve(pnt, m, q, vk, qw, .true., .false., flag)
    u = uk(p)
    do i = 1, npts_
        if (dproj_(i) .lt. cdist) then
            cpoints = points_(i, :)
            cparams = (/ u, params1d_(i) /)
            cdist = dproj_(i)
        end if
    end do
    
    ! In v-direction at u=1.
    qw = cpw(n, :, :)
    call project_point_to_curve(pnt, m, q, vk, qw, .true., .false., flag)
    u = uk(n + 1)
    do i = 1, npts_
        if (dproj_(i) .lt. cdist) then
            cpoints = points_(i, :)
            cparams = (/ u, params1d_(i) /)
            cdist = dproj_(i)
        end if
    end do
    
    ! Set results.
    call reset_results()
    allocate(points_(1, 0:2))
    allocate(params2d_(1, 0:1))
    allocate(dproj_(1))
    
    npts_ = cpts
    points_(1, :) = cpoints
    params2d_(1, :) = cparams
    dproj_(1) =cdist
    
    if (npts_ .gt. 0) then
        flag = 1
        return
    end if
    
    flag = 0
    call set_empty_results()
    
    ! Internal subroutines for recursive subdivision.
    contains
    
    subroutine obj(ni, x, fvec, fjac, ldfjac, iflag)
        ! Objective function for point refinement.
        ! Input/Output
        integer, intent(in) :: ni, ldfjac, iflag
        double precision, intent(in) :: x(2)
        double precision, intent(out) :: fvec(2), fjac(2, 2)
        
        ! Working
        double precision :: r(0:2)
        double precision :: sders(0:4, 0:4, 0:2)
        double precision :: pi(0:2), su(0:2), suu(0:2)
        double precision :: sv(0:2), svv(0:2), suv(0:2)
        
        ! Evaluate surface point and derivatives.
        call rat_surface_derivs(n, p, uk, m, q, vk, cpw, &
                                x(1), x(2), 4, sders)
        pi = sders(0, 0, :)
        su = sders(1, 0, :)
        suu = sders(2, 0, :)
        sv = sders(0, 1, :)
        svv = sders(0, 2, :)
        suv = sders(1, 1, :)
        
        r = pi - pnt
        ! Calculate function at x and return in fvec.
        if (iflag .eq. 1) then
            fvec(1) = dot_product(r, su)
            fvec(2) = dot_product(r, sv)
            return
        end if
        
        ! Calculate the jacobian at x and return in fjac.
        if (iflag .eq. 2) then
            fjac(1, 1) = dot_product(su, su) + dot_product(r, suu)
            fjac(1, 2) = dot_product(su, sv) + dot_product(r, suv)
            fjac(2, 1) = fjac(1, 2)
            fjac(2, 2) = dot_product(sv, sv) + dot_product(r, svv)
            return
        end if
    end subroutine obj
    
    recursive subroutine subdivide(n, m, cpw, ab)
        ! Recursive subdivsion
        integer, intent(in) :: n, m
        double precision, intent(in) :: cpw(0:n, 0:m, 0:3), ab(0:1, 0:1)
        
        logical :: flat, candidate
        double precision :: cp(0:n, 0:m, 0:2), w(0:n, 0:m), &
                            qw1(0:n, 0:m, 0:3), qw2(0:n, 0:m, 0:3), &
                            qw3(0:n, 0:m, 0:3), qw4(0:n, 0:m, 0:3), &
                            qw5(0:n, 0:m, 0:3), qw6(0:n, 0:m, 0:3), &
                            ab1(0:1, 0:1), ab2(0:1, 0:1), &
                            ab3(0:1, 0:1), ab4(0:1, 0:1), &
                            ab5(0:1, 0:1), ab6(0:1, 0:1)
        
        call dehomogenize_array2d(n, m, cpw, cp, w)
        
         call is_surface_flat(n, m, cp, ftol, flat)
        
        ! Force at least one subdivision since that seemed to help with
        ! Swept/skewed/tapered surfaces.
        if ((flat) .and. (nsub .gt. 0)) then
            call project(n, m, cp, ab)
        else
            nsub = nsub + 1
            call split_bezier_surface(cpw, n, m, 0.5d0, -1.0d0, &
                                      ab(0, 0), ab(0, 1), &
                                      ab(1, 0), ab(1, 1), &
                                      qw1, qw2, ab1, ab2)
            call split_bezier_surface(qw1, n, m, -1.0d0, 0.5d0, &
                                      ab1(0, 0), ab1(0, 1), &
                                      ab1(1, 0), ab1(1, 1), &
                                      qw3, qw4, ab3, ab4)
            call split_bezier_surface(qw2, n, m, -1.0d0, 0.5d0, &
                                      ab2(0, 0), ab2(0, 1), &
                                      ab2(1, 0), ab2(1, 1), &
                                      qw5, qw6, ab5, ab6)
            
            candidate = candidate_surface(n, m, qw3)
            if (candidate) then
                call subdivide(n, m, qw3, ab3)
            end if
            candidate = candidate_surface(n, m, qw4)
            if (candidate) then
                call subdivide(n, m, qw4, ab4)
            end if
            candidate = candidate_surface(n, m, qw5)
            if (candidate) then
                call subdivide(n, m, qw5, ab5)
            end if
            candidate = candidate_surface(n, m, qw6)
            if (candidate) then
                call subdivide(n, m, qw6, ab6)
            end if
        end if
        
    end subroutine subdivide
    
    subroutine project(n, m, cp, ab)
        ! Project point to a plane.
        integer, intent(in) :: n, m
        double precision, intent(in) :: cp(0:n, 0:m, 0:2), ab(0:1, 0:1)
        
        double precision :: ui, vi, vudot, vvdot, di, mag
        double precision :: vu(0:2), vv(0:2), vp(0:2), vn(0:2)
        double precision :: t3d(0:2, 0:2), t2d(0:2, 0:1)
                            
        vu = cp(n, 0, :) - cp(0, 0, :)
        vv = cp(0, m, :) - cp(0, 0, :)
        vudot = dot_product(vu, vu)
        vvdot = dot_product(vv, vv)
        if ((vudot .gt. 0.0d0) .and. (vvdot .gt. 0.0d0)) then        
            vp = pnt - cp(0, 0, :)
            t3d(0, :) = cp(0, 0, :)
            t3d(1, :) = cp(n, 0, :)
            t3d(2, :) = cp(0, m, :)
            t2d(0, :) = (/ ab(0, 0), ab(1, 0) /)
            t2d(1, :) = (/ ab(0, 1), ab(1, 0) /)
            t2d(2, :) = (/ ab(0, 0), ab(1, 1) /)
            call invert_point_in_triangle(pnt, t3d, t2d, ui, vi)
            if (all_pnts) then
                ! Store result.
                allocate(ps_ptr%proj)
                ps_ptr%proj%u0 = ui
                ps_ptr%proj%v0 = vi
                call list_insert(ps_list, transfer(ps_ptr, list_data))
                nsol = nsol + 1
            else
                vn = cross(vu, vv)
                mag = norm(vn)
                vn = vn / mag
                di = dabs(dot_product(vn, vp))
                if (nsol .eq. 0) then
                    ! Store result.
                    allocate(ps_ptr%proj)
                    ps_ptr%proj%u0 = ui
                    ps_ptr%proj%v0 = vi
                    call list_insert(ps_list, transfer(ps_ptr, list_data))
                    nsol = 1
                    dmin = di
                else
                    if (di .lt. dmin) then
                        ! Retrieve second node and store result.
                        ps_ptr = transfer(list_get(list_next(ps_list)), ps_ptr)
                        ps_ptr%proj%u0 = ui
                        ps_ptr%proj%v0 = vi
                        nsol = 1
                        dmin = di
                    end if
                end if
            end if
        end if
        
    end subroutine project
    
    function candidate_surface(n, m, cpw) result(candidate)
        ! Tets for potential projection to surface.
        integer, intent(in) :: n, m
        double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
        
        logical :: candidate
        
        double precision :: cp(0:n, 0:m, 0:2), w(0:n, 0:m)
        
        call dehomogenize_array2d(n, m, cpw, cp, w)
        candidate = point_nearest_bezier_patch(n, m, cp, pnt)
        
    end function candidate_surface

end subroutine project_point_to_surface

function point_nearest_bezier_curve(n, cp, p) result(is_near)
    !> Determine if the curve is a candidate for projection.
    !> n - Number of control points - 1.
    !> cp - Dehomogenized control points.
    !> p - Point to test.
    !> is_near - True if curve is candidate, False if not.
    
    !f2py intent(in) n, cp, p, atol
    !f2py intent(out) is_near
    !f2py depend(n) cp

    integer, intent(in) :: n
    double precision, intent(in) :: cp(0:n, 0:2), p(0:2)
    
    logical :: is_near
    
    integer :: i
    double precision :: r1, r2, d0, d1, p0(0:2)
    
    is_near = .true.
    r1 = dot_product(p - cp(0, :), cp(1, :) - cp(0, :))
    r2 = dot_product(cp(n, :) - p, cp(n, :) - cp(n - 1, :))
    if ((dabs(r1) .le. ptol) .or. (dabs(r2) .le. ptol)) then
        return
    end if
    
    d0 = dot_product(p - cp(0, :), p - cp(0, :))
    d1 = dot_product(p - cp(n, :), p - cp(n, :))
    if (d0 .le. d1) then
        p0 = cp(0, :)
        do i = 1, n
            if (dot_product(cp(i, :) - p0, p0 - p) .le. 0.0d0) then
                return
            end if
        end do
    else
        p0 = cp(n, :)
        do i = 0, n - 1
            if (dot_product(cp(i, :) - p0, p0 - p) .le. 0.0d0) then
                return
            end if
        end do
    end if
    is_near = .false.
    
end function point_nearest_bezier_curve

function point_nearest_bezier_patch(n, m, cp, p) result(is_near)
    !> Determine if the patch is a candidate for projection.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cp - Dehomogenized control points.
    !> p - Point to test.
    !> is_near - True is patch is candidate, False if not.
    
    !f2py intent(in) n, m, cp, p
    !f2py intent(out) is_near
    !f2py depend(n, m) cp

    integer, intent(in) :: n, m
    double precision, intent(in) :: cp(0:n, 0:m, 0:2), p(0:2)
    
    logical :: is_near

    integer :: i, j, i0, j0, indx
    double precision :: row(0:m, 0:2), col(0:n, 0:2), &
                        pc(0:2), dc, corners(4, 3), vc(0:2), &
                        dmin, di, p0(0:2), r, w(0:2), v0(0:2)
                        
    pc = cp(0, 0, :)
    vc = pc - p
    dc = norm(vc)
    corners(1, :) = (/ dc, 0.0d0, 0.0d0 /)
    
    pc = cp(n, 0, :)
    vc = pc - p
    dc = norm(vc)
    corners(2, :) = (/ dc, dble(n), 0.0d0 /)
    
    pc = cp(0, m, :)
    vc = pc - p
    dc = norm(vc)
    corners(3, :) = (/ dc, 0.0d0, dble(m) /)
    
    pc = cp(n, m, :)
    vc = pc - p
    dc = norm(vc)
    corners(4, :) = (/ dc, dble(n), dble(m) /)
    
    dmin = corners(1, 1)
    indx = 1
    do i = 2, 4
        di = corners(i, 1)
        if (di .lt. dmin) then
            dmin = di
            indx = i
        end if
    end do
    i0 = int(corners(indx, 2))
    j0 = int(corners(indx, 3))
    p0 = cp(i0, j0, :)
    
    is_near = .false.
    do i = 0, n
        do j = 0, m
            if ((i .eq. i0) .and. (j .eq. j0)) then
                ! continue
            else
                r = dot_product(cp(i, j, :) - p0, p0 - p)
                if (r .le. 0.0d0) then
                    is_near = .true.
                    exit
                end if
            end if
        end do
    end do
    if (is_near .eqv. .false.) then
        return
    end if
    
    ! Tangent cones.
    is_near = .false.
    do i = 0, n
        v0 = cp(i, 0, :) - p
        w = cp(i, 1, :) - cp(i, 0, :)
        r = dot_product(v0, w)
        if (r .le. 0.0d0) then
            is_near = .true.
            exit
        end if
    end do
    if (is_near .eqv. .false.) then
        return
    end if
    
    is_near = .false.
    do i = 0, n
        v0 = cp(i, m, :) - p
        w = cp(i, m - 1, :) - cp(i, m, :)
        r = dot_product(v0, w)
        if (r .le. 0.0d0) then
            is_near = .true.
            exit
        end if
    end do
    if (is_near .eqv. .false.) then
        return
    end if
    
    is_near = .false.
    do j = 0, m
        v0 = cp(0, j, :) - p
        w = cp(1, j, :) - cp(0, j, :)
        r = dot_product(v0, w)
        if (r .le. 0.0d0) then
            is_near = .true.
            exit
        end if
    end do
    if (is_near .eqv. .false.) then
        return
    end if
    
    is_near = .false.
    do j = 0, m
        v0 = cp(n, j, :) - p
        w = cp(n - 1, j, :) - cp(n, j, :)
        r = dot_product(v0, w)
        if (r .le. 0.0d0) then
            is_near = .true.
            exit
        end if
    end do
    if (is_near .eqv. .false.) then
        return
    end if
    
    is_near = .true.
    
end function point_nearest_bezier_patch

end module project