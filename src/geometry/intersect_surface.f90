module intersect_surface
use config, only: gtol, ptol, warnings
use math
use geom_utils
use evaluate
use modify
use divide
use calculate
use bounding_box
use intersect_bbox
use tessellate
use intersect_triangle
use invert, only: invert_points_on_plane
use nelder_mead, only: simplex
implicit none

private
public :: intersect_surface_plane, intersect_surface_surface
public :: refine_spi_point, refine_ssi_point

! Intersection data.
type si_data
    integer :: s1 = 0
    integer :: s2 = 0
end type si_data

type si_ptr
    type(si_data), pointer :: si
end type si_ptr

contains

subroutine intersect_surface_plane(n, p, uk, m, q, vk, cpw, p0, pnorm, vx, vy, ftol, flag)
    !> Find the intersection curve(s) between a surface and a plane.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> p0 - Origin of plane.
    !> pnorm - Unit vector normal to plane.
    !> vx - Vector defining x-axis of plane.
    !> vy - Vector defining y-axis of plane.
    !> ftol - Surface flatness tolerance.
    !> flag - Status flag.

    !f2py intent(in) n, p, uk, m, q, vk, cpw, ab, p0, pnorm, vx, vy, ftol
    !f2py intent(out) flag
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    
    ! Input
    integer, intent(in) :: n, m, p, q
    double precision, intent(in) :: ftol
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    double precision, intent(in) :: p0(0:2), pnorm(0:2), vx(0:2), vy(0:2)
    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    logical :: candidate, uclosed, vclosed
    integer :: i, j, indx, nmax, ndim
    integer :: etotal, umult, vmult
    integer, allocatable :: edges(:, :)
    integer, allocatable :: all_crv_ids(:, :)
    integer, allocatable :: no_filter(:)
    double precision :: ui, vi, u, v, up, vp, tol
    double precision :: pi(0:2)
    
    integer :: ncells, ncand, nsub
    integer, allocatable :: children(:, :)
    integer, allocatable :: parent(:)
    integer, allocatable :: acells(:, :)
    integer, allocatable :: position(:)
    integer, allocatable :: candidates(:)
    double precision :: tri(0:100, 0:2, 0:1)
    double precision, allocatable :: cell_params(:, :, :)
    double precision, allocatable :: points2d_s1(:, :), points2d_s2(:, :)
    
    ! Generic list for storing surface data.
    type(list_node_t), pointer :: cell_list => null()
    type(srf_ptr) :: cell_ptr
    
    ! Initialize.
    ncand = 0
    ncells = 0
    nsub = 0
    etotal = 0
    call list_init(cell_list)
    
    ! Step 1: See below for internal subroutines for recursive
    ! subdivision.
    
    ! Step 2: Use recursive subdivision to find potential intersection
    ! points.
    allocate(cell_ptr%ci)
    call candidate_intersect(n, m, cpw, candidate)
    if (candidate) then
        call subdivide(n, p, uk, m, q, vk, cpw, cell_ptr)
    end if
    
    ! Check if no intersections found.
    if (ncand .eq. 0) then
        flag = 0
        call list_free(cell_list)
        call set_empty_results()
        return
    end if
    
    ! Build arrays to tessellate the cells.
    call build_tess_data(ncells, cell_list, children, parent, acells, &
                         position, candidates, cell_params)

    ! Allocate based on number of potential intersection cells.
    call reset_results()
    nmax = 4 * ncells
    allocate(verts_(0:nmax - 1, 0:2))
    allocate(crv_size_(0:nmax - 1))
    allocate(crv_ids_(0:nmax - 1, 0:nmax - 1))
    allocate(all_crv_ids(0:nmax - 1, 0:nmax - 1))
    allocate(edges(0:nmax - 1, 0:1))
    allocate(ssi_params1_(0:nmax - 1, 0:1))
    allocate(ssi_params2_(0:nmax - 1, 0:1))
    allocate(points2d_s1(0:nmax - 1, 0:2))
    allocate(points2d_s2(0:nmax - 1, 0:2))
    allocate(no_filter(0:nmax - 1))
    
    etotal = 0
    edges(:, :) = 0
    ssi_params1_(:, :) = 0.0d0
    ssi_params2_(:, :) = 0.0d0
    verts_(:, :) = 0.0d0
    crv_size_(:) = 0
    crv_ids_(:, :) = 0
    all_crv_ids(:, :) = 0
    points2d_s1(:, :) = 0.0d0
    points2d_s2(:, :) = 0.0d0
    no_filter(:) = 0
    
    ! Step 3: Tessellate potential surfaces and intersect.
    tri(:, :, :) = 0.0d0
    ndim = size(children, 1)
    do i = 0, ncells
        if (candidates(i) .eq. 1) then
            call intersect(i)
        end if
    end do   
    
    ! Check for no intersection.    
    if (nverts_.le. 1) then
        flag = -1
        call list_free(cell_list)
        call set_empty_results()
        return
    end if
    
    ! Step 4: Find topology of the intersection curves.
    call trace_curves2(nverts_, etotal, points2d_s1, points2d_s2, edges, ptol, &
                       nmax, ncrvs_, crv_size_, all_crv_ids)
                      
    ! Step 5: Refine the intersection points by minimizing the distance
    ! between the initial intersection points and the true surface-plane
    ! intersection.
    tol = gtol / 100.0d0
    do i = 0, ncrvs_ - 1
        do j = 0, crv_size_(i) - 1
            indx = all_crv_ids(i, j)
            ui = ssi_params1_(indx, 0)
            vi = ssi_params1_(indx, 1)
            call refine_spi_point(n, p, uk, m, q, vk, cpw, p0, pnorm, vx, vy, &
                                  ui, vi, tol, u, v, up, vp, pi)
            ssi_params1_(indx, :) = (/ u, v /)
            ssi_params2_(indx, :) = (/ up, vp /)
            verts_(indx, :) = pi
        end do
    end do

    ! Build a list that specifies if a point should not be filtered.
    ! do i = 0, ncrvs_ - 1
    !     do j = 0, crv_size_(i) - 1
    !         indx = all_crv_ids(i, j)
    !         u = ssi_params1_(indx, 0)
    !         v = ssi_params1_(indx, 1)
    !         umult = find_mult(n, p, uk, u)
    !         vmult = find_mult(m, q, vk, v)
    !         if ((umult .ge. p) .or. (vmult .ge. q)) then
    !             no_filter(indx) = 1
    !         end if
    !     end do
    ! end do
    
    ! Step 6: Filter out points.
    ! call filter_points(ncrvs_, nmax, crv_size_, all_crv_ids, no_filter, verts_, 0.0d0, gtol, crv_ids_)
    call filter_points(ncrvs_, nmax, crv_size_, all_crv_ids, verts_, gtol, crv_ids_)
    
    ! Return results.
    call list_free(cell_list)
    flag = 1
    
    contains
    
    recursive subroutine subdivide(n, p, uk, m, q, vk, cpw, cell_ptr)
        ! Recursive subdivision.
        ! Input
        integer, intent(in) :: n, p, m, q
        double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
        double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
        
        type(srf_ptr), intent(inout) :: cell_ptr
        
        ! Working
        integer :: mid
        logical :: flat, t1, t2, t3, t4
        double precision :: cp(0:n, 0:m, 0:2), w(0:n, 0:m)
        
        ! Working arrays for surface splitting
        integer :: n1i, m1i, n2i, m2i, ukn, vkn
        double precision :: uv, ui, vi
        double precision, allocatable :: uk1i(:), uk2i(:)
        double precision, allocatable :: vk1i(:), vk2i(:)
        double precision, allocatable :: qw1i(:, :, :)
        double precision, allocatable :: qw2i(:, :, :)
        
        ! Working arrays for storing split surface data
        integer :: n1, n2, n3, n4, m1, m2, m3, m4
        double precision, allocatable :: uk1(:), uk2(:), uk3(:), uk4(:)
        double precision, allocatable :: vk1(:), vk2(:), vk3(:), vk4(:)
        double precision, allocatable :: qw1(:, :, :), qw2(:, : ,:)
        double precision, allocatable :: qw3(:, :, :), qw4(:, : ,:)
        
        ! New cells
        type(srf_ptr) :: ne_ptr, se_ptr, sw_ptr, nw_ptr
        
        nsub = nsub + 1
        
        ! Insert cell
        call list_insert(cell_list, transfer(cell_ptr, list_data))
        
        ! Store surface parameters
        cell_ptr%ci%u0 = uk(p)
        cell_ptr%ci%u1 = uk(n + 1)
        cell_ptr%ci%v0 = vk(q)
        cell_ptr%ci%v1 = vk(m + 1)

        ! Check flatness
        call dehomogenize_array2d(n, m, cpw, cp, w)
        call is_flat(n, m, cp, flat)
        if ((flat) .and. (nsub .gt. 1)) then
            ! Save candidate surfaces
            cell_ptr%ci%is_cand = .true.
            ncand = ncand + 1
        else
            ! Split surface into four regions at midpoint
            ui = 0.50d0 * (cell_ptr%ci%u0 + cell_ptr%ci%u1)
            vi = 0.50d0 * (cell_ptr%ci%v0 + cell_ptr%ci%v1)
                
            ! If linear, use middle interior knot
            ukn = size(uk, 1)
            vkn = size(vk, 1)
            if ((p .eq. 1) .and. (ukn .gt. 4)) then
                mid = ceiling((dble(ukn) - 1.0d0) / 2.0d0)
                ui = uk(mid)
            end if
            if ((q .eq. 1) .and. (vkn .gt. 4)) then
                mid = ceiling((dble(vkn) - 1.0d0) / 2.0d0)
                vi = vk(mid)
            end if
            ! Split surface at u-paramter
            call split_surface(n, p, uk, m, q, vk, cpw, ui, 'u', &
                               uk1i, vk1i, qw1i, uk2i, vk2i, qw2i)
            n1i = size(qw1i, 1) - 1
            m1i = size(qw1i, 2) - 1
            n2i = size(qw2i, 1) - 1
            m2i = size(qw2i, 2) - 1
            ! Split halves in v-direction and store
            call split_surface(n1i, p, uk1i, m1i, q, vk1i, qw1i, &
                               vi, 'v', uk1, vk1, qw1, &
                               uk2, vk2, qw2)
            n1 = size(qw1, 1) - 1
            m1 = size(qw1, 2) - 1
            n2 = size(qw2, 1) - 1
            m2 = size(qw2, 2) - 1
            call split_surface(n2i, p, uk2i, m2i, q, vk2i, qw2i, &
                               vi, 'v', uk3, vk3, qw3, &
                               uk4, vk4, qw4)
            n3 = size(qw3, 1) - 1
            m3 = size(qw3, 2) - 1
            n4 = size(qw4, 1) - 1
            m4 = size(qw4, 2) - 1
            
            ! Test for potential intersection
            call candidate_intersect(n1, m1, qw1, t1)
            call candidate_intersect(n2, m2, qw2, t2)
            call candidate_intersect(n3, m3, qw3, t3)
            call candidate_intersect(n4, m4, qw4, t4)
            
            ! Allocate new cells
            allocate(ne_ptr%ci, se_ptr%ci, sw_ptr%ci, nw_ptr%ci)
            
            if (t1 .or. t2 .or. t3 .or. t4) then
                cell_ptr%ci%has_child = .true.
            end if
            
            if (t1) then
                ! Cell 1 (SW)
                sw_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%sw => sw_ptr%ci
                sw_ptr%ci%parent => cell_ptr%ci
                sw_ptr%ci%position = 1
            end if
            
            if (t2) then
                ! Cell 2 (NW)
                nw_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%nw => nw_ptr%ci
                nw_ptr%ci%parent => cell_ptr%ci
                nw_ptr%ci%position = 2
            end if
            
            if (t3) then
                ! Cell 3 (SE)
                se_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%se => se_ptr%ci
                se_ptr%ci%parent => cell_ptr%ci
                se_ptr%ci%position = 3
            end if
            
            if (t4) then
                ! Cell 4 (NE)
                ne_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%ne => ne_ptr%ci
                ne_ptr%ci%parent => cell_ptr%ci
                ne_ptr%ci%position = 4
            end if
            
            if (t1) then
                ! Adjacent cells for cell 1 (SW)
                sw_ptr%ci%n => nw_ptr%ci
                sw_ptr%ci%e => se_ptr%ci
                sw_ptr%ci%s => cell_ptr%ci%s
                sw_ptr%ci%w => cell_ptr%ci%w
            end if
            
            if (t2) then
                ! Adjacent cells for cell 2 (NW)
                nw_ptr%ci%n => cell_ptr%ci%n
                nw_ptr%ci%e => ne_ptr%ci
                nw_ptr%ci%s => sw_ptr%ci
                nw_ptr%ci%w => cell_ptr%ci%w
            end if
            
            if (t3) then
                ! Adjacent cells for cell 3 (SE)
                se_ptr%ci%n => ne_ptr%ci
                se_ptr%ci%e => cell_ptr%ci%e
                se_ptr%ci%s => cell_ptr%ci%s
                se_ptr%ci%w => sw_ptr%ci
            end if
            
            if (t4) then
                ! Adjacent cells for cell 4 (NE)
                ne_ptr%ci%n => cell_ptr%ci%n
                ne_ptr%ci%e => cell_ptr%ci%e
                ne_ptr%ci%s => se_ptr%ci
                ne_ptr%ci%w => nw_ptr%ci
            end if
            
            ! Subdivide
            if (t1) then
                call subdivide(n1, p, uk1, m1, q, vk1, qw1, sw_ptr)
            end if
            if (t2) then
                call subdivide(n2, p, uk2, m2, q, vk2, qw2, nw_ptr)
            end if
            if (t3) then
                call subdivide(n3, p, uk3, m3, q, vk3, qw3, se_ptr)
            end if
            if (t4) then
                call subdivide(n4, p, uk4, m4, q, vk4, qw4, ne_ptr)
            end if
        end if
        
    end subroutine subdivide
    
    subroutine candidate_intersect(n, m, cpw, candidate)
        ! Check for possible intersection using bounding box.
        ! Input
        integer, intent(in) :: n, m
        double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
        
        ! Output
        logical, intent(out) :: candidate
        
        ! Working
        double precision :: bbox(0:2, 0:1)
        
        bbox = surface_bbox(n, m, cpw)
        candidate = bbox_intersects_plane(bbox, p0, pnorm, gtol)
        
    end subroutine candidate_intersect
    
    subroutine is_flat(n, m, cp, flat)
        ! Check surface flatness.
        integer, intent(in) :: n, m
        double precision, intent(in) :: cp(0:n, 0:m, 0:2)
        
        logical, intent(out) :: flat
        
        call is_surface_flat(n, m, cp, ftol, flat)
        
    end subroutine is_flat
    
    subroutine intersect(icell)
        ! Intersect cell using triangle-plane intersection.
        ! Input
        integer, intent(in) :: icell
        ! Working
        integer :: i, ni, ntri
        double precision :: ub, vb, area
        double precision :: uv0(0:1), uv1(0:1), uv2(0:1)
        double precision :: t0(0:2), t1(0:2), t2(0:2)
        double precision :: temp_tri(0:2, 0:2)
        double precision :: pi(0:1, 0:2)
        double precision :: ptemp(0:2)
        double precision :: plane_params(0:1, 0:1)
                
        ! Tessellate cell.
        call tessellate_cell(icell, ndim, children, acells, position, &
                             parent, cell_params, ntri, tri)
        ! Intersect each triangle with the plane.
        do i = 0, ntri - 1
            uv0 = tri(i, 0, :)
            uv1 = tri(i, 1, :)
            uv2 = tri(i, 2, :)
            call surface_point(n, p, uk, m, q, vk, cpw, &
                               uv0(0), uv0(1), t0)
            call surface_point(n, p, uk, m, q, vk, cpw, &
                               uv1(0), uv1(1), t1)
            call surface_point(n, p, uk, m, q, vk, cpw, &
                               uv2(0), uv2(1), t2)
            temp_tri(0, :) = t0
            temp_tri(1, :) = t1
            temp_tri(2, :) = t2
            ! Check for degenerate triangle.
            area = triangle_area(3, temp_tri)
            if (area .gt. 1.0d-12) then
                call intersect_triangle_plane(temp_tri, p0, pnorm, gtol, &
                                              ptol, ni, pi)
                if (ni .eq. 2) then
                    ! Add points and edges.
                    edges(etotal, :) = (/ nverts_, nverts_ + 1 /)
                    verts_(nverts_, :) = pi(0, :)
                    verts_(nverts_ + 1, :) = pi(1, :)
                    
                    ! Get barycentric coordinates and convert to surface
                    ! parameters.
                    ptemp = pi(0, :)
                    call barycentric_params(ptemp, temp_tri, ub, vb)
                    ssi_params1_(nverts_, :) = (1.0d0 - ub - vb) * uv0 + ub * uv1 &
                                                + vb * uv2
                    points2d_s1(nverts_, 0) = ssi_params1_(nverts_, 0)
                    points2d_s1(nverts_, 1) = ssi_params1_(nverts_, 1)
                    ptemp = pi(1, :)
                    call barycentric_params(ptemp, temp_tri, ub, vb)
                    ssi_params1_(nverts_ + 1, :) = (1.0d0 - ub - vb) * uv0 + &
                                                    ub * uv1 + vb * uv2
                    points2d_s1(nverts_ + 1, 0) = ssi_params1_(nverts_ + 1, 0)
                    points2d_s1(nverts_ + 1, 1) = ssi_params1_(nverts_ + 1, 1)
                    
                    ! Plane parameters.
                    call invert_points_on_plane(2, pi, p0, vx, vy, plane_params)
                    ssi_params2_(nverts_, :) = plane_params(0, :)
                    ssi_params2_(nverts_ + 1, :) = plane_params(1, :)
                    points2d_s2(nverts_, 0:1) = ssi_params2_(nverts_, :)
                    points2d_s2(nverts_ + 1, 0:1) = ssi_params2_(nverts_ + 1, :)
                    
                    ! Update totals.
                    etotal = etotal + 1
                    nverts_ = nverts_ + 2
                end if
            end if
        end do
        
    end subroutine intersect
    
end subroutine intersect_surface_plane

subroutine intersect_surface_surface(n1, p1, uk1, m1, q1, vk1, cpw1, &
                                      n2, p2, uk2, m2, q2, vk2, cpw2, ftol, flag)
    !> Find the intersection curve(s) between two surfaces.
    !> n1 - Number of control points - 1 in u-direction for surface 1.
    !> p1 - Degree in u-direction for surface 1.
    !> uk1 - Knot vector in u-direction for surface 1.
    !> m1 - Number of control points - 1 in v-direction for surface 1.
    !> q1 - Degree in v-direction for surface 1.
    !> vk1 - Knot vector in v-direction for surface 1.
    !> cpw1 - Control points for surface 1.
    !> n2 - Number of control points - 1 in u-direction for surface 2.
    !> p2 - Degree in u-direction for surface 2.
    !> uk2 - Knot vector in u-direction for surface 2.
    !> m2 - Number of control points - 1 in v-direction for surface 2.
    !> q2 - Degree in v-direction for surface 2.
    !> vk2 - Knot vector in v-direction for surface 2.
    !> cpw2 - Control points for surface 2.
    !> ftol - Surface flatness tolerance.
    !> flag - Status flag.
    
    !f2py intent(in) n1, p1, uk1, m1, q1, vk1, cpw1, n2, p2, uk2
    !f2py intent(in) m2, q2, vk2, cpw2, ftol
    !f2py intent(out) flag
    !f2py depend(n1, p1) uk1
    !f2py depend(m1, q1) vk1
    !f2py depend(n1, m1) cpw1
    !f2py depend(n2, p2) uk2
    !f2py depend(m2, q2) vk2
    !f2py depend(n2, m2) cpw2
    
    ! Input
    integer, intent(in) :: n1, p1, m1, q1, n2, p2, m2, q2
    double precision, intent(in) :: ftol
    double precision, intent(in) :: uk1(0:n1 + p1 + 1)
    double precision, intent(in) :: vk1(0:m1 + q1 + 1)
    double precision, intent(in) :: uk2(0:n2 + p2 + 1)
    double precision, intent(in) :: vk2(0:m2 + q2 + 1)
    double precision, intent(in) :: cpw1(0:n1, 0:m1, 0:3)
    double precision, intent(in) :: cpw2(0:n2, 0:m2, 0:3)
    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    logical :: candidate, uclosed1, vclosed1, uclosed2, vclosed2
    integer :: i, j, k, indx, nmax, nsub, ndim
    integer :: ncand, ncells1, ncells2, cell1, cell2, etotal
    
    integer, allocatable :: candidates(:, :)
    integer, allocatable :: candidates1(:), candidates2(:)
    integer, allocatable :: used_cells1(:), used_cells2(:)
    integer, allocatable :: children1(:, :), children2(:, :)
    integer, allocatable :: parent1(:), parent2(:)
    integer, allocatable :: acells1(:, :), acells2(:, :)
    integer, allocatable :: position1(:), position2(:)
    integer, allocatable :: edges(:, :), all_crv_ids(:, :)
    integer, allocatable :: no_filter(:)
    double precision :: u1, v1, u2, v2, tol
    double precision :: u01, v01, u02, v02
    double precision :: pi(0:2)
    double precision :: tri1(0:100, 0:2, 0:1), tri2(0:100, 0:2, 0:1)
    double precision, allocatable :: cell_params1(:, :, :)
    double precision, allocatable :: cell_params2(:, :, :)
    double precision, allocatable :: points2d_s1(:, :), points2d_s2(:, :)
    
    ! Generic list for storing surface data.
    type(list_node_t), pointer :: cell_list1 => null()
    type(list_node_t), pointer :: cell_list2 => null()
    type(list_node_t), pointer :: ssi_list => null()
    type(srf_ptr) :: cell_ptr1, cell_ptr2
    type(si_ptr) :: ssi_ptr
    
    ! Initialize
    ncand = 0
    ncells1 = 0
    ncells2 = 0
    nsub = 0
    call list_init(cell_list1)
    call list_init(cell_list2)
    call list_init(ssi_list)
    
    ! Step 1: See below for internal subroutines for recursive
    ! subdivision.
    
    ! Step 2: Use recursive subdivision to find candidate surfaces.
    allocate(cell_ptr1%ci, cell_ptr2%ci)
    call candidate_intersect(n1, m1, cpw1, n2, m2, cpw2, candidate)
    if (candidate) then
        call subdivide(n1, p1, uk1, m1, q1, vk1, cpw1, &
                       n2, p2, uk2, m2, q2, vk2, cpw2, &
                       cell_ptr1, cell_ptr2)
    end if

    ! Check if no intersections are found.
    if (ncand .eq. 0) then
        flag = 0
        call list_free(cell_list1)
        call list_free(cell_list2)
        call list_free(ssi_list)
        call set_empty_results()
        return
    end if
    
    ! Build arrays to tessellate data.
    call build_tess_data(nsub - 1, cell_list1, children1, parent1, acells1, &
                         position1, candidates1, cell_params1)
    call build_tess_data(nsub - 1, cell_list2, children2, parent2, acells2, &
                         position2, candidates2, cell_params2)
    
    ! Build candidate arrays for intersecting cells.
    allocate(candidates(0:nsub - 1, 0:nsub - 1))
    candidates(:, :) = 0
    do i = 1, ncand
        ! Get next cell.
        ssi_list = list_next(ssi_list)
        ssi_ptr = transfer(list_get(ssi_list), ssi_ptr)
        candidates(ssi_ptr%si%s1, ssi_ptr%si%s2) = 1
    end do

    ! It's possible that a parent surface and its children may be in the
    ! potential intersection list. This can sometimes cause issues in the
    ! tessellation algorithm leading to non-congruent edges. For this reason,
    ! if the parent of a surface is flat and also in the potential intersection
    ! list, replace the child surface with its parent.
    do i = 0, ncells1
        do j = 0, ncells2 
            if (candidates(i, j) .eq. 1) then
                ! Erase intersection.
                candidates(i, j) = 0
                cell1 = i
                cell2 = j
                ! Surface 1
                do while (candidates1(parent1(cell1)) .eq. 1)
                    cell1 = parent1(cell1)
                    children1(cell1, :) = 0
                    if (cell1 .eq. 0) then
                        exit
                    end if
                end do
                ! Surface 2
                do while (candidates2(parent2(cell2)) .eq. 1)
                    cell2 = parent2(cell2)
                    children2(cell2, :) = 0
                    if (cell2 .eq. 0) then
                        exit
                    end if
                end do
                ! Reset intersection.
                candidates(cell1, cell2) = 1
            end if
        end do
    end do
    
    ! Allocate based on number of potential intersection cells.
    call reset_results()
    nmax = 4 * max(ncells1, ncells2)
    allocate(edges(0:nmax - 1, 0:1))
    allocate(verts_(0:nmax - 1, 0:2))
    allocate(crv_size_(0:nmax - 1))
    allocate(crv_ids_(0:nmax - 1, 0:nmax - 1))
    allocate(all_crv_ids(0:nmax - 1, 0:nmax - 1))
    allocate(ssi_params1_(0:nmax - 1, 0:1))
    allocate(ssi_params2_(0:nmax - 1, 0:1))
    allocate(points2d_s1(0:nmax - 1, 0:2))
    allocate(points2d_s2(0:nmax - 1, 0:2))
    allocate(no_filter(0:nmax - 1))
    
    etotal = 0
    edges(:, :) = 0
    verts_(:, :) = 0.0d0
    crv_size_(:) = 0
    crv_ids_(:, :) = 0
    all_crv_ids(:, :) = 0
    ssi_params1_(:, :) = 0.0d0
    ssi_params2_(:, :) = 0.0d0
    points2d_s1(:, :) = 0.0d0
    points2d_s2(:, :) = 0.0d0
    no_filter(:) = 0
    
    ! Step 3: Tessellate potential surfaces and intersect.
    ndim = size(children1, 1)
    tri1(:, :, :) = 0.0d0
    tri2(:, :, :) = 0.0d0
    do i = 0, ncells1
        do j = 0, ncells2
            if (candidates(i, j) .eq. 1) then
                call intersect(i, j)
            end if
        end do
    end do
    
    ! Check for no intersection.
    if (nverts_ .le. 1) then
        flag = -1
        call list_free(cell_list1)
        call list_free(cell_list2)
        call list_free(ssi_list)
        call set_empty_results()
        return
    end if
    
    !Step 4: Find topology of the curves.
    call trace_curves2(nverts_, etotal, points2d_s1, points2d_s2, edges, ptol, &
                       nmax, ncrvs_, crv_size_, all_crv_ids)
    
    ! Step 5: Refine intersection points.
    tol = gtol / 100.0d0
    do i = 0, ncrvs_ - 1
        do j = 0, crv_size_(i) - 1
            indx = all_crv_ids(i, j)
            u01 = ssi_params1_(indx, 0)
            v01 = ssi_params1_(indx, 1)
            u02 = ssi_params2_(indx, 0)
            v02 = ssi_params2_(indx, 1)
            call refine_ssi_point(n1, p1, uk1, m1, q1, vk1, cpw1, &
                                  n2, p2, uk2, m2, q2, vk2, cpw2, &
                                  u01, v01, u02, v02, tol, &
                                  u1, v1, u2, v2, pi)
            ssi_params1_(indx, :) = (/ u1, v1 /)
            ssi_params2_(indx, :) = (/ u2, v2 /)
            verts_(indx, :) = pi
        end do
    end do
    
    ! Step 6: Filter out points.
    ! call filter_points(ncrvs_, nmax, crv_size_, all_crv_ids, no_filter, verts_, 0.0d0, gtol, crv_ids_)
    call filter_points(ncrvs_, nmax, crv_size_, all_crv_ids, verts_, gtol, crv_ids_)
    
    ! Return results.
    call list_free(cell_list1)
    call list_free(cell_list2)
    call list_free(ssi_list)
    flag = 1

    ! Internal subroutines for recursive subdivision.
    contains
    
    recursive subroutine subdivide(n1, p1, uk1, m1, q1, vk1, cpw1, &
                                   n2, p2, uk2, m2, q2, vk2, cpw2, &
                                   cell_ptr1, cell_ptr2)
        ! Recursive subdivision.
        ! Input
        integer, intent(in) :: n1, p1, m1, q1, n2, p2, m2, q2
        double precision, intent(in) :: uk1(0:n1 + p1 + 1)
        double precision, intent(in) :: vk1(0:m1 + q1 + 1)
        double precision, intent(in) :: uk2(0:n2 + p2 + 1)
        double precision, intent(in) :: vk2(0:m2 + q2 + 1)
        double precision, intent(in) :: cpw1(0:n1, 0:m1, 0:3)
        double precision, intent(in) :: cpw2(0:n2, 0:m2, 0:3)
        
        type(srf_ptr), intent(inout) :: cell_ptr1, cell_ptr2
        
        ! Working
        logical :: is_flat1, is_flat2
        integer :: mid
        double precision :: cp1(0:n1, 0:m1, 0:2), w1(0:n1, 0:m1)
        double precision :: cp2(0:n2, 0:m2, 0:2), w2(0:n2, 0:m2)
        
        logical :: t11_t21, t11_t22, t11_t23, t11_t24
        logical :: t12_t21, t12_t22, t12_t23, t12_t24
        logical :: t13_t21, t13_t22, t13_t23, t13_t24
        logical :: t14_t21, t14_t22, t14_t23, t14_t24
        integer :: n11, n12, n13, n14, m11, m12, m13, m14
        integer :: n21, n22, n23, n24, m21, m22, m23, m24
        integer :: csn11, csn12, csn13, csn14, csn21, csn22, csn23, csn24
        double precision, allocatable :: uk11(:), uk12(:)
        double precision, allocatable :: uk13(:), uk14(:)
        double precision, allocatable :: vk11(:), vk12(:)
        double precision, allocatable :: vk13(:), vk14(:)
        double precision, allocatable :: uk21(:), uk22(:)
        double precision, allocatable :: uk23(:), uk24(:)
        double precision, allocatable :: vk21(:), vk22(:)
        double precision, allocatable :: vk23(:), vk24(:)
        double precision, allocatable :: qw11(:, :, :), qw12(:, : ,:)
        double precision, allocatable :: qw13(:, :, :), qw14(:, : ,:)
        double precision, allocatable :: qw21(:, :, :), qw22(:, : ,:)
        double precision, allocatable :: qw23(:, :, :), qw24(:, : ,:)

        integer :: n1i, m1i, n2i, m2i, ukn, vkn
        double precision :: ui, vi
        double precision, allocatable :: uk1i(:), uk2i(:)
        double precision, allocatable :: vk1i(:), vk2i(:)
        double precision, allocatable :: qw1i(:, :, :)
        double precision, allocatable :: qw2i(:, :, :)
        
        ! New cells
        type(srf_ptr) :: ne_ptr1, se_ptr1, sw_ptr1, nw_ptr1
        type(srf_ptr) :: ne_ptr2, se_ptr2, sw_ptr2, nw_ptr2
        
        nsub = nsub + 1

        ! Insert cells
        call list_insert(cell_list1, transfer(cell_ptr1, list_data))        
        call list_insert(cell_list2, transfer(cell_ptr2, list_data))
        
        ! Store surface parameters.
        ! Surface 1
        cell_ptr1%ci%u0 = uk1(p1)
        cell_ptr1%ci%u1 = uk1(n1 + 1)
        cell_ptr1%ci%v0 = vk1(q1)
        cell_ptr1%ci%v1 = vk1(m1 + 1)
        
        ! Surface 2
        cell_ptr2%ci%u0 = uk2(p2)
        cell_ptr2%ci%u1 = uk2(n2 + 1)
        cell_ptr2%ci%v0 = vk2(q2)
        cell_ptr2%ci%v1 = vk2(m2 + 1)
        
        ! Check flatness.
        call dehomogenize_array2d(n1, m1, cpw1, cp1, w1)
        call dehomogenize_array2d(n2, m2, cpw2, cp2, w2)
        call is_flat(n1, m1, cp1, is_flat1)
        call is_flat(n2, m2, cp2, is_flat2)
        if (is_flat1 .and. is_flat2 .and. (nsub .gt. 1)) then
            ! Save candidates.
            ncand = ncand + 1
            cell_ptr1%ci%is_cand = .true.
            cell_ptr2%ci%is_cand = .true.
            allocate(ssi_ptr%si)
            ssi_ptr%si%s1 = cell_ptr1%ci%sid
            ssi_ptr%si%s2 = cell_ptr2%ci%sid
            call list_insert(ssi_list, transfer(ssi_ptr, list_data))
        else
            if ((is_flat1 .eqv. .false.) .and. is_flat2) then
                ! Subdivide suface 1.
                ui = 0.50d0 * (cell_ptr1%ci%u0 + cell_ptr1%ci%u1)
                vi = 0.50d0 * (cell_ptr1%ci%v0 + cell_ptr1%ci%v1)
                
                ! If linear, split at middle knot.
                ukn = size(uk1, 1)
                vkn = size(vk1, 1)
                if ((p1 .eq. 1) .and. (ukn .gt. 4)) then
                    mid = ceiling((dble(ukn) - 1.0d0) / 2.0d0)
                    ui = uk1(mid)
                end if
                if ((q1 .eq. 1) .and. (vkn .gt. 4)) then
                    mid = ceiling((dble(vkn) - 1.0d0) / 2.0d0)
                    vi = vk1(mid)
                end if
                ! Split surface at u-paramter.
                call split_surface(n1, p1, uk1, m1, q1, vk1, cpw1, ui, 'u', &
                                   uk1i, vk1i, qw1i, uk2i, vk2i, qw2i)
                n1i = size(qw1i, 1) - 1
                m1i = size(qw1i, 2) - 1
                n2i = size(qw2i, 1) - 1
                m2i = size(qw2i, 2) - 1
                ! Split halves in v-direction and store.
                call split_surface(n1i, p1, uk1i, m1i, q1, vk1i, qw1i, vi, 'v', &
                                   uk11, vk11, qw11, uk12, vk12, qw12)
                n11 = size(qw11, 1) - 1
                m11 = size(qw11, 2) - 1
                n12 = size(qw12, 1) - 1
                m12 = size(qw12, 2) - 1
                call split_surface(n2i, p1, uk2i, m2i, q1, vk2i, qw2i, vi, 'v', &
                                   uk13, vk13, qw13, uk14, vk14, qw14)
                n13 = size(qw13, 1) - 1
                m13 = size(qw13, 2) - 1
                n14 = size(qw14, 1) - 1
                m14 = size(qw14, 2) - 1
                ! Deallocate working arrays.
                deallocate(uk1i, uk2i, vk1i, vk2i, qw1i, qw2i)
                
                ! Check potential intersections.
                call candidate_intersect(n11, m11, qw11, n2, m2, cpw2, t11_t21)
                call candidate_intersect(n12, m12, qw12, n2, m2, cpw2, t12_t21)
                call candidate_intersect(n13, m13, qw13, n2, m2, cpw2, t13_t21)
                call candidate_intersect(n14, m14, qw14, n2, m2, cpw2, t14_t21)
                
                ! Allocate new cells
                allocate(ne_ptr1%ci, se_ptr1%ci, sw_ptr1%ci, nw_ptr1%ci)
                
                ! Assign children (if any) in case cells are duplicated.
                if (associated(cell_ptr1%ci%ne)) then
                    ne_ptr1%ci => cell_ptr1%ci%ne
                end if
                if (associated(cell_ptr1%ci%se)) then
                    se_ptr1%ci => cell_ptr1%ci%se
                end if
                if (associated(cell_ptr1%ci%sw)) then
                    sw_ptr1%ci => cell_ptr1%ci%sw
                end if
                if (associated(cell_ptr1%ci%nw)) then
                    nw_ptr1%ci => cell_ptr1%ci%nw
                end if
                
                ! Assign children.
                if (t11_t21) then
                    ! Surface 1: Cell 1 (SW)                
                    if (sw_ptr1%ci%sid .eq. 0) then
                        sw_ptr1%ci%sid = ncells1 + 1
                        ncells1 = ncells1 + 1
                    end if
                    cell_ptr1%ci%sw => sw_ptr1%ci
                    sw_ptr1%ci%parent => cell_ptr1%ci
                    sw_ptr1%ci%position = 1
                    cell_ptr1%ci%has_child = .true.
                end if
                
                if (t12_t21) then
                    ! Surface 1: Cell 2 (NW)
                    if (nw_ptr1%ci%sid .eq. 0) then
                        nw_ptr1%ci%sid = ncells1 + 1  
                        ncells1 = ncells1 + 1
                    end if
                    cell_ptr1%ci%nw => nw_ptr1%ci
                    nw_ptr1%ci%parent => cell_ptr1%ci
                    nw_ptr1%ci%position = 2
                    cell_ptr1%ci%has_child = .true.                
                end if
                
                if (t13_t21) then
                    ! Surface 1: Cell 3 (SE)
                    if (se_ptr1%ci%sid .eq. 0) then
                        se_ptr1%ci%sid = ncells1 + 1
                        ncells1 = ncells1 + 1       
                    end if
                    cell_ptr1%ci%se => se_ptr1%ci
                    se_ptr1%ci%parent => cell_ptr1%ci
                    se_ptr1%ci%position = 3
                    cell_ptr1%ci%has_child = .true.
                end if
                
                if (t14_t21) then
                    ! Surface 1: Cell 4 (NE)
                    if (ne_ptr1%ci%sid .eq. 0) then
                        ne_ptr1%ci%sid = ncells1 + 1
                        ncells1 = ncells1 + 1       
                    end if
                    cell_ptr1%ci%ne => ne_ptr1%ci
                    ne_ptr1%ci%parent => cell_ptr1%ci
                    ne_ptr1%ci%position = 4
                    cell_ptr1%ci%has_child = .true.
                end if
                
                ! Adjacent cells.
                if (t11_t21) then
                    ! Surface 1: Cell 1 (SW)
                    sw_ptr1%ci%n => nw_ptr1%ci
                    sw_ptr1%ci%e => se_ptr1%ci
                    sw_ptr1%ci%s => cell_ptr1%ci%s
                    sw_ptr1%ci%w => cell_ptr1%ci%w
                end if
                
                if (t12_t21) then
                    ! Surface 1: Cell 2 (NW)
                    nw_ptr1%ci%n => cell_ptr1%ci%n
                    nw_ptr1%ci%e => ne_ptr1%ci
                    nw_ptr1%ci%s => sw_ptr1%ci
                    nw_ptr1%ci%w => cell_ptr1%ci%w
                end if
                
                if (t13_t21) then
                    ! Surface 1: Cell 3 (SE)
                    se_ptr1%ci%n => ne_ptr1%ci
                    se_ptr1%ci%e => cell_ptr1%ci%e
                    se_ptr1%ci%s => cell_ptr1%ci%s
                    se_ptr1%ci%w => sw_ptr1%ci
                end if
                
                if (t14_t21) then
                    ! Surface 1: Cell 4 (NE)
                    ne_ptr1%ci%n => cell_ptr1%ci%n
                    ne_ptr1%ci%e => cell_ptr1%ci%e
                    ne_ptr1%ci%s => se_ptr1%ci
                    ne_ptr1%ci%w => nw_ptr1%ci
                end if
                
                ! Subdivide
                ! Surface 1 Cell 1
                if (t11_t21) then
                    call subdivide(n11, p1, uk11, m11, q1, vk11, qw11, &
                                   n2, p2, uk2, m2, q2, vk2, cpw2, &
                                   sw_ptr1, cell_ptr2)
                end if

                ! Surface 1 Cell 2
                if (t12_t21) then
                    call subdivide(n12, p1, uk12, m12, q1, vk12, qw12, &
                                   n2, p2, uk2, m2, q2, vk2, cpw2, &
                                   nw_ptr1, cell_ptr2)
                end if
                
                ! Surface 1 Cell 3
                if (t13_t21) then
                    call subdivide(n13, p1, uk13, m13, q1, vk13, qw13, &
                                   n2, p2, uk2, m2, q2, vk2, cpw2, &
                                   se_ptr1, cell_ptr2)
                end if
                
                ! Surface 1 Cell 4
                if (t14_t21) then
                    call subdivide(n14, p1, uk14, m14, q1, vk14, qw14, &
                                   n2, p2, uk2, m2, q2, vk2, cpw2, &
                                   ne_ptr1, cell_ptr2)
                end if
            
            elseif ((is_flat2 .eqv. .false.) .and. is_flat1) then
                ! Subdivide surface 2
                ui = 0.50d0 * (cell_ptr2%ci%u0 + cell_ptr2%ci%u1)
                vi = 0.50d0 * (cell_ptr2%ci%v0 + cell_ptr2%ci%v1)
                
                ! Split at first interior knot if linear.
                ukn = size(uk2, 1)
                vkn = size(vk2, 1)
                if ((p2 .eq. 1) .and. (ukn .gt. 4)) then
                    mid = ceiling((dble(ukn) - 1.0d0) / 2.0d0)
                    ui = uk2(mid)
                end if
                if ((q2 .eq. 1) .and. (vkn .gt. 4)) then
                    mid = ceiling((dble(vkn) - 1.0d0) / 2.0d0)
                    vi = vk2(mid)
                end if
                ! Split surface at u-paramter.
                call split_surface(n2, p2, uk2, m2, q2, vk2, cpw2, ui, 'u', &
                                   uk1i, vk1i, qw1i, uk2i, vk2i, qw2i)
                n1i = size(qw1i, 1) - 1
                m1i = size(qw1i, 2) - 1
                n2i = size(qw2i, 1) - 1
                m2i = size(qw2i, 2) - 1
                ! Split halves in v-direction and store.
                call split_surface(n1i, p2, uk1i, m1i, q2, vk1i, qw1i, vi, 'v', &
                                   uk21, vk21, qw21, uk22, vk22, qw22)
                n21 = size(qw21, 1) - 1
                m21 = size(qw21, 2) - 1
                n22 = size(qw22, 1) - 1
                m22 = size(qw22, 2) - 1
                call split_surface(n2i, p2, uk2i, m2i, q2, vk2i, qw2i, vi, 'v', &
                                   uk23, vk23, qw23, uk24, vk24, qw24)
                n23 = size(qw23, 1) - 1
                m23 = size(qw23, 2) - 1
                n24 = size(qw24, 1) - 1
                m24 = size(qw24, 2) - 1
                ! Deallocate working arrays.
                deallocate(uk1i, uk2i, vk1i, vk2i, qw1i, qw2i)
                
                ! Check potential intersections.
                call candidate_intersect(n1, m1, cpw1, n21, m21, qw21, t11_t21)
                call candidate_intersect(n1, m1, cpw1, n22, m22, qw22, t11_t22)
                call candidate_intersect(n1, m1, cpw1, n23, m23, qw23, t11_t23)
                call candidate_intersect(n1, m1, cpw1, n24, m24, qw24, t11_t24)
            
                ! Allocate new cells
                allocate(ne_ptr2%ci, se_ptr2%ci, sw_ptr2%ci, nw_ptr2%ci)
            
                ! Assign children (if any) in case cells are duplicated.
                if (associated(cell_ptr2%ci%ne)) then
                    ne_ptr2%ci => cell_ptr2%ci%ne
                end if
                if (associated(cell_ptr2%ci%se)) then
                    se_ptr2%ci => cell_ptr2%ci%se
                end if
                if (associated(cell_ptr2%ci%sw)) then
                    sw_ptr2%ci => cell_ptr2%ci%sw
                end if
                if (associated(cell_ptr2%ci%nw)) then
                    nw_ptr2%ci => cell_ptr2%ci%nw
                end if
                
                ! Assign children
                if (t11_t21) then
                    ! Surface 2: Cell 1 (SW)
                    if (sw_ptr2%ci%sid .eq. 0) then
                        sw_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%sw => sw_ptr2%ci
                    sw_ptr2%ci%parent => cell_ptr2%ci
                    sw_ptr2%ci%position = 1
                    cell_ptr2%ci%has_child = .true.         
                end if
                
                if (t11_t22) then
                    ! Surface 2: Cell 2 (NW)
                    if (nw_ptr2%ci%sid .eq. 0) then
                        nw_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%nw => nw_ptr2%ci
                    nw_ptr2%ci%parent => cell_ptr2%ci
                    nw_ptr2%ci%position = 2
                    cell_ptr2%ci%has_child = .true.     
                end if
                
                if (t11_t23) then
                    ! Surface 2: Cell 3 (SE)
                    if (se_ptr2%ci%sid .eq. 0) then
                        se_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%se => se_ptr2%ci
                    se_ptr2%ci%parent => cell_ptr2%ci
                    se_ptr2%ci%position = 3
                    cell_ptr2%ci%has_child = .true.           
                end if
                
                if (t11_t24) then
                    ! Surface 2: Cell 4 (NE)
                    if (ne_ptr2%ci%sid .eq. 0) then
                        ne_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%ne => ne_ptr2%ci
                    ne_ptr2%ci%parent => cell_ptr2%ci
                    ne_ptr2%ci%position = 4
                    cell_ptr2%ci%has_child = .true.
                end if
                
                ! Adjacent cells
                if (t11_t21) then
                    ! Surface 2: Cell 1 (SW)
                    sw_ptr2%ci%n => nw_ptr2%ci
                    sw_ptr2%ci%e => se_ptr2%ci
                    sw_ptr2%ci%s => cell_ptr2%ci%s
                    sw_ptr2%ci%w => cell_ptr2%ci%w
                end if
                
                if (t11_t22) then
                    ! Surface 2: Cell 2 (NW)
                    nw_ptr2%ci%n => cell_ptr2%ci%n
                    nw_ptr2%ci%e => ne_ptr2%ci
                    nw_ptr2%ci%s => sw_ptr2%ci
                    nw_ptr2%ci%w => cell_ptr2%ci%w
                end if
                
                if (t11_t23) then
                    ! Surface 2: Cell 3 (SE)
                    se_ptr2%ci%n => ne_ptr2%ci
                    se_ptr2%ci%e => cell_ptr2%ci%e
                    se_ptr2%ci%s => cell_ptr2%ci%s
                    se_ptr2%ci%w => sw_ptr2%ci
                end if
                
                if (t11_t24) then
                    ! Surface 2: Cell 4 (NE)
                    ne_ptr2%ci%n => cell_ptr2%ci%n
                    ne_ptr2%ci%e => cell_ptr2%ci%e
                    ne_ptr2%ci%s => se_ptr2%ci
                    ne_ptr2%ci%w => nw_ptr2%ci
                end if
                
                ! Subdivide
                ! Surface 2 Cell 1
                if (t11_t21) then
                    call subdivide(n1, p1, uk1, m1, q1, vk1, cpw1, &
                                   n21, p2, uk21, m21, q2, vk21, qw21, &
                                   cell_ptr1, sw_ptr2)
                end if
                
                ! Surface 2 Cell 2
                if (t11_t22) then
                    call subdivide(n1, p1, uk1, m1, q1, vk1, cpw1, &
                                   n22, p2, uk22, m22, q2, vk22, qw22, &
                                   cell_ptr1, nw_ptr2)
                end if
                
                ! Surface 2 Cell 3
                if (t11_t23) then
                    call subdivide(n1, p1, uk1, m1, q1, vk1, cpw1, &
                                   n23, p2, uk23, m23, q2, vk23, qw23, &
                                   cell_ptr1, se_ptr2)
                end if
                
                ! Surface 2 Cell 4
                if (t11_t24) then
                    call subdivide(n1, p1, uk1, m1, q1, vk1, cpw1, &
                                   n24, p2, uk24, m24, q2, vk24, qw24, &
                                   cell_ptr1, ne_ptr2)
                end if
            
            else
                ! Split surface 1.
                ui = 0.50d0 * (cell_ptr1%ci%u0 + cell_ptr1%ci%u1)
                vi = 0.50d0 * (cell_ptr1%ci%v0 + cell_ptr1%ci%v1)
                
                ! If linear, split at middle knot.
                ukn = size(uk1, 1)
                vkn = size(vk1, 1)
                if ((p1 .eq. 1) .and. (ukn .gt. 4)) then
                    mid = ceiling((dble(ukn) - 1.0d0) / 2.0d0)
                    ui = uk1(mid)
                end if
                if ((q1 .eq. 1) .and. (vkn .gt. 4)) then
                    mid = ceiling((dble(vkn) - 1.0d0) / 2.0d0)
                    vi = vk1(mid)
                end if
                ! Split surface at u-paramter.
                call split_surface(n1, p1, uk1, m1, q1, vk1, cpw1, ui, 'u', &
                                   uk1i, vk1i, qw1i, uk2i, vk2i, qw2i)
                n1i = size(qw1i, 1) - 1
                m1i = size(qw1i, 2) - 1
                n2i = size(qw2i, 1) - 1
                m2i = size(qw2i, 2) - 1
                ! Split halves in v-direction and store.
                call split_surface(n1i, p1, uk1i, m1i, q1, vk1i, qw1i, vi, 'v', &
                                   uk11, vk11, qw11, uk12, vk12, qw12)
                n11 = size(qw11, 1) - 1
                m11 = size(qw11, 2) - 1
                n12 = size(qw12, 1) - 1
                m12 = size(qw12, 2) - 1
                call split_surface(n2i, p1, uk2i, m2i, q1, vk2i, qw2i, vi, 'v', &
                                   uk13, vk13, qw13, uk14, vk14, qw14)
                n13 = size(qw13, 1) - 1
                m13 = size(qw13, 2) - 1
                n14 = size(qw14, 1) - 1
                m14 = size(qw14, 2) - 1
                ! Deallocate working arrays.
                deallocate(uk1i, uk2i, vk1i, vk2i, qw1i, qw2i)
                
                ! Split surface 2.
                ui = 0.50d0 * (cell_ptr2%ci%u0 + cell_ptr2%ci%u1)
                vi = 0.50d0 * (cell_ptr2%ci%v0 + cell_ptr2%ci%v1)
                
                ! Split at first interior knot if linear.
                ukn = size(uk2, 1)
                vkn = size(vk2, 1)
                if ((p2 .eq. 1) .and. (ukn .gt. 4)) then
                    mid = ceiling((dble(ukn) - 1.0d0) / 2.0d0)
                    ui = uk2(mid)
                end if
                if ((q2 .eq. 1) .and. (vkn .gt. 4)) then
                    mid = ceiling((dble(vkn) - 1.0d0) / 2.0d0)
                    vi = vk2(mid)
                end if
                ! Split surface at u-paramter.
                call split_surface(n2, p2, uk2, m2, q2, vk2, cpw2, ui, 'u', &
                                   uk1i, vk1i, qw1i, uk2i, vk2i, qw2i)
                n1i = size(qw1i, 1) - 1
                m1i = size(qw1i, 2) - 1
                n2i = size(qw2i, 1) - 1
                m2i = size(qw2i, 2) - 1
                ! Split halves in v-direction and store.
                call split_surface(n1i, p2, uk1i, m1i, q2, vk1i, qw1i, vi, 'v', &
                                   uk21, vk21, qw21, uk22, vk22, qw22)
                n21 = size(qw21, 1) - 1
                m21 = size(qw21, 2) - 1
                n22 = size(qw22, 1) - 1
                m22 = size(qw22, 2) - 1
                call split_surface(n2i, p2, uk2i, m2i, q2, vk2i, qw2i, vi, 'v', &
                                   uk23, vk23, qw23, uk24, vk24, qw24)
                n23 = size(qw23, 1) - 1
                m23 = size(qw23, 2) - 1
                n24 = size(qw24, 1) - 1
                m24 = size(qw24, 2) - 1
                ! Deallocate working arrays.
                deallocate(uk1i, uk2i, vk1i, vk2i, qw1i, qw2i)
                
                ! Check potential intersections.
                call candidate_intersect(n11, m11, qw11, n21, m21, qw21, t11_t21)
                call candidate_intersect(n11, m11, qw11, n22, m22, qw22, t11_t22)
                call candidate_intersect(n11, m11, qw11, n23, m23, qw23, t11_t23)
                call candidate_intersect(n11, m11, qw11, n24, m24, qw24, t11_t24)
                call candidate_intersect(n12, m12, qw12, n21, m21, qw21, t12_t21)
                call candidate_intersect(n12, m12, qw12, n22, m22, qw22, t12_t22)
                call candidate_intersect(n12, m12, qw12, n23, m23, qw23, t12_t23)
                call candidate_intersect(n12, m12, qw12, n24, m24, qw24, t12_t24)
                call candidate_intersect(n13, m13, qw13, n21, m21, qw21, t13_t21)
                call candidate_intersect(n13, m13, qw13, n22, m22, qw22, t13_t22)
                call candidate_intersect(n13, m13, qw13, n23, m23, qw23, t13_t23)
                call candidate_intersect(n13, m13, qw13, n24, m24, qw24, t13_t24)
                call candidate_intersect(n14, m14, qw14, n21, m21, qw21, t14_t21)
                call candidate_intersect(n14, m14, qw14, n22, m22, qw22, t14_t22)
                call candidate_intersect(n14, m14, qw14, n23, m23, qw23, t14_t23)
                call candidate_intersect(n14, m14, qw14, n24, m24, qw24, t14_t24)
                
                ! Allocate new cells
                allocate(ne_ptr1%ci, se_ptr1%ci, sw_ptr1%ci, nw_ptr1%ci)
                allocate(ne_ptr2%ci, se_ptr2%ci, sw_ptr2%ci, nw_ptr2%ci)
                
                ! Assign children (if any) in case cells are duplicated.
                if (associated(cell_ptr1%ci%ne)) then
                    ne_ptr1%ci => cell_ptr1%ci%ne
                end if
                if (associated(cell_ptr1%ci%se)) then
                    se_ptr1%ci => cell_ptr1%ci%se
                end if
                if (associated(cell_ptr1%ci%sw)) then
                    sw_ptr1%ci => cell_ptr1%ci%sw
                end if
                if (associated(cell_ptr1%ci%nw)) then
                    nw_ptr1%ci => cell_ptr1%ci%nw
                end if
                
                if (associated(cell_ptr2%ci%ne)) then
                    ne_ptr2%ci => cell_ptr2%ci%ne
                end if
                if (associated(cell_ptr2%ci%se)) then
                    se_ptr2%ci => cell_ptr2%ci%se
                end if
                if (associated(cell_ptr2%ci%sw)) then
                    sw_ptr2%ci => cell_ptr2%ci%sw
                end if
                if (associated(cell_ptr2%ci%nw)) then
                    nw_ptr2%ci => cell_ptr2%ci%nw
                end if
                
                ! Assign children.
                if (t11_t21 .or. t11_t22 .or. t11_t23 .or. t11_t24) then
                    ! Surface 1: Cell 1 (SW)                
                    if (sw_ptr1%ci%sid .eq. 0) then
                        sw_ptr1%ci%sid = ncells1 + 1
                        ncells1 = ncells1 + 1
                    end if
                    cell_ptr1%ci%sw => sw_ptr1%ci
                    sw_ptr1%ci%parent => cell_ptr1%ci
                    sw_ptr1%ci%position = 1
                    cell_ptr1%ci%has_child = .true.
                end if
                
                if (t12_t21 .or. t12_t22 .or. t12_t23 .or. t12_t24) then
                    ! Surface 1: Cell 2 (NW)
                    if (nw_ptr1%ci%sid .eq. 0) then
                        nw_ptr1%ci%sid = ncells1 + 1  
                        ncells1 = ncells1 + 1
                    end if
                    cell_ptr1%ci%nw => nw_ptr1%ci
                    nw_ptr1%ci%parent => cell_ptr1%ci
                    nw_ptr1%ci%position = 2
                    cell_ptr1%ci%has_child = .true.                
                end if
                
                if (t13_t21 .or. t13_t22 .or. t13_t23 .or. t13_t24) then
                    ! Surface 1: Cell 3 (SE)
                    if (se_ptr1%ci%sid .eq. 0) then
                        se_ptr1%ci%sid = ncells1 + 1
                        ncells1 = ncells1 + 1       
                    end if
                    cell_ptr1%ci%se => se_ptr1%ci
                    se_ptr1%ci%parent => cell_ptr1%ci
                    se_ptr1%ci%position = 3
                    cell_ptr1%ci%has_child = .true.
                end if
                
                if (t14_t21 .or. t14_t22 .or. t14_t23 .or. t14_t24) then
                    ! Surface 1: Cell 4 (NE)
                    if (ne_ptr1%ci%sid .eq. 0) then
                        ne_ptr1%ci%sid = ncells1 + 1
                        ncells1 = ncells1 + 1       
                    end if
                    cell_ptr1%ci%ne => ne_ptr1%ci
                    ne_ptr1%ci%parent => cell_ptr1%ci
                    ne_ptr1%ci%position = 4
                    cell_ptr1%ci%has_child = .true.
                end if
                
                if (t11_t21 .or. t12_t21 .or. t13_t21 .or. t14_t21) then
                    ! Surface 2: Cell 1 (SW)
                    if (sw_ptr2%ci%sid .eq. 0) then
                        sw_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%sw => sw_ptr2%ci
                    sw_ptr2%ci%parent => cell_ptr2%ci
                    sw_ptr2%ci%position = 1
                    cell_ptr2%ci%has_child = .true.         
                end if
                
                if (t11_t22 .or. t12_t22 .or. t13_t22 .or. t14_t22) then
                    ! Surface 2: Cell 2 (NW)
                    if (nw_ptr2%ci%sid .eq. 0) then
                        nw_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%nw => nw_ptr2%ci
                    nw_ptr2%ci%parent => cell_ptr2%ci
                    nw_ptr2%ci%position = 2
                    cell_ptr2%ci%has_child = .true.     
                end if
                
                if (t11_t23 .or. t12_t23 .or. t13_t23 .or. t14_t23) then
                    ! Surface 2: Cell 3 (SE)
                    if (se_ptr2%ci%sid .eq. 0) then
                        se_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%se => se_ptr2%ci
                    se_ptr2%ci%parent => cell_ptr2%ci
                    se_ptr2%ci%position = 3
                    cell_ptr2%ci%has_child = .true.           
                end if
                
                if (t11_t24 .or. t12_t24 .or. t13_t24 .or. t14_t24) then
                    ! Surface 2: Cell 4 (NE)
                    if (ne_ptr2%ci%sid .eq. 0) then
                        ne_ptr2%ci%sid = ncells2 + 1
                        ncells2 = ncells2 + 1   
                    end if
                    cell_ptr2%ci%ne => ne_ptr2%ci
                    ne_ptr2%ci%parent => cell_ptr2%ci
                    ne_ptr2%ci%position = 4
                    cell_ptr2%ci%has_child = .true.
                end if
                
                ! Adjacent cells.
                if (t11_t21 .or. t11_t22 .or. t11_t23 .or. t11_t24) then
                    ! Surface 1: Cell 1 (SW)
                    sw_ptr1%ci%n => nw_ptr1%ci
                    sw_ptr1%ci%e => se_ptr1%ci
                    sw_ptr1%ci%s => cell_ptr1%ci%s
                    sw_ptr1%ci%w => cell_ptr1%ci%w
                end if
                
                if (t12_t21 .or. t12_t22 .or. t12_t23 .or. t12_t24) then
                    ! Surface 1: Cell 2 (NW)
                    nw_ptr1%ci%n => cell_ptr1%ci%n
                    nw_ptr1%ci%e => ne_ptr1%ci
                    nw_ptr1%ci%s => sw_ptr1%ci
                    nw_ptr1%ci%w => cell_ptr1%ci%w
                end if
                
                if (t13_t21 .or. t13_t22 .or. t13_t23 .or. t13_t24) then
                    ! Surface 1: Cell 3 (SE)
                    se_ptr1%ci%n => ne_ptr1%ci
                    se_ptr1%ci%e => cell_ptr1%ci%e
                    se_ptr1%ci%s => cell_ptr1%ci%s
                    se_ptr1%ci%w => sw_ptr1%ci
                end if
                
                if (t14_t21 .or. t14_t22 .or. t14_t23 .or. t14_t24) then
                    ! Surface 1: Cell 4 (NE)
                    ne_ptr1%ci%n => cell_ptr1%ci%n
                    ne_ptr1%ci%e => cell_ptr1%ci%e
                    ne_ptr1%ci%s => se_ptr1%ci
                    ne_ptr1%ci%w => nw_ptr1%ci
                end if
                
                if (t11_t21 .or. t12_t21 .or. t13_t21 .or. t14_t21) then
                    ! Surface 2: Cell 1 (SW)
                    sw_ptr2%ci%n => nw_ptr2%ci
                    sw_ptr2%ci%e => se_ptr2%ci
                    sw_ptr2%ci%s => cell_ptr2%ci%s
                    sw_ptr2%ci%w => cell_ptr2%ci%w
                end if
                
                if (t11_t22 .or. t12_t22 .or. t13_t22 .or. t14_t22) then
                    ! Surface 2: Cell 2 (NW)
                    nw_ptr2%ci%n => cell_ptr2%ci%n
                    nw_ptr2%ci%e => ne_ptr2%ci
                    nw_ptr2%ci%s => sw_ptr2%ci
                    nw_ptr2%ci%w => cell_ptr2%ci%w
                end if
                
                if (t11_t23 .or. t12_t23 .or. t13_t23 .or. t14_t23) then
                    ! Surface 2: Cell 3 (SE)
                    se_ptr2%ci%n => ne_ptr2%ci
                    se_ptr2%ci%e => cell_ptr2%ci%e
                    se_ptr2%ci%s => cell_ptr2%ci%s
                    se_ptr2%ci%w => sw_ptr2%ci
                end if
                
                if (t11_t24 .or. t12_t24 .or. t13_t24 .or. t14_t24) then
                    ! Surface 2: Cell 4 (NE)
                    ne_ptr2%ci%n => cell_ptr2%ci%n
                    ne_ptr2%ci%e => cell_ptr2%ci%e
                    ne_ptr2%ci%s => se_ptr2%ci
                    ne_ptr2%ci%w => nw_ptr2%ci
                end if
                
                ! Subdivide
                ! Surface 1 Cell 1 with...
                if (t11_t21) then
                    call subdivide(n11, p1, uk11, m11, q1, vk11, qw11, &
                                   n21, p2, uk21, m21, q2, vk21, qw21, &
                                   sw_ptr1, sw_ptr2)
                end if
                if (t11_t22) then
                    call subdivide(n11, p1, uk11, m11, q1, vk11, qw11, &
                                   n22, p2, uk22, m22, q2, vk22, qw22, &
                                   sw_ptr1, nw_ptr2)
                end if
                if (t11_t23) then
                    call subdivide(n11, p1, uk11, m11, q1, vk11, qw11, &
                                   n23, p2, uk23, m23, q2, vk23, qw23, &
                                   sw_ptr1, se_ptr2)
                end if
                if (t11_t24) then
                    call subdivide(n11, p1, uk11, m11, q1, vk11, qw11, &
                                   n24, p2, uk24, m24, q2, vk24, qw24, &
                                   sw_ptr1, ne_ptr2)
                end if
                
                ! Surface 1 Cell 2 with...
                if (t12_t21) then
                    call subdivide(n12, p1, uk12, m12, q1, vk12, qw12, &
                                   n21, p2, uk21, m21, q2, vk21, qw21, &
                                   nw_ptr1, sw_ptr2)
                end if
                if (t12_t22) then
                    call subdivide(n12, p1, uk12, m12, q1, vk12, qw12, &
                                   n22, p2, uk22, m22, q2, vk22, qw22, &
                                   nw_ptr1, nw_ptr2)
                end if
                if (t12_t23) then
                    call subdivide(n12, p1, uk12, m12, q1, vk12, qw12, &
                                   n23, p2, uk23, m23, q2, vk23, qw23, &
                                   nw_ptr1, se_ptr2)
                end if
                if (t12_t24) then
                    call subdivide(n12, p1, uk12, m12, q1, vk12, qw12, &
                                   n24, p2, uk24, m24, q2, vk24, qw24, &
                                   nw_ptr1, ne_ptr2)
                end if
                
                ! Surface 1 Cell 3 with...
                if (t13_t21) then
                    call subdivide(n13, p1, uk13, m13, q1, vk13, qw13, &
                                   n21, p2, uk21, m21, q2, vk21, qw21, &
                                   se_ptr1, sw_ptr2)
                end if
                if (t13_t22) then
                    call subdivide(n13, p1, uk13, m13, q1, vk13, qw13, &
                                   n22, p2, uk22, m22, q2, vk22, qw22, &
                                   se_ptr1, nw_ptr2)
                end if
                if (t13_t23) then
                    call subdivide(n13, p1, uk13, m13, q1, vk13, qw13, &
                                   n23, p2, uk23, m23, q2, vk23, qw23, &
                                   se_ptr1, se_ptr2)
                end if
                if (t13_t24) then
                    call subdivide(n13, p1, uk13, m13, q1, vk13, qw13, &
                                   n24, p2, uk24, m24, q2, vk24, qw24, &
                                   se_ptr1, ne_ptr2)
                end if
                
                ! Surface 1 Cell 4 with...
                if (t14_t21) then
                    call subdivide(n14, p1, uk14, m14, q1, vk14, qw14, &
                                   n21, p2, uk21, m21, q2, vk21, qw21, &
                                   ne_ptr1, sw_ptr2)
                end if
                if (t14_t22) then
                    call subdivide(n14, p1, uk14, m14, q1, vk14, qw14, &
                                   n22, p2, uk22, m22, q2, vk22, qw22, &
                                   ne_ptr1, nw_ptr2)
                end if
                if (t14_t23) then
                    call subdivide(n14, p1, uk14, m14, q1, vk14, qw14, &
                                   n23, p2, uk23, m23, q2, vk23, qw23, &
                                   ne_ptr1, se_ptr2)
                end if
                if (t14_t24) then
                    call subdivide(n14, p1, uk14, m14, q1, vk14, qw14, &
                                   n24, p2, uk24, m24, q2, vk24, qw24, &
                                   ne_ptr1, ne_ptr2)
                end if
            end if
        end if
    end subroutine subdivide

    subroutine candidate_intersect(n1, m1, cpw1, n2, m2, cpw2, candidate)
        ! Check for possible intersection using bounding box.
        ! Input
        integer, intent(in) :: n1, m1, n2, m2
        double precision, intent(in) :: cpw1(0:n1, 0:m1, 0:3)
        double precision, intent(in) :: cpw2(0:n2, 0:m2, 0:3)
        
        ! Output
        logical, intent(out) :: candidate
        
        ! Working
        double precision :: bbox1(0:2, 0:1), bbox2(0:2, 0:1)
        
        bbox1 = surface_bbox(n1, m1, cpw1)
        bbox2 = surface_bbox(n2, m2, cpw2)
        
        candidate = bboxes_intersect(3, bbox1, bbox2, gtol)
        
    end subroutine candidate_intersect
    
    subroutine is_flat(n, m, cp, flat)
        ! Check surface flatness.
        ! Input
        integer, intent(in) :: n, m
        double precision, intent(in) :: cp(0:n, 0:m, 0:2)
        
        ! Output
        logical, intent(out) :: flat
        
        call is_surface_flat(n, m, cp, ftol, flat)
        
    end subroutine is_flat

    subroutine intersect(icell1, icell2)
        ! Intersect cells using triangle-triangle intersection.
        ! Input
        integer, intent(in) :: icell1, icell2
        
        ! Working
        integer :: ii, jj
        integer :: ntri1, ntri2, ni
        double precision :: ub, vb, area1, area2
        double precision :: uv10(0:1), uv11(0:1), uv12(0:1)
        double precision :: uv20(0:1), uv21(0:1), uv22(0:1)
        double precision :: t0(0:2), t1(0:2), t2(0:2)
        double precision :: ti1(0:2, 0:2)
        double precision :: ti2(0:2, 0:2)
        double precision :: pi(0:1, 0:2), pb(0:2)

        ! Tessellate cells.
        call tessellate_cell(icell1, ndim, children1, acells1, position1, &
                             parent1, cell_params1, ntri1, tri1)
        call tessellate_cell(icell2, ndim, children2, acells2, position2, &
                             parent2, cell_params2, ntri2, tri2)
                             
        ! Intersect triangles.
        do ii = 0, ntri1 - 1
            uv10 = tri1(ii, 0, :)
            uv11 = tri1(ii, 1, :)
            uv12 = tri1(ii, 2, :)
            call surface_point(n1, p1, uk1, m1, q1, vk1, cpw1, &
                               uv10(0), uv10(1), t0)
            call surface_point(n1, p1, uk1, m1, q1, vk1, cpw1, &
                               uv11(0), uv11(1), t1)
            call surface_point(n1, p1, uk1, m1, q1, vk1, cpw1, &
                               uv12(0), uv12(1), t2)
            ti1(0, :) = t0
            ti1(1, :) = t1
            ti1(2, :) = t2
            ! Check for degenerate triangle.
            area1 = triangle_area(3, ti1)
            if (area1 .gt. 1.0d-12) then            
                do jj = 0, ntri2 - 1
                    uv20 = tri2(jj, 0, :)
                    uv21 = tri2(jj, 1, :)
                    uv22 = tri2(jj, 2, :)
                    call surface_point(n2, p2, uk2, m2, q2, vk2, cpw2, &
                                       uv20(0), uv20(1), t0)
                    call surface_point(n2, p2, uk2, m2, q2, vk2, cpw2, &
                                       uv21(0), uv21(1), t1)
                    call surface_point(n2, p2, uk2, m2, q2, vk2, cpw2, &
                                       uv22(0), uv22(1), t2)
                    ti2(0, :) = t0
                    ti2(1, :) = t1
                    ti2(2, :) = t2
                    ! Check for degenerate triangle.
                    area2 = triangle_area(3, ti2)
                    if (area2 .gt. 1.0d-12) then
                        call intersect_triangles(ti1, ti2, gtol, ptol, ni, pi)                        
                        if (ni .eq. 2) then
                            ! Add points and edges.
                            edges(etotal, :) = (/ nverts_, nverts_ + 1 /)
                            verts_(nverts_, :) = pi(0, :)
                            verts_(nverts_ + 1, :) = pi(1, :)
                            ! Surface 1 parameters.
                            pb = verts_(nverts_, :)
                            call barycentric_params(pb, ti1, ub, vb)
                            ssi_params1_(nverts_, :) = (1.0d0 - ub - vb) * uv10 + ub * uv11 + &
                                                        vb * uv12
                            points2d_s1(nverts_, 0:1) = ssi_params1_(nverts_, :)
                            pb = verts_(nverts_ + 1, :)
                            call barycentric_params(pb, ti1, ub, vb)
                            ssi_params1_(nverts_ + 1, :) = (1.0d0 - ub - vb) * uv10 + ub * uv11 + &
                                                            vb * uv12
                            points2d_s1(nverts_ + 1, 0:1) = ssi_params1_(nverts_ + 1, :)
                            
                            ! Surface 2 parameters.
                            pb = verts_(nverts_, :)
                            call barycentric_params(pb, ti2, ub, vb)
                            ssi_params2_(nverts_, :) = (1.0d0 - ub - vb) * uv20 + ub * uv21 + &
                                                        vb * uv22
                            points2d_s2(nverts_, 0:1) = ssi_params2_(nverts_, :)
                            pb = verts_(nverts_ + 1, :)
                            call barycentric_params(pb, ti2, ub, vb)
                            ssi_params2_(nverts_ + 1, :) = (1.0d0 - ub - vb) * uv20 + ub * uv21 + &
                                                            vb * uv22
                            points2d_s2(nverts_ + 1, 0:1) = ssi_params2_(nverts_ + 1, :)
                            
                            ! Update totals.
                            etotal = etotal + 1
                            nverts_ = nverts_ + 2
                        end if
                    end if
                end do
            end if
        end do
        
    end subroutine intersect
    
end subroutine intersect_surface_surface

subroutine trace_curves(ptotal, etotal, points, edges, tol, nmax, &
                        ncrvs, crv_size, crv_ids)
    !> Trace topology of an unsorted set of edges.
    !> ptotal - Total number of points.
    !> etotal - Total number of edges.
    !> points - Points defining the vertices of the edges.
    !> edges - Unsorted edges where each row defines the index of its
    !> vertices.
    !> tol - Tolerance for equivalencing points.
    !> nmax - Maximum array size.
    !> ncrvs - Number of curves found.
    !> crv_size - Number of points in each curve.
    !> crv_ids - Defines the ordered vertices for each curve where each
    !> row is a curve and each column is an index in "points" array.
    
    ! Input
    integer, intent(in) :: ptotal, etotal, nmax, edges(0:nmax - 1, 0:1)
    double precision, intent(in) :: points(0:nmax - 1, 0:2), tol
    
    ! Output                                    
    integer, intent(out) :: ncrvs, crv_size(0:nmax - 1), &
                            crv_ids(0:nmax - 1, 0:nmax - 1)
                            
    ! Working
    logical :: unique, process_curves, point_found, is_closed
    integer :: i, j, k, nverts, point_use(0:ptotal - 1), &
               vert_to_point(0:ptotal - 1)
    double precision :: verts(0:ptotal - 1, 0:2)
    double precision :: v(0:2), p(0:2), d2
    
    integer :: new_edges(0:etotal - 1, 0:1), eid, &
               point_count(0:ptotal - 1), &
               visited(0:ptotal - 1), &
               adj_edges(0:ptotal - 1, 0:etotal - 1), &
               point_to_point(0:ptotal - 1, 0:ptotal - 1), pid1, pid2
               
    integer :: indx, nstack, stack(0:ptotal - 1), pid, adj_edge, vi
    
    ! Initialize
    ncrvs = 0
    crv_size(:) = 0
    crv_ids(:, :) = 0
    nverts = 0
    point_use(:) = 0
    vert_to_point(:) = 0
    verts(:, :) = 0.0d0
    new_edges(:, :) = 0
    point_count(:) = 0
    adj_edges(:, :) = 0
    visited(:) = 0
    point_to_point(:, :) = 0
    stack(:) = 0
    
    ! Equivalence array tracking point use id's.
    verts(0, :) = points(0, :)
    nverts = 1
    do i = 1, ptotal - 1
        p = points(i, :)
        unique = .true.
        do j = 0, nverts - 1
            v = verts(j, :)
            d2 = dsqrt((v(0) - p(0)) * (v(0) - p(0)) + &
                 (v(1) - p(1)) * (v(1) - p(1)) + &
                 (v(2) - p(2)) * (v(2) - p(2)))
            if (d2 .le. tol) then
                unique = .false.
                point_use(i) = j
                vert_to_point(j) = i
                exit
            end if
        end do
        if (unique) then
            verts(nverts, :) = points(i, :)
            point_use(i) = nverts
            vert_to_point(nverts) = i
            nverts = nverts + 1
        end if
    end do
    
    ! Build new edges.
    eid = 0
    do i = 0, etotal - 1
        pid1 = point_use(edges(i, 0))
        pid2 = point_use(edges(i, 1))
        ! Attempt to avoid duplicate edges.
        if ((point_to_point(pid1, pid2) .eq. 0) .and. &
                (point_to_point(pid2, pid1) .eq. 0) .and. &
                (pid1 .ne. pid2)) then
            point_to_point(pid1, pid2) = 1
            point_to_point(pid2, pid1) = 1
            new_edges(eid, :) = (/ pid1, pid2 /)
            adj_edges(pid1, point_count(pid1)) = eid
            adj_edges(pid2, point_count(pid2)) = eid
            point_count(pid1) = point_count(pid1) + 1
            point_count(pid2) = point_count(pid2) + 1
            eid = eid + 1
        end if
    end do
    
    ! Process curves until all points are visited.
    process_curves = .true.
    do while (process_curves .eqv. .true.)
        ! Try to find a point with only one adjacent edge. If all have
        ! more than one adjacent edge it may imply a closed curve. In
        ! that case start anywhere.
        point_found = .false.
        is_closed = .false.
        pid1 = 0
        do i = 0, nverts - 1
            if ((point_count(i) .eq. 1) .and. (visited(i) .eq. 0)) then
                pid1 = i
                point_found = .true.
                exit
            end if
        end do
        ! Select first unvisited point if no single point was found.
        if (point_found .eqv. .false.) then
            do i = 0, nverts - 1
                if ((point_count(i) .gt. 0) .and. (visited(i) .eq. 0)) then
                    pid1 = i
                    point_found = .true.
                    is_closed = .true.
                    exit
                end if
            end do
        end if
        
        ! Trace the topology of the curve using a depth-first search.
        if (point_found) then
            stack(0) = pid1
            nstack = 1
            indx = 0
            do while (nstack .gt. 0)
                pid = stack(nstack - 1)
                nstack = nstack - 1
                visited(pid) = 1
                crv_ids(ncrvs, indx) = vert_to_point(pid)
                indx = indx + 1
                k = point_count(pid)
                do i = 0, k - 1
                    adj_edge = adj_edges(pid, i)
                    do j = 0, 1
                        vi = new_edges(adj_edge, j)
                        if (visited(vi) .eq. 0) then
                            stack(nstack) = vi
                            nstack = nstack + 1
                            visited(vi) = 1
                        end if
                    end do
                end do
            end do
            if (is_closed) then
                crv_ids(ncrvs, indx) = crv_ids(ncrvs, 0)
                indx = indx + 1
            end if
            crv_size(ncrvs) = indx
            ncrvs = ncrvs + 1
        else
            process_curves = .false.
        end if
    end do
    
end subroutine trace_curves

subroutine trace_curves2(ptotal, etotal, points2d_s1, points2d_s2, edges, tol, nmax, &
                         ncrvs, crv_size, crv_ids)
    !> Trace topology of an unsorted set of edges.
    !> ptotal - Total number of points.
    !> etotal - Total number of edges.
    !> points2d_s1 - Points defining the vertices of the edges in the first surface.
    !> points2d_s2 - Points defining the vertices of the edges in the second surface.
    !> edges - Unsorted edges where each row defines the index of its
    !> vertices.
    !> tol - Tolerance for equivalencing points.
    !> nmax - Maximum array size.
    !> ncrvs - Number of curves found.
    !> crv_size - Number of points in each curve.
    !> crv_ids - Defines the ordered vertices for each curve where each
    !> row is a curve and each column is an index in "points" array.
    
    ! Input
    integer, intent(in) :: ptotal, etotal, nmax, edges(0:nmax - 1, 0:1)
    double precision, intent(in) :: tol
    double precision, intent(in) :: points2d_s1(0:nmax - 1, 0:2), points2d_s2(0:nmax - 1, 0:2)
    
    ! Output                                    
    integer, intent(out) :: ncrvs, crv_size(0:nmax - 1), &
                            crv_ids(0:nmax - 1, 0:nmax - 1)
                            
    ! Working
    logical :: unique, process_curves, point_found, is_closed
    integer :: i, j, k, nverts, point_use(0:ptotal - 1), &
               vert_to_point(0:ptotal - 1)
    double precision :: d1, d2
    double precision :: verts1(0:ptotal - 1, 0:2), verts2(0:ptotal - 1, 0:2)
    double precision :: v1(0:2), p1(0:2), v2(0:2), p2(0:2)
    
    integer :: new_edges(0:etotal - 1, 0:1), eid, &
               point_count(0:ptotal - 1), &
               visited(0:ptotal - 1), &
               adj_edges(0:ptotal - 1, 0:etotal - 1), &
               point_to_point(0:ptotal - 1, 0:ptotal - 1), pid1, pid2
               
    integer :: indx, nstack, stack(0:ptotal - 1), pid, adj_edge, vi
    
    ! Initialize
    ncrvs = 0
    crv_size(:) = 0
    crv_ids(:, :) = 0
    nverts = 0
    point_use(:) = 0
    vert_to_point(:) = 0
    verts1(:, :) = 0.0d0
    verts2(:, :) = 0.0d0
    new_edges(:, :) = 0
    point_count(:) = 0
    adj_edges(:, :) = 0
    visited(:) = 0
    point_to_point(:, :) = 0
    stack(:) = 0
    
    ! Equivalence array tracking point use id's.
    verts1(0, :) = points2d_s1(0, :)
    verts2(0, :) = points2d_s2(0, :)
    nverts = 1
    do i = 1, ptotal - 1
        p1 = points2d_s1(i, :)
        p2 = points2d_s2(i, :)
        unique = .true.
        do j = 0, nverts - 1
            v1 = verts1(j, :)
            d1 = dsqrt((v1(0) - p1(0)) * (v1(0) - p1(0)) + &
                 (v1(1) - p1(1)) * (v1(1) - p1(1)) + &
                 (v1(2) - p1(2)) * (v1(2) - p1(2)))
            v2 = verts2(j, :)
            d2 = dsqrt((v2(0) - p2(0)) * (v2(0) - p2(0)) + &
                 (v2(1) - p2(1)) * (v2(1) - p2(1)) + &
                 (v2(2) - p2(2)) * (v2(2) - p2(2)))
            if ((d1 .le. tol) .and. (d2 .le. tol)) then
                unique = .false.
                point_use(i) = j
                vert_to_point(j) = i
                exit
            end if
        end do
        if (unique) then
            verts1(nverts, :) = points2d_s1(i, :)
            verts2(nverts, :) = points2d_s2(i, :)
            point_use(i) = nverts
            vert_to_point(nverts) = i
            nverts = nverts + 1
        end if
    end do
    
    ! Build new edges.
    eid = 0
    do i = 0, etotal - 1
        pid1 = point_use(edges(i, 0))
        pid2 = point_use(edges(i, 1))
        ! Attempt to avoid duplicate edges.
        if ((point_to_point(pid1, pid2) .eq. 0) .and. &
                (point_to_point(pid2, pid1) .eq. 0) .and. &
                (pid1 .ne. pid2)) then
            point_to_point(pid1, pid2) = 1
            point_to_point(pid2, pid1) = 1
            new_edges(eid, :) = (/ pid1, pid2 /)
            adj_edges(pid1, point_count(pid1)) = eid
            adj_edges(pid2, point_count(pid2)) = eid
            point_count(pid1) = point_count(pid1) + 1
            point_count(pid2) = point_count(pid2) + 1
            eid = eid + 1
        end if
    end do
    
    ! Process curves until all points are visited.
    process_curves = .true.
    do while (process_curves .eqv. .true.)
        ! Try to find a point with only one adjacent edge. If all have
        ! more than one adjacent edge it may imply a closed curve. In
        ! that case start anywhere.
        point_found = .false.
        is_closed = .false.
        pid1 = 0
        do i = 0, nverts - 1
            if ((point_count(i) .eq. 1) .and. (visited(i) .eq. 0)) then
                pid1 = i
                point_found = .true.
                exit
            end if
        end do
        ! Select first unvisited point if no single point was found.
        if (point_found .eqv. .false.) then
            do i = 0, nverts - 1
                if ((point_count(i) .gt. 0) .and. (visited(i) .eq. 0)) then
                    pid1 = i
                    point_found = .true.
                    is_closed = .true.
                    exit
                end if
            end do
        end if
        
        ! Trace the topology of the curve using a depth-first search.
        if (point_found) then
            stack(0) = pid1
            nstack = 1
            indx = 0
            do while (nstack .gt. 0)
                pid = stack(nstack - 1)
                nstack = nstack - 1
                visited(pid) = 1
                crv_ids(ncrvs, indx) = vert_to_point(pid)
                indx = indx + 1
                k = point_count(pid)
                do i = 0, k - 1
                    adj_edge = adj_edges(pid, i)
                    do j = 0, 1
                        vi = new_edges(adj_edge, j)
                        if (visited(vi) .eq. 0) then
                            stack(nstack) = vi
                            nstack = nstack + 1
                            visited(vi) = 1
                        end if
                    end do
                end do
            end do
            if (is_closed) then
                crv_ids(ncrvs, indx) = crv_ids(ncrvs, 0)
                indx = indx + 1
            end if
            crv_size(ncrvs) = indx
            ncrvs = ncrvs + 1
        else
            process_curves = .false.
        end if
    end do
    
end subroutine trace_curves2

subroutine filter_points(ncrvs, nmax, crv_size, all_crv_ids, points, gtol, crv_ids)
    !> Filter out points based on tolerance.
    
    ! Input
    integer, intent(in) :: ncrvs, nmax
    integer, intent(in) :: all_crv_ids(0:nmax - 1, 0:nmax - 1)
    double precision, intent(in) :: gtol
    double precision, intent(in) :: points(0:nmax - 1, 0:2)
    
    ! Input/Output
    integer, intent(inout) :: crv_size(0:nmax - 1)
    
    ! Output
    integer, intent(out) :: crv_ids(0:nmax - 1, 0:nmax - 1)
    
    ! Working
    logical :: flat
    integer :: i, nids, i0, i1, i2, npts
    double precision :: d, dp, angle
    double precision :: v0(0:2), v1(0:2)

    do i = 0, ncrvs - 1
        nids = crv_size(i)
        if (nids .eq. 2) then
            crv_ids(i, 0) = all_crv_ids(i, 0)
            crv_ids(i, 1) = all_crv_ids(i, 1)
        else
            i0 = 0
            i1 = 1
            i2 = 2
            crv_ids(i, 0) = all_crv_ids(i, i0)
            npts = 1
            do while (i2 .lt. nids)
                flat = .false.
                v0 = points(all_crv_ids(i, i1), :) - points(all_crv_ids(i, i0), :)
                v1 = points(all_crv_ids(i, i2), :) - points(all_crv_ids(i, i1), :)
                ! Check minimum distance.
                d = norm(v0)
                if (d .le. gtol) then
                    flat = .true.
                end if
                ! Check for reversed points using the dot product and angle
                ! between the vectors.
                if (flat .eqv. .false.) then
                    dp = dot_product(v0, v1)
                    if (dp .lt. 0.0d0) then
                        angle = angle_between_vecs(v0, v1)
                        if (angle .gt. 170.0d0) then
                            flat = .true.
                        end if
                    end if
                end if
                ! Adjust indices and/or add point to curve.
                if (flat) then
                    i1 = i1 + 1
                    i2 = i1 + 1
                else
                    crv_ids(i, npts) = all_crv_ids(i, i1)
                    i0 = i1
                    i1 = i0 + 1
                    i2 = i1 + 1
                    npts = npts + 1
                end if
            end do
            ! Add last point or replace if previous point is coincident.
            v0 = points(crv_ids(i, npts - 1), :) - points(all_crv_ids(i, nids - 1), :)
            d = norm(v0)
            if (d .gt. gtol) then
                crv_ids(i, npts) = all_crv_ids(i, nids - 1)
                crv_size(i) = npts + 1
            else
                crv_ids(i, npts - 1) = all_crv_ids(i, nids - 1)
                crv_size(i) = npts
            end if
        end if
    end do

end subroutine filter_points

subroutine refine_spi_point(n, p, uk, m, q, vk, cpw, origin, pnorm, vx, vy, &
                            u0, v0, tol, u, v, up, vp, pi)
    !> Refine surface-plane intersection point.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> origin - Origin of plane.
    !> pnorm - Unit vector normal to plane.
    !> vx - Vector defining plane major axis.
    !> vy - Vector defining plane minor axis.
    !> u0 - Initial parameter on surface.
    !> v0 - Initial parameter on surface.
    !> tol - Refinement tolerance.
    !> u - Refined parameter on surface.
    !> v - Refined parameter on surface.
    !> up - Refined parameter on plane.
    !> vp - Refined parameter on plane.
    !> pi - Refined intersection point.
    
    !f2py intent(in) n, p, uk, m, q, vk, cpw, origin, pnorm, vx, vy, u0, v0, tol
    !f2py intent(out) u, v, up, vp, pi
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    
    ! Input
    integer, intent(in) :: n, m, p, q
    double precision, intent(in) :: u0, v0, tol
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    double precision, intent(in) :: origin(0:2), pnorm(0:2), vx(0:2), vy(0:2)
    
    ! Output
    double precision, intent(out) :: u, v, up, vp
    double precision, intent(out) :: pi(0:2)
    
    ! Working
    integer :: k, umult, vmult
    double precision :: du, dv, umin, umax, vmin, vmax
    double precision :: dd, dist, dp, dq, dn, dline, numer, vnorm, dpq, mag
    double precision :: p0(0:2), vec(0:2), q0(0:2), su(0:2), pq0(0:2)
    double precision :: sv(0:2), skl(0:2, 0:2, 0:2), np(0:2)
    double precision :: nn(0:2), xi(0:2), dp0(0:2), ru(0:2), rv(0:2)
    double precision :: vcross1(0:2), vcross2(0:2), vcross3(0:2)
    double precision :: temp(0:2), vx_norm(3), vy_norm(3)
    double precision :: plane_params(1, 2), temp_pnts(1, 3)
    
    umin = uk(0)
    umax = uk(n + p + 1)
    vmin = vk(0)
    vmax = vk(m + q + 1)
    u = u0
    v = v0
    k = 0
    do while ( k .lt. 100)
        ! Evaluate point on surface and derivatives.
        call rat_surface_derivs(n, p, uk, m, q, vk, cpw, u, v, 2, skl)
        p0 = skl(0, 0, :)
        su = skl(1, 0, :)
        sv = skl(0, 1, :)
        ! Project point to plane.
        vec = p0 - origin
        dd = dot_product(pnorm, vec)
        q0 = p0 - dd * pnorm
        temp = q0 - p0
        ! Check for convergence criteria.
        dist = norm(temp)
        if (dist .le. tol) then
            exit
        end if
        ! Surface and plane normals.
        np = cross(su, sv)
        vnorm = norm(np)
        if (vnorm .le. 1.0d-12) then
            exit
        end if
        np = np / vnorm
        ! Intersection of all planes.
        dp = dot_product(np, p0)
        dq = dot_product(pnorm, q0)
        nn = cross(np, pnorm)
        vnorm = norm(nn)
        if (vnorm .le. 1.0d-12) then
            exit
        end if
        nn = nn / vnorm
        pq0 = 0.50d0 * (p0 + q0)
        dn = dot_product(nn, pq0)
        vcross1 = cross(pnorm, nn)
        vcross2 = cross(nn, np)
        vcross3 = cross(np, pnorm)
        xi = (dp * vcross1 + dq * vcross2 + dn * vcross3) / &
             dot_product(vcross3, nn)
        ! Update parameters.
        dp0 = xi - p0
        ru = cross(su, np)
        rv = cross(sv, np)
        ! Check to see if current parameter is on a isoparameter of
        ! the surface. If it is and its multiplicity is equal to the
        ! degree, constrain the refinement process along the
        ! isoparameter direction.
        umult = find_mult(n, p, uk, u)
        vmult = find_mult(m, q, vk, v)
        dpq = dot_product(pnorm, p0 - q0)
        if (p .le. umult) then
            ! Adjust v only.
            du = 0.0d0
            if (dot_product(pnorm, sv) * dpq .ge. 0.d0) then
                dv = -dabs(dot_product(ru, dp0) / dot_product(ru, sv))
            else
                dv = dabs(dot_product(ru, dp0) / dot_product(ru, sv))
            
            end if
        elseif (q .le. vmult) then
            ! Adjust u only.
            dv = 0.0d0
            if (dot_product(pnorm, su) * dpq .ge. 0.0d0) then
                du = -dabs(dot_product(rv, dp0) / dot_product(rv, su))
            else
                du = dabs(dot_product(rv, dp0) / dot_product(rv, su))
            end if
        else
            du = dot_product(rv, dp0) / dot_product(rv, su)
            dv = dot_product(ru, dp0) / dot_product(ru, sv)
        end if
        u = u + du
        v = v + dv
        ! Check parameters.
        if (u .lt. umin) then
            u = umin
        elseif (u .gt. umax) then
            u = umax
        end if
        if (v .lt. vmin) then
            v = vmin
        elseif (v .gt. vmax) then
            v = vmax
        end if
        k = k + 1
    end do
    if ((k .ge. 100) .or. (dist .gt. tol)) then
        ! Attemp Nelder-Mead.
        umult = find_mult(n, p, uk, u0)
        vmult = find_mult(m, q, vk, v0)
        if (p .le. umult) then
            call refine_spi_nm(n, p, uk, m, q, vk, cpw, origin, pnorm, &
                               u0, v0, tol, 'v', u, v)
        elseif (q .le. vmult) then
            call refine_spi_nm(n, p, uk, m, q, vk, cpw, origin, pnorm, &
                               u0, v0, tol, 'u', u, v)
        else
            call refine_spi_nm(n, p, uk, m, q, vk, cpw, origin, pnorm, &
                               u0, v0, tol, 'b', u, v)
        end if
        ! Check parameters.
        if (u .lt. umin) then
            u = umin
        elseif (u .gt. umax) then
            u = umax
        end if
        if (v .lt. vmin) then
            v = vmin
        elseif (v .gt. vmax) then
            v = vmax
        end if
        call surface_point(n, p, uk, m, q, vk, cpw, u, v, p0)
        ! Invert points on plane.
        temp_pnts(1, :) = p0
        call invert_points_on_plane(1, temp_pnts, origin, vx, vy, plane_params)
        up = plane_params(1, 1)
        vp = plane_params(1, 2)
        mag = norm(vx)
        vx_norm = vx / mag
        mag = norm(vy)
        vy_norm = vy / mag
        q0 = origin + up * vx_norm + vp * vy_norm
        temp = p0 - q0
        dist = norm(temp)
        if ((dist .gt. tol) .and. (warnings)) then
            print *, "WARNING: Distance in SPI refinement exceeds tolerance. Distance=", dist
        end if
    end if
    ! Invert points on plane.
    temp_pnts(1, :) = p0
    call invert_points_on_plane(1, temp_pnts, origin, vx, vy, plane_params)
    up = plane_params(1, 1)
    vp = plane_params(1, 2)
    
    pi = 0.50d0 * (p0 + q0)
    
end subroutine refine_spi_point

subroutine refine_spi_nm(n, p, uk, m, q, vk, cpw, origin, pnorm, &
                         u0, v0, tol, d, u, v)
    !> Refine SPI point using Nelder-Mead optimization.
    
    ! Input
    character, intent(in) :: d
    integer, intent(in) :: n, m, p, q
    double precision, intent(in) :: u0, v0, tol
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    double precision, intent(in) :: origin(3), pnorm(3)

    ! Output
    double precision, intent(out) :: u, v
    
    ! Working
    logical :: success
    double precision, allocatable :: xin(:)
    
    if (d .eq. 'u') then
        allocate(xin(1))
        xin(1) = u0
        v = v0
        call simplex(obj, 1, xin, tol, success)
        u = xin(1)
    elseif (d .eq. 'v') then
        allocate(xin(1))
        xin(1) = v0
        u = u0
        call simplex(obj, 1, xin, tol, success)
        v = xin(1)
    else
        allocate(xin(2))
        xin(:) = (/ u0, v0 /)
        call simplex(obj, 2, xin, tol, success)
        u = xin(1)
        v = xin(2)
    end if
    
    return

    contains
    
    subroutine obj(nx, x, fx)
        integer, intent(in) :: nx
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: fx
        
        ! Working
        double precision :: dd
        double precision :: p0(3), q0(3), vp(3), dpq(3)
        
        if (d .eq. 'u') then
            call surface_point(n, p, uk, m, q, vk, cpw, x(1), v, p0)
        elseif (d .eq. 'v') then
            call surface_point(n, p, uk, m, q, vk, cpw, u, x(1), p0)
        else
            call surface_point(n, p, uk, m, q, vk, cpw, x(1), x(2), p0)
        end if
        vp = p0 - origin
        dd = dot_product(vp, pnorm)
        q0 = p0 - dd * pnorm
        dpq = p0 - q0
        fx = norm(dpq)
    end subroutine obj

end subroutine refine_spi_nm

subroutine refine_ssi_point(n1, p1, uk1, m1, q1, vk1, cpw1, &
                            n2, p2, uk2, m2, q2, vk2, cpw2, &
                            u01, v01, u02, v02, tol, &
                            u1, v1, u2, v2, pi)
    !> Refine surface-surface intersection point.
    !> n1 - Number of control points - 1 in u-direction for surface 1.
    !> p1 - Degree in u-direction for surface 1.
    !> uk1 - Knot vector in u-direction for surface 1.
    !> m1 - Number of control points - 1 in v-direction for surface 1.
    !> q1 - Degree in v-direction for surface 1.
    !> vk1 - Knot vector in v-direction for surface 1.
    !> cpw1 - Control points for surface 1.
    !> n2 - Number of control points - 1 in u-direction for surface 2.
    !> p2 - Degree in u-direction for surface 2.
    !> uk2 - Knot vector in u-direction for surface 2.
    !> m2 - Number of control points - 1 in v-direction for surface 2.
    !> q2 - Degree in v-direction for surface 2.
    !> vk2 - Knot vector in v-direction for surface 2.
    !> cpw2 - Control points for surface 2.
    !> u01 - Initial parameter for surface 1.
    !> v01 - Initial parameter for surface 1.
    !> u02 - Initial parameter for surface 2.
    !> v02 - Initial parameter for surface 2.
    !> tol - Refinement tolerance.
    !> u1 - Refined parameter for surface 1.
    !> v1 - Refined parameter for surface 1.
    !> u2 - Refined parameter for surface 2.
    !> v2 - Refined parameter for surface 2.
    !> pi - Refined intersection point.
    
    !f2py intent(in) n1, p1, uk1, m1, q1, vk1, cpw1
    !f2py intent(in) n2, p2, uk2, m2, q2, vk2, cpw2
    !f2py intent(in) u01, v01, u02, v02, tol
    !f2py intent(out) u1, v1, u2, v2, pi
    !f2py depend(n1, p1) uk1
    !f2py depend(m1, q1) vk1
    !f2py depend(n1, m1) cpw1
    !f2py depend(n2, p2) uk2
    !f2py depend(m2, q2) vk2
    !f2py depend(n2, m2) cpw2
    
    ! Input
    integer, intent(in) :: n1, m1, p1, q1, n2, m2, p2, q2
    double precision, intent(in) :: u01, v01, u02, v02, tol
    double precision, intent(in) :: uk1(0:n1 + p1 + 1), vk1(0:m1 + q1 + 1)
    double precision, intent(in) :: cpw1(0:n1, 0:m1, 0:3)
    double precision, intent(in) :: uk2(0:n2 + p2 + 1), vk2(0:m2 + q2 + 1)
    double precision, intent(in) :: cpw2(0:n2, 0:m2, 0:3)
    
    ! Output
    double precision, intent(out) :: u1, v1, u2, v2
    double precision, intent(out) :: pi(0:2)
    
    ! Working
    integer :: k
    double precision :: vnorm, dp, dq, dn, d0
    double precision :: umin1, umin2, umax1, umax2
    double precision :: vmin1, vmin2, vmax1, vmax2
    double precision :: skl(0:2, 0:2, 0:2)
    double precision :: p0(0:2), su1(0:2), sv1(0:2)
    double precision :: q0(0:2), su2(0:2), sv2(0:2)
    double precision :: dpq(0:2), np(0:2), nq(0:2), nn(0:2), pq0(0:2)
    double precision :: vcross1(0:2), vcross2(0:2), vcross3(0:2)
    double precision :: xi(0:2), dp0(0:2), dq0(0:2)
    double precision :: ru1(0:2), rv1(0:2), ru2(0:2), rv2(0:2)
    
    ! Initial values.
    u1 = u01
    v1 = v01
    u2 = u02
    v2 = v02
    k = 0
    umin1 = uk1(0)
    umax1 = uk1(n1 + p1 + 1)
    umin2 = uk2(0)
    umax2 = uk2(n2 + p2 + 1)
    vmin1 = vk1(0)
    vmax1 = vk1(m1 + q1 + 1)
    vmin2 = vk2(0)
    vmax2 = vk2(m2 + q2 + 1)
    call rat_surface_derivs(n1, p1, uk1, m1, q1, vk1, cpw1, u1, v1, 2, skl)
    p0 = skl(0, 0, :)
    su1 = skl(1, 0, :)
    sv1 = skl(0, 1, :)
    call rat_surface_derivs(n2, p2, uk2, m2, q2, vk2, cpw2, u2, v2, 2, skl)
    q0 = skl(0, 0, :)
    su2 = skl(1, 0, :)
    sv2 = skl(0, 1, :)
    dpq = p0 - q0
    d0 = norm(dpq)
    do while (k .lt. 100)
        ! Check for convergence criteria.
        if (d0 .le. tol) then
            exit
        end if
        ! Surface unit normals.
        np = cross(su1, sv1)
        vnorm = norm(np)
        if (vnorm .le. 1.0d-12) then
            exit
        end if
        np = np / vnorm
        nq = cross(su2, sv2)
        vnorm = norm(nq)
        if (vnorm .le. 1.0d-12) then
            exit
        end if
        nq = nq / vnorm
        ! Intersection of all three planes.
        dp = dot_product(np, p0)
        dq = dot_product(nq, q0)
        nn = cross(np, nq)
        vnorm = norm(nn)
        if (vnorm .le. 1.0d-12) then
            exit
        end if
        nn = nn / vnorm
        pq0 = 0.50d0 * (p0 + q0)
        dn = dot_product(nn, pq0)
        vcross1 = cross(nq, nn)
        vcross2 = cross(nn, np)
        vcross3 = cross(np, nq)
        xi = (dp * vcross1 + dq * vcross2 + dn * vcross3) / &
             dot_product(vcross3, nn)
        ! New parameters for surface 1.
        dp0 = xi - p0
        ru1 = cross(su1, np)
        rv1 = cross(sv1, np)
        u1 = u1 + dot_product(rv1, dp0) / dot_product(rv1, su1)
        v1 = v1 + dot_product(ru1, dp0) / dot_product(ru1, sv1)
        ! Check parameters.
        if (u1 .lt. umin1) then
            u1 = umin1
        elseif (u1 .gt. umax1) then
            u1 = umax1
        end if
        if (v1 .lt. vmin1) then
            v1 = vmin1
        elseif (v1 .gt. vmax1) then
            v1 = vmax1
        end if
        ! New parameters for surface 2.
        dq0 = xi - q0
        ru2 = cross(su2, nq)
        rv2 = cross(sv2, nq)
        u2 = u2 + dot_product(rv2, dq0) / dot_product(rv2, su2)
        v2 = v2 + dot_product(ru2, dq0) / dot_product(ru2, sv2)
        ! Check parameters.
        if (u2 .lt. umin2) then
            u2 = umin2
        elseif (u2 .gt. umax2) then
            u2 = umax2
        end if
        if (v2 .lt. vmin2) then
            v2 = vmin2
        elseif (v2 .gt. vmax2) then
            v2 = vmax2
        end if
        ! New location.
        call rat_surface_derivs(n1, p1, uk1, m1, q1, vk1, cpw1, u1, v1, 2, skl)
        p0 = skl(0, 0, :)
        su1 = skl(1, 0, :)
        sv1 = skl(0, 1, :)
        call rat_surface_derivs(n2, p2, uk2, m2, q2, vk2, cpw2, u2, v2, 2, skl)
        q0 = skl(0, 0, :)
        su2 = skl(1, 0, :)
        sv2 = skl(0, 1, :)
        dpq = p0 - q0
        d0 = norm(dpq)
        k = k + 1
    end do
    if ((k .ge. 100) .or. (d0 .gt. tol)) then
        ! Attempt Nelder-Mead.
        call refine_ssi_nm(n1, p1, uk1, m1, q1, vk1, cpw1, &
                           n2, p2, uk2, m2, q2, vk2, cpw2, &
                           u01, v01, u02, v02, tol, u1, v1, u2, v2)
        ! Check parameters.
        if (u1 .lt. umin1) then
            u1 = umin1
        elseif (u1 .gt. umax1) then
            u1 = umax1
        end if
        if (v1 .lt. vmin1) then
            v1 = vmin1
        elseif (v1 .gt. vmax1) then
            v1 = vmax1
        end if
        if (u2 .lt. umin2) then
            u2 = umin2
        elseif (u2 .gt. umax2) then
            u2 = umax2
        end if
        if (v2 .lt. vmin2) then
            v2 = vmin2
        elseif (v2 .gt. vmax2) then
            v2 = vmax2
        end if
        call surface_point(n1, p1, uk1, m1, q1, vk1, cpw1, u1, v1, p0)
        call surface_point(n2, p2, uk2, m2, q2, vk2, cpw2, u2, v2, q0)
        dpq = p0 - q0
        d0 = norm(dpq)
        if ((d0 .gt. tol) .and. (warnings)) then
            print *, "WARNING: Distance in SSI refinement exceeds tolerance. Distance=", d0
        end if
    end if
    
    pi = 0.50d0 * (p0 + q0)

    end subroutine refine_ssi_point

subroutine refine_ssi_nm(n1, p1, uk1, m1, q1, vk1, cpw1, &
                         n2, p2, uk2, m2, q2, vk2, cpw2, &
                         u10, v10, u20, v20, tol, u1, v1, u2, v2)
    !> Refine SSI point using Nelder-Mead optimization.
    
    ! Input
    integer, intent(in) :: n1, m1, p1, q1, n2, m2, p2, q2
    double precision, intent(in) :: u10, v10, u20, v20, tol
    double precision, intent(in) :: uk1(0:n1 + p1 + 1), vk1(0:m1 + q1 + 1)
    double precision, intent(in) :: uk2(0:n2 + p2 + 1), vk2(0:m2 + q2 + 1)
    double precision, intent(in) :: cpw1(0:n1, 0:m1, 0:3)
    double precision, intent(in) :: cpw2(0:n2, 0:m2, 0:3)

    ! Output
    double precision, intent(out) :: u1, v1, u2, v2
    
    ! Working
    logical :: success
    ! double precision :: au1, bu1, av1, bv1
    ! double precision :: au2, bu2, av2, bv2
    double precision :: xin(4)
    
    ! au1 = uk1(0)
    ! bu1 = uk1(n1 + p1 + 1)
    ! av1 = vk1(0)
    ! bv1 = vk1(m1 + q1 + 1)
    
    ! au2 = uk2(0)
    ! bu2 = uk2(n2 + p2 + 1)
    ! av2 = vk2(0)
    ! bv2 = vk2(m2 + q2 + 1)
    
    xin(:) = (/ u10, v10, u20, v20 /)
    call simplex(obj, 4, xin, tol, success)
    u1 = xin(1)
    v1 = xin(2)
    u2 = xin(3)
    v2 = xin(4)
    
    return

    contains
    
    subroutine obj(nx, x, fx)
        integer, intent(in) :: nx
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: fx
        
        ! Working
        ! double precision :: factor
        double precision :: p0(3), q0(3), dpq(3)
        
        ! Penalize objective function if outside domain.
        ! factor = 1.0d0
        ! if ((x(1) .lt. au1) .or. (x(1) .gt. bu1)) then
        !     factor = 1000.0d0
        ! elseif ((x(2) .lt. av1) .or. (x(2) .gt. bv1)) then
        !     factor = 1000.0d0
        ! elseif ((x(3) .lt. au2) .or. (x(3) .gt. bu2)) then
        !     factor = 1000.0d0
        ! elseif ((x(4) .lt. av2) .or. (x(4) .gt. bv2)) then
        !     factor = 1000.0d0
        ! end if
        
        call surface_point(n1, p1, uk1, m1, q1, vk1, cpw1, x(1), x(2), p0)
        call surface_point(n2, p2, uk2, m2, q2, vk2, cpw2, x(3), x(4), q0)
        dpq = p0 - q0
        fx = norm(dpq)
    end subroutine obj

end subroutine refine_ssi_nm


end module intersect_surface