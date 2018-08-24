module tessellate
use math
use generic_list
use geom_results
use geom_utils
use evaluate
use divide
implicit none

private :: crv_data, crv_ptr, srf_data
private :: find_neighbors, find_edge_params

! Tessellation data.
type crv_data
    double precision :: u
end type crv_data

type crv_ptr
    type(crv_data), pointer :: c
end type crv_ptr

type srf_data
    integer :: sid = 0
    integer :: position = 0
    logical :: is_cand = .false.
    logical :: has_child = .false.
    double precision :: u0 = 0.0d0
    double precision :: u1 = 0.0d0
    double precision :: v0 = 0.0d0
    double precision :: v1 = 0.0d0
    type(srf_data), pointer :: parent => null()
    type(srf_data), pointer :: ne => null()
    type(srf_data), pointer :: se => null()
    type(srf_data), pointer :: sw => null()
    type(srf_data), pointer :: nw => null()
    type(srf_data), pointer :: n => null()
    type(srf_data), pointer :: e => null()
    type(srf_data), pointer :: s => null()
    type(srf_data), pointer :: w => null()
end type srf_data

type srf_ptr
    type(srf_data), pointer :: ci
end type srf_ptr

contains

subroutine adaptive_curve_tessellate(n, p, uk, cpw, tol, flag)
    !> Adaptive curve tessellation.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> tol - Tolerance for checking curve flatness.
    !> flag - Status flag.
    
    !f2py intent(in) n, p, uk, cpw, tol
    !f2py intent(out) flag
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: tol
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    integer :: i
    double precision :: u
    double precision :: pi(0:2)
    double precision :: cp(0:n, 0:2), w(0:n)
    
    ! Generic list for storing results.
    type(list_node_t), pointer :: tess_list => null()
    type(crv_ptr) :: ptr
    
    ! Use control points if linear.
    ! if (p .eq. 1) then
    !     call reset_results()
    !     nverts_ = n + 1
    !     allocate(verts_(0:nverts_ - 1, 0:2))
    !     allocate(params1d_(0:nverts_ - 1))
    !     call dehomogenize_array1d(n, cpw, cp, w)
    !     verts_ = cp
    !     params1d_ = uk(1:n + 1)
    !     return
    ! end if

    ! Initialize
    call list_init(tess_list)
    call reset_results()
    flag = 0
    
    ! Step 1: Method for recursive subdivision are below.

    ! Step 2: Use recursvie subdivision to generate flat curves.
    call subdivide(n, p, uk, cpw)
    
    ! Store last point.
    allocate(ptr%c)
    ptr%c%u = uk(n + 1)
    call list_insert(tess_list, transfer(ptr, list_data))
    nverts_ = nverts_ + 1
    
    ! Allocate result arrays.
    allocate(verts_(0:nverts_ - 1, 0:2))
    allocate(params1d_(0:nverts_ - 1))
    verts_(:, :) = 0.0d0
    params1d_(:) = 0.0d0
    
    ! Step 3: Generate vertices and edges. Work backwards due to linked list.
    do i = nverts_ - 1, 0, -1
        ! Get next
        tess_list = list_next(tess_list)
        ptr = transfer(list_get(tess_list), ptr)
        u = ptr%c%u
        call curve_point(n, p, uk, cpw, u, pi)
        verts_(i, :) = pi
        params1d_(i) = u
    end do

    ! Return.
    flag = 1

    ! Internal subroutines.
    contains
    
    recursive subroutine subdivide(n, p, uk, cpw)
        ! Recursive subdivision.
        ! Input
        integer, intent(in) :: n, p
        double precision, intent(in) :: uk(0:n + p + 1)
        double precision, intent(in) :: cpw(0:n, 0:3)
        
        ! Working
        logical :: flat
        integer :: n1, n2, ukn, mid
        double precision :: ui
        double precision :: cp(0:n, 0:2), w(0:n)
        double precision, allocatable :: uk1(:), uk2(:)
        double precision, allocatable :: qw1(:, :)
        double precision, allocatable :: qw2(:, :)

        ! Check flatness.
        call dehomogenize_array1d(n, cpw, cp, w)
        call is_curve_flat(n, cp, tol, flat)
        if (flat .eqv. .false.) then
            ! Split curve.
            if (p .eq. 1) then
                ukn = size(uk, 1)
                mid = ceiling((dble(ukn) - 1.0d0) / 2.0d0)
                ui = uk(mid)
            else
                ui = 0.50d0 * (uk(p) + uk(n + 1))
            end if
            call split_curve(n, p, uk, cpw, ui, uk1, qw1, uk2, qw2)
            n1 = size(qw1, 1) - 1
            n2 = size(qw2, 1) - 1
            ! Subdivide
            call subdivide(n1, p, uk1, qw1)
            call subdivide(n2, p, uk2, qw2)
        else
            ! Store results.
            allocate(ptr%c)
            ptr%c%u = uk(p)
            call list_insert(tess_list, transfer(ptr, list_data))
            nverts_ = nverts_ + 1
        end if
    end subroutine subdivide

end subroutine adaptive_curve_tessellate

subroutine adaptive_surface_tessellate(n, p, uk, m, q, vk, cpw, tol, flag)
    !> Adaptive surface tessellation.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> tol - Tolerance for checking surface flatness.
    !> flag - Status flag.
    
    !f2py intent(in) n, p, uk, m, q, vk, cpw, tol
    !f2py intent(out) flag
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    
    ! Input
    integer, intent(in) :: n, p, m, q
    double precision, intent(in) :: tol
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    integer, intent(out) :: flag
    
    ! Working
    integer :: i, j, ncells, nt, nflat
    integer, allocatable :: children(:, :)
    integer, allocatable :: parent(:)
    integer, allocatable :: acells(:, :)
    integer, allocatable :: position(:)
    integer, allocatable :: candidates(:)
    double precision :: cp(0:n, 0:m, 0:2), w(0:n, 0:m)
    double precision :: uv0(0:1), uv1(0:1), uv2(0:1)
    double precision :: t0(0:2), t1(0:2), t2(0:2)
    double precision :: tri(0:100, 0:2, 0:1)
    double precision, allocatable :: cell_params(:, :, :)
    
    ! Generic list for storing flat cells.
    type(list_node_t), pointer :: cell_list => null()
    type(srf_ptr) :: cell_ptr
    
    ! Tessellate control net if linear in both directions.
    if ((p .eq. 1) .and. (q .eq. 1)) then
        call reset_results()
        nverts_ = (n + 1) * (m + 1)
        ntri_ = 2 * n * m
        allocate(verts_(0:nverts_ - 1, 0:2))
        allocate(triangles_(0:ntri_ - 1, 0:2))
        verts_(:, :) = 0.0d0
        triangles_(:, :) = 0
        call dehomogenize_array2d(n, m, cpw, cp, w)
        call tessellate_cp_net(n, m, cp, verts_, triangles_)
        return
    end if
    
    ! Initialize.
    ncells = 0
    nflat = 0
    call list_init(cell_list)
    
    ! Step 1: Method for recursive subdivision are below.
    
    ! Step 2: Use recursvie subdivision to generate flat surfaces.
    allocate(cell_ptr%ci)
    call subdivide(n, p, uk, m, q, vk, cpw, cell_ptr)
    
    ! Build arrays to tessellate the cells.
    call build_tess_data(ncells, cell_list, children, parent, acells, &
                         position, candidates, cell_params)
    
    ! Step 3: Tessellate flat cells. Allocate arrays based on number
    ! of flat cells. This may need adjusted in the future based on experience.
    ! Vertices are not equivalenced.
    call reset_results()
    allocate(verts_(0:8 * nflat, 0:2), triangles_(0:4 * nflat, 0:2))
    verts_(:, :) = 0.0d0
    triangles_(:, :) = 0
    tri(:, :, :) = 0.0d0
    do i = 0, ncells
        if (candidates(i) .eq. 1) then
            ! Get parameters of tessellated cell.
            call tessellate_cell(i, ncells + 1, children, acells, position, &
                                 parent, cell_params, nt, tri)
            ! Determine and save triangle points.
            do j = 0, nt - 1, 1
                uv0 = tri(j, 0, :)
                uv1 = tri(j, 1, :)
                uv2 = tri(j, 2, :)
                call surface_point(n, p, uk, m, q, vk, cpw, &
                                   uv0(0), uv0(1), t0)
                call surface_point(n, p, uk, m, q, vk, cpw, &
                                   uv1(0), uv1(1), t1)
                call surface_point(n, p, uk, m, q, vk, cpw, &
                                   uv2(0), uv2(1), t2)         
                verts_(nverts_, :) = t0
                verts_(nverts_ + 1, :) = t1
                verts_(nverts_ + 2, :) = t2
                triangles_(ntri_, :) = (/ nverts_, nverts_ + 1, nverts_ + 2 /)
                nverts_ = nverts_ + 3
                ntri_ = ntri_ + 1
            end do
        end if
    end do
    
    ! Normal completion.
    call list_free(cell_list)
    flag = 1
    
    ! ! Step 4: Equivalence vertices and build connectivity array.
    ! allocate(point_use(0:ptotal - 1))
    ! point_use(:) = 0
    ! verts(0, :) = points(0, :)
    ! nverts = 1
    ! gtol2 = gtol * gtol
    ! do i = 1, ptotal - 1
        ! pi = points(i, :)
        ! unique = .true.
        ! do j = 0, nverts - 1
            ! vi = verts(j, :)
            ! d2 = (vi(0) - pi(0)) * (vi(0) - pi(0)) + &
                 ! (vi(1) - pi(1)) * (vi(1) - pi(1)) + &
                 ! (vi(2) - pi(2)) * (vi(2) - pi(2))
            ! if (d2 .le. gtol2) then
                ! unique = .false.
                ! point_use(i) = j
                ! exit
            ! end if
        ! end do
        ! if (unique) then    
            ! verts(nverts, :) = points(i, :)
            ! point_use(i) = nverts
            ! nverts = nverts + 1
            ! if (nverts .gt. nmax - 1) then
                ! success = .false.
                ! resize = resize + 1
                ! nverts = nmax - 1
            ! end if
        ! end if
    ! end do
    ! ! Check for failure.
    ! if (success .eqv. .false.) then
        ! flag = -2
        ! resize = resize + nmax
        ! return
    ! end if
    
    ! ! Build connectivity matrix.    
    ! ntri = 0
    ! do i = 0, tri_total - 1
        ! triangles(ntri, 0) = point_use(conn(i, 0))
        ! triangles(ntri, 1) = point_use(conn(i, 1))
        ! triangles(ntri, 2) = point_use(conn(i, 2))
        ! ntri = ntri + 1
        ! if (ntri .gt. nmax - 1) then
            ! success = .false.
            ! resize = resize + 1
            ! ntri = nmax - 1
        ! end if
    ! end do
    ! ! Check for failure.
    ! if (success .eqv. .false.) then
        ! flag = -3
        ! resize = resize + nmax
        ! return
    ! end if
    
    ! Internal subroutines.
    contains
    
    recursive subroutine subdivide(n, p, uk, m, q, vk, cpw, cell_ptr)
        ! Recursive subdivision.
        ! Input
        integer, intent(in) :: n, p, m, q
        double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
        double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
        
        type(srf_ptr), intent(inout) :: cell_ptr
        
        ! Working
        logical :: flat
        integer :: mid
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
        
        ! Insert cell
        call list_insert(cell_list, transfer(cell_ptr, list_data))
        
        ! Store surface parameters
        cell_ptr%ci%u0 = uk(p)
        cell_ptr%ci%u1 = uk(n + 1)
        cell_ptr%ci%v0 = vk(q)
        cell_ptr%ci%v1 = vk(m + 1)
        
        ! Check flatness
        call dehomogenize_array2d(n, m, cpw, cp, w)
        call is_surface_flat(n, m, cp, tol, flat)
        if (flat .and. (ncells .gt. 1)) then
            ! Save candidate surfaces.
            cell_ptr%ci%is_cand = .true.
            nflat = nflat + 1
        else
            ! Split surface
            ukn = size(uk, 1)
            vkn = size(vk, 1)
            if ((q .eq. 1) .and. (ukn .eq. 4)) then
                ! Split at v-parameter only if bilinear surface and
                ! u-direction has no interior knots. Split at middle
                ! interior knot.
                mid = ceiling((dble(vkn) - 1.0d0) / 2.0d0)
                uv = vk(mid)
                call split_surface(n, p, uk, m, q, vk, cpw, uv, 'v', &
                                   uk1, vk1, qw1, uk2, vk2, qw2)
                n1 = size(qw1, 1) - 1
                m1 = size(qw1, 2) - 1
                n2 = size(qw2, 1) - 1
                m2 = size(qw2, 2) - 1
                
                ! Allocate new cells.
                allocate(sw_ptr%ci, nw_ptr%ci)
                
                ! Cell 1 (SW).
                sw_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%has_child = .true.
                cell_ptr%ci%sw => sw_ptr%ci
                sw_ptr%ci%parent => cell_ptr%ci
                sw_ptr%ci%position = 1
                
                ! Cell 2 (NW).
                nw_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%nw => nw_ptr%ci
                nw_ptr%ci%parent => cell_ptr%ci
                nw_ptr%ci%position = 2
                
                ! Adjacent cells for cell 1 (SW).
                sw_ptr%ci%n => nw_ptr%ci
                sw_ptr%ci%e => cell_ptr%ci%e
                sw_ptr%ci%s => cell_ptr%ci%s
                sw_ptr%ci%w => cell_ptr%ci%w
                
                ! Adjacent cells for cell 2 (NW).
                nw_ptr%ci%n => cell_ptr%ci%n
                nw_ptr%ci%e => cell_ptr%ci%e
                nw_ptr%ci%s => sw_ptr%ci
                nw_ptr%ci%w => cell_ptr%ci%w
                
                ! Subdivide.
                call subdivide(n1, p, uk1, m1, q, vk1, qw1, sw_ptr)
                call subdivide(n2, p, uk2, m2, q, vk2, qw2, nw_ptr)
            elseif ((p .eq. 1) .and. (vkn .eq. 4)) then
                ! Split at u-parameter only if bilinear surface and
                ! v-direction has no interior knots. Split at middle
                ! interior knot.
                mid = ceiling((dble(ukn) - 1.0d0) / 2.0d0)
                uv = uk(mid)
                call split_surface(n, p, uk, m, q, vk, cpw, uv, 'u', &
                                   uk1, vk1, qw1, uk3, vk3, qw3)
                n1 = size(qw1, 1) - 1
                m1 = size(qw1, 2) - 1
                n3 = size(qw3, 1) - 1
                m3 = size(qw3, 2) - 1
                
                ! Allocate new cells.
                allocate(sw_ptr%ci, se_ptr%ci)
                
                ! Cell 1 (SW).
                sw_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%has_child = .true.
                cell_ptr%ci%sw => sw_ptr%ci
                sw_ptr%ci%parent => cell_ptr%ci
                sw_ptr%ci%position = 1
                
                ! Cell 3 (SE).
                se_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%se => se_ptr%ci
                se_ptr%ci%parent => cell_ptr%ci
                se_ptr%ci%position = 3
               
               ! Adjacent cells for cell 1 (SW).
                sw_ptr%ci%n => cell_ptr%ci%n
                sw_ptr%ci%e => se_ptr%ci
                sw_ptr%ci%s => cell_ptr%ci%s
                sw_ptr%ci%w => cell_ptr%ci%w
                
                ! Adjacent cells for cell 3 (SE).
                se_ptr%ci%n => cell_ptr%ci%n
                se_ptr%ci%e => cell_ptr%ci%e
                se_ptr%ci%s => cell_ptr%ci%s
                se_ptr%ci%w => sw_ptr%ci
               
                call subdivide(n1, p, uk1, m1, q, vk1, qw1, sw_ptr)
                call subdivide(n3, p, uk3, m3, q, vk3, qw3, se_ptr)
            else
                ! Split surface into four regions at midpoint
                ui = 0.50d0 * (cell_ptr%ci%u0 + cell_ptr%ci%u1)
                vi = 0.50d0 * (cell_ptr%ci%v0 + cell_ptr%ci%v1)
                
                ! Split at middle interior knot if linear
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
                
                ! Allocate new cells
                allocate(ne_ptr%ci, se_ptr%ci, sw_ptr%ci, nw_ptr%ci)
                
                ! Cell 1 (SW)
                sw_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%has_child = .true.
                cell_ptr%ci%sw => sw_ptr%ci
                sw_ptr%ci%parent => cell_ptr%ci
                sw_ptr%ci%position = 1
                
                ! Cell 2 (NW)
                nw_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%nw => nw_ptr%ci
                nw_ptr%ci%parent => cell_ptr%ci
                nw_ptr%ci%position = 2
                
                ! Cell 3 (SE)
                se_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%se => se_ptr%ci
                se_ptr%ci%parent => cell_ptr%ci
                se_ptr%ci%position = 3
               
                ! Cell 4 (NE)
                ne_ptr%ci%sid = ncells + 1
                ncells = ncells + 1
                cell_ptr%ci%ne => ne_ptr%ci
                ne_ptr%ci%parent => cell_ptr%ci
                ne_ptr%ci%position = 4
                
                ! Adjacent cells for cell 1 (SW)
                sw_ptr%ci%n => nw_ptr%ci
                sw_ptr%ci%e => se_ptr%ci
                sw_ptr%ci%s => cell_ptr%ci%s
                sw_ptr%ci%w => cell_ptr%ci%w
                
                ! Adjacent cells for cell 2 (NW)
                nw_ptr%ci%n => cell_ptr%ci%n
                nw_ptr%ci%e => ne_ptr%ci
                nw_ptr%ci%s => sw_ptr%ci
                nw_ptr%ci%w => cell_ptr%ci%w
                
                ! Adjacent cells for cell 3 (SE)
                se_ptr%ci%n => ne_ptr%ci
                se_ptr%ci%e => cell_ptr%ci%e
                se_ptr%ci%s => cell_ptr%ci%s
                se_ptr%ci%w => sw_ptr%ci
               
                ! Adjacent cells for cell 4 (NE)
                ne_ptr%ci%n => cell_ptr%ci%n
                ne_ptr%ci%e => cell_ptr%ci%e
                ne_ptr%ci%s => se_ptr%ci
                ne_ptr%ci%w => nw_ptr%ci
                
                ! Subdivide
                call subdivide(n1, p, uk1, m1, q, vk1, qw1, sw_ptr)
                call subdivide(n2, p, uk2, m2, q, vk2, qw2, nw_ptr)
                call subdivide(n3, p, uk3, m3, q, vk3, qw3, se_ptr)
                call subdivide(n4, p, uk4, m4, q, vk4, qw4, ne_ptr)
            end if
        end if
        
    end subroutine subdivide
    
end subroutine adaptive_surface_tessellate

subroutine build_tess_data(ncells, cell_list, children, parent, acells, &
                           position, candidates, cell_params)
    !> Build arrays for tessellation methods given linked list of
    !> surface cells.
    
    ! Input
    integer, intent(in) :: ncells
    type(list_node_t), pointer, intent(in) :: cell_list
    
    ! Output
    integer, allocatable, intent(out) :: children(:, :)
    integer, allocatable, intent(out) :: parent(:)
    integer, allocatable, intent(out) :: acells(:, :)
    integer, allocatable, intent(out) :: position(:)
    integer, allocatable, intent(out) :: candidates(:)
    double precision, allocatable, intent(out) :: cell_params(:, :, :)
    
    ! Working
    integer :: i, sid
    type(srf_ptr) :: cell_ptr
    
    ! Allocate
    if (allocated(children)) then
        deallocate(children)
    end if
    if (allocated(parent)) then
        deallocate(parent)
    end if
    if (allocated(acells)) then
        deallocate(acells)
    end if
    if (allocated(position)) then
        deallocate(position)
    end if
    if (allocated(candidates)) then
        deallocate(candidates)
    end if
    if (allocated(cell_params)) then
        deallocate(cell_params)
    end if

    allocate(children(0:ncells, 0:4))
    allocate(parent(0:ncells))
    allocate(acells(0:ncells, 0:4))
    allocate(position(0:ncells))
    allocate(candidates(0:ncells))
    allocate(cell_params(0:ncells, 0:3, 0:1))
    
    children(:, :) = 0
    parent(:) = 0
    acells(:, :) = 0
    position(:) = 0
    candidates(:) = 0
    cell_params(:, :, :) = 0.0d0
    
    do i = 0, ncells
        ! Get next cell.
        cell_list = list_next(cell_list)
        cell_ptr = transfer(list_get(cell_list), cell_ptr)
        sid = cell_ptr%ci%sid
        
        ! Candidate cell.
        if (cell_ptr%ci%is_cand) then
            candidates(sid) = 1
        end if
        
        ! Parameters.
        cell_params(sid, 0, :) = (/ cell_ptr%ci%u0, cell_ptr%ci%v0 /)
        cell_params(sid, 1, :) = (/ cell_ptr%ci%u1, cell_ptr%ci%v0 /)
        cell_params(sid, 2, :) = (/ cell_ptr%ci%u1, cell_ptr%ci%v1 /)
        cell_params(sid, 3, :) = (/ cell_ptr%ci%u0, cell_ptr%ci%v1 /)
        
        ! Parent
        if (associated(cell_ptr%ci%parent)) then
            parent(sid) = cell_ptr%ci%parent%sid
        end if
        
        ! Position
        if (cell_ptr%ci%position .gt. 0) then
            position(sid) = cell_ptr%ci%position
        end if
        
        ! Children
        if (associated(cell_ptr%ci%sw)) then
            children(sid, 0) = 1
            children(sid, 1) = cell_ptr%ci%sw%sid
        end if
        
        if (associated(cell_ptr%ci%nw)) then
            children(sid, 0) = 1
            children(sid, 2) = cell_ptr%ci%nw%sid
        end if
        
        if (associated(cell_ptr%ci%se)) then
            children(sid, 0) = 1
            children(sid, 3) = cell_ptr%ci%se%sid
        end if
        
        if (associated(cell_ptr%ci%ne)) then
            children(sid, 0) = 1
            children(sid, 4) = cell_ptr%ci%ne%sid
        end if
        
        ! Adjacent cells.
        if (associated(cell_ptr%ci%n)) then
            acells(sid, 0) = 4
            acells(sid, 4) = cell_ptr%ci%n%sid
        end if
        if (associated(cell_ptr%ci%e)) then
            acells(sid, 0) = 4
            acells(sid, 3) = cell_ptr%ci%e%sid
        end if
        if (associated(cell_ptr%ci%s)) then
            acells(sid, 0) = 4
            acells(sid, 1) = cell_ptr%ci%s%sid
        end if
        if (associated(cell_ptr%ci%w)) then
            acells(sid, 0) = 4
            acells(sid, 2) = cell_ptr%ci%w%sid
        end if
    end do
    
end subroutine build_tess_data

subroutine tessellate_cell(csn, n, children, acells, position, &
                           parent, cell_params, ntri, triangles)
    !> Tessellate a cell.
    !> csn - Cell number.
    !> n - Maximum array size.
    !> children - Children of each cell.
    !> acells - Adjacency array for each cell.
    !> position - Position of each cell.
    !> parent - Parent of each cell.
    !> cell_params - Corner parameters of each cell.
    !> ntri - Number of triangles.
    !> triangles - Parameters of each vertex of each triangle.
    
    ! Input
    integer, intent(in) :: csn, n
    integer, intent(in) :: children(0:n - 1, 0:4), parent(0:n - 1)
    integer, intent(in) :: acells(0:n - 1, 0:4), position(0:n - 1)
    double precision, intent(in) :: cell_params(0:n - 1, 0:3, 0:1)
    
    ! Output
    integer, intent(out) :: ntri
    double precision, intent(out) :: triangles(0:100, 0:2, 0:1)
    
    ! Working
    logical :: simple
    integer :: i, nprms_i, nadj, esn, j, cid, ndim
    integer :: adj_cells(0:100), nprms(0:100), edge_order(0:3)
    double precision :: edge_prms(0:3, 0:100, 0:1)
    double precision :: prms_i(0:100, 0:1)
    double precision ::all_params(0:100, 0:1)
    double precision :: uv0(0:1), uv1(0:1), uvc(0:1)
    double precision :: eprms(0:3, 0:1)
    
    simple = .false.
    nprms(:) = 0
    edge_prms(:, :, :) = 0.0d0
    all_params(:, :) = 0.0d0
    ntri = 0
    triangles(:, :, :) = 0.0d0
    
    ! Determine number of interior points and parameters on each edge.
    edge_order(:) = (/ 1, 3, 4, 2 /)
    do i = 0, 3
        esn = edge_order(i)
        call find_neighbors(csn, esn, n, children, acells, &
                            position, parent, nadj, adj_cells)
        call find_edge_params(esn, nadj, n, adj_cells, cell_params, &
                              nprms_i, prms_i)
        nprms(i) = nprms_i
        do j = 0, nprms_i - 1
            edge_prms(i, j, :) = prms_i(j, :)
        end do
        if (nprms_i .gt. 1) then
            simple = .true.
        end if
    end do
    
    ! Use simple triangulation if any edge has more than one interior point.
    if (simple) then
        ! Make a single list of parameters in counter-clockwise order.
        i = 0
        ! Edge 1
        all_params(i, :) = cell_params(csn, 0, :)
        i = i + 1
        do j = 0, nprms(0) - 1
            all_params(i, :) = edge_prms(0, j, :)
            i = i + 1
        end do
        all_params(i, :) = cell_params(csn, 1, :)
        i = i + 1
        ! Edge 3
        do j = 0, nprms(1) - 1
            all_params(i, :) = edge_prms(1, j, :)
            i = i + 1
        end do
        all_params(i, :) = cell_params(csn, 2, :)
        i = i + 1
        ! Edge 4
        do j = 0, nprms(2) - 1
            all_params(i, :) = edge_prms(2, j, :)
            i = i + 1
        end do
        all_params(i, :) = cell_params(csn, 3, :)
        i = i + 1
        ! Edge 2
        do j = 0, nprms(3) - 1
            all_params(i, :) = edge_prms(3, j, :)
            i = i + 1
        end do
        all_params(i, :) = cell_params(csn, 0, :)
        i = i + 1
        
        ! Middle parameter.
        uv0 = cell_params(csn, 0, :)
        uv1 = cell_params(csn, 2, :)
        uvc = 0.50d0 * (uv0 + uv1)
        ! Generate triangles.
        do j = 0, i - 2
            uv0 = all_params(j, :)
            uv1 = all_params(j + 1, :)
            triangles(j, 0, :) = uv0
            triangles(j, 1, :) = uv1
            triangles(j, 2, :) = uvc
            ntri = ntri + 1
        end do
        return
    end if
    
    ! Use predefined triangles.
    ! Determine triangulation case by the number of interior points on
    ! each edge.
    cid = nprms(0) + nprms(1) * 10 + nprms(2) * 100 + nprms(3) * 1000
    ! Put edge parameters into array.
    eprms(:, :) = 0.0d0
    j = 0
    do i = 0, 3
        if (nprms(i) .eq. 1) then
            eprms(j, :) = edge_prms(i, 0, :)
            j = j + 1
        end if
    end do
    
    ! Case 0
    if (cid .eq. 0) then
        ntri = 2
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = cell_params(csn, 1, :)
        triangles(0, 2, :) = cell_params(csn, 2, :)
        
        triangles(1, 0, :) = cell_params(csn, 2, :)
        triangles(1, 1, :) = cell_params(csn, 3, :)
        triangles(1, 2, :) = cell_params(csn, 0, :)
        return
    end if
    
    ! Case 1
    if (cid .eq. 1) then
        ntri = 3
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = cell_params(csn, 3, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 1, :)
        triangles(1, 2, :) = cell_params(csn, 2, :)
        
        triangles(2, 0, :) = cell_params(csn, 2, :)
        triangles(2, 1, :) = cell_params(csn, 3, :)
        triangles(2, 2, :) = eprms(0, :)
        return
    end if
    
    ! Case 2
    if (cid .eq. 10) then
        ntri = 3
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = cell_params(csn, 1, :)
        triangles(0, 2, :) = eprms(0, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 2, :)
        triangles(1, 2, :) = cell_params(csn, 3, :)
        
        triangles(2, 0, :) = cell_params(csn, 3, :)
        triangles(2, 1, :) = cell_params(csn, 0, :)
        triangles(2, 2, :) = eprms(0, :)
        return
    end if
    
    ! Case 3
    if (cid .eq. 11) then
        ntri = 4
        triangles(0, 0, :) = cell_params(csn, 0 , :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = cell_params(csn, 3 , :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 1 , :)
        triangles(1, 2, :) = eprms(1, :)
        
        triangles(2, 0, :) = eprms(1, :)
        triangles(2, 1, :) = cell_params(csn, 2 , :)
        triangles(2, 2, :) = cell_params(csn, 3 , :)
        
        triangles(3, 0, :) = cell_params(csn, 3 , :)
        triangles(3, 1, :) = eprms(0, :)
        triangles(3, 2, :) = eprms(1, :)
        return
    end if
    
    ! Case 4
    if (cid .eq. 100) then
        ntri = 3
        triangles(0, 0, :) = cell_params(csn, 0 ,:)
        triangles(0, 1, :) = cell_params(csn, 1 ,:)
        triangles(0, 2, :) = eprms(0, :)
        
        triangles(1, 0, :) = cell_params(csn, 1 ,:)
        triangles(1, 1, :) = cell_params(csn, 2 ,:)
        triangles(1, 2, :) = eprms(0, :)
        
        triangles(2, 0, :) = eprms(0, :)
        triangles(2, 1, :) = cell_params(csn, 3 ,:)
        triangles(2, 2, :) = cell_params(csn, 0 ,:)
        return
    end if
    
    ! Case 5
    if (cid .eq. 101) then
        ntri = 4
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = cell_params(csn, 3, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = eprms(1, :)
        triangles(1, 2, :) = cell_params(csn, 3, :)
        
        triangles(2, 0, :) = eprms(0, :)
        triangles(2, 1, :) = cell_params(csn, 1, :)
        triangles(2, 2, :) = eprms(1, :)
        
        triangles(3, 0, :) = cell_params(csn, 1, :)
        triangles(3, 1, :) = cell_params(csn, 2, :)
        triangles(3, 2, :) = eprms(1, :)
        return
    end if
    
    ! Case 6
    if (cid .eq. 110) then
        ntri = 4
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = cell_params(csn, 1, :)
        triangles(0, 2, :) = eprms(0, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 2, :)
        triangles(1, 2, :) = eprms(1, :)
        
        triangles(2, 0, :) = eprms(1, :)
        triangles(2, 1, :) = cell_params(csn, 0, :)
        triangles(2, 2, :) = eprms(0, :)
        
        triangles(3, 0, :) = eprms(1, :)
        triangles(3, 1, :) = cell_params(csn, 3, :)
        triangles(3, 2, :) = cell_params(csn, 0, :)
        return
    end if
    
    ! Case 7
    if (cid .eq. 111) then
        ntri = 5
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = eprms(2, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 1, :)
        triangles(1, 2, :) = eprms(1, :)
        
        triangles(2, 0, :) = eprms(1, :)
        triangles(2, 1, :) = cell_params(csn, 2, :)
        triangles(2, 2, :) = eprms(2, :)
        
        triangles(3, 0, :) = eprms(2, :)
        triangles(3, 1, :) = eprms(0, :)
        triangles(3, 2, :) = eprms(1, :)
    
        triangles(4, 0, :) = eprms(2, :)
        triangles(4, 1, :) = cell_params(csn, 3, :)
        triangles(4, 2, :) = cell_params(csn, 0, :)
        return
    end if
    
    ! Case 8
    if (cid .eq. 1000) then
        ntri = 3
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = cell_params(csn, 1, :)
        triangles(0, 2, :) = eprms(0, :)
        
        triangles(1, 0, :) = cell_params(csn, 1, :)
        triangles(1, 1, :) = cell_params(csn, 2, :)
        triangles(1, 2, :) = eprms(0, :)
        
        triangles(2, 0, :) = cell_params(csn, 2, :)
        triangles(2, 1, :) = cell_params(csn, 3, :)
        triangles(2, 2, :) = eprms(0, :)
        return
    end if
    
    ! Case 9
    if (cid .eq. 1001) then
        ntri = 4
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = eprms(1, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 1, :)
        triangles(1, 2, :) = cell_params(csn, 2, :)
        
        triangles(2, 0, :) = cell_params(csn, 2, :)
        triangles(2, 1, :) = eprms(1, :)
        triangles(2, 2, :) = eprms(0, :)
        
        triangles(3, 0, :) = cell_params(csn, 2, :)
        triangles(3, 1, :) = cell_params(csn, 3, :)
        triangles(3, 2, :) = eprms(1, :)
        return
    end if
    
    ! Case 10
    if (cid .eq. 1010) then
        ntri = 4
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = cell_params(csn, 1, :)
        triangles(0, 2, :) = eprms(0, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 2, :)
        triangles(1, 2, :) = eprms(1, :)
        
        triangles(2, 0, :) = cell_params(csn, 2, :)
        triangles(2, 1, :) = cell_params(csn, 3, :)
        triangles(2, 2, :) = eprms(1, :)
        
        triangles(3, 0, :) = eprms(1, :)
        triangles(3, 1, :) = cell_params(csn, 0, :)
        triangles(3, 2, :) = eprms(0, :)
        return
    end if
    
    ! Case 11
    if (cid .eq. 1011) then
        ntri = 5
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = eprms(2, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = eprms(1, :)
        triangles(1, 2, :) = eprms(2, :)
        
        triangles(2, 0, :) = eprms(0, :)
        triangles(2, 1, :) = cell_params(csn, 1, :)
        triangles(2, 2, :) = eprms(1, :)
        
        triangles(3, 0, :) = eprms(1, :)
        triangles(3, 1, :) = cell_params(csn, 2, :)
        triangles(3, 2, :) = eprms(2, :)
    
        triangles(4, 0, :) = cell_params(csn, 2, :)
        triangles(4, 1, :) = cell_params(csn, 3, :)
        triangles(4, 2, :) = eprms(2, :)
        return
    end if
    
    ! Case 12
    if (cid .eq. 1100) then
        ntri = 4
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = cell_params(csn, 1, :)
        triangles(0, 2, :) = eprms(1, :)
        
        triangles(1, 0, :) = cell_params(csn, 1, :)
        triangles(1, 1, :) = cell_params(csn, 2, :)
        triangles(1, 2, :) = eprms(0, :)
        
        triangles(2, 0, :) = eprms(0, :)
        triangles(2, 1, :) = eprms(1, :)
        triangles(2, 2, :) = cell_params(csn, 1, :)
        
        triangles(3, 0, :) = eprms(0, :)
        triangles(3, 1, :) = cell_params(csn, 3, :)
        triangles(3, 2, :) = eprms(1, :)
        return
    end if
    
    ! Case 13
    if (cid .eq. 1101) then
        ntri = 5
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = eprms(2, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = cell_params(csn, 1, :)
        triangles(1, 2, :) = cell_params(csn, 2, :)
        
        triangles(2, 0, :) = cell_params(csn, 2, :)
        triangles(2, 1, :) = eprms(1, :)
        triangles(2, 2, :) = eprms(0, :)
        
        triangles(3, 0, :) = eprms(1, :)
        triangles(3, 1, :) = eprms(2, :)
        triangles(3, 2, :) = eprms(0, :)
    
        triangles(4, 0, :) = eprms(1, :)
        triangles(4, 1, :) = cell_params(csn, 3, :)
        triangles(4, 2, :) = eprms(2, :)
        return
    end if
    
    ! Case 14
    if (cid .eq. 1110) then
        ntri = 5
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = cell_params(csn, 1, :)
        triangles(0, 2, :) = eprms(0, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = eprms(2, :)
        triangles(1, 2, :) = cell_params(csn, 0, :)
        
        triangles(2, 0, :) = eprms(0, :)
        triangles(2, 1, :) = cell_params(csn, 2, :)
        triangles(2, 2, :) = eprms(1, :)
        
        triangles(3, 0, :) = eprms(1, :)
        triangles(3, 1, :) = eprms(2, :)
        triangles(3, 2, :) = eprms(0, :)
    
        triangles(4, 0, :) = eprms(1, :)
        triangles(4, 1, :) = cell_params(csn, 3, :)
        triangles(4, 2, :) = eprms(2, :)
        return
    end if
    
    ! Case 15
    if (cid .eq. 1111) then
        ntri = 6
        triangles(0, 0, :) = cell_params(csn, 0, :)
        triangles(0, 1, :) = eprms(0, :)
        triangles(0, 2, :) = eprms(3, :)
        
        triangles(1, 0, :) = eprms(0, :)
        triangles(1, 1, :) = eprms(1, :)
        triangles(1, 2, :) = eprms(3, :)
        
        triangles(2, 0, :) = eprms(0, :)
        triangles(2, 1, :) = cell_params(csn, 1, :)
        triangles(2, 2, :) = eprms(1, :)
        
        triangles(3, 0, :) = eprms(1, :)
        triangles(3, 1, :) = cell_params(csn, 2, :)
        triangles(3, 2, :) = eprms(2, :)
    
        triangles(4, 0, :) = eprms(2, :)
        triangles(4, 1, :) = eprms(3, :)
        triangles(4, 2, :) = eprms(1, :)
    
        triangles(5, 0, :) = eprms(2, :)
        triangles(5, 1, :) = cell_params(csn, 3, :)
        triangles(5, 2, :) = eprms(3, :)
        return
    end if
end subroutine tessellate_cell

subroutine find_neighbors(csn, esn, nmax, children, acells, position, &
                          parent, nadj, adj_cells)
    !> Find adjacent cells of a given cell.
    !> csn - Cell number.
    !> esn - Edge number.
    !> nmax - Maximum array size.
    !> children - Children of each cell.
    !> acells - Adjacency array for each cell.
    !> position - Position of each cell.
    !> parent - Parent of each cell.
    !> nadj - Number of adjacent cells.
    !> adj_cells - Adjacent cells.
    
    !f2py intent(in) csn, esn, nmax, children, acells, position, parent
    !f2py intent(out) nadj, adj_cells
    !f2py depend(nmax) children, acells, position, parent, adj_cells
    
    ! Input
    integer, intent(in) :: csn, esn, nmax
    integer, intent(in) :: children(0:nmax - 1, 0:4)
    integer, intent(in) :: acells(0:nmax - 1, 0:4)
    integer, intent(in) :: position(0:nmax - 1)
    integer, intent(in) :: parent(0:nmax - 1)
    
    ! Output
    integer, intent(out) :: nadj
    integer, intent(out) :: adj_cells(0:100)
    
    ! Working
    integer :: acell, pcell, adj_pcell, icell, pos, nstack, pos1, pos2
    integer :: i, n, s1, s2
    integer :: stack(0:100), parents(0:100)
    integer :: adj_cell_pos1(0:4, 0:4)
    integer :: adj_cell_pos2(0:4, 0:1)
    
    ! Array for looking up adjacent cell based on cell position and edge.
    adj_cell_pos1(:, :) = 0
    adj_cell_pos1(1, :) = (/ 0, 2, 3, 1, 1 /)
    adj_cell_pos1(2, :) = (/ 0, 2, 4, 2, 1 /)
    adj_cell_pos1(3, :) = (/ 0, 4, 3, 1, 3 /)
    adj_cell_pos1(4, :) = (/ 0, 4, 4, 2, 3 /)

    ! Array for looking up two adjacent cell positions.
    adj_cell_pos2(:, :) = 0
    adj_cell_pos2(1, :) = (/ 2, 4 /)
    adj_cell_pos2(2, :) = (/ 4, 3 /)
    adj_cell_pos2(3, :) = (/ 1, 2 /)
    adj_cell_pos2(4, :) = (/ 3, 1 /)
    
    ! Initialize.
    adj_cells(:) = 0
    
    ! Find adjacent cell.
    acell = acells(csn, esn)
    
    ! If adjacent cell has no children or is at a subdivision level
    ! equal to or higher than the current cell, there are no interior
    ! points.
    if ((acell .eq. 0) .or. (children(acell, 0) .eq. 0)) then
        nadj = 0
        return
    end if
    
    ! Traverse up the tree from the current cell until you find a
    ! common parent. Track the parent cell numbers along the way.
    pcell = csn
    adj_pcell = parent(acell)
    parents(:) = 0
    parents(0) = pcell
    n = 1
    do while (parent(pcell) .ne. adj_pcell)
        parents(n) = parent(pcell)
        pcell = parent(pcell)
        n = n + 1
    end do
    
    ! Traverse back down the tree from the highest parent cell to the
    ! current cell, tracking adjacent cells along the way. If a cell is
    ! reached before you get back to current cell, then the current
    ! cell has no interior points.
    do i = n - 2, 0, -1
        icell = parents(i)
        pos = position(icell)
        acell = children(acell, adj_cell_pos1(pos, esn))
        if (acell .eq. 0) then
            nadj = 0
            return
        end if
    end do
    
    ! For each subdivision level beyond current cell, gather the two
    ! adjacent cells next to specified edge. These will have interior
    ! points on the current edge.
    stack(:) = 0
    stack(0) = acell
    nstack = 1
    pos1 = adj_cell_pos2(esn, 0)
    pos2 = adj_cell_pos2(esn, 1)
    nadj = 0
    do while (nstack > 0)
        acell = stack(nstack - 1)
        nstack = nstack - 1
        if ((acell .eq. 0) .or. (children(acell, 0) .eq. 0)) then
            adj_cells(nadj) = acell
            nadj = nadj + 1
        else
            stack(nstack) = children(acell, pos2)
            nstack = nstack + 1
            stack(nstack) = children(acell, pos1)
            nstack = nstack + 1
        end if
    end do
end subroutine find_neighbors

subroutine find_edge_params(esn, nadj, nmax, adj_cells, cell_params, &
                            nparams, params)
    !> Find interior parameters along an edge in a cell.
    !> esn - Edge number.
    !> nadj - Number of adjacent cells.
    !> adj_cells - Adjacent cells.
    !> cell_params - Corner parameters of each cell.
    !> nparams - Number of parameters.
    !> params - Edge parameters.
    
    !f2py intent(in) esn, nadj, nmax, adj_cells, cell_params
    !f2py intent(out) nparams, params
    !f2py depend(nmax) adj_cells, cell_params, params
    
    ! Input
    integer, intent(in) :: esn, nadj, nmax
    integer, intent(in) :: adj_cells(0:100)
    double precision, intent(in) :: cell_params(0:nmax - 1, 0:3, 0:1)
    
    ! Output
    integer, intent(out) :: nparams
    double precision, intent(out) :: params(0:100, 0:1)
    
    ! Working
    integer :: i, icell, indx
    integer :: adj_edge_param(0:4, 0:1)
    
    ! Initialize
    nparams = 0
    params(:, :) = 0.0d0
    
    ! Array for looking up adjacent edge parameters.
    adj_edge_param(:, :) = 0
    adj_edge_param(1, :) = (/ 3, 2 /)
    adj_edge_param(2, :) = (/ 2, 1 /)
    adj_edge_param(3, :) = (/ 0, 3 /)
    adj_edge_param(4, :) = (/ 1, 0 /)
    
    if (nadj .eq. 1) then
        return
    end if
    
    do i = 0, nadj - 1
        icell = adj_cells(i)
        if (icell .ne. 0) then
            if (i .eq. 0) then
                ! Get single corner.
                indx = adj_edge_param(esn, 1)
                params(nparams, :) = cell_params(icell, indx, :)
                nparams = nparams + 1
            elseif ((i .gt. 0) .and. (i .lt. nadj - 1)) then
                if (adj_cells(i - 1) .eq. 0) then
                    ! Get both corners if previous cell was empty.
                    indx = adj_edge_param(esn, 0)
                    params(nparams, :) = cell_params(icell, indx, :)
                    nparams = nparams + 1
                    indx = adj_edge_param(esn, 1)
                    params(nparams, :) = cell_params(icell, indx, :)
                    nparams = nparams + 1
                else
                    ! Get single corner.
                    indx = adj_edge_param(esn, 1)
                    params(nparams, :) = cell_params(icell, indx, :)
                    nparams = nparams + 1
                end if
            elseif (i .eq. nadj - 1) then
                ! Get parameter if previous cell was empty.
                if (adj_cells(i - 1) .eq. 0) then
                    ! Get single corner.
                    indx = adj_edge_param(esn, 0)
                    params(nparams, :) = cell_params(icell, indx, :)
                    nparams = nparams + 1
                end if
            end if
        end if
    end do
    
end subroutine find_edge_params

end module tessellate