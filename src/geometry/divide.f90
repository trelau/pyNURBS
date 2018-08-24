module divide
use evaluate
use modify
implicit none

private
public :: split_bezier_curve, extract_bezier_curve, split_bezier_surface
public :: split_curve, split_surface

contains

subroutine split_bezier_curve(cpw, n, u, a, b, qw1, a1, b1, qw2, a2, b2)
    !> Split Bezier curve into two segments.
    !> cpw - Control points.
    !> n - Number of control points - 1.
    !> u - Parameter to split curve at (0 <= u <= 1).
    !> a - Lower domain of curve (0 <= a < b).
    !> b - Upper domain of curve (a < b <= 1).
    !> qw1 - Control points for curve 1.
    !> a1 - Lower domain for curve 1.
    !> b1 - Upper domain for curve 1.
    !> qw2 - Control points for curve 2.
    !> a2 - Lower domain for curve 2.
    !> b2 - Upper domain for curve 2.
    
    !f2py intent(in) cpw, n, u, a, b
    !f2py intent(out) qw1, a1, b1, qw2, a2, b2
    !f2py depend(n) cpw, qw1, qw2
    
    ! Input    
    integer, intent(in) :: n
    double precision, intent(in) :: u, a, b
    double precision, intent(in) :: cpw(0:n,0:3)
    
    ! Output
    double precision, intent(out) :: a1, b1, a2, b2
    double precision, intent(out) :: qw1(0:n, 0:3), qw2(0:n, 0:3)
    
    ! Working
    integer :: k, i
    
    qw1(:, :) = 0.0d0
    qw2(:, :) = 0.0d0
    
    ! New control points.
    qw2 = cpw
    qw1(0, :) = qw2(0, :)    
    do k = 1, n
        do i = 0, n - k
            qw2(i, :) = (1.d0 - u) * qw2(i, :) + u * qw2(i + 1, :)
            qw1(k, :) = qw2(0, :)
        end do
    end do
    ! New domain.
    a1 = a
    b1 = a + u * (b - a)
    a2 = b1
    b2 = b
    
end subroutine split_bezier_curve

subroutine extract_bezier_curve(cpw, n, u0, u1, a, b, qw, ae, be)
    !> Extract a Bezier curve between parameters.
    !> cpw - Control points.
    !> n - Number of control points - 1.
    !> u0 - Starting parameter (0 <= u0 < u1).
    !> u1 - Ending parameter (u0 < u1 <= 1).
    !> a - Lower domain of curve (0 <= a < b).
    !> b - Upper domain of curve (a < b <= 1).
    !> qw - Control points of extracted Bezier curve.
    !> ae - Lower domain of extracted curve.
    !> be - Upper domain of extracted curve.
    
    !f2py intent(in) cpw, n, u0, u1, a, b
    !f2py intent(out) qw, ae, be
    !f2py depend(n) cpw, qw
    
    ! Input
    integer, intent(in) :: n
    double precision, intent(in) :: u0, u1, a, b
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: ae, be
    double precision, intent(out) :: qw(0:n, 0:3)
    
    ! Working
    double precision :: u1b, a1, b1, a2, b2, a3, b3
    double precision :: qw1(0:n, 0:3), qw2(0:n, 0:3), qw3(0:n, 0:3)
    
    call split_bezier_curve(cpw, n, u0, a, b, qw1, a1, b1, qw2, a2, b2)
    u1b = (u1 - u0) / (1.0d0 - u0)
    call split_bezier_curve(qw2, n, u1b, a2, b2, qw, ae, be, qw3, a3, b3)
    
end subroutine extract_bezier_curve

subroutine split_bezier_surface(cpw, n, m, u, v, au, bu, av, bv, &
                                qw1, qw2, ab1, ab2)
    !> Split Bezier surface into two patches.
    !> cpw - Number of control points - 1.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> u - Parameter to split at in v-direction (v should equal 0).
    !> v - Parameter to split at in u-direction (u should equal 0).
    !> au - Lower domain of surface in u-direction.
    !> bu - Upper domain of surface in u-direction.
    !> av - Lower domain of surface in v-direction.
    !> bv - Upper domain of surface in v-direction.
    !> qw1 - Control points of surface 1.
    !> qw2 - Control points of surface 2.
    !> ab1 - Bounds of surface 1 ((au1, bu1), (av1, bv1))
    !> ab2 - Bounds of surface 2 ((au2, bu2), (av2, bv2))
    
    !f2py intent(in) cpw, n, m, u, v, au, bu, av, bv
    !f2py intent(out) qw1, qw2, ab1, ab2
    !f2py depend(n, m) cpw, qw1, qw2
    
    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: u, v, au, bu, av, bv
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: qw1(0:n, 0:m, 0:3)
    double precision, intent(out) :: qw2(0:n, 0:m, 0:3)
    double precision, intent(out) :: ab1(0:1, 0:1)
    double precision, intent(out) :: ab2(0:1, 0:1)
                                     
    ! Working
    integer :: i, j, k
    
    qw1(:, :, :) = 0.0d0
    qw2(:, :, :) = 0.0d0
    ab1(:, :) = 0.0d0
    ab2(:, :) = 0.0d0
    
    qw2 = cpw
    if (v .lt. 0.0d0) then
        ! New control points.
        qw1(0, :, :) = qw2(0, :, :)
        do j = 0, m
            do k = 1, n
                do i = 0, n - k
                    qw2(i, j, :) = (1.0d0 - u) * qw2(i, j, :) + &
                                    u * qw2(i + 1, j, :)
                end do
                qw1(k, j, :) = qw2(0, j, :)
            end do
        end do
        ! New domain.
        ab1(0, 0) = au
        ab1(0, 1) = au + u * (bu - au)
        ab1(1, 0) = av
        ab1(1, 1) = bv
        ab2(0, 0) = ab1(0, 1)
        ab2(0, 1) = bu
        ab2(1, 0) = av
        ab2(1, 1) = bv
        return
    elseif (u .lt. 0.0d0) then
        ! New control points.
        qw1(:, 0, :) = qw2(:, 0, :)
        do i = 0, n
            do k = 1, m
                do j = 0, m - k
                    qw2(i, j, :) = (1.0d0 - v) * qw2(i, j, :) + &
                                    v * qw2(i, j + 1, :)
                end do
                qw1(i, k, :) = qw2(i, 0, :)
            end do
        end do
        ! New domain.
        ab1(0, 0) = au
        ab1(0, 1) = bu
        ab1(1, 0) = av
        ab1(1, 1) = av + v * (bv - av)
        ab2(0, 0) = au
        ab2(0, 1) = bu
        ab2(1, 0) = ab1(1, 1)
        ab2(1, 1) = bv
        return
    end if
    
end subroutine split_bezier_surface

subroutine split_curve(n, p, uk, cpw, u, uk1, qw1, uk2, qw2)
    !> Split a NURBS curve.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> u - Parameter to split at.
    !> uk1 - Knot vector of curve 1.
    !> qw1 - Control points of curve 1.
    !> uk2 - Knot vector of curve 2.
    !> qw2 - Control points of curve 2.
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, allocatable, intent(out) :: uk1(:), uk2(:)
    double precision, allocatable, intent(out) :: qw1(:, :)
    double precision, allocatable, intent(out) :: qw2(:, :)
    
    ! Working
    integer :: i, k, s, r, nq
    double precision, allocatable :: uq(:)
    double precision, allocatable :: qw(:, :)
    
    ! Find multiplicity of knot.
    call find_span_mult(n, p, u, uk, k, s)
    
    ! Insert knot p - s times and update knot vector and control points.
    if (s .ge. p) then
        r = 0
        allocate(uq(0:n + p + 1))
        allocate(qw(0:n, 0:3))
        do i = 0, n + p + 1
            uq(i) = uk(i)
        end do
        do i = 0, n
            qw(i, :) = cpw(i, :)
        end do
    else
        r = p - s
        allocate(uq(0:n + r + p + 1))
        allocate(qw(0:n + r, 0:3))
        call curve_knot_ins(n, p, uk, cpw, u, r, uq, qw)
    end if
    
    ! Control points.
    nq = n + r
    allocate(qw1(0:k - s, 0:3))
    allocate(qw2(0:nq - (k - s), 0:3))
    qw1(:, :) = 0.0d0
    qw2(:, :) = 0.0d0
    ! Curve 1
    do i = 0, k - s
        qw1(i, :) = qw(i, :)
    end do
    ! Curve 2
    do i = k - s, nq
        qw2(i - (k - s), :) = qw(i, :)
    end do

    ! Knot vectors.
    allocate(uk1(0:k + r + 1))
    allocate(uk2(0:(nq + p + 1) - (k - s + 1) + 1))
    uk1(:) = 0.0d0
    uk2(:) = 0.0d0
    ! Curve 1
    do i = 0, k + r
        uk1(i) = uq(i)
    end do
    uk1(i) = u  
    ! Curve 2
    uk2(0) = u
    do i = k - s + 1, nq + p + 1
        uk2((i + 1) - (k - s + 1)) = uq(i)
    end do
    
end subroutine split_curve

subroutine split_surface(n, p, uk, m, q, vk, cpw, uv, dir, &
                         uk1, vk1, qw1, uk2, vk2, qw2)
    !> Split NURBS surface.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> uv - Parameter to split at.
    !> dir - Parameter direction ("u" or "v"). For example, 
    !> dir="u", then the surface is split at u=uv and along
    !> the v-direction.
    !> uk1 - Knot vector of surface 1 in u-direction.
    !> vk1 - Knot vector of surface 1 in v-direction.
    !> qw1 - Control points of surface 1.
    !> uk2 - Knot vector of surface 2 in u-direction.
    !> vk2 - Knot vector of surface 2 in v-direction.
    !> qw2 - Control points of surface 2.

    ! Input
    character, intent(in) :: dir
    integer, intent(in) :: n, p, m, q
    double precision, intent(in) :: uv
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, allocatable, intent(out) :: uk1(:), uk2(:)
    double precision, allocatable, intent(out) :: vk1(:), vk2(:)
    double precision, allocatable, intent(out) :: qw1(:, :, :)
    double precision, allocatable, intent(out) :: qw2(:, :, :)
    
    ! Working
    integer :: i, j, nq, mq, k, s, r
    double precision, allocatable :: uq(:)
    double precision, allocatable :: vq(:)
    double precision, allocatable :: qw(:, :, :)
    
    if (dir .eq. 'u') then
        ! Find multiplicity of knot.
        call find_span_mult(n, p, uv, uk, k, s)
        ! Insert knot p - s times and update knot vector and control points.
        if (s .ge. p) then
            r = 0
            allocate(uq(0:n + p + 1))
            allocate(qw(0:n, 0:m, 0:3))
            do i = 0, n + p + 1
                uq(i) = uk(i)
            end do
            do i = 0, n
                do j = 0, m
                    qw(i, j, :) = cpw(i, j, :)
                end do
            end do
        else
            r = p - s
            allocate(uq(0:n + r + p + 1))
            allocate(vq(0:m + q + 1))
            allocate(qw(0:n + r, 0:m, 0:3))
            call surface_knot_ins(n, p, uk, m, q, vk, cpw, 'u', uv, r, &
                                  n + r + p + 1, m + q + 1, n + r, m, &
                                  uq, vq, qw)
        end if
        nq = n + r
        ! Control points.
        allocate(qw1(0:k - s, 0:m, 0:3))
        allocate(qw2(0:nq - (k - s), 0:m, 0:3))
        qw1(:, :, :) = 0.0d0
        qw2(:, :, :) = 0.0d0
        do j = 0, m
            ! Surface 1.
            do i = 0, k - s
                qw1(i, j, :) = qw(i, j, :)
            end do
            ! Surface 2.
            do i = k - s, nq
                qw2(i - (k - s), j, :) = qw(i, j, :)
            end do
        end do
        ! Knot vector for surface 1.
        allocate(uk1(0:k + r + 1))
        uk1(:) = 0.0d0
        do i = 0, k + r
            uk1(i) = uq(i)
        end do
        uk1(i) = uv
        ! Knot vector for surface 2.
        allocate(uk2(0:(nq + p + 1) - (k - s + 1) + 1))
        uk2(:) = 0.0d0
        uk2(0) = uv
        do i = k - s + 1, nq + p + 1
            uk2((i + 1) - (k - s + 1)) = uq(i)
        end do
        ! Knot vectors in v-direction.
        allocate(vk1(0:m + q + 1))
        allocate(vk2(0:m + q + 1))
        vk1 = vk
        vk2 = vk
        return
    elseif (dir .eq. 'v') then
        ! Find multiplicity of knot.
        call find_span_mult(m, q, uv, vk, k, s)
        ! Insert knot q - s times and update knot vector and
        ! control points.
        if (s .ge. q) then
            r = 0
            allocate(vq(0:m + q + 1))
            allocate(qw(0:n, 0:m, 0:3))
            do i = 0, m + q + 1
                vq(i) = vk(i)
            end do
            do i = 0, n
                do j = 0, m
                    qw(i, j, :) = cpw(i, j, :)
                end do
            end do
        else
            r = q - s
            allocate(uq(0:n + p + 1))
            allocate(vq(0:m + r + q + 1))
            allocate(qw(0:n, 0:m + r, 0:3))
            call surface_knot_ins(n, p, uk, m, q, vk, cpw, 'v', uv, r, &
                                  n + p + 1, m + r + q + 1, n, m + r, &
                                  uq, vq, qw)
        end if
        mq = m + r
        ! Control points.
        allocate(qw1(0:n, 0:k - s, 0:3))
        allocate(qw2(0:n, 0:mq - (k - s), 0:3))
        qw1(:, :, :) = 0.0d0
        qw2(:, :, :) = 0.0d0
        do i = 0, n
            ! Surface 1.
            do j = 0, k - s
                qw1(i, j, :) = qw(i, j, :)
            end do
            ! Surface 2.
            do j = k - s, mq
                qw2(i, j - (k - s), :) = qw(i, j, :)
            end do
        end do
        ! Knot vector for surface 1.
        allocate(vk1(0:k + r + 1))
        vk1(:) = 0.0d0
        do i = 0, k + r
            vk1(i) = vq(i)
        end do
        vk1(i) = uv
        ! Knot vector for surface 2.
        allocate(vk2(0:(mq + q + 1) - (k - s + 1) + 1))
        vk2(:) = 0.0d0
        vk2(0) = uv
        do i = k - s + 1, mq + q + 1
            vk2((i + 1) - (k - s + 1)) = vq(i)
        end do
        ! Knot vectors in u-direction.
        allocate(uk1(0:n + p + 1))
        allocate(uk2(0:n + p + 1))
        uk1 = uk
        uk2 = uk
        return
    end if
    
end subroutine split_surface

end module divide