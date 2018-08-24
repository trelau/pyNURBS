module evaluate
use config, only: ptol
use compare
use math
implicit none

private
public :: basis_funs, ders_basis_funs
public :: decasteljau1, decasteljau2, bezier_curve_points, bezier_surface_points
public :: curve_point, curve_points, rat_curve_derivs, curve_derivs_alg1
public :: surface_point, surface_points, rat_surface_derivs, surface_derivs_alg1
public :: find_span, find_span_mult, find_mult_knots, find_mult
public :: check_param

contains

subroutine decasteljau1(cpw, n, u, pnt)
    !> Compute a point on a Bezier curve using deCasteljau algorithm.
    !> n - Number of control points - 1.
    !> cpw - Control points.
    !> u - Parameter (0. <= u <= 1.).
    !> pnt - Point on Bezier curve.
    
    !f2py inent(in) cpw, n, u
    !f2py intent(out) pnt
    !f2py depend(n) cpw
    
    ! Input
    integer, intent(in) :: n
    double precision, intent(in) :: u
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: pnt(0:3)
    
    ! Working
    integer :: i, k
    double precision :: u_
    double precision :: temp(0:n, 0:3)
    
    u_ = check_bounds(u, 0.0d0, 1.0d0)
    temp = cpw
    do k = 1, n
        do i = 0, n - k
            temp(i, :) = (1.0d0 - u_) * temp(i, :) + u_ * temp(i + 1, :)
        end do
    end do
    pnt = temp(0, :)
    
end subroutine decasteljau1

subroutine decasteljau2(cpw, n, m, u, v, pnt)
    !> Compute a point on a Bezier surface using deCasteljau algorithm.
    !> cpw - Control points.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> u - Parameter in u-direction (0. <= u <= 1.).
    !> v - Parameter in v-direction (0. <= v <= 1.).
    !> pnt - Point on Bezier surface.
    
    !f2py intent(in) cpw, n, m, u, v
    !f2py intent(out) pnt
    !f2py depend(n, m) cpw
    
    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: u, v
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: pnt(0:3)
    
    ! Working
    integer :: i, j, k
    double precision :: temp1(0:min(n, m), 0:3)
    double precision :: temp2(0:max(n, m), 0:3)
    
    temp1(:, :) = 0.0d0
    temp2(:, :) = 0.0d0
    if (n .le. m) then
        do j = 0, m
            temp1 = cpw(0:n, j, :)
            call decasteljau1(temp1, n, u, pnt)
            temp2(j, :) = pnt
        end do
        call decasteljau1(temp2, m, v, pnt)
    else
        do i = 0, n
            temp1 = cpw(i, 0:m, :)
            call decasteljau1(temp1, m, v, pnt)
            temp2(i, :) = pnt
        end do
        call decasteljau1(temp2, n, u, pnt)
    end if
    
end subroutine decasteljau2

subroutine find_span(n, p, u, uk, mid)
    !> Determine the knot span index.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> u - Parameter .
    !> uk - Knot vector.
    !> mid - Knot span.
    
    !f2py inent(in) n, p, u, uk
    !f2py intent(out) mid
    !f2py depend(n, p) uk
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:n + p + 1)
    
    ! Output
    integer, intent(out) :: mid
    
    ! Working
    integer :: low, high
    
    ! Special case
    if (GreaterThanEq(u, uk(n + 1), ptol)) then
        mid = n
        return
    end if
    if (LessThanEq(u, uk(p), ptol)) then
        mid = p
        return
    end if
    
    low = p
    high = n + 1
    mid = (low + high) / 2
    do while ((LessThan(u, uk(mid), ptol)) .or. (GreaterThanEq(u, uk(mid + 1), ptol)))
        if (LessThan(u, uk(mid), ptol)) then
            high = mid
        else
            low = mid
        end if
        mid = (low + high) / 2
    end do
    
end subroutine find_span

subroutine basis_funs(i, u, p, uk, m, bf)
    !> Compute the non-vanishing basis functions.
    !> i - Knot span index.
    !> u - Parameter.
    !> p - Degree.
    !> uk - Knot vector.
    !> m - Size of knot vector  - 1.
    !> bf - Non-vanishing basis functions.
    
    !f2py intent(in) i, u, p, uk, m
    !f2py intent(out) bf
    !f2py depend(m) uk
    !f2py depend(p) bf
    
    ! Input
    integer, intent(in) :: i, p, m
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:m)
    
    ! Output
    double precision, intent(out) :: bf(0:p)
    
    ! Working
    integer :: j, r
    double precision :: temp, saved, u_
    double precision :: left(0:p), right(0:p)
    
    bf(:) = 0.0d0
    left(:) = 0.0d0
    right(:) = 0.0d0
    
    u_ = check_bounds(u, uk(0), uk(m))
    bf(0) = 1.0d0
    do j = 1, p
        left(j) = u_ - uk(i + 1 - j)
        right(j) = uk(i + j) - u_
        saved = 0.0d0
        do r = 0, j - 1
            temp = bf(r) / (right(r + 1) + left(j - r))
            bf(r) = saved + right(r + 1) * temp
            saved = left(j - r) * temp
        end do
        bf(j) = saved
    end do
    
end subroutine basis_funs

subroutine ders_basis_funs(i, u, p, n, uk, m, ders)
    !> Compute nonzero basis functions and their derivatives.
    !> i - Knot span index.
    !> u - Parameter.
    !> p - Degree.
    !> n - Number of derivatives to compute (n <= p).
    !> uk - Knot vector.
    !> m - Size of knot vector - 1.
    !> ders - Basis functions and derivatives ders(k, j) where k is
    !> the k-th derivative.
    
    !f2py intent(in) i, u, p, n, uk, m
    !f2py intent(out) ders
    !f2py depend(m) uk
    !f2py depend(n, p) ders
    
    ! Input
    integer, intent(in) :: i, p, n, m
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:m)
    
    ! Output
    double precision, intent(out) :: ders(0:n, 0:p)
    
    ! Working
    integer :: j, k, r, j1, j2, rk, pk, s1, s2
    double precision :: saved, temp, d, u_
    double precision :: left(0:p), right(0:p)
    double precision :: ndu(0:p, 0:p), a(0:1, 0:p)
    
    ders(:, :) = 0.0d0
    left(:) = 0.0d0
    right(:) = 0.0d0
    ndu(:, :) = 0.0d0
    a(:, :) = 0.0d0
    
    u_ = check_bounds(u, uk(0), uk(m))
    ndu(0, 0) = 1.0d0
    do j = 1, p
        left(j) = u_ - uk(i + 1 - j)
        right(j) = uk(i + j) - u_
        saved = 0.d0
        do r = 0, j - 1
            ndu(j, r) = right(r + 1) + left(j - r)
            temp = ndu(r, j - 1) / ndu(j, r)
            ndu(r, j) = saved + right(r + 1) * temp
            saved = left(j - r) * temp
        end do
        ndu(j, j) = saved
    end do
    do j = 0, p
        ders(0, j) = ndu(j, p)
    end do
    do r = 0, p
        s1 = 0
        s2 = 1
        a(0, 0) = 1.d0
        do k = 1, n
            d = 0.d0
            rk = r - k
            pk = p - k
            if (r .ge. k) then
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk)
                d = a(s2, 0) * ndu(rk, pk)
            end if
            if (rk .ge. -1) then
                j1 = 1
            else
                j1 = -rk
            end if
            if (r - 1 .le. pk) then
                j2 = k - 1
            else
                j2 = p - r
            end if
            do j = j1, j2
                a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j)
                d = d + a(s2, j) * ndu(rk + j, pk)
            end do
            if (r .le. pk) then
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r)
                d = d + a(s2, k) * ndu(r, pk)
            end if
            ders(k, r) = d
            j = s1
            s1 = s2
            s2 = j
        end do
    end do
    r = p
    do k = 1, n
        do j = 0, p
            ders(k, j) = ders(k, j) * r
        end do
        r = r * (p - k)
    end do
end subroutine ders_basis_funs

subroutine curve_point(n, p, uk, cpw, u, pnt)
    !> Compute curve point.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> u - Parameter.
    !> pnt - Point on curve.
    
    !f2py intent(in) n, p, uk, cpw, u
    !f2py intent(out) pnt
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: pnt(0:2)
    
    ! Working
    integer :: i, span, m
    double precision :: bf(0:p), temp(0:3)
    
    m = n + p + 1
    call find_span(n, p, u, uk, span)
    call basis_funs(span, u, p, uk, m, bf)
    temp(:) = 0.0d0
    do i = 0, p
        temp = temp + bf(i) * cpw(span - p + i, :)
    end do
    pnt = temp(0:2) / temp(3)
    
end subroutine curve_point

subroutine curve_derivs_alg1(n, p, uk, cpw, u, d, ck)
    !> Compute curve derivatives.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> u - Parameter.
    !> d - Maximum derivative to compute.
    !> ck - Curve derivatives ck(k, :) where k is the k-th derivative and
    !> 0 <= k <= d.
    
    !f2py intent(in) n, p, uk, cpw, u, d
    !f2py intent(out) ck
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    !f2py depend(d) ck
    
    ! Input
    integer, intent(in) :: n, p, d
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: ck(0:d, 0:3)
    
    ! Working
    integer :: du, k, span, j, m
    double precision :: nders(0:min(d, p), 0:p)
    
    ck(:, :) = 0.0d0
    du = min(d, p)
    m = n + p + 1
    call find_span(n, p, u, uk, span)
    call ders_basis_funs(span, u, p, du, uk, m, nders)
    do k = 0, du
        do j = 0, p
            ck(k, :) = ck(k, :) + nders(k, j) * cpw(span - p + j, :)
        end do
    end do
    
end subroutine curve_derivs_alg1

subroutine rat_curve_derivs(n, p, uk, cpw, u, d, ck)
    !> Compute rational curve derivative.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> u - Parameter.
    !> d - Maximum derivative to compute.
    !> ck - Curve derivatives ck(k, :) where k is the k-th derivative and
    !> 0 <= k <= d.

    !f2py intent(in) n, p, uk, cpw, u, d
    !f2py intent(out) ck
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    !f2py depend(d) ck

    ! Input
    integer, intent(in) :: n, p, d
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: ck(0:d, 0:2)
    
    integer :: i, k
    double precision :: bc
    double precision :: cders(0:d, 0:3), aders(0:d, 0:2), wders(0:d)
    double precision :: v(0:2)
                        
    cders(:, :) = 0.0d0
    aders(:, :) = 0.0d0
    wders(:) = 0.0d0
    
    call curve_derivs_alg1(n, p, uk, cpw, u, d, cders)
    aders = cders(:, 0:2)
    wders = cders(:, 3)
    ck(:, :) = 0.0d0
    do k = 0, d
        v = aders(k, :)
        do i = 1, k
            bc = binomial(k, i)
            v = v - bc * wders(i) * ck(k - i, :)
        end do
        ck(k, :) = v / wders(0)
    end do
    
end subroutine rat_curve_derivs

subroutine surface_point(n, p, uk, m, q, vk, cpw, u, v, pnt)
    !> Compute surface point.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> u - Parameter in u-direction.
    !> v - Parameter in v-direction.
    !> pnt - Point on surface.
    
    !f2py intent(in) n, p, uk, m, q, vk, cpw, u, v
    !f2py intent(out) pnt
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw

    ! Input
    integer, intent(in) :: n, p, m, q
    double precision, intent(in) :: u, v
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: pnt(0:2)
    
    ! Working
    integer :: uspan, vspan, uind, vind, l, k, r, s
    double precision :: nu(0:p), nv(0:q)
    double precision :: temp1(0:3), temp2(0:3)
    
    call find_span(n, p, u, uk, uspan)
    r = n + p + 1
    call basis_funs(uspan, u, p, uk, r, nu)
    call find_span(m, q, v, vk, vspan)
    s = m + q + 1
    call basis_funs(vspan, v, q, vk, s, nv)
    uind = uspan - p
    temp2(:) = 0.0d0
    do l = 0, q
        temp1(:) = 0.0d0
        vind = vspan - q + l
        do k = 0, p
            temp1 = temp1 + nu(k) * cpw(uind + k, vind, :)
        end do
        temp2 = temp2 + nv(l) * temp1
    end do
    pnt = temp2(0:2) / temp2(3)
    
end subroutine surface_point

subroutine surface_derivs_alg1(n, p, uk, m, q, vk, cpw, u, v, d, skl)
    !> Compute surface derivatives.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> u - Parameter in u-direction.
    !> v - Parameter in v-direction.
    !> d - Maximum derivative to compute.
    !> skl - Surface derivatives skl(k, l) with respect to u k times
    !> and v l times and 0 <= k, l <= d.

    !f2py intent(in) n, p, uk, m, q, vk, cpw, u, v, d
    !f2py intent(out) skl
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    !f2py depend(d) skl
    
    ! Input
    integer, intent(in) :: n, p, m, q, d
    double precision, intent(in) :: u, v
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: skl(0:d, 0:d, 0:3)
    
    ! Working
    integer :: du, dv, uspan, vspan, k, s, r, dd, l
    double precision :: temp(0:q, 0:3)
    double precision, allocatable :: nu(:, :), nv(:, :)
    
    du = min(d, p)
    dv = min(d, q)
    allocate(nu(0:du, 0:p), nv(0:dv, 0:q))
    
    call find_span(n, p, u, uk, uspan)
    r = n + p + 1
    call ders_basis_funs(uspan, u, p, du, uk, r, nu)
    call find_span(m, q, v, vk, vspan)
    s = m + q + 1
    call ders_basis_funs(vspan, v, q, dv, vk, s, nv)

    skl(:, :, :) = 0.0d0
    temp(:, :) = 0.0d0
    do k = 0, du
        do s = 0, q
            temp(s, :) = 0.0d0
            do r = 0, p
                temp(s, :) = temp(s, :) + nu(k, r) * cpw(uspan - p + r, &
                                                         vspan - q + s, :)
            end do
        end do
        dd = min(d - k, dv)
        do l = 0, dd
            do s = 0, q
                skl(k, l, :) = skl(k, l, :) + nv(l, s) * temp(s, :)
            end do
        end do
    end do
    
end subroutine surface_derivs_alg1

subroutine rat_surface_derivs(n, p, uk, m, q, vk, cpw, u, v, d, skl)
    !> Compute rational surface derivatives.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> u - Parameter in u-direction.
    !> v - Parameter in v-direction.
    !> d - Maximum derivative to compute.
    !> skl - Surface derivatives skl(k, l) with respect to u k times
    !> and v l times and 0 <= k, l <= d.

    !f2py intent(in) n, p, uk, m, q, vk, cpw, u, v, d
    !f2py intent(out) skl
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    !f2py depend(d) skl

    ! Input
    integer, intent(in) :: n, p, m, q, d
    double precision, intent(in) :: u, v
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: skl(0:d, 0:d, 0:2)
    
    ! Working
    integer :: k, l, i, j
    double precision :: bc
    double precision :: v1(0:2), v2(0:2)
    double precision :: sders(0:d, 0:d, 0:3)
    double precision :: aders(0:d, 0:d, 0:2)
    double precision :: wders(0:d, 0:d)
    
    skl(:, : ,:) = 0.0d0
    sders(:, :, :) = 0.0d0
    aders(:, :, :) = 0.0d0
    wders(:, :) = 0.0d0
    
    call surface_derivs_alg1(n, p, uk, m, q, vk, cpw, u, v, d, sders)
    do i = 0, d
        do j = 0, d
            aders(i, j, :) = sders(i, j, 0:2)
            wders(i, j) = sders(i, j, 3)
        end do
    end do

    do k = 0, d
        do l = 0, d - k
            v1 = aders(k, l, :)
            do j = 1, l
                bc = binomial(l, j)
                v1 = v1 - bc * wders(0, j) * skl(k, l - j, :)
            end do
            do i = 1, k
                bc = binomial(k, i)
                v1 = v1 - bc * wders(i, 0) * skl(k - i, l, :)
                v2(:) = 0.0d0
                do j = 1, l
                    bc = binomial(l, j)
                    v2 = v2 + bc * wders(i, j) * skl(k - i, l - j, :)
                end do
                bc = binomial(k, i)
                v1 = v1 - bc * v2
            end do
            skl(k, l, :) = v1 / wders(0, 0)
        end do
    end do
    
end subroutine rat_surface_derivs

subroutine find_span_mult(n, p, u, uk, span, s)
    !> Determine the knot span index and multiplicity.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> u - Parameter.
    !> uk - Knot vector.
    !> span - Knot span index.
    !> s - Multiplicity of knot.
    
    !f2py intent(in) n, p, u, uk
    !f2py intent(out) span, s
    !f2py depend(n, p) uk

    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:n + p + 1)
    
    ! Output
    integer, intent(out) :: span, s
    
    ! Working
    integer :: i, first, last, step

    call find_span(n, p, u, uk, span)
    
    if ((GreaterThan(u, uk(span), ptol)) .and. (LessThan(u, uk(span + 1), ptol))) then
        s = 0
        return
    end if
    
    if ((GreaterThan(u, uk(span), ptol)) .and. (LessThanEq(u, uk(span + 1), ptol))) then
        first = span + 1
        last = n + p + 1
        step = 1
    else
        first = span
        last = 0
        step = -1
    end if
    
    s = 0
    do i = first, last, step
        ! if (u .eq. uk(i)) then
        if (EqualTo(u, uk(i), ptol)) then
            s = s + 1
        else
            return
        end if
    end do
    
end subroutine find_span_mult

function find_mult(n, p, uk, u) result(mult)
    !> Determine the multiplicity of the parameter.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> u - Parameter.
    !> mult - Multiplicity of parameter.
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0: n + p + 1)
    
    ! Output
    integer :: mult
    
    ! Working
    integer :: k
    
    call find_span_mult(n, p, u, uk, k, mult)
    
end function find_mult

subroutine find_mult_knots(n, p, uk, nu, um, uq)
    !> Find the multiplicities and the unique knots of the knot vector.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> nu - Number of unique knots.
    !> um - Multiplicity of each unique knot.
    !> uq - Unique knot vector.
    
    !f2py intent(in) n, p, uk
    !f2py intent(out) nu, um, uq
    !f2py depend(n, p) uk, um, uq
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: uk(0:n + p + 1)
    
    ! Output
    integer, intent(out) :: nu
    integer, intent(out) :: um(0:n + p + 1)
    double precision, intent(out) :: uq(0:n + p + 1)
    
    ! Working
    integer :: i, m, mult
    
    um(:) = 0
    uq(:) = 0.0d0
    
    i = 0
    m = n + p + 1
    nu = 0
    do while (i .le. m)
        uq(nu) = uk(i)
        mult = 0
        ! do while ((i .le. m) .and. (dabs(uk(i) .eq. uq(nu))))
        do while ((i .le. m) .and. (EqualTo(uk(i), uq(nu), ptol)))
            i = i + 1
            mult = mult + 1
            if ( i .gt. m) then
                exit
            end if
        end do
        um(nu) = mult
        nu = nu + 1
    end do
    
end subroutine find_mult_knots

subroutine check_param(n, p, uk, u, is_closed)
    !> Check that the parameter is within the knot vector or is within
    !> tolerance of a unique knot value.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> u - Parameter.
    !> Option to specify a closed knot vector.
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: uk(0:n + p + 1)
    
    ! In/Out
    double precision, intent(inout) :: u
    
    ! Optional
    logical, intent(in), optional :: is_closed
    
    ! Working
    logical :: is_closed_
    integer :: k
    
    if (present(is_closed)) then
        is_closed_ = is_closed
    else
        is_closed_ = .false.
    end if
    
    if (is_closed_) then
        if (u .gt. uk(n + 1)) then
            u = u - uk(n + 1) + uk(p)
        elseif (u .lt. uk(p)) then
            u = uk(n + 1) + u - uk(p)
        end if
    end if
    
    call find_span(n, p, u, uk, k)
    
    if (LessThanEq(u, uk(k), ptol)) then
        u = uk(k)
        return
    end if
    
    if (GreaterThanEq(u, uk(k + 1), ptol)) then
        u = uk(k + 1)
        return
    end if
    
end subroutine check_param

subroutine bezier_curve_points(n, cpw, nu, arr_u, pnts)
    !> Compute Bezier curve points.
    !> n - Number of control points - 1.
    !> cpw - Control points.
    !> nu - Number of parameters.
    !> arr_u - Parameters.
    !> pnts - Points on curve.
    
    !f2py intent(in) n, cpw, nu, arr_u
    !f2py intent(out) pnts
    !f2py depend(n) cpw
    !f2py depend(nu) arr_u
    !f2py depend(nu) pnts
    
    ! Input
    integer, intent(in) :: n, nu
    double precision, intent(in) :: arr_u(nu)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: pnts(nu, 3)
    
    ! Working
    integer :: i
    double precision :: pi(4)
    
    pnts(:, :) = 0.0d0
    do i = 1, nu
        call decasteljau1(cpw, n, arr_u(i), pi)
        pnts(i, :) = pi(1:3) / pi(4)
    end do

end subroutine bezier_curve_points

subroutine curve_points(n, p, uk, cpw, nu, arr_u, pnts)
    !> Compute curve points.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> nu - Number of parameters.
    !> arr_u - Parameters.
    !> pnts - Points on curve.
    
    !f2py intent(in) n, p, uk, cpw, nu, arr_u
    !f2py intent(out) pnts
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    !f2py depend(nu) arr_u
    !f2py depend(nu) pnts
    
    ! Input
    integer, intent(in) :: n, p, nu
    double precision, intent(in) :: arr_u(nu)
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: pnts(nu, 3)
    
    ! Working
    integer :: i
    double precision :: pi(3)
    
    pnts(:, :) = 0.0d0
    do i = 1, nu
        call curve_point(n, p, uk, cpw, arr_u(i), pi)
        pnts(i, :) = pi
    end do

end subroutine curve_points

subroutine bezier_surface_points(n, m, cpw, np, arr_u, arr_v, pnts)
    !> Compute points on a Bezier surface.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cpw - Control points.
    !> np - Number of parameters.
    !> arr_u - Parameters in u-direction.
    !> arr_v - Parameters in v-direction.
    !> pnts - Points on Bezier surface.
    
    !f2py intent(in) cpw, n, m, arr_u, arr_v, np
    !f2py intent(out) pnts
    !f2py depend(n, m) cpw
    !f2py depend(np) arr_u, arr_v, pnts
    
    ! Input
    integer, intent(in) :: n, m, np
    double precision, intent(in) :: arr_u(np), arr_v(np)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: pnts(np, 3)
    
    ! Working
    integer :: i
    double precision :: pi(4)
    
    pnts(:, :) = 0.0d0
    do i = 1, np
        call decasteljau2(cpw, n, m, arr_u(i), arr_v(i), pi)
        pnts(i, :) = pi(1:3) / pi(4)
    end do

end subroutine bezier_surface_points

subroutine surface_points(n, p, uk, m, q, vk, cpw, np, arr_u, arr_v, pnts)
    !> Compute surface points.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> np - Number of parameters.
    !> arr_u - Parameters in u-direction.
    !> arr_v - Parameters in v-direction.
    !> pnt - Points on surface.
    
    !f2py intent(in) n, p, uk, m, q, vk, cpw, np, arr_u, arr_v
    !f2py intent(out) pnts
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    !f2py depend(np) arr_u, arr_v, pnts

    ! Input
    integer, intent(in) :: n, p, m, q, np
    double precision, intent(in) :: arr_u(np), arr_v(np)
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: pnts(np, 3)
    
    ! Working
    integer :: i
    double precision :: pi(3)
    
    pnts(:, :) = 0.0d0
    do i = 1, np
        call surface_point(n, p, uk, m, q, vk, cpw, arr_u(i), arr_v(i), pi)
        pnts(i, :) = pi
    end do

end subroutine surface_points

end module evaluate