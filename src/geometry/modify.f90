module modify
use config, only: ptol
use compare
use evaluate
implicit none

private
public :: curve_knot_ins, surface_knot_ins
public :: decompose_curve, decompose_surface

contains

subroutine curve_knot_ins(n, p, uk, cpw, u, r, uq, qw)
    !> Curve knot insertion.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> u - Parameter to insert.
    !> r - Number of times to insert parameter.
    !> uq - New knot vector.
    !> qw - New control points.
    
    !f2py intent(in) n, p, uk, cpw, u, r
    !f2py intent(out) uq, qw
    !f2py depend(n, p) uk
    !f2py depend(n) cpw
    !f2py depend(n, r, p) uq
    !f2py depend(n, r) qw
    
    ! Input
    integer, intent(in) :: n, p, r
    double precision, intent(in) :: u
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision, intent(out) :: uq(0:n + p + r + 1)
    double precision, intent(out) :: qw(0:n + r, 0:3)
    
    ! Working
    integer :: i, j, t, k, s
    double precision :: alpha
    double precision :: rw(0:p, 0:3)
    
    call find_span_mult(n, p, u, uk, k, s)
    ! generate knot vector
    uq(:) = 0.0d0
    do i = 0, k
        uq(i) = uk(i)
    end do
    do i = 1, r
        uq(k + i) = u
    end do
    do i = k + 1, n + p + 1
        uq(i + r) = uk(i)
    end do
    ! save unchanged control points
    qw(:, :) = 0.0d0
    do i = 0, k - p
        qw(i, :) = cpw(i, :)
    end do
    do i = k - s, n
        qw(i + r, :) = cpw(i, :)
    end do
    rw(:, :) = 0.0d0
    do i = 0, p - s
        rw(i, :) = cpw(k - p + i, :)
    end do
    ! insert knot r times
    t = k - p
    do j = 1, r
        t = k - p + j
        do i = 0, p - j - s
            alpha = (u - uk(t + i)) / (uk(i + k + 1) - uk(t + i))
            rw(i, :) =alpha * rw(i + 1, :) + (1.0d0 - alpha) * rw(i, :)
        end do
        qw(t, :) = rw(0, :)
        qw(k + r - j - s, :) = rw(p - j - s, :)
    end do
    ! load remaining control points
    do i = t + 1, k - s - 1
        qw(i, :) = rw(i - t, :)
    end do
    
end subroutine curve_knot_ins

subroutine surface_knot_ins(n, p, uk, m, q, vk, cpw, dir, uv, r, &
                            d1, d2, d3, d4, uq, vq, qw)
    !> Surface knot insertion.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> dir - Direction to insert knot ('u' or 'v').
    !> uv - Parameter to insert.
    !> r - Number of times to insert parameter.
    !> d1 - Dimension for u-direction knot vector.
    !> d2 - Dimension for v-direction knot vector.
    !> d3 - Dimension 1 for control point output.
    !> d4 - Dimension 2 for control point output.
    !> uq - New knot vector in u-direction.
    !> vq - New knot vector in v-direction.
    !> qw - New control points.
    
    !f2py intent(in) n, p, uk, m, q, vk, cpw, dir, uv, r, d1, d2, d3, d4
    !f2py intent(out) uq, vq, qw
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    !f2py depend(d1) uq
    !f2py depend(d2) vq
    !f2py depend(d3, d4) qw
    
    ! Input
    character, intent(in) :: dir
    integer, intent(in) :: n, p, m, q, r, d1, d2, d3, d4
    double precision, intent(in) :: uv
    double precision, intent(in) :: uk(0:n + p + 1), vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision, intent(out) :: uq(0:d1), vq(0:d2)
    double precision, intent(out) :: qw(0:d3, 0:d4, 0:3)
    
    ! Working
    integer :: i, j, t, k, s, row, col
    double precision :: alpha_u(0:p, 0:r), alpha_v(0:q, 0:r)
    double precision :: rw_u(0:p, 0:3), rw_v(0:q, 0:3)
    
    if (dir .eq. 'u') then
        ! Find span and multiplicity.
        call find_span_mult(n, p, uv, uk, k, s)
        uq(:) = 0.0d0
        ! Generate knot vector.
        do i = 0, k
            uq(i) = uk(i)
        end do
        do i = 1, r
            uq(k + i) = uv
        end do
        do i = k + 1, n + p + 1
            uq(i + r) = uk(i)
        end do
        ! Copy v-direction knot vector.
        vq = vk
        ! Save alphas.
        alpha_u(:, :) = 0.0d0
        do j = 1, r
            t = k - p + j
            do i = 0, p - j - s
                alpha_u(i, j) = (uv - uk(t + i)) / (uk(i + k + 1) - uk(t + i))
            end do
        end do
        ! Insert knot for each row at each column.
        qw(:, :, :) = 0.0d0
        do col = 0, m
            ! Save unchanged control points.
            do i = 0, k - p
                qw(i, col, :) = cpw(i, col, :)
            end do
            do i = k - s, n
                qw(i + r, col, :) = cpw(i, col, :)
            end do
            ! Load auxiliary control points.
            rw_u(:, :) = 0.0d0
            do i = 0, p - s
                rw_u(i, :) = cpw(k - p + i, col, :)
            end do
            ! Insert the knot r times.
            t = k - p
            do j = 1, r
                t = k - p + j
                do i = 0, p - j - s
                    rw_u(i, :) = alpha_u(i, j) * rw_u(i + 1, :) + &
                                 (1.0d0 - alpha_u(i, j)) * rw_u(i, :)
                end do
                qw(t, col, :) = rw_u(0, :)
                qw(k + r - j - s, col, :) = rw_u(p - j - s, :)
            end do
            ! Load remaining control points.
            do i = t + 1, k - s - 1
                qw(i, col, :) = rw_u(i - t, :)
            end do
        end do
        return
    elseif (dir .eq. 'v') then
        ! Find span and multiplicity.
        call find_span_mult(m, q, uv, vk, k, s)
        vq(:) = 0.0d0
        ! Generate knot vector.
        do i = 0, k
            vq(i) = vk(i)
        end do
        do i = 1, r
            vq(k + i) = uv
        end do
        do i = k + 1, m + q + 1
            vq(i + r) = vk(i)
        end do
        ! Copy u-direction knot vector.
        uq = uk
        ! Save alphas.
        alpha_v(:, :) = 0.0d0
        do j = 1, r
            t = k - q + j
            do i = 0, q - j - s
                alpha_v(i, j) = (uv - vk(t + i)) / (vk(i + k + 1) - vk(t + i))
            end do
        end do
        ! Insert knot for each row at each column.
        qw(:, :, :) = 0.0d0
        do row = 0, n
            ! Save unchanged control points.
            do i = 0, k - q
                qw(row, i, :) = cpw(row, i, :)
            end do
            do i = k - s, m
                qw(row, i + r, :) = cpw(row, i, :)
            end do
            ! Load auxiliary control points.
            rw_v(:, :) = 0.0d0
            do i = 0, q - s
                rw_v(i, :) = cpw(row, k - q + i, :)
            end do
            ! Insert the knot r times.
            t = k - q
            do j = 1, r
                t = k - q + j
                do i = 0, q - j - s
                    rw_v(i, :) = alpha_v(i, j) * rw_v(i + 1, :) + &
                                 (1.0d0 - alpha_v(i, j)) * rw_v(i, :)
                end do
                qw(row, t, :) = rw_v(0, :)
                qw(row, k + r - j - s, :) = rw_v(q - j - s, :)
            end do
            ! Load remaining control points.
            do i = t + 1, k - s - 1
                qw(row, i, :) = rw_v(i - t, :)
            end do
        end do
        return
    endif    
    
end subroutine surface_knot_ins
    
subroutine decompose_curve(n, p, uk, cpw, nb, qw, ab)
    !> Decompose NURBS curve into Bezier segments.
    !> n - Number of control points - 1.
    !> p - Degree.
    !> uk - Knot vector.
    !> cpw - Control points.
    !> nb - Number of Bezier segments.
    !> qw - Control points for each Bezier segment.
    !> ab - Domain for each Bezier segment.
    
    !f2py intent(in) n, p, uk, cpw
    !f2py intent(out) nb, qw, ab
    !f2py depend(n, p) uk, qw, ab
    !f2py depend(n) cpw
    
    ! Input
    integer, intent(in) :: n, p
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    integer, intent(out) :: nb
    double precision, intent(out) :: qw(0:n + p + 1, 0:p, 0:3)
    double precision, intent(out) :: ab(0:n + p + 1, 0:1)
    
    ! Working
    integer :: a, b, i, j, k, m, mult, s1, s2
    double precision ::numer, alpha
    double precision :: alphas(0:p)
    
    qw(:, :, :) = 0.0d0
    ab(:, :) = 0.0d0
    
    m = n + p + 1
    a = p
    b = p + 1
    nb = 0
    do i = 0, p
        qw(nb, i, :) = cpw(i, :)
    end do
    do while (b .lt. m)
        i = b
        ab(nb, 0) = uk(b - 1)
        do while ((b .lt. m) .and. (EqualTo(uk(b + 1), uk(b), ptol)))
            b = b + 1
            if (b + 1 .gt. m) then
                exit
            end if
        end do
        ab(nb, 1) = uk(b)
        mult = b - i + 1
        if (mult .lt. p) then
            numer = uk(b) - uk(a)
            alphas(:) = 0.0d0
            do j = p, mult + 1, -1
                alphas(j - mult - 1) = numer / (uk(a + j) - uk(a))
            end do
            do j = 1, p - mult
                s1 = p - mult - j
                s2 = mult + j
                do k = p, s2, -1
                    alpha = alphas(k - s2)
                    qw(nb, k, :) = alpha * qw(nb, k, :) + &
                                   (1.0d0 - alpha) * qw(nb, k - 1, :)
                end do
                if (b .lt. m) then
                    qw(nb + 1, s1, :) = qw(nb, p, :)
                end if
            end do
        end if
        nb = nb + 1
        if (b .lt. m) then
            do i = p - mult, p
                qw(nb, i, :) = cpw(b - p + i, :)
            end do
            a = b
            b = b + 1
        end if
    end do
    
end subroutine decompose_curve

subroutine decompose_surface(n, p, uk, m, q, vk, cpw, dir, d1, d2, d3, &
                             nb, qw, ab)
    !> Decompose a NURBS surface into Bezier patches in a specified direction.
    !> n - Number of control points - 1 in u-direction.
    !> p - Degree in u-direction.
    !> uk - Knot vector in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> q - Degree in v-direction.
    !> vk - Knot vector in v-direction.
    !> cpw - Control points.
    !> dir - Direction to decompose surface knot ('u' or 'v').
    !> d1 - Dimension 1 for output array.
    !> d2 - Dimension 2 for output array.
    !> d3 - Dimension 3 for output array.
    !> nb - Number of Bezier patches.
    !> qw - Control points for each Bezier patch.
    !> ab - Surface domain for each Bezier patch.
    
    !f2py intent(in) n, p, uk, m, q, vk, cpw, dir, d1, d2, d3
    !f2py intent(out) nb, qw, ab
    !f2py depend(n, p) uk
    !f2py depend(m, q) vk
    !f2py depend(n, m) cpw
    !f2py depend(d1, d2, d3) qw
    !f2py depend(d1) ab
    
    ! Input
    character, intent(in) :: dir
    integer, intent(in) :: n, p, m, q, d1, d2, d3
    double precision, intent(in) :: uk(0:n + p + 1)
    double precision, intent(in) :: vk(0:m + q + 1)
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    integer, intent(out) :: nb 
    double precision, intent(out) :: qw(0:d1, 0:d2, 0:d3, 0:3)
    double precision, intent(out) :: ab(0:d1, 0:1)
                                        
    ! Working
    integer :: a, b, i, j, k, mult, s1, s2, row, col, r
    double precision :: numer, alpha
    double precision :: alphas(0:max(p, q))
    
    qw(:, :, :, :) = 0.0d0
    ab(:, :) = 0.0d0
    
    if (dir .eq. 'u') then
        a = p
        b = p + 1
        nb = 0
        r = n + p + 1
        do i = 0, p
            do col = 0, m
                qw(nb, i, col, :) = cpw(i, col, :)
            end do
        end do
        do while (b .lt. r)
            i = b
            ab(nb, 0) = uk(b - 1)
            do while ((b .lt. r) .and. (EqualTo(uk(b + 1), uk(b), ptol)))
                b = b + 1
                if (b + 1 .gt. r) then
                    exit
                end if
            end do
            ab(nb, 1) = uk(b)
            mult = b - i + 1
            if (mult .lt. p) then
                numer = uk(b) - uk(a)
                alphas(:) = 0.0d0
                do j = p, mult + 1, -1
                    alphas(j - mult - 1) = numer / (uk(a + j) - uk(a))
                end do
                do j = 1, p - mult
                    s1 = p - mult - j
                    s2 = mult + j
                    do k = p, s2, -1
                        alpha = alphas(k - s2)
                        do col = 0, m
                            qw(nb, k, col, :) = alpha * &
                                                qw(nb, k, col, :) + &
                                                (1.0d0 - alpha) * &
                                                qw(nb, k - 1, col, :)
                        end do
                    end do
                    if (b .lt. r) then
                        do col = 0, m
                            qw(nb + 1, s1, col, :) = qw(nb, p, col, :)
                        end do
                    end if
                end do
            end if
            nb = nb + 1
            if (b .lt. r) then
                do i = p - mult, p
                    do col = 0, m
                        qw(nb, i, col, :) = cpw(b - p + i, col, :)
                    end do
                end do
                a = b
                b = b + 1
            end if
        end do
    elseif (dir .eq. 'v') then
        a = q
        b = q + 1
        nb = 0
        r = m + q + 1
        do i = 0, q
            do row = 0, n
                qw(nb, row, i, :) = cpw(row, i, :)
            end do
        end do
        do while (b .lt. r)
            i = b
            ab(nb, 0) = vk(b - 1)
            do while ((b .lt. r) .and. (EqualTo(vk(b + 1), vk(b), ptol)))
                b = b + 1
                if (b + 1 .gt. r) then
                    exit
                end if
            end do
            ab(nb, 1) = vk(b)
            mult = b - i + 1
            if (mult .lt. q) then
                numer = vk(b) - vk(a)
                alphas(:) = 0.0d0
                do j = q, mult + 1, -1
                    alphas(j - mult - 1) = numer / (vk(a + j) - vk(a))
                end do
                do j = 1, q - mult
                    s1 = q - mult - j
                    s2 = mult + j
                    do k = q, s2, -1
                        alpha = alphas(k - s2)
                        do row = 0, n
                            qw(nb, row, k, :) = alpha * &
                                                qw(nb, row, k, :) + &
                                                (1.0d0 - alpha) * &
                                                qw(nb, row, k - 1, :)
                        end do
                    end do
                    if (b .lt. r) then
                        do row = 0, n
                            qw(nb + 1, row, s1, :) = qw(nb, row, q, :)
                        end do
                    end if
                end do
            end if
            nb = nb + 1
            if (b .lt. r) then
                do i = q - mult, q
                    do row = 0, n
                        qw(nb, row, i, :) = cpw(row, b - q + i, :)
                    end do
                end do
                a = b
                b = b + 1
            end if
        end do
    end if
end subroutine decompose_surface

end module modify