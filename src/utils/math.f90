module math
implicit none

private
public :: dummy_math
public :: norm_sq, norm, dot, cross, binomial

contains

subroutine dummy_math()
    !> Dummy routine for f2py.
end subroutine dummy_math

function norm_sq(v) result(vn_sq)
    !> Squared vector norm.
    !> v - Vector.
    !> vn_sq - Squared vector norm.
    
    ! Input
    double precision, intent(in) :: v(:)
    
    ! Output
    double precision :: vn_sq
    
    ! Working
    integer :: n, i
    
    vn_sq = 0.0d0
    n = size(v, 1)
    do i = 1, n
        vn_sq = vn_sq + v(i) * v(i)
    end do
    
end function norm_sq

function norm(v) result(vn)
    !> Vector norm.
    !> v - Vector.
    !> vn - Vector norm.

    ! Input
    double precision, intent(in) :: v(:)
    
    ! Output
    double precision :: vn
    
    ! Working
    double precision :: vn_sq
    
    vn_sq = norm_sq(v)
    vn = dsqrt(vn_sq)

end function norm

function dot(v1, v2) result(dp)
    !> Dot product of two vector.
    !> v1 - Vector 1.
    !> v2 - Vector 2.
    
    ! Input
    double precision, intent(in) :: v1(:), v2(:)
    
    ! Output
    double precision :: dp
    
    ! Working
    integer :: n1, n2, i
    
    dp = 0.0d0
    n1 = size(v1, 1)
    n2 = size(v2, 1)
    do i = 1, min(n1, n2)
        dp = dp + v1(i) * v2(i)
    end do

end function dot

function cross(v1, v2) result(vc)
    !> Cross product of two vectors.
    !> v1 - Vector 1.
    !> v2 - Vector 2.
    !> vcross - Vector cross product.
    
    ! Input
    double precision, intent(in) :: v1(3), v2(3)
    
    ! Output
    double precision :: vc(3)
    
    vc(1) = v1(2) * v2(3) - v1(3) * v2(2)
    vc(2) = -(v1(1) * v2(3) - v1(3) * v2(1))
    vc(3) = v1(1) * v2(2) - v1(2) * v2(1)
end function cross

function binomial(n, k) result(bc)
    !> Compute binomial coefficients.
    !> n - Non-negative integer (n > 0).
    !> k - Non-negative integer (k < n).
    !> bc - Binomial coefficient.
    
    ! Input
    integer, intent(in) :: n, k
    
    ! Output
    double precision :: bc
    
    ! Working
    integer :: i, kk
    
    if ((k .lt. 0) .or. ( k .gt. n)) then
        bc = 0.0d0
        return
    end if
    if ((k .eq. 0) .or. ( k .eq. n)) then
        bc = 1.0d0
        return
    end if
    kk = min(k, n - k)
    bc = 1.0d0
    do i = 0, kk - 1
        bc = bc * dble(n - i) / dble(i + 1)
    end do
    
end function binomial

end module math