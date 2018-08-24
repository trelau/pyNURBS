module invert
use config, only: ptol
use math
implicit none

private
public :: invert_points_on_plane

contains

subroutine invert_points_on_plane(n, pnts, p0, vu, vv, params)
    !> Invert the array of points on a plane.
    !> n - Number of points to invert.
    !> pnts - Points to invert.
    !> p0 - Origin of plane.
    !> vu - Vector defining major axis of plane.
    !> vv - Vector defining minor axis of plane (should be orthognal to vu).
    !> params - Parameters of points on the plane.
    
    ! f2py intent(in) n, pnts, origin, pnorm
    ! f2py intent(out) params
    ! f2py depend(n) pnts, params
    
    ! Input
    integer, intent(in) :: n
    double precision, intent(in) :: p0(3), vu(3), vv(3)
    double precision, intent(in) :: pnts(n, 3)
    
    ! Output
    double precision, intent(out) :: params(n, 2)
    
    ! Working
    integer :: i
    double precision :: mag, vu_norm2, vv_norm2, u, v
    double precision :: vu_norm(3), vv_norm(3), pi(3)
    
    
    ! Get unit vectors and dot products.
    mag = norm(vu)
    vu_norm = vu / mag
    mag = norm(vv)
    vv_norm = vv / mag
    vu_norm2 = dot_product(vu_norm, vu_norm)
    vv_norm2 = dot_product(vv_norm, vv_norm)
    
    
    ! Calculate parameters on plane for each point.
    params(:, :) = 0.0d0
    do i = 1, n
        pi = pnts(i, :)
        u = dot_product(pi - p0, vu_norm) / vu_norm2
        v = dot_product(pi - p0, vv_norm) / vv_norm2
        if (dabs(u) .le. ptol) then
            u = 0.0d0
        end if
        if (dabs(v) .le. ptol) then
            v = 0.0d0
        end if
        params(i, :) = (/ u, v /)
    end do

end subroutine invert_points_on_plane

end module invert