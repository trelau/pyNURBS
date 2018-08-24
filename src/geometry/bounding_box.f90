module bounding_box
use geom_utils
implicit none

private
public :: dummy_bbox
public :: curve_bbox, surface_bbox

contains

subroutine dummy_bbox()
    !> Dummy routine for f2py.
end subroutine dummy_bbox

function curve_bbox(n, cpw) result(bbox)
    !> Generate a 3-D axis-aligned bounding box from a curve.
    !> n - Number of control points - 1.
    !> cpw - Control points.
    !> bbox - Bouning box array where bbox(d, 0) is the minimum value
    !> for the d-dimension and bbox(d, 1) is the maximum value.
    
    ! Input
    integer, intent(in) :: n
    double precision, intent(in) :: cpw(0:n, 0:3)
    
    ! Output
    double precision :: bbox(0:2, 0:1)
    
    ! Working
    integer :: i, j
    double precision :: cp(0:n, 0:2), w(0:n)
    
    call dehomogenize_array1d(n, cpw, cp, w)
    bbox(0, :) = cp(0, 0)
    bbox(1, :) = cp(0, 1)
    bbox(2, :) = cp(0, 2)
    do i = 0, n
        do j = 0, 2
            if (cp(i, j) .lt. bbox(j, 0)) then
                bbox(j, 0) = cp(i, j)
            end if
            if (cp(i, j) .gt. bbox(j, 1)) then
                bbox(j, 1) = cp(i, j)
            end if
        end do
    end do
    
end function curve_bbox

function surface_bbox(n, m, cpw) result(bbox)
    !> Generate a 3-D axis-aligned bounding box from a surface.
    !> n - Number of control points - 1 in u-direction.
    !> m - Number of control points - 1 in v-direction.
    !> cpw - Control points.
    !> bbox - Bouning box array where bbox(d, 0) is the minimum value
    !> for the d-dimension and bbox(d, 1) is the maximum value.
    
    ! Input
    integer, intent(in) :: n, m
    double precision, intent(in) :: cpw(0:n, 0:m, 0:3)
    
    ! Output
    double precision :: bbox(0:2, 0:1)
    
    ! Working
    integer :: i, j, k
    double precision :: cp(0:n, 0:m, 0:2)
    double precision :: w(0:n, 0:m)
    
    call dehomogenize_array2d(n, m, cpw, cp, w)
    do k = 0, 2
    bbox(k, :) = cp(0, 0, k)
    end do
    do i = 0, n
        do j = 0, m    
            do k = 0, 2
                if (cp(i, j, k) .lt. bbox(k, 0)) then
                    bbox(k, 0) = cp(i, j, k)
                end if
                if (cp(i, j, k) .gt. bbox(k, 1)) then
                    bbox(k, 1) = cp(i, j, k)
                end if
            end do
        end do
    end do
    
end function surface_bbox

end module bounding_box