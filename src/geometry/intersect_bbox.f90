module intersect_bbox
use math, only: norm
implicit none

private
public :: dummy_bi
public :: bboxes_intersect, bbox_intersects_plane

contains

subroutine dummy_bi()
    !> Dummy routine for f2py.
end subroutine dummy_bi

function bboxes_intersect(k, bbox1, bbox2, tol) result(intersects)
    !> Test if the bounding boxes intersect.
    !> k - Dimensions of bounding box.
    !> bbox1 - Bounding box 1.
    !> bbox2 - Bounding box 2.
    !> tol - Tolerance for bounding boxes.
    !> intersects - True is boxes intersect, False if not.
    
    ! Input
    integer, intent(in) :: k
    double precision, intent(in) :: tol
    double precision, intent(in) :: bbox1(k, 2), bbox2(k, 2)
    
    ! Output
    logical :: intersects
    
    ! Working
    integer :: i
    double precision :: bbox1_(k, 2), bbox2_(k, 2)
    
    ! Adjust bounding boxes by the tolerance.
    do i = 1, k
        bbox1_(i, 1) = bbox1(i, 1) - tol
        bbox1_(i, 2) = bbox1(i, 2) + tol
        bbox2_(i, 1) = bbox2(i, 1) - tol
        bbox2_(i, 2) = bbox2(i, 2) + tol
    end do

    intersects = .true.
    do i = 1, k
        if ((bbox1_(i, 1) .gt. bbox2_(i, 2)) .or. (bbox2_(i, 1) .gt. bbox1_(i, 2))) then
            intersects = .false.
            return
        end if
    end do
    
end function bboxes_intersect

function bbox_intersects_plane(bbox, p0, pnorm, tol) result(intersects)
    !> Test if the bounding box intersects a plane.
    !> bbox - Bounding box.
    !> p0 - Origin point of plane.
    !> pnorm - Unit normal vector of plane.
    !> tol - Tolerance for bounding boxes.
    !> intersects - True if box is intersected by plane, False if not.
    
    ! Input
    double precision, intent(in) :: tol
    double precision, intent(in) :: bbox(0:2, 0:1)
    double precision, intent(in) :: p0(0:2)
    double precision, intent(in) :: pnorm(0:2)
    
    ! Output
    logical :: intersects
    
    ! Working
    integer :: i
    double precision :: d, r
    double precision :: cg(0:2), v(0:2)
    double precision :: bbox_(0:2, 0:1)
    
    ! Adjust bounding boxes by the tolerance.
    do i = 0, 2
        bbox_(i, 0) = bbox(i, 0) - tol
        bbox_(i, 1) = bbox(i, 1) + tol
        ! Centroid and diameter vector.
        cg(i) = (bbox_(i, 0) + bbox_(i, 1)) / 2.0d0
        v(i) = bbox_(i, 1) - bbox_(i, 0)
    end do
    
    ! Radius of bounding box.
    r = norm(v) / 2.0d0
    
    ! Distance from cg to plane.
    d = dabs(dot_product(pnorm, cg - p0))
    
    intersects = .false.
    if (d .le. r) then
        intersects = .true.
    end if

end function bbox_intersects_plane

end module intersect_bbox