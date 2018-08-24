module geom_results
implicit none

! Global result variables.
integer :: npts_, nverts_, ntri_, ncrvs_
integer, allocatable :: triangles_(:, :)
integer, allocatable :: crv_size_(:)
integer, allocatable :: crv_ids_(:, :)
double precision, allocatable :: verts_(:, :)
double precision, allocatable :: points_(:, :)
double precision, allocatable :: params1d_(:)
double precision, allocatable :: params2d_(:, :)
double precision, allocatable :: dproj_(:)
double precision, allocatable :: ssi_params1_(:, :)
double precision, allocatable :: ssi_params2_(:, :)
double precision, allocatable :: cinterp_(:, :)

contains

subroutine reset_results()
    !> Reset and deallocate all geometry results arrays.
    
    npts_ = 0
    nverts_ = 0
    ntri_ = 0
    ncrvs_ = 0
    
    if (allocated(triangles_)) then
        deallocate(triangles_)
    end if
    
    if (allocated(crv_size_)) then
        deallocate(crv_size_)
    end if
    
    if (allocated(crv_ids_)) then
        deallocate(crv_ids_)
    end if
    
    if (allocated(verts_)) then
        deallocate(verts_)
    end if
    
    if (allocated(points_)) then
        deallocate(points_)
    end if
    
    if (allocated(params1d_)) then
        deallocate(params1d_)
    end if
    
    if (allocated(params2d_)) then
        deallocate(params2d_)
    end if
    
    if (allocated(dproj_)) then
        deallocate(dproj_)
    end if
    
    if (allocated(ssi_params1_)) then
        deallocate(ssi_params1_)
    end if
    
    if (allocated(ssi_params2_)) then
        deallocate(ssi_params2_)
    end if
    
    if (allocated(cinterp_)) then
        deallocate(cinterp_)
    end if

end subroutine reset_results

subroutine set_empty_results()
    !> Allocate result arrays but set empty.
    
    call reset_results()
    
    allocate(triangles_(1, 3))
    triangles_(:, :) = 0
    
    allocate(crv_size_(1))
    crv_size_(:) = 0
    
    allocate(crv_ids_(1, 1))
    crv_ids_(:, :) = 0
    
    allocate(verts_(1, 3))
    verts_(:, :) = 0.0d0
    
    allocate(points_(1, 3))
    points_(:, :) = 0.0d0
    
    allocate(params1d_(1))
    params1d_(:) = 0.0d0
    
    allocate(params2d_(1, 1))
    params2d_(:, :) = 0.0d0
    
    allocate(dproj_(1))
    dproj_(:) = 0.0d0
    
    allocate(ssi_params1_(1, 1))
    ssi_params1_(:, :) = 0.0d0
    
    allocate(ssi_params2_(1, 1))
    ssi_params2_(:, :) = 0.0d0
    
    allocate(cinterp_(1, 2))
    cinterp_(:, :) = 0.0d0
    
end subroutine set_empty_results

end module geom_results