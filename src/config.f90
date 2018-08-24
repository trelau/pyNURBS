module config
implicit none

! Print warnings.
logical :: warnings

! Tolerance values.
double precision :: gtol, ptol, atol, stol, ftol

contains

subroutine dummy_config()
    !> Dummy routine for module so f2py works.
end subroutine dummy_config

end module config