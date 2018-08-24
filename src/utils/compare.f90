module compare
implicit none

private
public :: dummy_compare
public :: EqualTo, GreaterThan, LessThan, GreaterThanEq, LessThanEq
public :: check_bounds

contains

subroutine dummy_compare()
    !> Dummy routine for module so f2py works.
end subroutine dummy_compare

function EqualTo(a, b, abs_tol) result(equal_to)
    !> Check if two floats are almost equal.
    !> a - First float.
    !> b - Other float.
    !> abs_tol - Absolute toleranace.
    
    ! Input
    double precision, intent(in) :: a, b
    double precision, intent(in), optional :: abs_tol
    
    ! Output
    logical :: equal_to
    
    ! Working
    double precision :: atol, diff
    
    ! In case they are exact.
    if (a .eq. b) then
        equal_to = .true.
        return
    end if
    
    ! Set defaults.
    atol = 1.0d-8
    if (present(abs_tol)) then
        atol = abs_tol
    end if
    
    diff = dabs(b - a)
    ! Absolute comparison.
    if (diff .le. atol) then
        equal_to = .true.
        return
    end if
    equal_to = .false.
end function EqualTo

function GreaterThan(a, b, abs_tol) result(greater_than)
    !> Check if float is almost greater than other. 
    !> a - First float.
    !> b - Other float.
    !> abs_tol - Absolute toleranace.
    
    ! Input
    double precision, intent(in) :: a, b
    double precision, intent(in), optional :: abs_tol
    
    ! Output
    logical :: greater_than
    
    ! Working
    double precision :: atol, diff
    
    ! Set defaults.
    atol = 1.0d-8
    if (present(abs_tol)) then
        atol = abs_tol
    end if
    
    ! Check almost equal first.
    if (EqualTo(a, b, atol)) then
        greater_than = .false.
        return
    end if
    if (a .gt. b) then
        greater_than = .true.
    else
        greater_than = .false.
    end if
end function GreaterThan

function LessThan(a, b, abs_tol) result(less_than)
    !> Check if float is almost less than other. 
    !> a - First float.
    !> b - Other float.
    !> abs_tol - Absolute toleranace.
    
    ! Input
    double precision, intent(in) :: a, b
    double precision, intent(in), optional :: abs_tol
    
    ! Output
    logical :: less_than
    
    ! Working
    double precision :: atol, diff
    
    ! Set defaults.
    atol = 1.0d-8
    if (present(abs_tol)) then
        atol = abs_tol
    end if
    
    ! Check almost equal first.
    if (EqualTo(a, b, atol)) then
        less_than = .false.
        return
    end if
    if (a .lt. b) then
        less_than = .true.
    else
        less_than = .false.
    end if
end function LessThan

function GreaterThanEq(a, b, abs_tol) result(greater_than_eq)
    !> Check if float is almost greater than or almost equal to other. 
    !> a - First float.
    !> b - Other float.
    !> rel_tol - Relative tolerance.
    !> abs_tol - Absolute toleranace.
    
    ! Input
    double precision, intent(in) :: a, b
    double precision, intent(in), optional :: abs_tol
    
    ! Output
    logical :: greater_than_eq
    
    ! Working
    double precision :: atol, diff
    
    ! Set defaults.
    atol = 1.0d-8
    if (present(abs_tol)) then
        atol = abs_tol
    end if
    
    ! Check almost equal first.
    if (EqualTo(a, b, atol)) then
        greater_than_eq = .true.
        return
    end if
    if (a .gt. b) then
        greater_than_eq = .true.
    else
        greater_than_eq = .false.
    end if
end function GreaterThanEq

function LessThanEq(a, b, abs_tol) result(less_than_eq)
    !> Check if float is almost less than or equal to other. 
    !> a - First float.
    !> b - Other float.
    !> abs_tol - Absolute toleranace.
    
    ! Input
    double precision, intent(in) :: a, b
    double precision, intent(in), optional :: abs_tol
    
    ! Output
    logical :: less_than_eq
    
    ! Working
    double precision :: atol, diff
    
    ! Set defaults.
    atol = 1.0d-8
    if (present(abs_tol)) then
        atol = abs_tol
    end if
    
    ! Check almost equal first.
    if (EqualTo(a, b, atol)) then
        less_than_eq = .true.
        return
    end if
    if (a .lt. b) then
        less_than_eq = .true.
    else
        less_than_eq = .false.
    end if
end function LessThanEq

function check_bounds(x, a, b) result(xout)
    !> Check that a float is between the bounds and return an updated value
    !> if not.
    !> x - Float to check.
    !> a - Lower bound.
    !> b - Upper bound.
    
    ! Input
    double precision, intent(in) :: x, a, b
    
    ! Output
    double precision :: xout
    
    xout = x
    
    if ( x .lt. a) then
        xout = a
        return
    end if
    
    if ( x .gt. b) then
        xout = b
        return
    end if

end function check_bounds


end module compare