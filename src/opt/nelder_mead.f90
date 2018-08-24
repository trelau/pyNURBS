! Michael F. Hutt
! http://www.mikehutt.com
! Mar. 31, 1998
! $Id: frosen.f90,v 1.4 2007/07/10 12:45:32 mike Exp $
!
! * Copyright (c) 1998-2004 <Michael F. Hutt>
! *
! * Permission is hereby granted, free of charge, to any person obtaining
! * a copy of this software and associated documentation files (the
! * "Software"), to deal in the Software without restriction, including
! * without limitation the rights to use, copy, modify, merge, publish,
! * distribute, sublicense, and/or sell copies of the Software, and to
! * permit persons to whom the Software is furnished to do so, subject to
! * the following conditions:
! *
! * The above copyright notice and this permission notice shall be
! * included in all copies or substantial portions of the Software.
! *
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
! * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
! * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
! Converted program from Fortran 77 to Fortran 90
! This program will attempt to minimize Rosenbrock's function using the 
! Nelder-Mead simplex method. The program was originally developed in C. 
! To be consistent with the way arrays are handled in C, all arrays will
! start from 0.
! compiles with ELF90

! Modifications by trelau - 2016/06/30
! - Converted to double precision.
! - Made interface for external subroutine for objective function.
! - General code cleanup.
! - Removed write statements and print options.
! - Rearranged interface.
! - Made initial variables in/out to retrieve results.
! - Use absolute tolerance for convergence criteria since this is only used
!   to drive the distance between geometry to zero.
!
module nelder_mead
implicit none

private
public :: dummy_simplex
public :: simplex

contains

subroutine dummy_simplex()
    !> Dummy routine for module so f2py works.
end subroutine dummy_simplex

subroutine simplex(fcn, n, x, tol, success)
    ! This is the simplex routine
    
    integer, intent(in) :: n
    double precision, intent(inout), dimension(0:n - 1) :: x
    double precision, intent(in) :: tol
    logical, intent(out) :: success

    ! External subroutine to minimize.
    interface
        subroutine fcn(n, x, f)
            integer, intent(in) :: n
            double precision, intent(in) :: x(n)
            double precision, intent(out) :: f
        end subroutine
    end interface

    ! Define Constants
    integer, parameter :: max_it = 1000
    double precision, parameter :: factor = 1.0d0
    double precision, parameter :: alpha = 1.0d0
    double precision, parameter :: beta = 0.50d0
    double precision, parameter :: gamma = 2.0d0

    ! ======================================================================
    ! Variable Definitions
    ! 
    ! Integer vs = vertex with the smallest value
    ! Integer vh = vertex with next smallest value 
    ! Integer vg = vertex with largest value 
    ! Integer i,j,m,row
    ! Integer k = track the number of function evaluations 
    ! Integer itr = track the number of iterations
    ! real v = holds vertices of simplex 
    ! real pn,qn = values used to create initial simplex 
    ! real f = value of function at each vertex 
    ! real fr = value of function at reflection point 
    ! real fe = value of function at expansion point 
    ! real fc = value of function at contraction point 
    ! real vr = reflection - coordinates 
    ! real ve = expansion - coordinates 
    ! real vc = contraction - coordinates 
    ! real vm = centroid - coordinates 
    ! real xmin
    ! real fsum,favg,s,cent
    ! real vtmp = temporary array passed to fcn
    ! ======================================================================

    integer :: vs, vh, vg
    integer :: i, j, k, itr, m, row
    double precision :: pn, qn, xi
    double precision :: fr, fe, fc
    double precision :: xmin, fsum, favg, cent, s, ds
    double precision, dimension(:,:), allocatable :: v
    double precision, dimension(:), allocatable  :: f
    double precision, dimension(:), allocatable :: vr
    double precision, dimension(:), allocatable :: ve
    double precision, dimension(:), allocatable :: vc
    double precision, dimension(:), allocatable :: vm
    double precision, dimension(:), allocatable :: vtmp

    allocate(v(0:n,0:n - 1))
    allocate(f(0:n))
    allocate(vr(0:n - 1))
    allocate(ve(0:n - 1))
    allocate(vc(0:n - 1))
    allocate(vm(0:n - 1))
    allocate(vtmp(0:n - 1))

    success = .false.
    
    ! create the initial simplex
    ! assume one of the vertices is 0.0

    pn = factor * (dsqrt(dble(n) + 1.0d0) - 1.0d0 + dble(n)) / (dble(n) * dsqrt(2.0d0))
    qn = factor * (dsqrt(dble(n) + 1.0d0) - 1.0d0) / (dble(n) * dsqrt(2.0d0))

    do i = 0, n - 1
        v(0,i) = x(i)
    end do

    do i = 1, n
        do j = 0, n - 1
            if (i - 1 .eq. j) then
                v(i,j) = pn + x(j)
            else
                v(i,j) = qn + x(j)
            end if
        end do
    end do


    ! find the initial function values

    do j = 0, n
        ! put coordinates into single dimension array
        ! to pass it to fcn
        do m = 0, n - 1
            vtmp(m) = v(j,m)
        end do
        call fcn(n, vtmp, xi)
        f(j) = xi
    end do

    k = n + 1

    ! begin main loop of the minimization

    do itr = 1, max_it
        ! find the index of the largest value
        vg = 0
        do j=0, n
            if (f(j) .gt. f(vg)) then
                vg = j
            end if
        end do

        ! find the index of the smallest value
        vs = 0
        do j=0, n
            if (f(j) .lt. f(vs)) then
                vs = j
            end if
        end do

        ! find the index of the second largest value
        vh = vs
        do j=0, n
            if ((f(j) .gt. f(vh)) .and. (f(j) .lt. f(vg))) then
                vh = j
            end if
        end do

        ! calculate the centroid
        do j = 0, n - 1
            cent = 0.0d0
            do m = 0, n
                if (m .ne. vg) then
                    cent = cent + v(m,j)
                end if
            end do
            vm(j) = cent / dble(n)
        end do

        ! reflect vg to new vertex vr
        do j = 0, n - 1
            vr(j) = (1.0d0 + alpha) * vm(j) - alpha * v(vg,j)
        end do
        call fcn(n, vr, fr)
        k = k + 1

        if ((fr .le. f(vh)) .and. (fr .gt. f(vs))) then
            do j = 0, n - 1
                v(vg,j) = vr(j)
            end do
            f(vg) = fr
        end if

        ! investigate a step further in this direction
        if (fr .le. f(vs)) then
            do j = 0, n - 1
                ve(j) = gamma * vr(j) + (1.0d0 - gamma) * vm(j)
            end do
            call fcn(n, ve, fe)
            k = k + 1

            ! by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
            ! takes 62 iterations as opposed to 64. 

            if (fe .lt. fr) then
                do j = 0, n - 1
                    v(vg,j) = ve(j)
                end do
                f(vg) = fe
            else
                do j = 0, n - 1
                    v(vg,j) = vr(j)
                end do
                f(vg) = fr
            end if
        end if

        ! check to see if a contraction is necessary
        if (fr .gt. f(vh)) then
            do j = 0, n - 1
                vc(j) = beta * v(vg,j) + (1.0d0 - beta) * vm(j)
            end do
            call fcn(n, vc, fc)
            k = k + 1
            if (fc .lt. f(vg)) then
                do j = 0, n - 1
                    v(vg,j) = vc(j)
                end do
                f(vg) = fc

            ! at this point the contraction is not successful,
            ! we must halve the distance from vs to all the
            ! vertices of the simplex and then continue.
            ! 10/31/97 - modified C program to account for 
            ! all vertices.
            else
                do row = 0, n
                    if (row .ne. vs) then
                        do j = 0, n - 1
                            v(row,j) = v(vs,j) + (v(row,j) - v(vs,j)) / 2.0d0
                        end do
                    end if
                end do
    
                do m = 0, n - 1
                    vtmp(m) = v(vg,m)
                end do
                call fcn(n, vtmp, xi)
                f(vg) = xi
                k = k + 1

                do m = 0, n - 1
                    vtmp(m) = v(vh,m)
                end do
                call fcn(n, vtmp, xi)
                f(vh) = xi
                k = k + 1
            end if
        end if

        ! test for convergence
        ! fsum = 0.0d0
        ! do j = 0, n
        !     fsum = fsum + f(j)
        ! end do
        ! favg = fsum / (n + 1.0d0)
        ! s = 0.0d0
        ! do j = 0, n
        !     s = s + ((f(j) - favg) ** 2.0d0) / dble(n)
        ! end do
        ! s = dsqrt(s)
        
        s = 0.0d0
        do j = 0, n
            ds = dabs(f(j))
            if (ds .gt. s) then
                s = ds
            end if
        end do
        
        if (s .le. tol) then
            success = .true.
            exit ! Nelder Mead has converged - exit main loop
        end if        
    end do
    
    ! end main loop of the minimization
    ! find the index of the smallest value

    vs = 0
    do j = 0, n
        if (f(j) .lt. f(vs)) then
            vs = j
        end if
    end do

    ! Put results back in x.
    do j = 0, n - 1
        x(j) = v(vs, j)
    end do

end subroutine simplex

end module nelder_mead