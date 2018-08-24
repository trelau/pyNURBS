! generic_list.f90 -- A Generic Linked List Implementation in Fortran 95
!
! Copyright (C) 2009, 2012 Jason R. Blevins
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.

! Revision History:
!
! 1. July 21, 2012: In the list_free subroutine, line 11 should read
! nullify(current%data) instead of nullify(self%data).  Thanks to
! Michael Quinlan.

module generic_list
  implicit none

  private
  public :: dummy_list
  public :: list_node_t, list_data
  public :: list_init, list_free
  public :: list_insert, list_put, list_get, list_next

  ! A public variable used as a MOLD for transfer()
  integer, dimension(:), allocatable :: list_data

  ! Linked list node
  type :: list_node_t
     private
     integer, dimension(:), pointer :: data => null()
     type(list_node_t), pointer :: next => null()
  end type list_node_t

contains

    subroutine dummy_list()
        !> Dummy routine for module so f2py works.
    end subroutine dummy_list

  ! Initialize a head node SELF and optionally store the provided DATA.
  subroutine list_init(self, data)
    type(list_node_t), pointer :: self
    integer, dimension(:), intent(in), optional :: data

    allocate(self)
    nullify(self%next)

    if (present(data)) then
       allocate(self%data(size(data)))
       self%data = data
    else
       nullify(self%data)
    end if
  end subroutine list_init

  ! Free the entire list and all data, beginning at SELF
  subroutine list_free(self)
    type(list_node_t), pointer :: self
    type(list_node_t), pointer :: current
    type(list_node_t), pointer :: next

    current => self
    do while (associated(current))
       next => current%next
       if (associated(current%data)) then
          deallocate(current%data)
          nullify(current%data)
       end if
       deallocate(current)
       nullify(current)
       current => next
    end do
  end subroutine list_free

  ! Insert a list node after SELF containing DATA (optional)
  subroutine list_insert(self, data)
    type(list_node_t), pointer :: self
    integer, dimension(:), intent(in), optional :: data
    type(list_node_t), pointer :: next

    allocate(next)

    if (present(data)) then
       allocate(next%data(size(data)))
       next%data = data
    else
       nullify(next%data)
    end if

    next%next => self%next
    self%next => next
  end subroutine list_insert

  ! Store the encoded DATA in list node SELF
  subroutine list_put(self, data)
    type(list_node_t), pointer :: self
    integer, dimension(:), intent(in) :: data

    if (associated(self%data)) then
       deallocate(self%data)
       nullify(self%data)
    end if
    self%data = data
  end subroutine list_put

  ! Return the DATA stored in the node SELF
  function list_get(self) result(data)
    type(list_node_t), pointer :: self
    integer, dimension(:), pointer :: data
    data => self%data
  end function list_get

  ! Return the next node after SELF
  function list_next(self)
    type(list_node_t), pointer :: self
    type(list_node_t), pointer :: list_next
    list_next => self%next
  end function list_next

end module generic_list