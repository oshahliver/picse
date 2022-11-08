MODULE class_dict

implicit none

type dict

integer :: n
character(len=30), dimension(:), allocatable :: names
real(kind=8), dimension(:), allocatable :: real_vals
integer, dimension(:), allocatable :: int_vals
real(kind=8), dimension(:,:), allocatable :: real_arr
integer, dimension(:,:), allocatable :: int_arr

end type dict

contains

!#######################################################################
SUBROUTINE init_dict(self, n, n1, n2)

type(dict), intent(inout) :: self
integer, intent(in) :: n, n1, n2

if(.not.allocated(self%real_arr))then
	allocate(self%real_vals(n))
	allocate(self%int_vals(n))
	allocate(self%int_arr(n,n1))
	allocate(self%real_arr(n,n2))
endif

END SUBROUTINE init_dict

!#######################################################################
SUBROUTINE access_key(self, key, real_out, int_out)

type(dict), intent(inout) :: self
character(len=30), intent(in) :: key
integer, intent(out), optional :: int_out
real(kind=8), intent(out), optional :: real_out
integer :: i

do i=1, self%n
  if(key==self%names(i))then
    int_out = self%int_vals(i)
    real_out = self%real_vals(i)
  endif
enddo

END SUBROUTINE access_key

END MODULE class_dict

!#######################################################################
!#######################################################################
MODULE class_weird_array

implicit none

type alloc_array

real(kind=8), dimension(:), allocatable :: real_array
integer, dimension(:), allocatable :: int_array

end type alloc_array


type weird_array

integer :: ndim
integer, dimension(:), allocatable :: axes_lengths
type(alloc_array), dimension(:), allocatable :: axes

end type weird_array

contains

SUBROUTINE init_alloc_array(self, length)

type(alloc_array), intent(inout) :: self
integer, intent(in) :: length

allocate(self%int_array(length))
allocate(self%real_array(length))

END SUBROUTINE init_alloc_array

!#######################################################################
SUBROUTINE init_weird_array(self, ndim, axes_lengths)

type(weird_array), intent(inout) :: self
integer, intent(in) :: ndim
integer, dimension(ndim), intent(in) :: axes_lengths
type(alloc_array), dimension(ndim) :: axes

integer :: i, j

self%ndim = ndim
self%axes_lengths = axes_lengths
self%axes = axes


do i=1, ndim
  allocate(self%axes(i)%int_array(self%axes_lengths(i)))
  allocate(self%axes(i)%real_array(self%axes_lengths(i)))
enddo

END SUBROUTINE init_weird_array

END MODULE class_weird_array
