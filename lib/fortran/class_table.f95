MODULE class_table

   use linalg

   implicit none

   type table

      integer :: label
      integer :: inner_count = 1
      integer :: outer_count = 1
      integer :: order = 1
      real(kind=8), dimension(:, :), allocatable :: table_data
      integer, dimension(:), allocatable :: axes_lengths, alphas
      real(kind=8), dimension(:, :), allocatable :: matrix, vec_b, vec_a
      integer, dimension(:, :), allocatable :: index_matrix, indices
      character(len=3), dimension(:), allocatable :: scalings
      real(kind=8), dimension(:), allocatable :: starts, ends
      integer :: n_in, n_out
      integer, dimension(:), allocatable :: which_out
      real(kind=8), dimension(:), allocatable :: results

   end type table

   type table_list

      type(table), dimension(:), allocatable :: objects
      integer, dimension(:), allocatable :: n_in_list

   end type table_list

contains

!#######################################################################
   SUBROUTINE init(self, label)

      implicit none

      type(table), intent(inout) :: self
      integer, intent(in) :: label

      self%label = label

   END SUBROUTINE init

   SUBROUTINE load(self, file_name)

      implicit none

      type(table), intent(inout) :: self
      character(len=*), intent(in) :: file_name
      character(len=30) :: str
      integer :: n, u, i, j, ipos
      integer :: n_input_params
      integer :: a, n_data_points, n_output_params
      character(len=2) :: c

!write(*,*) 'Reading data table file ', file_name
      open (20, file=file_name, status="old", action="read")

!Extract meta data

!axes dimesions
      read (20, *) str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) n_input_params

!For the axes_lengths also the number of output parameters must be
!stored. The array hence requires one additional element
      allocate (self%axes_lengths(n_input_params + 1))

!For these things only the values for the input parameters are required
      allocate (self%alphas(n_input_params))
      allocate (self%starts(n_input_params))
      allocate (self%ends(n_input_params))
      allocate (self%scalings(n_input_params))

!once the axes lengths are known, allocate them to the arrays
      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) self%axes_lengths

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) self%alphas

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) self%scalings

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) self%starts

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) self%ends

!The total number of parameters minus the input parameters gives the
!number of output parameters
      n_output_params = self%axes_lengths(size(self%axes_lengths)) - n_input_params

!print *, n_input_params, ' input params found'
!print *, 'the dimensions of the axes are:', axes_lengths
!print *, 'the axis alphas are:', alphas
!print *, 'starts/ends are:', starts,'/', ends
!print *, 'axis scalings are:', scalings

!Compute the total number of grid points
      n_data_points = 1
      do i = 1, size(self%axes_lengths) - 1
         n_data_points = n_data_points*self%axes_lengths(i)
      end do

!The data array used for table interpolation is 2d regardless of the
!dimensions of the input and output vector. The first axis contains
!all grid points. The second axis contains all parameter values
!for the given gridpoint. A parameter value q at given grid point i can
!then be accessed by table_data(i,q)
      allocate (self%table_data(n_data_points, self%axes_lengths(n_input_params + 1)))

      self%n_in = n_input_params
      self%n_out = n_output_params

!print *, 'total number of grid points:', n_data_points
!print *, 'total number of parameters:', axes_lengths(n_input_params+1)

!Assing the parameter values for each grid point from the input data
!This depends on the dimension of the input data
      do i = 1, self%axes_lengths(n_input_params + 1)
         read (20, *) self%table_data(:, i)
      end do

1     close (20)

   END SUBROUTINE load

   SUBROUTINE convert_indices(self, inds, n, s)

      type(table), intent(inout) :: self
      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: inds
      integer i, j, tot, ss
      integer, intent(out) :: s

      tot = 1
      do i = 1, n
         tot = tot*self%axes_lengths(i)
      end do

      s = 0
      do i = 1, n - 1
         !print *, 'i=', i, inds(i)-1
         ss = (inds(i) - 1)
         do j = 1, n - i
            !print *, 'j=', j, dims(j+i)
            ss = ss*self%axes_lengths(j + i)
         end do
         s = s + ss
      end do

      s = s + inds(n)

   END SUBROUTINE convert_indices

   RECURSIVE SUBROUTINE for_recursive(self, n, c, iter_list, ranges, &
                                      sub, vals)

      type(table), intent(inout) :: self
      integer, intent(in) :: n, c
      external :: sub
      real(kind=8), dimension(n), intent(inout), optional :: vals
      integer, dimension(n), intent(inout), optional :: iter_list
      integer, dimension(n) :: iter_list_dummy
      integer, dimension(n), intent(in) :: ranges
      integer :: i, j

!Check if iter_list is given and if not
!assign zeros to each element
      if (present(iter_list)) then
         iter_list_dummy = iter_list

      else
!print *, 'No iter_list given'
         do i = 1, n
            iter_list_dummy(i) = 1
         end do
      end if

      if (c == n) then
         do i = 1, ranges(c)
            iter_list_dummy(c) = i
            !print *, 'iter_list after=', iter_list_dummy
            if (present(vals)) then
               !print *, 'additional args found'
               call sub(self, n, iter_list_dummy, vals)
            else
               !print *, 'no additional args found'
               call sub(self, n, iter_list_dummy)
            end if
         end do

      else
         do i = 1, ranges(c)
            iter_list_dummy(c) = i
            !print *, 'iter_list after=', iter_list_dummy
            if (present(vals)) then
               !print *, 'additional args found'
               call for_recursive(self=self, c=c + 1, iter_list=iter_list_dummy, &
                                  ranges=ranges, n=n, sub=sub, &
                                  vals=vals)
            else
               !print *, 'no additional args found'
               call for_recursive(self=self, c=c + 1, iter_list=iter_list_dummy, &
                                  ranges=ranges, n=n, sub=sub)
            end if

         end do
      end if

   END SUBROUTINE for_recursive

   SUBROUTINE get_indices(self, x, which, n, o, inds)

      type(table), intent(inout) :: self
      real(kind=8), intent(in) :: x
      real(kind=8) :: a, left_val, right_val, delta
      integer :: left_ind, right_ind
      integer, intent(in) :: which, n, o
      integer, dimension(self%order + 1) :: relative_coordinates
      integer, dimension(o + 1), intent(out) :: inds
      integer :: i, j, start, exponent, offset

      start = -self%order
      do i = 1, self%order + 1
         relative_coordinates(i) = start + i
      end do

      if (self%scalings(which) == 'log') then
         exponent = int(log10(x))
         a = 10.0**(exponent - self%alphas(which))
         left_val = int(x/a)*a
         right_val = int(x/a + 1.0d0)*a

         left_ind = int(left_val/a - 10.0**self%alphas(which))

         offset = int(left_val/a - 10**self%alphas(which)) + 1

         do i = 1, self%order + 1
            inds(i) = 9*10**self%alphas(which)*(exponent - self%starts(which)) + &
                      offset + relative_coordinates(i)
         end do

      elseif (self%scalings(which) == 'lin') then
         delta = (self%ends(which) - self%starts(which))/(self%axes_lengths(which) - 1)

         do i = 1, self%order + 1
            inds(i) = int((x - self%starts(which))/delta) + relative_coordinates(i) + 1
         end do

      end if

   END SUBROUTINE get_indices

   SUBROUTINE index_row(self, n, iter_list)

      type(table), intent(inout) :: self
      integer, dimension(n), intent(in) :: iter_list
      integer, intent(in) :: n
      integer :: i

      do i = 1, n
         self%index_matrix(self%inner_count, i) = self%indices(i, iter_list(i))
      end do

      self%inner_count = self%inner_count + 1

   END SUBROUTINE index_row

   SUBROUTINE construct_index_matrix(self, n, vals)

      type(table), intent(inout) :: self
      integer, intent(in) :: n
      integer, dimension(n) :: iter_ranges
      real(kind=8), dimension(n), intent(in) :: vals
      integer :: i

      do i = 1, n
         call get_indices(self=self, inds=self%indices(i, :), x=vals(i), which=i, n=n, o=self%order)
         print *, 'ind =', self%indices(i, :)
      end do

      do i = 1, n
         iter_ranges(i) = self%order + 1
      end do

      self%inner_count = 1
      call for_recursive(self=self, n=n, c=1, ranges=iter_ranges, sub=index_row)
      self%inner_count = 1

   END SUBROUTINE construct_index_matrix

   SUBROUTINE partial_row(self, n, iter_list, vals)

      type(table), intent(inout) :: self
      integer, dimension(n), intent(in) :: iter_list
      real(kind=8), dimension(n), intent(in) :: vals
      real(kind=8) :: res
      integer, intent(in) :: n
      integer :: i

      res = 1.0
      do i = 1, n
         res = res*vals(i)**(iter_list(i) - 1)
      end do

!print *, 'row_index =', outer_count, 'inner_count =', inner_count
!print *, 'iter_list', iter_list,' vals', vals,' res =', res
      self%matrix(self%outer_count, self%inner_count) = res

      self%inner_count = self%inner_count + 1

   END SUBROUTINE partial_row

   SUBROUTINE row(self, n, iter_list)

      type(table), intent(inout) :: self
      integer, intent(in) :: n
      real(kind=8), dimension(n) :: vals_in
      integer, dimension(n), intent(in) :: iter_list
      integer, dimension(n) :: iter_ranges
      integer :: i, j, ind

      call convert_indices(self=self, inds=self%index_matrix(self%inner_count, :), &
                           n=n, s=ind)

!call mat_print('index matrix', self%index_matrix)

      if (ind < 1) then
         ind = 1
      end if

      do i = 1, n
         iter_ranges(i) = self%order + 1
         vals_in(i) = self%table_data(ind, i)
      end do

      do i = 1, self%n_out
         print *, 'ind/which =', ind, self%which_out(i)
         self%vec_b(self%inner_count, i) = self%table_data(ind, self%which_out(i))
      end do
      print *, 'vals_in =', vals_in(:)
      self%outer_count = self%inner_count
      self%inner_count = 1

      call for_recursive(self=self, n=n, c=1, ranges=iter_ranges, sub=partial_row, &
                         vals=vals_in)

      self%inner_count = self%outer_count

      self%inner_count = self%inner_count + 1

   END SUBROUTINE row

   SUBROUTINE construct_matrix(self, n)

      type(table), intent(inout) :: self
      integer, intent(in) :: n
      integer, dimension(n) :: iter_ranges
!integer, dimension(n_out) :: which
      integer :: i, j

      do i = 1, (self%order + 1)**n
         do j = 1, self%n_out
            self%vec_a(i, j) = 0.0d0
         end do
      end do
!allocate(which_out(n_out))

!do i=1, n_out
      ! which_out(i) = which(i)
!enddo
      do i = 1, n
         iter_ranges(i) = self%order + 1
      end do

      self%inner_count = 1

      call for_recursive(self=self, sub=row, ranges=iter_ranges, n=n, c=1)

      self%inner_count = 1

   END SUBROUTINE construct_matrix

   SUBROUTINE compute_coeffs(self, solver, n)

      type(table), intent(inout) :: self
      integer, intent(in) :: n
      integer, intent(inout), optional :: solver
      integer :: i

      if (.not. present(solver)) then
         solver = 0
      end if

!Solve system via inverse matrix
      if (solver == 0) then
         call solve_linear_system(n=(self%order + 1)**n, l=self%n_out, &
                                  matrix=self%matrix, b=self%vec_b, x=self%vec_a)

!Solve system via LU decomposition
      elseif (solver == 1) then
         print *, 'WARNING: LU decomposition is currently not available!'
         !call solve_LU(a=self%matrix, b=self%vec_b, x=self%vec_a)

!Solve system via classical Gauss algorithm
      elseif (solver == 2) then

         call gauss_elimination(a=self%matrix, sol=self%vec_a, b=self%vec_b)

      end if

   END SUBROUTINE compute_coeffs

   SUBROUTINE compute_result(self, n, iter_list, vals)

      type(table), intent(inout) :: self
      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: iter_list
      real(kind=8), dimension(self%n_out) :: prod
      real(kind=8) :: delta
      real(kind=8), dimension(n), intent(in) :: vals
      integer :: i, j

      do i = 1, self%n_out
         prod(i) = 1.0d0
      end do

      do i = 1, self%n_out
         do j = 1, n
            prod(i) = prod(i)*vals(j)**(iter_list(j) - 1)
         end do

         delta = self%vec_a(self%inner_count, i)*prod(i)
         self%results(i) = self%results(i) + delta
      end do

      self%inner_count = self%inner_count + 1

   END SUBROUTINE compute_result

   SUBROUTINE interpolate(self, n_params, ord, which, vals, ll)

      implicit none

      type(table), intent(inout) :: self
      integer, intent(in) :: n_params
      integer, intent(in), optional :: ll
      integer, intent(inout), optional :: ord
      integer, dimension(n_params), intent(in) :: which
      integer, dimension(self%n_in) :: iter_ranges
      real(kind=8), dimension(self%n_in), intent(inout) :: vals
      integer :: i, j, solver

      solver = 2

!self%n_out = n_params
!self%n_in = n

!Temp and Pres are restricted to 10 K and 1000 Pa
      vals(1) = max(2d0, vals(1))
      vals(2) = max(2d3, vals(2))

      if (present(ord)) then
         self%order = ord
      end if

      do i = 1, self%n_out
         self%which_out(i) = which(i) + self%n_in
         self%results(i) = 0.0d0
      end do

      do i = 1, self%n_in
         iter_ranges(i) = self%order + 1
      end do

      call construct_index_matrix(self, self%n_in, vals)

      call construct_matrix(self, n=self%n_in)

      call compute_coeffs(self, n=self%n_in, solver=solver)

      self%inner_count = 1

      call for_recursive(self=self, sub=compute_result, c=1, n=self%n_in, &
                         vals=vals, ranges=iter_ranges)

      self%inner_count = 1

   END SUBROUTINE interpolate

END MODULE class_table

