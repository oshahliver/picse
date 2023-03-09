
MODULE table_tools

   use linalg

   implicit none

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
   character(len=30), dimension(:), allocatable :: table_file_names

contains

   SUBROUTINE read_data(file_name)

      implicit none

      character(len=*), intent(in) :: file_name
      character(len=30) :: str
      integer :: n, u, i, j, ipos
      integer :: n_input_params
      integer :: a, n_data_points, n_output_params
      character(len=2) :: c

      write (*, *) 'Reading data table file ', file_name
      open (20, file=file_name, status="old", action="read")

!Extract meta data

!axes dimesions
      read (20, *) str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) n_input_params

!For the axes_lengths also the number of output parameters must be
!stored. The array hence requires one additional element
      allocate (axes_lengths(n_input_params + 1))

!For these things only the values for the input parameters are required
      allocate (alphas(n_input_params))
      allocate (starts(n_input_params))
      allocate (ends(n_input_params))
      allocate (scalings(n_input_params))

!once the axes lengths are known, allocate them to the arrays
      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) axes_lengths

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) alphas

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) scalings

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) starts

      read (20, '(a)') str
      ipos = scan(str, '=', back=.true.)
      read (str(ipos + 1:), *) ends

!print *, n_input_params, ' input params found'
!print *, 'the dimensions of the axes are:', axes_lengths
!print *, 'the axis alphas are:', alphas
!print *, 'starts/ends are:', starts,'/', ends
!print *, 'axis scalings are:', scalings

!Compute the total number of grid points
      n_data_points = 1
      do i = 1, size(axes_lengths) - 1
         n_data_points = n_data_points*axes_lengths(i)
      end do

!The data array used for table interpolation is 2d regardless of the
!dimensions of the input and output vector. The first axis contains
!all grid points. The second axis contains all parameter values
!for the given gridpoint. A parameter value q at given grid point i can
!then be accessed by table_data(i,q)
      allocate (table_data(n_data_points, axes_lengths(n_input_params + 1)))

!print *, 'total number of grid points:', n_data_points
!print *, 'total number of parameters:', axes_lengths(n_input_params+1)

!Assing the parameter values for each grid point from the input data
!This depends on the dimension of the input data
      do i = 1, axes_lengths(n_input_params + 1)
         read (20, *) table_data(:, i)
      end do

1     close (20)

   END SUBROUTINE read_data

   SUBROUTINE convert_indices(inds, dims, n, s)

      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: dims, inds
      integer i, j, tot, ss
      integer, intent(out) :: s

      tot = 1
      do i = 1, n
         tot = tot*axes_lengths(i)
      end do

      s = 0
      do i = 1, n - 1
         !print *, 'i=', i, inds(i)-1
         ss = (inds(i) - 1)
         do j = 1, n - i
            !print *, 'j=', j, dims(j+i)
            ss = ss*axes_lengths(j + i)
         end do
         s = s + ss
      end do

      s = s + inds(n)

   END SUBROUTINE convert_indices

   RECURSIVE SUBROUTINE for_recursive(n, c, iter_list, ranges, &
                                      sub, vals)

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
               call sub(n, iter_list_dummy, vals)
            else
               !print *, 'no additional args found'
               call sub(n, iter_list_dummy)
            end if
         end do

      else
         do i = 1, ranges(c)
            iter_list_dummy(c) = i
            !print *, 'iter_list after=', iter_list_dummy
            if (present(vals)) then
               !print *, 'additional args found'
               call for_recursive(c=c + 1, iter_list=iter_list_dummy, &
                                  ranges=ranges, n=n, sub=sub, &
                                  vals=vals)
            else
               !print *, 'no additional args found'
               call for_recursive(c=c + 1, iter_list=iter_list_dummy, &
                                  ranges=ranges, n=n, sub=sub)
            end if

         end do
      end if

   END SUBROUTINE for_recursive

   SUBROUTINE get_indices(x, which, n, o, inds)

      real(kind=8), intent(in) :: x
      real(kind=8) :: a, left_val, right_val, delta
      integer :: left_ind, right_ind
      integer, intent(in) :: which, n, o
      integer, dimension(order + 1) :: relative_coordinates
      integer, dimension(o + 1), intent(out) :: inds
      integer :: i, j, start, exponent

      start = -order
      do i = 1, order + 1
         relative_coordinates(i) = start + i
      end do

      if (scalings(which) == 'log') then
         exponent = int(log10(x))
         a = 10.0**(exponent - alphas(which))
         left_val = int(x/a)*a
         right_val = int(x/a + 1.0d0)*a

         left_ind = int(left_val/a - 10.0**alphas(which))

         do i = 1, order + 1
            inds(i) = left_ind + relative_coordinates(i) + 1
         end do

      elseif (scalings(which) == 'lin') then
         delta = (ends(which) - starts(which))/(axes_lengths(which) - 1)

         do i = 1, order + 1
            inds(i) = int((x - starts(which))/delta) + relative_coordinates(i) + 1
         end do

      end if

   END SUBROUTINE get_indices

   SUBROUTINE index_row(n, iter_list)

      integer, dimension(n), intent(in) :: iter_list
      integer, intent(in) :: n
      integer :: i

      do i = 1, n
         index_matrix(inner_count, i) = indices(i, iter_list(i))
      end do

      inner_count = inner_count + 1

   END SUBROUTINE index_row

   SUBROUTINE construct_index_matrix(n, vals)

      integer, intent(in) :: n
      integer, dimension(n) :: iter_ranges
      real(kind=8), dimension(n), intent(in) :: vals
      integer :: i

      allocate (index_matrix((order + 1)**n, n))
      allocate (indices(n, order + 1))

      do i = 1, n
         call get_indices(inds=indices(i, :), x=vals(i), which=i, n=n, o=order)
         !print *, 'ind =', indices(i,:)
      end do

      do i = 1, n
         iter_ranges(i) = order + 1
      end do

      inner_count = 1
      call for_recursive(n=n, c=1, ranges=iter_ranges, sub=index_row)
      inner_count = 1

   END SUBROUTINE construct_index_matrix

   SUBROUTINE partial_row(n, iter_list, vals)

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
      matrix(outer_count, inner_count) = res

      inner_count = inner_count + 1

   END SUBROUTINE partial_row

   SUBROUTINE row(n, iter_list)

      integer, intent(in) :: n
      real(kind=8), dimension(n) :: vals_in
      integer, dimension(n), intent(in) :: iter_list
      integer, dimension(n) :: iter_ranges
      integer :: i, j, ind

      call convert_indices(inds=index_matrix(inner_count, :), &
                           dims=axes_lengths, n=n, s=ind)

      do i = 1, n
         iter_ranges(i) = order + 1
         vals_in(i) = table_data(ind, i)
      end do

      do i = 1, n_out
         vec_b(inner_count, i) = table_data(ind, which_out(i))
      end do

      outer_count = inner_count
      inner_count = 1

      call for_recursive(n=n, c=1, ranges=iter_ranges, sub=partial_row, &
                         vals=vals_in)

      inner_count = outer_count

      inner_count = inner_count + 1

   END SUBROUTINE row

   SUBROUTINE construct_matrix(n)

      integer, intent(in) :: n
      integer, dimension(n) :: iter_ranges
!integer, dimension(n_out) :: which
      integer :: i

      allocate (matrix((order + 1)**n, (order + 1)**n))
      allocate (vec_b((order + 1)**n, n_out))
      allocate (vec_a((order + 1)**n, n_out))
!allocate(which_out(n_out))

!do i=1, n_out
      ! which_out(i) = which(i)
!enddo

      do i = 1, n
         iter_ranges(i) = order + 1
      end do

      inner_count = 1
      call for_recursive(sub=row, ranges=iter_ranges, n=n, c=1)
      inner_count = 1

   END SUBROUTINE construct_matrix

   SUBROUTINE compute_coeffs(n)

      integer, intent(in) :: n

      call solve_linear_system(n=(order + 1)**n, l=n_out, matrix=matrix, &
                               b=vec_b, x=vec_a)

   END SUBROUTINE compute_coeffs

   SUBROUTINE compute_result(n, iter_list, vals)

      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: iter_list
      real(kind=8), dimension(n_out) :: prod
      real(kind=8) :: delta
      real(kind=8), dimension(n), intent(in) :: vals
      integer :: i, j

      do i = 1, n_out
         prod(i) = 1.0d0
      end do

      do i = 1, n_out
         do j = 1, n
            prod(i) = prod(i)*vals(j)**(iter_list(j) - 1)
         end do

         delta = vec_a(inner_count, i)*prod(i)
         results(i) = results(i) + delta
      end do

      inner_count = inner_count + 1

   END SUBROUTINE compute_result

   SUBROUTINE interpolate(output_for_python, n, n_params, ord, which, vals, &
                          ll)

!use table_tools

      implicit none

      integer, intent(in) :: n, n_params
      integer, intent(in), optional :: ll
      integer, intent(inout), optional :: ord
      integer, dimension(n_params), intent(in) :: which
      integer, dimension(n) :: iter_ranges
      real(kind=8), dimension(n), intent(inout) :: vals
      real(kind=8), dimension(n_params), intent(out) :: output_for_python
      integer :: i, j

      n_out = n_params
      n_in = n

      if (present(ord)) then
         order = ord
      end if

      allocate (which_out(n_out))
      allocate (results(n_out))

      do i = 1, n_out
         which_out(i) = which(i)
      end do

      do i = 1, n
         iter_ranges(i) = order + 1
      end do

      call construct_index_matrix(n, vals)

      call construct_matrix(n=n)
      call compute_coeffs(n=n)

      inner_count = 1

      call for_recursive(sub=compute_result, c=1, n=n, vals=vals, &
                         ranges=iter_ranges)

      inner_count = 1

      do i = 1, n_out
         output_for_python(i) = results(i)
      end do

!print *, 'results=', results(:)

   END SUBROUTINE interpolate

END MODULE table_tools

