subroutine alloc(a, x, b, n)

   implicit none

   integer, intent(in) :: n
   real(8), intent(inout) :: a(:, :)
   real(8), intent(inout) :: x(:), b(:)
   integer :: i, j

   do i = 1, n
      do j = 1, n
         a(i, j) = i + sqrt(real(i))
      end do

      b(i) = i
      x(i) = 0d0
   end do

end subroutine alloc

program test

   use class_weird_array
   use linalg

   implicit none

   integer :: i, j, N_run
   real(kind=8) :: z
   real(kind=8), dimension(3, 3) :: array1
   real(8), dimension(3, 3) :: matrix = reshape([1, 1, 2, 1, -2, 3, -1, 3, 1], [3, 3])
   real(8), dimension(3) :: b = (/4, -6, 7/), x = (/0.0d0, 0.0d0, 0.0d0/)
   real(8), dimension(3, 4) :: gauss_matrix
   real(kind=8), dimension(3) :: array2
   character(len=35) :: file_name, str
   type(weird_array) :: test_array
   integer, dimension(3) :: axlgth
   real(8) :: t_start, t_end

   N_run = 1
   axlgth = (/2, 3, 4/)

   array2 = (/1./7., 2., 3./)

   array1 = transpose(reshape((/1, 1, 1, 2, 2, 2, 7, 8, 9/), shape(array1)))

   print *, "array1 = ", array1
   print *, "array2 =", array2

   print *, array2(1)
   do i = 1, 3
      file_name = "Output"//trim(str(i))//".txt"
      print *, file_name
   end do

   call cpu_time(t_start)

   x = (/0.0d0, 0.0d0, 0.0d0/)
   b = (/4, -6, 7/)
   matrix = reshape([1, 1, 2, 1, -2, 3, -1, 3, 1], [3, 3])
   do i = 1, N_run
      call solve_LU(matrix, b, x)
   end do

   call cpu_time(t_end)

   print *, ''
   print *, 'time =', (t_end - t_start)*1.0d3, 'ms'
   print *, 'x =', x(:)

   call cpu_time(t_start)

   x = (/0.0d0, 0.0d0, 0.0d0/)
   b = (/4, -6, 7/)
   matrix = reshape([1, 1, 2, 1, -2, 3, -1, 3, 1], [3, 3])
   do i = 1, N_run
      call solve_linear_system(3, 1, matrix, b, x)
   end do

   call cpu_time(t_end)

   print *, ''
   print *, 'time =', (t_end - t_start)*1.0d3, 'ms'

!call mat_print('matrix', matrix)

   print *, 'x =', x(:)

   call cpu_time(t_start)

   x = (/0.0d0, 0.0d0, 0.0d0/)
   b = (/4, -6, 7/)
   matrix = reshape([1, 1, 2, 1, -2, 3, -1, 3, 1], [3, 3])

   do i = 1, 3
      do j = 1, 3
         gauss_matrix(i, j) = matrix(i, j)
      end do
   end do

   do i = 1, 3
      gauss_matrix(i, 4) = b(i)
   end do

   do i = 1, N_run
      call gauss_elimination(gauss_matrix, x)
   end do

   call cpu_time(t_end)
   print *, ''
   print *, 'time =', (t_end - t_start)*1.0d3, 'ms'
   print *, 'x =', x(:)

end program test

character(len=35) function str(k)
   integer, intent(in) :: k
   write (str, *) k
   str = adjustl(str)
end function str
