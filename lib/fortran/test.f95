module mods

   implicit none

contains

   subroutine init(a, a_gauss, x, b, n)

      implicit none

      integer, intent(in) :: n
      real(8), intent(inout) :: a(:, :), a_gauss(:, :)
      real(8), intent(inout) :: x(:, :), b(:, :)
      integer :: i, j

      do i = 1, n
         do j = 1, n
            a(i, j) = i**j
         end do

         b(i, 1) = j
         b(i, 2) = j*i

         x(i, :) = 0d0

!call random_number(a)

      end do

   end subroutine init

end module mods

program test

   use class_weird_array
   use linalg
   use mods

   implicit none

   integer :: i, j, N_run
   integer, parameter :: dim = 8
   real(kind=8) :: z
   real(kind=8), dimension(3, 3) :: array1
   real(8), allocatable :: matrix(:, :)
   real(8), allocatable :: b(:, :), x(:, :)
   real(8), allocatable :: gauss_matrix(:, :)
   real(kind=8), dimension(3) :: array2
   character(len=35) :: file_name, str
   type(weird_array) :: test_array
   integer, dimension(3) :: axlgth
   real(8) :: t_start, t_end

   allocate (matrix(dim, dim), b(dim, 2), x(dim, 2), gauss_matrix(dim, dim + 1))
   call init(matrix, gauss_matrix, x, b, dim)

!matrix = reshape((/1,1,2,1,-2,3,-1,3,1/), [dim,dim])
!b = (/4, -6, 7/)

   N_run = 100000
   axlgth = (/2, 3, 4/)

   array2 = (/1./7., 2., 3./)

   array1 = transpose(reshape((/1, 1, 1, 2, 2, 2, 7, 8, 9/), shape(array1)))

!print*, "array1 = ", array1
!print *, "array2 =", array2

   do i = 1, 3
      file_name = "Output"//trim(str(i))//".txt"
      print *, file_name
   end do

!call cpu_time(t_start)

!do i=1, N_run
!  call solve_LU(matrix, b, x)
!enddo

!call cpu_time(t_end)

!print *, ''
!print *, 'Time LU =', (t_end-t_start)*1.0d3, 'ms'
!print *, 'x =', x(:)

   call init(matrix, gauss_matrix, x, b, dim)

!matrix = reshape((/1,1,2,1,-2,3,-1,3,1/), [dim,dim])

   call cpu_time(t_start)

   do i = 1, N_run
      call solve_linear_system(dim, 2, matrix, b, x)
   end do

   call cpu_time(t_end)

   print *, ''
   print *, 'Time Inverse =', (t_end - t_start)*1.0d3, 'ms'

!do i=1, size(b, 1)
!  print *, 'b =', b(i,:)
!enddo

!do i=1, size(x, 1)
!  print *, 'x =', x(i,:)
!enddo

   call init(matrix, gauss_matrix, x, b, dim)

!matrix = reshape((/1,1,2,1,-2,3,-1,3,1/), [dim,dim])

   call init(matrix, gauss_matrix, x, b, dim)

   call cpu_time(t_start)

   do i = 1, N_run
      call gauss_elimination(matrix, x, b)
   end do

   call cpu_time(t_end)
   print *, ''
   print *, 'Time Gauss =', (t_end - t_start)*1.0d3, 'ms'

end program test

character(len=35) function str(k)
   integer, intent(in) :: k
   write (str, *) k
   str = adjustl(str)
end function str
