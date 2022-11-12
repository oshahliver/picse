MODULE linalg

   implicit none

contains

   SUBROUTINE matrix_multiplication(n, m, l, a1, a2, a)
!a1xa2 = a, nxm x mxl = nxl matrix multiplication

      integer, intent(in) :: n, m, l
      real(kind=8), dimension(n, m), intent(in) :: a1
      real(kind=8), dimension(m, l), intent(in) :: a2
      real(kind=8), dimension(m, l), intent(out) :: a
      integer :: i, j, k

      do i = 1, n
         do j = 1, l
            do k = 1, m
               a(i, j) = a(i, j) + a1(i, k)*a2(k, j)
            end do
         end do
      end do

   END SUBROUTINE matrix_multiplication

   SUBROUTINE gauss_elimination(a, sol, b)

      implicit none

      real(8), intent(in) :: a(:, :), b(:, :)
      real(8), intent(out) :: sol(:, :)
      real(8), allocatable :: aa(:, :)
      integer :: i, j, k, n, nn

!Size of the system n= (order+1)**n_input_params
      n = size(a, 1)

!Number of parameters for which system should be solvec
      nn = size(b, 2)

      allocate (aa(n, n + nn))

      do i = 1, n
         aa(i, :n) = a(i, :)
         aa(i, n + 1:) = b(i, :)
      end do

      do i = 1, n
         sol(i, :) = 0.0d0
      end do

!Forward Elimination
      do k = 1, n - 1
         do i = k + 1, n
            aa(i, k) = aa(i, k)/aa(k, k)
            do j = k + 1, n + nn
               aa(i, j) = aa(i, j) - aa(i, k)*aa(k, j)
            end do
         end do
      end do

!Backward Substitution
      do i = 1, n
         sol(n - i + 1, :) = aa(n - i + 1, n + 1:)
         do j = 1, i - 1
            sol(n - i + 1, :) = sol(n - i + 1, :) - aa(n - i + 1, n - j + 1)*sol(n - j + 1, :)
         end do
         sol(n - i + 1, :) = sol(n - i + 1, :)/aa(n - i + 1, n - i + 1)
      end do

   END SUBROUTINE gauss_elimination

   subroutine solve_LU(a, b, x)
      real(8), intent(in)   :: a(:, :)
      real(8), intent(inout) :: b(:), x(:)
      integer               :: i, j, n
      real(8), allocatable  :: aa(:, :), l(:, :), u(:, :)
      real(8), allocatable  :: p(:, :), bb(:)
      integer, allocatable  :: ipiv(:)

      n = size(a, 1)
      allocate (aa(n, n), l(n, n), u(n, n), p(n, n), ipiv(n), bb(n))

      forall (j=1:n, i=1:n)
         aa(i, j) = a(i, j)
         u(i, j) = 0d0
         p(i, j) = merge(1, 0, i .eq. j)
         l(i, j) = merge(1d0, 0d0, i .eq. j)
      end forall

      call lu(aa, ipiv)

      do i = 1, n
         l(i, :i - 1) = aa(ipiv(i), :i - 1)
         u(i, i:) = aa(ipiv(i), i:)
         bb(i) = b(ipiv(i))
      end do

      p(ipiv, :) = p

      call backward_substitution(l, b, bb, 'lower')
      call backward_substitution(u, x, b, 'upper')
      !call matrix_multiplication(n=n, m=n, l=1, a1=p, a2=b, a=x)
      !call matrix_multiplication(n=n, m=n, l=n, a1=l, a2=u, a=test)
      !call matrix_multiplication(n=n, m=n, l=n, a1=p, a2=test, a=testtest)

      !call mat_print('test',test)
      !call mat_print('testtest',testtest)
!    call mat_print('a',a)
!    call mat_print('p',p)
!    call mat_print('l',l)
!    call mat_print('u',u)

!    print *, 'bb'
!    print *, bb(:)
!    print *, 'x'
!    print *, x(:)

      !print *, "residual"
      !print *, "|| P.A - L.U || =  ", maxval(abs(matmul(p,a)-matmul(l,u)))
   end subroutine solve_LU

   subroutine lu(a, p)
!   in situ decomposition, corresponds to LAPACK's dgebtrf
      real(8), intent(inout) :: a(:, :)
      integer, intent(out) :: p(:)
      integer                :: n, i, j, k, kmax
      n = size(a, 1)
      p = [(i, i=1, n)]
      do k = 1, n - 1
         kmax = maxloc(abs(a(p(k:), k)), 1) + k - 1
         if (kmax /= k) p([k, kmax]) = p([kmax, k])
         a(p(k + 1:), k) = a(p(k + 1:), k)/a(p(k), k)
         forall (j=k + 1:n) a(p(k + 1:), j) = a(p(k + 1:), j) - a(p(k + 1:), k)*a(p(k), j)
      end do
   end subroutine

   subroutine backward_substitution(matrix, sol, b, type)
      real(8), intent(in) :: matrix(:, :), b(:)
      real(8), intent(inout) :: sol(:)
      character(len=5), intent(in) :: type
      integer :: n, i, j
      n = size(matrix, 1)

      if (type == 'lower') then
      do i = 1, n

         sol(i) = b(i)
         do j = 1, i - 1
            sol(i) = sol(i) - matrix(i, j)*sol(j)
         end do
         sol(i) = sol(i)/matrix(i, i)
      end do

      elseif (type == 'upper') then
      do i = 1, n
         sol(n - i + 1) = b(n - i + 1)
         do j = 1, i - 1
            sol(n - i + 1) = sol(n - i + 1) - matrix(n - i + 1, n - j + 1)*sol(n - j + 1)
         end do
         sol(n - i + 1) = sol(n - i + 1)/matrix(n - i + 1, n - i + 1)
      end do
      end if
   end subroutine backward_substitution

   subroutine mat_print(amsg, a)
      character(*), intent(in) :: amsg
      class(*), intent(in) :: a(:, :)
      integer                  :: i
      print *, ' '
      print *, amsg
      do i = 1, size(a, 1)
         select type (a)
         type is (real(8)); print'(100f8.3)', a(i, :)
         type is (integer); print'(100i8  )', a(i, :)
         end select
      end do
      print *, ' '
   end subroutine

!#######################################################################
   SUBROUTINE compute_inverse_matrix(n, a, a1)

      integer, intent(in) :: n
      real(kind=8), intent(inout) :: a(:, :)
      real(kind=8), dimension(n, n), intent(inout) :: a1
      integer :: i, j, k, l
      real(kind=8) :: z

      do i = 1, n
         do j = 1, n
            if (i == j) then
               a1(i, j) = 1
            else
               a1(i, j) = 0
            end if

         end do
      end do

      do i = 1, n
         z = a(i, i)
         do j = 1, n
            a(i, j) = a(i, j)/z
            a1(i, j) = a1(i, j)/z
         end do

         do j = i + 1, n
            z = a(j, i)
            do k = 1, n
               a(j, k) = a(j, k) - z*a(i, k)
               a1(j, k) = a1(j, k) - z*a1(i, k)
            end do
         end do
      end do

      do i = 1, n - 1
         do j = i + 1, n
            z = a(i, j)
            do l = 1, n
               a(i, l) = a(i, l) - z*a(j, l)
               a1(i, l) = a1(i, l) - z*a1(j, l)
            end do
         end do
      end do

   END SUBROUTINE compute_inverse_matrix

!#######################################################################
   SUBROUTINE solve_linear_system(n, l, matrix, b, x)

      integer, intent(in) :: n, l
      real(8), intent(in) :: matrix(:, :)
      real(kind=8), dimension(n, n) :: matrix_dummy
      real(kind=8), dimension(n, n) :: inverse
      real(kind=8), dimension(n, l), intent(in) :: b
      real(kind=8), dimension(n, l), intent(out) :: x
      integer :: i, j

      matrix_dummy = matrix

      call compute_inverse_matrix(n=n, a=matrix_dummy, a1=inverse)
      call matrix_multiplication(n=n, m=n, l=l, a1=inverse, a2=b, a=x)

   END SUBROUTINE solve_linear_system

END MODULE linalg

