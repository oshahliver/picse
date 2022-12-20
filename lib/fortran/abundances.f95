MODULE abund

   use linalg
   use constants

contains

   subroutine construct_abundance_matrix(SiMg, FeMg, n_mats, &
                                         YMgi, YSii, xiH2Oi, matrix, b, xiFei, xiAlSii, xiAlMgi, additional)

      implicit none

      integer, intent(in) :: n_mats
      integer, dimension(n_mats), intent(in) :: YMgi, YSii
      real(kind=8), dimension(n_mats), intent(in) :: xiH2Oi, xiAlSii, xiAlMgi, &
                                                     xiFei
      real(8), intent(in) :: additional(:)
      real(kind=8), intent(in) :: SiMg, FeMg
      real(kind=8), dimension(n_mats, n_mats), intent(out) :: matrix
      real(8), intent(out) :: b(:, :)
      integer :: i, j, k

      do j = 1, n_mats
         matrix(1, j) = (1.0d0 - xiH2Oi(j))*(SiMg*(1.0d0 - xiAlMgi(j))* &
                                             (1.0d0 - xiFei(j))*YMgi(j) - (1.0d0 - xiAlSii(j))*YSii(j))

         !matrix(2,j) = (1.0d0-xiH2Oi(j))*(FeMg*(1.0d0-xiAlMgi(j))*&
         ! (1.0d0-xiFei(j))*YMgi(j) - (1.0d0-xiAlMgi(j))*xiFei(j)*YMgi(j))

         matrix(2, j) = 1.0d0

      end do

      if (n_mats > 2) then
         do i = 3, n_mats
            do j = 1, n_mats
               matrix(i, j) = 0d0

            end do
            b(i, 1) = additional(i - 2)
            matrix(i, i) = 1d0
         end do
      end if

      b(1, 1) = 0d0
      b(2, 1) = 1d0

!b(1,1)=-xiStv*(1.0d0-xiH2Oi(n_mats))*(SiMg*&
!(1.0d0-xiAlMgi(n_mats))*(1.0d0-xiFei(n_mats))*YMgi(n_mats) - &
!(1.0d0-xiAlSii(n_mats))*YSii(n_mats))

!-xiStv

!print *, 'matrix'
!do i=1, n_mats - 1
      ! print*, matrix(i,:)
!enddo

   end subroutine construct_abundance_matrix

   subroutine compute_abundance_vector(SiMg, FeMg, n_mats, abundances, &
                                       contents, additional, xiFei)

      implicit none

      integer, intent(in) :: n_mats
      integer, dimension(n_mats), intent(in) :: contents
      integer, dimension(n_mats) :: YMgi, YSii
      real(kind=8), dimension(n_mats) :: xiH2Oi, xiAlSii, xiAlMgi
      real(8), dimension(n_mats), intent(in) :: xiFei
      real(kind=8), intent(in) :: SiMg, FeMg
      real(kind=8), dimension(n_mats, n_mats) :: matrix
      real(kind=8), dimension(n_mats), intent(out) :: abundances
      real(kind=8), dimension(n_mats, 1) :: abundances_dummy
      real(8), dimension(n_mats, 1) :: b_vec
      real(8), intent(in), optional :: additional(:)
      real(8), allocatable :: additional_dummy(:)
      integer :: i, j

      if (n_mats > 2) then
         allocate (additional_dummy(n_mats - 2))
         if (present(additional)) then
            additional_dummy(:) = additional

         else
            additional_dummy(:) = 0d0
         end if

      else
         allocate (additional_dummy(n_mats))
         additional_dummy(1:) = 0d0

      end if

      do i = 1, n_mats
         abundances_dummy(i, 1) = 0.0d0
         YMgi(i) = material_YMg(contents(i))
         YSii(i) = material_YSi(contents(i))
         xiH2Oi(i) = 0d0
         xiAlSii(i) = 0d0
         xiAlMgi(i) = 0d0
      end do

!print *, 'contents =', contents(:)
!print *, 'xiAlSii =', xiAlSii(:)
!print *, 'xiAlMgi =', xiAlMgi(:)
!print *, 'FeMg =', FeMg
!print *, 'SiMg =', SiMg
!print *, 'xiH2Oi =', xiH2Oi
!print *, 'xiFei =', xiFei
!print *, 'additional =', additional_dummy

!In the core the first material is pure iron. There no compositional
!gradients are currently allowed and the fractions need not be computed
      if (contents(1) == 2 .or. contents(1) == 1) then
!pass

      elseif (size(contents) == 1) then
         abundances(1) = 1d0

      else

         call construct_abundance_matrix(SiMg=SiMg, FeMg=FeMg, n_mats=n_mats, &
                                         YMgi=YMgi, YSii=YSii, xiH2Oi=xiH2Oi, xiFei=xiFei, xiAlSii=xiAlSii, &
                                         xiAlMgi=xiAlMgi, matrix=matrix, b=b_vec, &
                                         additional=additional_dummy)

!call mat_print('matrix', matrix)
!print *, 'b_vec =', b_vec

         call gauss_elimination(a=matrix, sol=abundances_dummy, b=b_vec)

!print *, 'abundances =', abundances_dummy(:,1)

         do i = 1, n_mats
            abundances(i) = abundances_dummy(i, 1)
         end do
      end if

   end subroutine compute_abundance_vector

END MODULE abund
