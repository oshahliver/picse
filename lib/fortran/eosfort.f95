MODULE eosfort

   use class_table
   use run_params
   use phase

   implicit none

!Why do I need to do it like this? Why can I not just
!Define type(table_list) :: tables ??
   private
   public :: initiate, compute
   type(table_list), public :: tables

contains

   character(len=20) function str(k)
!Convert integer to string

      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)

   end function str

   SUBROUTINE initiate(n_tables, n_tables_hyd, auto, table_dir)

      implicit none

      integer, intent(inout) :: n_tables, n_tables_hyd
!integer :: n_tables_var
      logical, intent(in), optional :: auto
      character(*), intent(in) :: table_dir
      character(34), dimension(:), allocatable :: file_names
      character(8) :: fmt
      logical :: auto_dummy
      integer :: ord = 1
      integer :: i, j

      if (.not. present(auto)) then
         auto_dummy = .true.
      else
         auto_dummy = auto
      end if

!Create table file names to read in from files
      if (auto_dummy) then
         n_tables = n_tables
         n_tables_hyd = 5
      end if

!n_tables_var = size(eos_tables_variable_Fe)
      
      allocate (file_names(n_tables + n_tables_hyd))
      allocate (tables%objects(n_tables + n_tables_hyd))
      allocate (tables%n_in_list(n_tables + n_tables_hyd))

!Load tables with variable iron contents
      do i = 1, n_tables
         fmt = str(i)
         if (eos_tables_variable_Fe(i) == 0) then
            file_names(i) = 'eos_'//trim(fmt)//'_Fe_0_0_H2O_0_0_Al_0_0.tab'

         elseif (eos_tables_variable_Fe(i) == 1) then
            file_names(i) = 'eos_'//trim(fmt)//'_Fe_var_H2O_0_0_Al_0_0.tab'

         end if
      end do

!Load tables for hydrated materials
      j = 0
      do i = 1, n_tables
         fmt = str(i)
         if (eos_tables_hydrated(i) == 1) then
            j = j + 1
            if (eos_tables_variable_Fe(i) == 0) then
               file_names(j + n_tables) = 'eos_'//trim(fmt)//'_Fe_0_0_H2O_1_0_Al_0_0.tab'

            elseif (eos_tables_variable_Fe(i) == 1) then
               file_names(j + n_tables) = 'eos_'//trim(fmt)//'_Fe_var_H2O_1_0_Al_0_0.tab'
            end if
         end if
      end do

      do i = 1, n_tables + n_tables_hyd
         call init(self=tables%objects(i), label=i)
      end do

!print *, 'We have the following objects:'
!do i=1, n_tables
!  write(*,*) tables%objects(i)%label
!enddo

      do i = 1, n_tables + n_tables_hyd
         print *, 'loading EoS table file: ', table_dir//file_names(i)
         call load(tables%objects(i), file_name=table_dir//file_names(i))
!  print *, 'The axes_lengths are:'
!  print *, tables%objects(i)%axes_lengths
      end do

   END SUBROUTINE initiate

   SUBROUTINE compute(which, n_out, vals, order, alloc, ll, res, eps_H2O)

      implicit none

      logical, intent(in) :: alloc
      logical :: alloc_dummy
      integer, intent(inout), optional :: order
      integer, intent(in) :: n_out
      integer, intent(in) :: ll
      integer :: ll_dummy
      integer, dimension(n_out), intent(in) :: which
      real(8), intent(inout) :: vals(:)
      real(8), intent(in) :: eps_H2O
      real(8), dimension(n_out), intent(inout) :: res
      real(8) :: pres, temp
      integer :: i, j, n_in, phase

      temp = vals(1)
      pres = vals(2)

      if (.not. present(order)) then
         order = 1
      end if

      if (alloc) then
      do i = 1, size(tables%objects)
         if (allocated(tables%objects(i)%which_out)) then
            deallocate (tables%objects(i)%which_out)
            deallocate (tables%objects(i)%results)
            deallocate (tables%objects(i)%indices)
            deallocate (tables%objects(i)%index_matrix)
            deallocate (tables%objects(i)%matrix)
            deallocate (tables%objects(i)%vec_b)
            deallocate (tables%objects(i)%vec_a)
         end if

         n_in = tables%objects(i)%n_in

         allocate (tables%objects(i)%which_out(n_out))
         allocate (tables%objects(i)%results(n_out))
         allocate (tables%objects(i)%index_matrix((order + 1)**n_in, n_in))
         allocate (tables%objects(i)%indices(n_in, order + 1))
         allocate (tables%objects(i)%matrix((order + 1)**n_in, (order + 1)**n_in))
         allocate (tables%objects(i)%vec_b((order + 1)**n_in, n_out))
         allocate (tables%objects(i)%vec_a((order + 1)**n_in, n_out))
      end do
      end if
!If the material is in hydrated state (eps_H2O = 1), the hydrated
!EoS table has to be taken.
      if (eps_H2O .gt. 9.9d-1) then
         ll_dummy = which_eos_tables_hydrated(ll)
      else
         ll_dummy = ll
      end if
!ll_dummy = 10

!Here the phase must be updated according to the given input parameters
!and the material at hand.
      call get_phase(T=temp, P=pres, ll=ll_dummy, ph=phase, xiFe=vals(3))
!~ print *, 'phase in compute =', phase
      call interpolate(self=tables%objects(ll_dummy), n_params=n_out, ord=order, &
                       which=which, vals=vals, phase=phase)
      res(:) = tables%objects(ll_dummy)%results(:)
!~ print *, 'll =', ll_dummy
!~ print *, 'phase =', phase
!~ print *, 'vals =', vals
!~ print *, 'which =', which
!~ print *, 'res =', res
!res(5) = 1.0d-2
   END SUBROUTINE compute

END MODULE eosfort

