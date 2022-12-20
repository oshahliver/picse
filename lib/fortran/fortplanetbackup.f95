MODULE class_shell

   use class_mixture
   use functions
   use constants
   use class_dict

   implicit none

   type shell

      type(mixture) :: mixture
      integer, dimension(:), allocatable :: contents
      real(kind=8), dimension(:), allocatable :: fractions
      real(kind=8) :: temp, pres, dens, Si_number, Fe_number, mass, radius
      real(kind=8) :: indigenous_mass, eps_H2O, gravity, v_esc, X_H2O
      real(kind=8) :: dPdr, eps_T_zero, gammaG0, rho0, q, dPdrho
      real(kind=8), dimension(4) :: gradients
      character(len=30) :: status = 'bare'
      type(dict) :: initials
      integer :: lay, adiabatType, tempType
      logical :: saturation

   end type shell

contains

   SUBROUTINE init_shell(self, contents, fractions, n_mats, T, P, eps_H2O, &
                         Fe_number, lay, m, r, tempType, gammaG0, alloc, eps_T_zero, adiabatType, &
                         q)

      logical, optional, intent(in) :: alloc
      logical :: alloc_dummy
      type(shell), intent(inout) :: self
      integer, intent(in) :: n_mats
      integer, intent(in), optional :: lay, tempType, adiabatType
      integer, dimension(n_mats), intent(in), optional :: contents
      real(kind=8), dimension(n_mats), intent(in), optional :: fractions
      real(kind=8), intent(in), optional :: gammaG0, q
      real(kind=8), intent(in), optional :: T, P, eps_H2O, Fe_number, m, r
      real(kind=8), intent(in), optional :: eps_T_zero

      if (present(alloc)) then
         alloc_dummy = alloc

      else

         alloc_dummy = .true.
      end if

      print *, 'alloc dummy =', alloc_dummy

      if (alloc_dummy) then

         allocate (self%contents(n_mats))
         allocate (self%fractions(n_mats))

         call init_mixture(self=self%mixture, contents=contents, &
                           fractions=fractions, n_mats=n_mats, T=T, P=P, eps_H2O=eps_H2O, &
                           Fe_number=Fe_number)

!Compute ambient density
         call compute_mixture(self=self%mixture, T=300.0d0, P=1.0d4)

         self%rho0 = self%mixture%dens

         self%gammaG0 = gammaG0
         self%q = q
         self%tempType = tempType
         self%mass = m
         self%radius = r
         self%temp = T
         self%pres = P
         self%contents = contents
         self%fractions = fractions
         self%lay = lay
         self%eps_T_zero = eps_T_zero
         self%adiabatType = adiabatType
!self%Fe_number = Fe_number
!self%eps_H2O = eps_H2O

         call init_dict(self=self%initials, n=11, n1=4, n2=4)

         call update_shell(self=self)
         print *, 'initial grads =', self%gradients(:)
         self%initials%real_vals(1) = self%radius
         self%initials%real_vals(2) = self%temp
         self%initials%real_vals(3) = self%mass
         self%initials%real_vals(4) = self%pres
         self%initials%real_vals(5) = self%dens
         self%initials%int_vals(6) = self%lay
         self%initials%int_vals(7) = self%tempType
         self%initials%real_vals(8) = self%dPdrho
         self%initials%real_vals(9) = self%gammaG0
         self%initials%real_vals(10) = self%rho0
         self%initials%real_arr(11, :) = self%gradients(:)

!If shell is being reset, don-t allocate everything and don-t
!update all parameters. Just reset all parameters to their original
!values.
      else

         self%radius = self%initials%real_vals(1)
         self%temp = self%initials%real_vals(2)
         self%mass = self%initials%real_vals(3)
         self%pres = self%initials%real_vals(4)
         self%dens = self%initials%real_vals(5)
         self%lay = self%initials%int_vals(6)
         self%tempType = self%initials%int_vals(7)
         self%dPdrho = self%initials%real_vals(8)
         self%gammaG0 = self%initials%real_vals(9)
         self%gradients(:) = self%initials%real_arr(11, :)

      end if

   END SUBROUTINE init_shell

   SUBROUTINE get_shell_abundances(self)

      type(mixture), intent(inout) :: self
      integer :: i

!Do stuff here to update the material fractions

   END SUBROUTINE get_shell_abundances

   SUBROUTINE get_shell_contents(self)

      type(shell), intent(inout) :: self

   END SUBROUTINE get_shell_contents

   SUBROUTINE update_shell_gradients(self)

      type(shell), intent(inout) :: self
      real(kind=8), dimension(4) :: y
      integer :: i

      y(1) = self%pres
      y(2) = self%mass
      y(3) = self%temp
      y(4) = self%dens

      call gradients(grads=self%gradients, &
                     r=self%radius, &
                     y=y, &
                     fractions=self%fractions, &
                     nmat=size(self%fractions), &
                     ll=self%contents, &
                     gammaG0=self%gammaG0, &
                     tempType=self%tempType, &
                     q=self%q, &
                     d0=self%rho0, &
                     adiabatType=self%adiabatType, &
                     eps_T_zero=self%eps_T_zero)

   END SUBROUTINE update_shell_gradients

   SUBROUTINE update_shell(self, T, P)

      type(shell), intent(inout) :: self
      real(kind=8), intent(in), optional :: P, T

      if (present(T)) then
         self%temp = T
      end if

      if (present(P)) then
         self%pres = P
      end if

      call compute_mixture(self=self%mixture, T=self%temp, P=self%pres)

      self%dens = self%mixture%dens
      self%dPdrho = self%mixture%dPdrho

      call update_shell_gradients(self=self)

   END SUBROUTINE update_shell

   SUBROUTINE print_shell(self)

      type(shell), intent(inout) :: self

      print *, ''
      print *, 'These are the shell properties:'
      print *, 'r/m:', self%radius, self%mass
      print *, 'P/T/rho/alpha:', self%pres, self%temp, self%dens, self%mixture%alpha_th
      print *, 'gradients:', self%gradients

   END SUBROUTINE print_shell

   SUBROUTINE reset_shell(self)

      type(shell), intent(inout) :: self

      call init_shell(self=self, n_mats=size(self%contents), alloc=.false.)

      self%status = 'bare'

!call update_shell_gradients(self=self)

   END SUBROUTINE reset_shell

   SUBROUTINE construct_shell(self, overconstruct, dr)

      type(shell), intent(inout) :: self
      logical, intent(in), optional :: overconstruct
      logical :: overconstruct_dummy
      real(kind=8), dimension(4) :: params
      real(kind=8), intent(in) :: dr
      real(kind=8) :: r_dummy

      if (present(overconstruct)) then
         overconstruct_dummy = overconstruct
      else
         overconstruct_dummy = .false.
      end if

      if (.not. self%status == 'constructed' .or. overconstruct_dummy) then

         params(1) = self%pres
         params(2) = self%mass
         params(3) = self%temp
         params(4) = self%dens

         r_dummy = self%radius

         call integrateRK(r_start=r_dummy, &
                          h=dr, &
                          y_in=params, &
                          gammaG0=self%gammaG0, &
                          fractions=self%fractions, &
                          tempType=self%tempType, &
                          nmat=size(self%contents), &
                          ll=self%contents, &
                          q=self%q, &
                          d0=self%rho0, &
                          adiabatType=self%adiabatType, &
                          eps_T_zero=self%eps_T_zero, &
                          r=self%radius, &
                          grads=self%gradients, &
                          y_out=params)

         self%pres = params(1)
         self%mass = params(2)
         self%temp = params(3)
         self%dens = params(4)

      else
         print *, 'WARNING: This shell has already been constructed.'
         print *, 'Pass overconstruct=.true. to ignore this message.'
      end if

   END SUBROUTINE construct_shell

END MODULE class_shell

!#######################################################################

MODULE class_layer

   use class_shell

   implicit none

   type layer

      integer :: lay, tempType, rhoType, adiabatType
      real(kind=8) :: r_in, r_out, eps_r, radius, FeMg, SiMg, Fe_number, X_H2O
      real(kind=8) :: eps_H2O, eps_Al, N_tot, N_Mg, N_Si, N_Al, N_H2O, N_Fe
      real(kind=8) :: rho0, q, eps_T_zero, pres, temp, mass, dens
      real(kind=8) :: mass_should, indigenous_mass, dr, gammaG0
      real(kind=8), dimension(4) :: params, gradients
      integer, dimension(:), allocatable :: contents
      real(kind=8), dimension(:), allocatable :: fractions
      logical, dimension(:), allocatable :: saturation
      logical :: bisec = .false.
      logical :: shell_iteration = .true.
      integer :: shell_iteration_count = 0
      logical :: change_bisec = .true.
      type(shell), dimension(:), allocatable :: shells
      logical :: force_bisection = .false.

   end type layer

contains

   SUBROUTINE init_layer(self, contents, fractions, n_mats, r_in, m, &
                         T, P, tempType, rhoType, adiabatType, q, gammaG0, eps_r, rho0, &
                         eps_T_zero, n_shells, eps_H2O, Fe_number, lay)

      type(layer), intent(inout) :: self
      integer, intent(in) :: n_mats, n_shells, lay
      real(kind=8), dimension(n_mats), intent(in) :: fractions
      integer, dimension(n_mats), intent(in) :: contents
      real(kind=8), intent(in) :: r_in, m, T, P, q, gammaG0, eps_r, rho0, eps_T_zero
      real(kind=8), intent(in) :: Fe_number, eps_H2O
      integer, intent(in) :: tempType, rhoType, adiabatType

      allocate (self%contents(n_mats))
      allocate (self%fractions(n_mats))
      allocate (self%shells(n_shells))

      self%lay = lay
      self%contents = contents
      self%fractions = fractions
      self%r_in = r_in
      self%mass = m
      self%temp = T
      self%pres = P
      self%q = q
      self%rho0 = rho0
      self%gammaG0 = gammaG0
      self%eps_r = eps_r
      self%eps_T_zero = eps_T_zero
      self%rhoType = rhoType
      self%tempType = tempType
      self%adiabatType = adiabatType
      self%Fe_number = Fe_number
      self%eps_H2O = eps_H2O

!call init_shell(self=self%shells(1), T=T, P=P, contents=contents, &
!fractions = fractions, tempType=tempType, adiabatType=adiabatType, &
!q=q, gammaG0=gammaG0, eps_T_zero=eps_T_zero, alloc=.true., &
!eps_H2O=eps_H2O, Fe_number=Fe_number, n_mats=n_mats, lay=lay)

!call print_shell(self=self%shells(1))

   END SUBROUTINE init_layer

END MODULE class_layer

!#######################################################################

MODULE class_planet

   use class_layer

   implicit none

   type planet

      real(kind=8) :: T_surface_should, P_surface_should, M_surface_should
      real(kind=8) :: R_surface_should, Mg_number_should, Si_number_should

   end type planet

contains

   SUBROUTINE init_planet(self)

      type(planet), intent(inout) :: self

   END SUBROUTINE init_planet

   SUBROUTINE construct_planet(self)

      type(planet), intent(inout) :: self

   END SUBROUTINE construct_planet

END MODULE class_planet
