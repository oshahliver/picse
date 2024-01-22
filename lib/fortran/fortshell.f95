MODULE class_shell

   use class_mixture
   use functions
   use constants
   use class_dict
   use functionspy

   implicit none

   type shell

      type(mixture) :: mixture
      integer, dimension(:), allocatable :: contents, YMg, YSi, YO, YS, YFe
      real(4) :: timer, t_start, t_end
      real(8), dimension(:), allocatable :: fractions, xi_Fe, xi_H2O, &
                                            X_H2O, xi_Si, weight_fractions, composition_gradients, &
                                            mean_fractions, mean_weight_fractions
      real(8) :: Si_number, Fe_number, radius
      real(8) :: indigenous_mass, gravity, v_esc, volume, &
                 N_Mg, N_H2O, N_Fe, N_Si, N_Al, N_O, N_S, N_H, N_tot, FeMg, SiMg, eps_Al, &
                 eps_H2O
      real(8) :: dPdr, eps_T_zero, gammaG0, rho0, q, dPdrho, xi_H, xi_Stv
      real(8), dimension(n_params_integration) :: integration_parameters
      real(8) :: omega, dr
      real(8), dimension(n_params_integration) :: gradients
      character(len=30) :: status = 'bare'
      type(dict) :: initials
      integer :: lay, adiabatType, tempType, n_mats
      logical :: saturation
      logical :: force_bisection = .false., constant_Fe = .true.
   end type shell

contains

   SUBROUTINE init_shell(self, contents, fractions, n_mats, eps_H2O, &
                         Fe_number, lay, r, tempType, gammaG0, alloc, eps_T_zero, adiabatType, &
                         q, rho0, eps_Al, Si_number, omega, xi_H, xi_Stv, &
                         composition_gradients, integration_parameters)
                         

      logical, optional, intent(in) :: alloc
      logical :: alloc_dummy
      type(shell), intent(inout) :: self
      integer, intent(in) :: n_mats
      integer, intent(in), optional :: lay, tempType, adiabatType
      integer, dimension(n_mats), intent(in), optional :: contents
      real(8), dimension(n_mats), intent(in), optional :: fractions, &
                                                          composition_gradients
      real(8), intent(in), optional :: gammaG0, q, rho0, omega, xi_H, xi_Stv
      real(8), intent(in), optional :: eps_H2O, Fe_number, r, &
                                       Si_number
      real(8), intent(in), optional :: eps_T_zero, eps_Al
      real(8), intent(in), dimension(n_params_integration), optional :: integration_parameters
      integer :: i

      if (present(alloc)) then
         alloc_dummy = alloc
      else
         alloc_dummy = .true.
      end if

      if (present(xi_H)) then
         self%xi_H = xi_H
      else
         self%xi_H = 0.0d0
      end if

      if (present(xi_Stv)) then
         self%xi_Stv = xi_Stv
      else
         self%xi_Stv = 0.0d0
      end if

      if (alloc_dummy) then
         if (.not. allocated(self%contents)) then
            allocate (self%contents(n_mats))
         end if

         if (.not. allocated(self%fractions)) then
            allocate (self%fractions(n_mats))
         end if

         if (.not. allocated(self%weight_fractions)) then
            allocate (self%weight_fractions(n_mats))
         end if

         if (.not. allocated(self%mean_weight_fractions)) then
            allocate (self%mean_weight_fractions(n_mats))
         end if

         if (.not. allocated(self%mean_fractions)) then
            allocate (self%mean_fractions(n_mats))
         end if

         if (.not. allocated(self%composition_gradients)) then
            allocate (self%composition_gradients(n_mats))
         end if

         if (.not. allocated(self%YMg)) then
            allocate (self%YMg(n_mats))
            allocate (self%YSi(n_mats))
            allocate (self%YO(n_mats))
            allocate (self%YS(n_mats))
            allocate (self%YFe(n_mats))
            allocate (self%xi_Fe(n_mats))
            allocate (self%xi_H2O(n_mats))
            allocate (self%X_H2O(n_mats))
         end if

         if (present(omega)) then
            self%omega = omega
         else
            self%omega = 0d0
         end if

         if (present(composition_gradients)) then
            self%composition_gradients = composition_gradients
         else
            self%composition_gradients(:) = 0d0
         end if
         
         self%integration_parameters = integration_parameters
         self%n_mats = n_mats
         self%gammaG0 = gammaG0
         self%q = q
         self%tempType = tempType
         self%radius = r
         self%contents = contents
         self%fractions = fractions
         self%lay = lay
         self%eps_T_zero = eps_T_zero
         self%adiabatType = adiabatType
         self%indigenous_mass = 0.0d0
         self%eps_Al = eps_Al
         self%eps_H2O = eps_H2O
         self%Fe_number = Fe_number
         self%Si_number = Si_number
         self%volume = 0d0
         self%timer = 0.0
         self%xi_Stv = xi_Stv
         self%dr = 0d0

         do i = 1, n_mats
            self%YO(i) = material_YO(self%contents(i))
            self%YMg(i) = material_YMg(self%contents(i))
            self%YSi(i) = material_YSi(self%contents(i))
            self%YS(i) = material_YS(self%contents(i))
            self%YFe(i) = material_YFe(self%contents(i))
         end do

         if (.not. Si_number == 1.0d0) then
            self%SiMg = self%Si_number/(1.0d0 - self%Si_number)
         else
            self%SiMg = 1.0d10
         end if

         if (.not. Fe_number == 1.0d0) then
            self%FeMg = self%Fe_number/(1.0d0 - self%Fe_number)
         else
            self%FeMg = 1.0d10
         end if

         call init_mixture(self=self%mixture, contents=contents, &
                           fractions=fractions, n_mats=n_mats, T=self%integration_parameters(3), &
                           P=self%integration_parameters(1), eps_H2O=eps_H2O, &
                           eps_Al=eps_Al, Si_number=self%Si_number, Fe_number=self%Fe_number, &
                           xi_H=self%xi_H, xi_Stv=self%xi_Stv, lay=self%lay)

         self%weight_fractions = self%mixture%weight_fractions
         self%fractions = self%mixture%fractions
         
         self%xi_H2O = self%mixture%xi_H2O
         self%X_H2O = self%mixture%X_H2O
         self%xi_Fe = self%mixture%xi_Fe

!Compute ambient density
!Note that the fractions at T=300 K and T=1.0d4 Pa would be in general
!different but we want to know the ambient conditions of the given shell
!at the given fractions, so the fractions must not be updated to compute
!the ambient parameters here
         call compute_mixture(self=self%mixture, T=300.0d0, P=1.0d4)
         self%rho0 = self%mixture%dens

         call init_dict(self=self%initials, n=21, n1=n_params_integration, &
                        n2=n_params_integration)

         self%initials%real_arr(19, 1:self%n_mats) = self%fractions(1:self%n_mats)
         self%initials%real_arr(20, 1:self%n_mats) = self%weight_fractions(1:self%n_mats)

!Here also the mixture must be computed and the shell gradients must
!be updated
         call update_shell(self=self)

         self%initials%real_vals(1) = self%radius
         self%initials%real_vals(2) = self%integration_parameters(3)
         self%initials%real_vals(3) = self%integration_parameters(2)
         self%initials%real_vals(4) = self%integration_parameters(1)
         self%initials%real_vals(5) = self%integration_parameters(4)
         self%initials%int_vals(6) = self%lay
         self%initials%int_vals(7) = self%tempType
         self%initials%real_vals(8) = self%dPdrho
         self%initials%real_vals(9) = self%gammaG0
         self%initials%real_vals(10) = self%rho0
         self%initials%real_arr(11, :) = self%gradients(:)
         self%initials%real_vals(12) = self%Fe_number
         self%initials%real_vals(13) = self%Si_number
         self%initials%real_vals(14) = self%eps_H2O
         self%initials%real_vals(15) = self%eps_Al
         self%initials%real_vals(16) = self%integration_parameters(5)
         self%initials%real_vals(17) = self%omega
         self%initials%real_vals(18) = self%xi_H
         self%initials%real_vals(19) = self%integration_parameters(7)
         self%initials%real_vals(20) = self%integration_parameters(8)

!If shell is being reset, don-t allocate everything and don-t
!update all parameters. Just reset all parameters to their original
!values.
      else
!~ print *, 'not allocating shell'
         self%radius = self%initials%real_vals(1)
         self%lay = self%initials%int_vals(6)
         self%tempType = self%initials%int_vals(7)
         self%dPdrho = self%initials%real_vals(8)
         self%gammaG0 = self%initials%real_vals(9)
         self%gradients(:) = self%initials%real_arr(11, :)
         self%FeMg = self%initials%real_vals(12)
         self%SiMg = self%initials%real_vals(13)
         self%eps_H2O = self%initials%real_vals(14)
         self%eps_Al = self%initials%real_vals(15)
         self%omega = self%initials%real_vals(17)
         self%xi_H = self%initials%real_vals(18)
         self%fractions(1:self%n_mats) = self%initials%real_arr(19, 1:self%n_mats)
         self%weight_fractions(1:self%n_mats) = self%initials%real_arr(20, 1:self%n_mats)
         self%indigenous_mass = 0.0d0

         self%N_tot = 0.0d0
         self%N_Mg = 0.0d0
         self%N_Si = 0.0d0
         self%N_Fe = 0.0d0
         self%N_O = 0.0d0
         self%N_H2O = 0.0d0
         self%N_Al = 0.0d0
         self%N_S = 0.0d0
         self%N_H = 0d0 !This counts the hydrogen of hydrated substances (e.g. FeH)

         self%integration_parameters(1) = self%initials%real_vals(4) ! Pressure
         self%integration_parameters(2) = self%initials%real_vals(3) ! Enclosed mass
         self%integration_parameters(3) = self%initials%real_vals(2) ! Temperature
         self%integration_parameters(4) = self%initials%real_vals(5) ! Density
         self%integration_parameters(5) = self%initials%real_vals(16) ! Weight fraction
         self%integration_parameters(6) = self%initials%real_vals(18) ! Moment of inertia factor
         self%integration_parameters(7) = self%initials%real_vals(19) ! Gravitational energy
         self%integration_parameters(8) = self%initials%real_vals(20) ! Internal energy
      end if

   END SUBROUTINE init_shell

!#######################################################################
   SUBROUTINE get_shell_contents(self)

      type(shell), intent(inout) :: self
      real(8), dimension(size(self%contents)) :: abundances, &
                                                 xiAlMgi, xiAlSii
      integer :: i, j

      do i = 1, size(self%contents)
         xiAlSii(i) = 0.0d0
         xiAlMgi(i) = 0.0d0
      end do

      self%N_tot = 0.0d0
      self%N_Mg = 0.0d0
      self%N_Si = 0.0d0
      self%N_Fe = 0.0d0
      self%N_O = 0.0d0
      self%N_H2O = 0.0d0
      self%N_Al = 0.0d0
      self%N_S = 0.0d0
      self%N_H = 0d0

!Compute number of atomic species per mole
      do i = 1, size(self%contents)

         self%N_Mg = self%N_Mg + self%mean_fractions(i)*(1.0d0 - self%xi_H2O(i))* &
                     (1.0d0 - self%xi_Fe(i))*(1.0d0 - xiAlMgi(i))*self%YMg(i)

         !Count iron from Mg-Fe substitution
         self%N_Fe = self%N_Fe + self%mean_fractions(i)*(1.0d0 - self%xi_H2O(i))* &
                     self%xi_Fe(i)*(1d0 - self%xi_H)*self%YMg(i)

         !Count iron from pure species
         self%N_Fe = self%N_Fe + self%mean_fractions(i)*self%YFe(i)
         self%N_Si = self%N_Si + self%mean_fractions(i)*(1.0d0 - self%xi_H2O(i))* &
                     (1.0d0 - xiAlSii(i))*self%YSi(i)
         self%N_Al = self%N_Al + self%mean_fractions(i)*(1.0d0 - self%xi_H2O(i))* &
                     (xiAlSii(i)*self%YSi(i) + (1.0d0 - self%xi_Fe(i))*xiAlMgi(i)*self%YMg(i))
         self%N_O = self%N_O + self%mean_fractions(i)*(1.0d0 - self%xi_H2O(i))*self%YO(i)

         self%N_S = self%N_S + self%mean_fractions(i)*self%YS(i)

         self%N_H2O = self%N_H2O + self%mean_fractions(i)*self%xi_H2O(i)

         !Note: xi_H is only a scalar at this point. This can be generalized later
         !By construction only the first material in the core (i.e. Fe) can
         !contain H at this point
         if (i == 1) then
            self%N_H = self%N_H + self%mean_fractions(i)*self%xi_H
         end if
      end do

!Compute total number of moles in the shell
      self%N_tot = self%indigenous_mass/(self%N_Mg*mMg + self%N_Fe*mFe + &
                                         self%N_Si*mSi + self%N_Al*mAl + self%N_H2O*mH2O + self%N_O*mO + &
                                         self%N_S*mS + self%N_H*mH)
      self%N_Mg = self%N_Mg*self%N_tot
      self%N_Si = self%N_Si*self%N_tot
      self%N_Fe = self%N_Fe*self%N_tot
      self%N_Al = self%N_Al*self%N_tot
      self%N_O = self%N_O*self%N_tot
      self%N_S = self%N_S*self%N_tot
      self%N_H = self%N_H*self%N_tot
      self%N_H2O = self%N_H2O*self%N_tot
  
!Note that N_tot referes to the number of moles of molecular substances
!while the atomic abundances (N_Mg etc.) refer to the number of moles
!of each species given the number of moles of molecules which are made
!of these species. Hence the sum of the atomic abundances does not equal
!the total number of moles of molecules: N_Mg + N_Si + ... is not = N_tot
   END SUBROUTINE get_shell_contents

!#######################################################################
   SUBROUTINE update_shell_gradients(self)

      type(shell), intent(inout) :: self
      integer :: i, n_mats

      n_mats = size(self%contents)

      call gradients(grads=self%gradients, &
                     r=self%radius, &
                     y=self%integration_parameters, &
                     fractions=self%fractions, &
                     nmat=size(self%fractions), &
                     ll=self%contents, &
                     gammaG0=self%gammaG0, &
                     tempType=self%tempType, &
                     q=self%q, &
                     d0=self%rho0, &
                     adiabatType=self%adiabatType, &
                     eps_T_zero=self%eps_T_zero, &
                     xi_Fe=self%Fe_number, &
                     eps_H2O=self%eps_H2O, &
                     eps_Al=self%eps_Al, &
                     omega=self%omega, &
                     phi=0d0, &
                     weight_fractions=self%weight_fractions, &
                     composition_gradients=self%composition_gradients, &
                     molar_masses=self%mixture%molar_masses, &
                     lay=self%lay, &
                     xi_H=self%xi_H)

   END SUBROUTINE update_shell_gradients

!#######################################################################
   SUBROUTINE update_shell(self, T, P, update_grads, compute_mix)

      type(shell), intent(inout) :: self
      real(8), intent(in), optional :: P, T
      logical, intent(inout), optional :: update_grads, compute_mix
      logical :: update_grads_dummy, compute_mix_dummy
      integer :: i

      if (.not. present(update_grads)) then
         update_grads_dummy = .true.
      else
         update_grads_dummy = update_grads
      end if

      if (.not. present(compute_mix)) then
         compute_mix_dummy = .true.
      else
         compute_mix_dummy = compute_mix
      end if

      if (present(T)) then
         self%integration_parameters(3) = T
      end if

      if (present(P)) then
         self%integration_parameters(1) = P
      end if

!Update the mean mole fractions and weight fractions
      do i = 1, self%n_mats
         self%mean_fractions(i) = self%fractions(i) + self%initials%real_arr(19, i)
         self%mean_fractions(i) = self%mean_fractions(i)/2d0
         self%mean_weight_fractions(i) = self%weight_fractions(i) + self%initials%real_arr(20, i)
         self%mean_weight_fractions(i) = self%mean_weight_fractions(i)/2d0
      end do

      self%mixture%fractions = self%mean_fractions
      self%mixture%weight_fractions = self%weight_fractions

      if (compute_mix_dummy) then
         call compute_mixture(self=self%mixture, T=self%integration_parameters(3), P=self%integration_parameters(1))

         self%dPdrho = self%mixture%dPdrho
         self%integration_parameters(4) = self%mixture%dens
      end if

      if (update_grads_dummy) then
         call update_shell_gradients(self=self)
      end if

      call get_shell_contents(self=self)

      self%gravity = G*self%integration_parameters(2)/self%radius**2

   END SUBROUTINE update_shell

!#######################################################################
   SUBROUTINE reset_shell(self)

      type(shell), intent(inout) :: self

      call init_shell(self=self, n_mats=size(self%contents), alloc=.false.)
      self%status = 'bare'

!call update_shell_gradients(self=self)

   END SUBROUTINE reset_shell

!#######################################################################
   SUBROUTINE construct_shell(self, overconstruct, dr)

      type(shell), intent(inout) :: self
      logical, intent(in), optional :: overconstruct
      logical :: overconstruct_dummy, update_grads, compute_mix
      real(8), dimension(n_params_integration) :: params_dummy
      real(8), intent(in) :: dr
      real(8) :: r_dummy
      integer :: n_mats

      n_mats = size(self%contents)

      if (present(overconstruct)) then
         overconstruct_dummy = overconstruct
      else
         overconstruct_dummy = .false.
      end if

      if (.not. self%status == 'constructed' .or. overconstruct_dummy) then

!For more than two materials the last one can have composition gradient
!Else just juse zero as a place holder during the integration
         if (n_mats .gt. 2) then
            self%integration_parameters(6) = self%weight_fractions(n_mats)
         else
            self%integration_parameters(6) = 0d0
         end if

         r_dummy = self%radius
         ! print *, "constructing shell with dr:", dr
         ! print *, "constructing shell from mass:", self%integration_parameters(2)
         call integrateRK(r_start=r_dummy, &
                          h=dr, &
                          y_in=self%integration_parameters, &
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
                          y_out=params_dummy, &
                          xi_Fe=self%Fe_number, &
                          eps_H2O=self%eps_H2O, &
                          eps_Al=self%eps_Al, &
                          omega=self%omega, &
                          phi=0d0, &
                          weight_fractions=self%weight_fractions, &
                          composition_gradients=self%composition_gradients, &
                          molar_masses=self%mixture%molar_masses, &
                          lay=self%lay, &
                          xi_H=self%xi_H)

         self%indigenous_mass = params_dummy(2) - self%integration_parameters(2)
         self%integration_parameters = params_dummy

         self%dr = dr

         self%volume = 4d0/3d0*PI*(self%radius**3 - (self%radius - dr)**3)

         update_grads = .false.
         compute_mix = .false.
         ! print *, "constructed shell to mass:", self%integration_parameters(2)
         call update_shell(self=self, update_grads=update_grads, &
                           compute_mix=compute_mix)

         self%status = 'constructed'
         ! print *, "constructed shell to mass:", self%integration_parameters(2)

      else
         print *, 'WARNING: This shell has already been constructed.'
         print *, 'Pass overconstruct=.true. to ignore this message.'
      end if

   END SUBROUTINE construct_shell

!#######################################################################
   SUBROUTINE merge_shells(self, other)

      type(shell), intent(inout) :: self
      type(shell), intent(in) :: other

      self%initials%real_vals(1) = other%radius
      self%initials%real_vals(2) = other%integration_parameters(3)
      self%initials%real_vals(3) = other%integration_parameters(2)
      self%initials%real_vals(4) = other%integration_parameters(1)
      self%initials%real_vals(5) = other%integration_parameters(4)
      self%initials%int_vals(6) = other%lay
      self%initials%int_vals(7) = other%tempType
      self%initials%real_vals(8) = other%dPdrho
      self%initials%real_vals(9) = other%gammaG0
      self%initials%real_vals(10) = other%rho0
      self%initials%real_arr(11, :) = other%gradients(:)
      self%initials%real_vals(12) = other%Fe_number
      self%initials%real_vals(13) = other%Si_number
      self%initials%real_vals(14) = other%eps_H2O
      self%initials%real_vals(15) = other%eps_Al
      self%initials%real_vals(16) = other%integration_parameters(5)
      self%initials%real_vals(17) = other%omega
      self%initials%real_vals(18) = other%xi_H
      self%initials%real_vals(19) = other%integration_parameters(7)
      self%initials%real_vals(20) = other%integration_parameters(8)
      self%initials%real_arr(19, 1:self%n_mats) = other%fractions(1:self%n_mats)
      self%initials%real_arr(20, 1:self%n_mats) = other%weight_fractions(1:self%n_mats)

      call init_shell(self=self, alloc=.false., fractions=other%fractions, &
                      contents=other%contents, n_mats=size(other%contents), &
                      integration_parameters=other%integration_parameters)

   END SUBROUTINE merge_shells

END MODULE class_shell

!#######################################################################

