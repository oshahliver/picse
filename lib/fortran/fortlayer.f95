
MODULE class_layer

   use class_shell

   implicit none

   type layer

      integer :: lay, tempType, rhoType, adiabatType
      integer :: shell_count = 0
      real(8) :: r_in, r_out, eps_r, radius, FeMg, SiMg, Fe_number, X_H2O, &
                 Si_number, omega, xi_H, xi_Stv, X_impurity_slope, X_impurity_0
      real(8) :: eps_H2O, eps_Al, N_tot, N_Mg, N_Si, N_Al, N_H2O, N_Fe, N_S
      real(8) :: rho0, q, eps_T_zero
      real(8) :: mass_should, indigenous_mass, dr, gammaG0
      real(8), dimension(n_params_integration) :: gradients
      integer, dimension(:), allocatable :: contents, YMg, YSi, YO, YS
      real(8), dimension(:), allocatable :: fractions, composition_gradients
      logical, dimension(:), allocatable :: saturation
      logical :: bisec = .false.
      logical :: shell_iteration = .true.
      integer :: shell_iteration_count = 0
      logical :: change_bisec = .true., major_complete = .false.
      type(shell), dimension(:), allocatable :: shells
      logical :: force_bisection = .false., overshoot = .false., &
                 major_overshoot = .false., minor_overshoot = .false., constant_Fe = .true.

   end type layer

contains

   SUBROUTINE init_layer(self, contents, fractions, n_mats, r_in, &
                         tempType, rhoType, adiabatType, q, gammaG0, eps_r, rho0, &
                         eps_T_zero, n_shells, eps_H2O, lay, eps_Al, Si_number, Fe_number, &
                         omega, xi_H, xi_Stv, X_impurity_0, X_impurity_slope, integration_parameters)

      type(layer), intent(inout) :: self
      integer, intent(in) :: n_mats, n_shells, lay
      real(8), dimension(n_mats), intent(in) :: fractions
      integer, dimension(n_mats), intent(in) :: contents
      real(8), intent(in) :: r_in, q, gammaG0, eps_r, rho0, &
                             eps_T_zero, eps_Al
      real(8), intent(in) :: Fe_number, eps_H2O, Si_number
      integer, intent(in) :: tempType, rhoType, adiabatType
      real(8), intent(in), optional :: omega, xi_H, xi_Stv, X_impurity_0, &
                                       X_impurity_slope
      real(8), dimension(n_params_integration) :: params
      real(8), dimension(n_params_integration), intent(in) :: integration_parameters
      integer :: i

      allocate (self%contents(n_mats))
      allocate (self%fractions(n_mats))
      allocate (self%composition_gradients(n_mats))
      allocate (self%shells(n_shells))
      allocate (self%YMg(n_mats))
      allocate (self%YSi(n_mats))
      allocate (self%YO(n_mats))
      allocate (self%YS(n_mats))

      self%shell_count = 1
      self%lay = lay
      self%contents = contents
      self%fractions = fractions
      self%X_impurity_0 = 0d0
      self%X_impurity_slope = 0d0

      if (present(X_impurity_0)) then
         self%X_impurity_0 = X_impurity_0
      end if

      if (present(X_impurity_slope)) then
         self%X_impurity_slope = X_impurity_slope
      end if

      self%composition_gradients(:) = 0d0
      self%composition_gradients(n_mats) = self%X_impurity_slope

      do i = 1, n_mats
         self%YO(i) = material_YO(self%contents(i))
         self%YMg(i) = material_YMg(self%contents(i))
         self%YSi(i) = material_YSi(self%contents(i))
         self%YS(i) = material_YS(self%contents(i))
      end do

      if (present(omega)) then
         self%omega = omega
      else
         self%omega = 0d0
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

      self%r_in = r_in
      self%radius = r_in
      self%q = q
      self%rho0 = rho0
      self%gammaG0 = gammaG0
      self%eps_r = eps_r
      self%eps_T_zero = eps_T_zero
      self%rhoType = rhoType
      self%tempType = tempType
      self%adiabatType = adiabatType
      self%Fe_number = Fe_number
      self%Si_number = Si_number
      self%eps_H2O = eps_H2O
      self%eps_Al = eps_Al

      if (.not. Si_number == 1.0d0) then
         self%SiMg = 1.0d0/(1.0d0 - self%Si_number)
      else
         self%SiMg = 1.0d10
      end if

      if (.not. Fe_number == 1.0d0) then
         self%FeMg = 1.0d0/(1.0d0 - self%Fe_number)
      else
         self%FeMg = 1.0d10
      end if

      call init_shell(self=self%shells(1), contents=contents, &
                      fractions=fractions, tempType=tempType, adiabatType=adiabatType, &
                      q=q, gammaG0=gammaG0, eps_T_zero=eps_T_zero, alloc=.true., &
                      eps_H2O=eps_H2O, eps_Al=self%eps_Al, Fe_number=self%Fe_number, &
                      n_mats=n_mats, lay=lay, r=r_in, Si_number=self%Si_number, &
                      omega=self%omega, xi_H=self%xi_H, xi_Stv=self%xi_Stv, &
                      composition_gradients=self%composition_gradients, integration_parameters = integration_parameters)

      call update_layer(self=self)

   END SUBROUTINE init_layer

!#######################################################################
   SUBROUTINE update_layer(self)

      type(layer), intent(inout) :: self
      real(8), dimension(n_params_integration) :: length_scales, params, gradients
      integer :: i
      ! print *, "Updating layer with shell:", self%shell_count
      self%radius = self%shells(self%shell_count)%radius
      if (self%shells(self%shell_count)%force_bisection) then
         self%force_bisection = .true.
      end if

      ! Extract integration parameters from the last shell
      do i=1, n_params_integration
         params(i) = self%shells(self%shell_count)%integration_parameters(i)
      enddo

      self%fractions = self%shells(self%shell_count)%fractions
      gradients = self%shells(self%shell_count)%gradients
      self%indigenous_mass = 0.0d0
      
      do i = 1, self%shell_count
         self%indigenous_mass = self%indigenous_mass + &
                                self%shells(i)%indigenous_mass
      end do

      if (.not. self%bisec) then
         do i = 1, size(gradients)
            length_scales(i) = abs(params(i) / gradients(i))
         end do

         !Ommit MOI for calculation of dr
         self%dr = minval(length_scales(1:4))*self%eps_r
         ! print *, "New dr in layer =", self%dr
      end if

      ! print *, "Updating successful."
   END SUBROUTINE update_layer

!#######################################################################
   SUBROUTINE construct_layer(self)

      type(layer), intent(inout) :: self
      real(8) :: X0, X1
      real(8), dimension(n_params_integration) :: params
      logical :: alloc
      integer :: i
      ! print *, "constructing shell in layer:", self%shell_count
!If layer bisection is ongoing and the constraint is currently not
!overshot that means that the integration step size has been reduces
!sufficiently to actually add one more useful shell before the layer
!transition. In this case, the next shell should start the integration
!from exactly this point. For this reason it is initially set to have
!the same properties as the previously integrated shell by merging them.
      if (self%bisec .and. .not. self%overshoot) then
         ! print *, "merging shells:", self%shell_count, self%shell_count -1 
         call merge_shells(self=self%shells(self%shell_count), &
                           other=self%shells(self%shell_count - 1))
      end if

      ! print *, "constructing shell in layer with idx:", self%shell_count
      ! print *, "Constructing layer with dr =", self%dr
      ! print *, "The shell that will be constructed has mass:", self%shells(self%shell_count)%integration_parameters(2)
      call construct_shell(self=self%shells(self%shell_count), dr=self%dr)
      ! print *, "The shell that was constructed has mass:", self%shells(self%shell_count)%integration_parameters(2)
      if (self%shells(self%shell_count)%force_bisection) then
         self%force_bisection = .true.
      end if

      call update_shell(self=self%shells(self%shell_count))
      ! print *, "The shell that was constructed has mass:", self%shells(self%shell_count)%integration_parameters(2)
      call update_layer(self=self)
      
      ! print *, "The shell that was constructed has mass:", self%shells(self%shell_count)%integration_parameters(2)
      self%r_out = self%radius

      ! TODO. Remove this if possible (think it through again!)
      ! Since in overshoot, the pre-initiated shell is now completely
      ! overwritten with an empty shell, it should never be required
      ! to use alloc = .false.
      if (.not. self%force_bisection) then
         if (self%overshoot) then
            alloc = .true.
         else
            alloc = .true.
         end if
         
         ! The integration parameters need to be extracted after the shell is updated
         ! Extract integration parameters from the last shell
         do i=1, n_params_integration
               params(i) = self%shells(self%shell_count)%integration_parameters(i)
         enddo
       
         !Here the fractions are not updated yet. However, the fraction of the
         !third component must already be given here. Thus, if the impurity
         !fraction is not constant the new value must be evaluated after
         !the last shell has been constructed but prior to initiating the new
         !shell. This is done in the update_layer call above.
         call init_shell(self=self%shells(self%shell_count + 1), &
                         contents=self%contents, &
                         fractions=self%fractions, &
                         r=self%radius, &
                         tempType=self%tempType, &
                         lay=self%lay, &
                         gammaG0=self%gammaG0, &
                         rho0=self%rho0, &
                         adiabatType=self%adiabatType, &
                         q=self%q, &
                         eps_T_zero=self%eps_T_zero, &
                         n_mats=size(self%contents), &
                         eps_H2O=self%eps_H2O, &
                         eps_Al=self%eps_Al, &
                         alloc=alloc, &
                         Fe_number=self%Fe_number, &
                         Si_number=self%Si_number, &
                         omega=self%omega, &
                         xi_H=self%xi_H, &
                         xi_Stv=self%xi_Stv, &
                         composition_gradients=self%composition_gradients, &
                         integration_parameters = params)
      end if

      ! Increment shell count by one
      self%shell_count = self%shell_count + 1
      ! print *, "Constructing shell successfull"
   END SUBROUTINE construct_layer

END MODULE class_layer

!#######################################################################
