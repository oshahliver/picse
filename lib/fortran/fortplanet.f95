MODULE class_planet

   use run_params
   use constants
   use phase
   use class_layer
   use class_weird_array

   implicit none

   type planet

      integer :: n_layers, n_mats, M, n_shells
      integer :: lay = 1, direction
      type(mixture), dimension(:), allocatable :: ambient_mixtures
      type(weird_array) :: contents, fractions, YMg_layers, YSi_layers, &
                           YO_layers, YS_layers, YH_layers
      real(8), dimension(:, :), allocatable :: profiles
      type(layer), dimension(:), allocatable :: layers
      real(8) :: R_surface_is, R_surface_should, M_surface_is, &
                 M_surface_should, P_surface_is, P_surface_should, T_surface_is, &
                 T_surface_should, T_center, P_center, Mg_number_should, Mg_number_is, &
                 Si_number_is, Si_number_should, R_seed, eps_r, M_ocean_is, M_ocean_should, &
                 M_H2O_is, M_H2O_should, M_DWR_is, constraint_value_is, N_Fe, N_Mg, N_Si, &
                 N_H2O, N_Al, N_O, constraint_value_should, eps_layer, eps_Al, eps_H2O, &
                 N_H, N_S, T_center_is, P_center_is, xi_H_core_predicted
      integer :: tempType, rhoType, adiabatType, layerType, subphase_res
      real(8) :: inner_core_fraction, MOI_is, omega, P_H2, xi_H_core, &
                 eps_T_zero, ocean_frac_should, ocean_frac_is, gravity, &
                 M_core_should, rho_mean, Si_number_mantle, Fe_number_mantle, xi_Stv, &
                 T_ICB, T_CS, P_CS, T_MTZ, P_MTZ, E_grav, E_int
      real(8), dimension(:), allocatable :: xi_all_core, X_all_core
      real(8), dimension(:), allocatable :: layer_masses, GammaG0_layers, &
                                            temp_jumps, layer_radii, layer_temps, layer_pres, q_layers, rho0_layers, &
                                            Si_number_layers, Fe_number_layers, fraction_procs_layers, &
                                            xi_H_layers, xi_Stv_layers, X_impurity_slope_layers, X_impurity_0_layers
      character(len=20) :: status, major_constraint
      integer, dimension(:), allocatable :: layer_constraints, n_shells_layers, &
                                            layer_dims
      real(8), dimension(5, 3) :: external_temp_profile
      logical :: silence = .false.
      logical :: constant_Fe = .true.
      integer :: shell_iteration_count = 1
      integer :: layer_iteration_count = 1
      logical :: change_bisec = .true., bisec = .false.
      logical :: shell_iteration = .true.
      logical :: layer_iteration = .true.
      logical :: layer_complete = .false.
      logical :: major_complete = .false.
      logical :: echo = .false.
      logical :: subphase_refine_inhibit = .false.
      logical :: subphase_refine_done = .false.
      logical :: skip_bisec = .false.
      logical :: mantle_exists = .true.
      logical :: inner_core_exists = .true.
      logical :: outer_core_exists = .true.
      logical :: core_segregation_model = .true.
      logical :: inner_core_segregation_model = .false.

   end type planet

contains

   SUBROUTINE init_planet(self, T_center, R_seed, contents, tempType, &
                          rhoType, adiabatType, fractions, layerType, P_surface_should, &
                          T_surface_should, R_surface_should, M_surface_should, eps_r, &
                          layer_masses, layer_pres, X_H2O, echo, Mg_number_should, &
                          Si_number_should, temp_jumps, q_layers, GammaG0_layers, &
                          M_core_should, n_layers, layer_radii, P_center, eps_T_zero, &
                          major_constraint, layer_constraints, Si_number_layers, &
                          Fe_number_layers, omega, xi_H_core_predicted, &
                          subphase_res, xi_Stv, X_impurity_0_layers, X_impurity_slope_layers, &
                          xi_all_core, X_all_core, P_CS, core_segregation_model, &
                          M_ocean_should, inner_core_segregation_model, external_temp_profile)

      type(planet), intent(inout) :: self
      integer, intent(in) :: n_layers
      type(weird_array), intent(in) :: contents, fractions
      type(weird_array) :: weight_fractions
      real(8), intent(in) :: R_seed, P_center, T_center, eps_r, &
                             eps_T_zero, P_CS
      real(8) :: eps_H2O, eps_Al
      real(8), intent(in), optional :: X_H2O
      real(8), intent(in), optional :: Si_number_should, omega, &
                                       xi_H_core_predicted, xi_Stv
      real(8), intent(in), dimension(5) :: xi_all_core, X_all_core
      real(8), intent(in) :: M_core_should, R_surface_should, &
                             P_surface_should, T_surface_should
      real(8), intent(in) :: M_surface_should, M_ocean_should, &
                             Mg_number_should
      real(8), dimension(n_layers), intent(in) :: layer_masses, &
                                                  layer_pres, layer_radii, temp_jumps, GammaG0_layers, q_layers
      real(8), dimension(n_layers), intent(in) :: Si_number_layers, &
                                                  Fe_number_layers
      real(8), dimension(n_layers), intent(in), optional :: X_impurity_0_layers, &
                                                            X_impurity_slope_layers
      integer, intent(in) :: tempType, rhoType, adiabatType, layerType
      logical, intent(in), optional :: echo
      real(8) :: M_seed
      integer, intent(in), optional :: subphase_res
      type(mixture) :: seed_mixture
      integer, dimension(n_layers), intent(in) :: layer_constraints
      character(len=20), intent(inout), optional :: major_constraint
      integer :: i, j, k, c
      logical :: alloc = .true.
      integer :: layer_dims(n_layers)
      integer :: core_segregation_model, inner_core_segregation_model
      real(8), dimension(5, 3), intent(in), optional :: external_temp_profile

      allocate (self%layer_dims(n_layers))
      allocate (self%ambient_mixtures(n_layers))
      allocate (self%layers(n_layers))
      allocate (self%layer_constraints(n_layers))
      allocate (self%layer_masses(n_layers))
      allocate (self%layer_radii(n_layers))
      allocate (self%layer_pres(n_layers))
      allocate (self%rho0_layers(n_layers))
      allocate (self%q_layers(n_layers))
      allocate (self%gammaG0_layers(n_layers))
      allocate (self%layer_temps(n_layers))
      allocate (self%temp_jumps(n_layers))
      allocate (self%Si_number_layers(n_layers))
      allocate (self%Fe_number_layers(n_layers))
      allocate (self%fraction_procs_layers(n_layers))
      allocate (self%xi_H_layers(n_layers))
      allocate (self%xi_Stv_layers(n_layers))
      allocate (self%n_shells_layers(n_layers))
      allocate (self%X_impurity_0_layers(n_layers))
      allocate (self%X_impurity_slope_layers(n_layers))

      if (present(external_temp_profile)) then
         self%external_temp_profile = external_temp_profile

      else
         self%external_temp_profile = &
            reshape((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                      0d0, 0d0, 0d0, 0d0/), shape(external_temp_profile))
      end if

      if (present(subphase_res)) then
         self%subphase_res = subphase_res
      else
         self%subphase_res = 16
      end if

      if (present(major_constraint)) then
         self%major_constraint = major_constraint
      else
         self%major_constraint = 'P_surface'
      end if

! if(present(P_surface_should))then
!   self%P_surface_should = P_surface_should
! else
!   self%P_surface_should = 1.0d5
! endif

! if(present(T_surface_should))then
!   self%T_surface_should = T_surface_should
! else
!   self%T_surface_should = 300.0d0
! endif

      self%M_surface_should = M_surface_should
      self%M_ocean_should = M_ocean_should
      self%Mg_number_should = Mg_number_should
      self%P_surface_should = P_surface_should
      self%T_surface_should = T_surface_should

      if (present(omega)) then
         self%omega = omega
      else
         self%omega = 0d0
      end if

      if (present(xi_H_core_predicted)) then
         self%xi_H_core_predicted = xi_H_core_predicted
      else
         self%xi_H_core_predicted = 0d0
      end if

      if (present(xi_Stv)) then
         self%xi_Stv = xi_Stv
      else
         self%xi_Stv = 0d0
      end if

      if (present(X_impurity_0_layers)) then
         self%X_impurity_0_layers = X_impurity_0_layers
      else
         do i = 1, n_layers
            self%X_impurity_0_layers(i) = 0d0
         end do
      end if

      if (present(X_impurity_slope_layers)) then
         self%X_impurity_slope_layers = X_impurity_slope_layers
      else
         do i = 1, n_layers
            self%X_impurity_slope_layers(i) = 0d0
         end do
      end if

!currently hydrogen is only dissolved in the iron core
!Note: The iron content in any given material is meant with respect
!to Mg, i.e. Mg_(1-xi_Fe) + Fe_(xi_Fe). In the cas of FeH there is no
!Mg and hence xi_Fe must be zero from a compositional viewpoint.
!However, in order to compute the EoS parameters correctly at given
!hydrogen content, the table interpolation needs the Fe number of FeH_x
!with respect to hydrogen, i.e. Fe_(xi_Fe) + H_(1-xi_Fe). It is given
!by (1-xi_H). Hence when the table interpolation is called xi_Fe = 1-xi_H
!needs to be passed. But to compute the contents in the shell xi_Fe = 1 in
!order to avoid Mg to be counted in the core.
      do i = 1, n_layers
         if (i == 2) then
            self%xi_H_layers(i) = self%xi_H_core_predicted
            self%xi_Stv_layers(i) = 0d0

            !In the solid inner core there are currently no lighter elements
         elseif (i == 1) then
            self%xi_H_layers(i) = 0d0
            self%xi_Stv_layers(i) = 0d0
         else
            self%xi_H_layers(i) = 0d0
            self%xi_Stv_layers(i) = self%xi_Stv
         end if
      end do

      eps_H2O = 0e0
      eps_Al = 0e0

      self%P_center_is = P_center
      self%T_center_is = T_center
      self%eps_r = eps_r
      self%adiabatType = adiabatType
      self%rhoType = rhoType
      self%tempType = tempType
      self%eps_T_zero = eps_T_zero
      self%temp_jumps = temp_jumps
      self%layer_masses = layer_masses
      self%layer_radii = layer_radii
      self%layer_pres = layer_pres
      self%q_layers = q_layers
      self%gammaG0_layers = gammaG0_layers
      self%Fe_number_layers = Fe_number_layers
      self%Si_number_layers = Si_number_layers
      self%eps_H2O = eps_H2O
      self%eps_Al = eps_Al
      self%N_Mg = 0.0d0
      self%N_Fe = 0.0d0
      self%N_Si = 0.0d0
      self%N_H2O = 0.0d0
      self%N_Al = 0.0d0
      self%N_O = 0.0d0
      self%N_S = 0.0d0
      self%N_H = 0d0
      self%n_shells = 0
      self%ocean_frac_is = -10d0
      self%core_segregation_model = core_segregation_model
      self%inner_core_segregation_model = inner_core_segregation_model
      self%contents = contents
      self%fractions = fractions
      self%layer_constraints = layer_constraints
!Individual mol and wt fractions of the elements in the liquid part
!of the core

      self%xi_all_core = xi_all_core
      self%X_all_core = X_all_core

!~ print *, '######################################'
!~ print *, 'fractions in planet init =' ,self%fractions%axes(2)%real_array
      self%status = 'shadow of the future'

!Compute the melting pressure of iron at current temperature
      self%T_ICB = T_melt_iron(P=P_center, &
                               X_Si=self%X_all_core(4), &
                               X_O=self%X_all_core(5), &
                               X_S=self%X_all_core(3))

!~ print *, 'T_melt, P_melt (GPa) =', self%T_ICB, P_center*1d-9
!~ print *, 'T_center =', T_center
!If ICB already exceeded directly jump to outer core
      if (self%T_ICB < T_center) then
!~ print *, 'Go directly to outer core'
         self%lay = 2
         self%inner_core_exists = .false.
      end if

!Compute pressure at IM-OM transition
      self%P_MTZ = P_Ol_Pv(T_center)

      do i = 1, n_layers
         self%layer_dims(i) = size(contents%axes(i)%int_array)
      end do

      call init_weird_array(self=self%YMg_layers, ndim=n_layers, &
                            axes_lengths=self%layer_dims)

      call init_weird_array(self=self%YSi_layers, ndim=n_layers, &
                            axes_lengths=self%layer_dims)

      call init_weird_array(self=self%YO_layers, ndim=n_layers, &
                            axes_lengths=self%layer_dims)

      call init_weird_array(self=self%YS_layers, ndim=n_layers, &
                            axes_lengths=self%layer_dims)

!Core segregation pressure
      self%P_CS = P_CS
!Compute core segregation temperature assuming chemical equilibration
!to have taken place at the pyrolite liquidus
      self%T_CS = T_liquidus_pyrolite(self%P_CS)
!~ print *, 'T_CS (K) / P_CS (GPa)=', self%T_CS, self%P_CS*1d-9
!Update contents to ensure the right eos tables are used in case
!the Fe content in the mantle is chosen to be variable
!!if (.not.self%constant_Fe)then
!!        do i=1, n_layers
!!                do j=1, layer_dims(i)
!!                        c=self%contents%axes(i)%int_array(j)
!!                        self%contents%axes(i)%int_array(j)=eos_translation(c)
!!                enddo
!!        enddo
!!endif

!~ do i=1, n_layers
!~         print *, 'fractions ',i, ' =', self%fractions%axes(i)%real_array
!~ enddo
      do i = 1, n_layers
!~ print *, "init mixture for layer", i
!~ print *, "conts =", self%contents%axes(i)%int_array
!~ print *, "fracs =", self%fractions%axes(i)%real_array
!~ print *, "stuff =", eps_H2O, eps_Al, self%Si_number_layers(i), self%Fe_number_layers(i)
         call init_mixture(self=self%ambient_mixtures(i), &
                           n_mats=size(self%contents%axes(i)%int_array), &
                           contents=self%contents%axes(i)%int_array, &
                           fractions=self%fractions%axes(i)%real_array, &
                           P=1.0d4, &
                           T=300.0d0, &
                           eps_H2O=eps_H2O, &
                           eps_Al=eps_Al, &
                           Si_number=self%Si_number_layers(i), &
                           Fe_number=self%Fe_number_layers(i), &
                           alloc=alloc, &
                           xi_H=self%xi_H_layers(i), &
                           xi_Stv=self%xi_Stv_layers(i))
         call compute_mixture(self=self%ambient_mixtures(i))
         call update_mixture_fractions(self=self%ambient_mixtures(i))
         self%rho0_layers(i) = self%ambient_mixtures(i)%dens

         do j = 1, self%layer_dims(i)
            self%YMg_layers%axes(i)%int_array(j) = material_YMg(self%contents%axes(i)%int_array(j))
            self%YO_layers%axes(i)%int_array(j) = material_YO(self%contents%axes(i)%int_array(j))
            self%YSi_layers%axes(i)%int_array(j) = material_YSi(self%contents%axes(i)%int_array(j))
            self%YS_layers%axes(i)%int_array(j) = material_YS(self%contents%axes(i)%int_array(j))
         end do

      end do

      call init_mixture(self=seed_mixture, &
                        n_mats=size(self%contents%axes(1)%int_array), &
                        contents=self%contents%axes(1)%int_array, &
                        fractions=self%fractions%axes(1)%real_array, &
                        P=P_center, &
                        T=T_center, &
                        eps_H2O=eps_H2O, &
                        eps_Al=eps_Al, &
                        Si_number=self%Si_number_layers(1), &
                        Fe_number=self%Fe_number_layers(1), &
                        xi_H=self%xi_H_core_predicted, &
                        xi_Stv=self%xi_Stv_layers(1))

      call compute_mixture(self=seed_mixture)
!call print_mixture(self=seed_mixture)
!initialize the first layer

      M_seed = 4.0d0/3.0d0*PI*R_seed**3*seed_mixture%dens
      self%MOI_is = 2d0/5d0*M_seed*R_seed**2

      do i = 1, self%lay

         call init_layer(self=self%layers(i), &
                         n_mats=size(self%contents%axes(i)%int_array), &
                         contents=self%contents%axes(i)%int_array, &
                         fractions=self%fractions%axes(i)%real_array, &
                         r_in=R_seed, &
                         m=M_seed, &
                         T=T_center, &
                         P=P_center, &
                         tempType=self%tempType, &
                         rhoType=self%rhoType, &
                         adiabatType=self%adiabatType, &
                         q=q_layers(i), &
                         gammaG0=self%gammaG0_layers(i), &
                         eps_r=self%eps_r, &
                         rho0=self%rho0_layers(i), &
                         eps_T_zero=self%eps_T_zero, &
                         n_shells=n_shells_max_layer, &
                         eps_H2O=self%eps_H2O, &
                         eps_Al=self%eps_Al, &
                         lay=i, &
                         Si_number=self%Si_number_layers(self%lay), &
                         Fe_number=self%Fe_number_layers(self%lay), &
                         MOI=self%MOI_is, &
                         omega=self%omega, &
                         xi_H=self%xi_H_layers(self%lay), &
                         xi_Stv=self%xi_Stv_layers(self%lay), &
                         X_impurity_0=self%X_impurity_0_layers(self%lay), &
                         X_impurity_slope=self%X_impurity_slope_layers(self%lay), &
                         external_temp_profile=self%external_temp_profile(self%lay, :))

         call update_layer(self=self%layers(i))

      end do

      self%status = 'seed'

   END SUBROUTINE init_planet

!#######################################################################
   SUBROUTINE get_integration_time(self)

      type(planet), intent(inout) :: self
      integer :: i, j, n_shells
      real(4) :: timer

      do i = 1, size(self%layers)
         do j = 1, size(self%layers(i)%shells)
            n_shells = n_shells + 1
            timer = timer + self%layers(i)%shells(j)%timer
         end do
      end do

      print *, 'Total time spent in integrate:', timer, ' sec'
      print *, 'Average time spent per fct call:', timer/real(n_shells), ' sec'

   END SUBROUTINE get_integration_time

!#######################################################################
   SUBROUTINE compute_E_grav(self)
!Computes gravitational energy contribution of each shell approximating the
!dnesity as a linear function between the bottom and the top of each shell.
!The total gravitational energy of the planet is obtained by summing up
!all contributions of the individual shells.

      type(planet), intent(inout) :: self
      real(8) :: dU, rho1, r, dr, ks
      integer :: i, j

      self%E_grav = 0d0

!Iterate over all layers
      do i = 1, self%lay
         !Iterate over individual shells within a layer
         do j = 1, self%layers(i)%shell_count
            self%E_grav = self%E_grav + self%layers(i)%shells(j)%dE_grav
         end do
      end do

   END SUBROUTINE compute_E_grav

!#######################################################################
   SUBROUTINE compute_E_int(self)

      type(planet), intent(inout) :: self
      real(8) :: dU, rho1, r, dr, ks
      integer :: i, j

      self%E_int = 0d0

!Iterate over all layers
      do i = 1, self%lay
         !Iterate over individual shells within a layer
         do j = 1, self%layers(i)%shell_count
            self%E_int = self%E_int + self%layers(i)%shells(j)%dE_int
         end do
      end do

   END SUBROUTINE compute_E_int

!#######################################################################
   ! SUBROUTINE compute_MOI(self)

   !    type(planet), intent(inout) :: self
   !    real(8) :: dr, r1, r2, rho1, rho2, moi
   !    integer :: i, j

   !    self%MOI_is = 2d0/5d0*self%layers(1)%shells(1)%mass* &
   !                  self%layers(1)%shells(1)%radius**2

   !    do i = 1, self%lay
   !       do j = 2, self%layers(i)%shell_count

   !          r1 = self%layers(i)%shells(j - 1)%radius
   !          rho1 = self%layers(i)%shells(j - 1)%dens

   !          r2 = self%layers(i)%shells(j)%radius
   !          rho2 = self%layers(i)%shells(j)%dens

   !          moi = MOI_integrand_linear(r1=r1, r2=r2, rho1=rho1, rho2=rho2)
   !          self%MOI_is = self%MOI_is + 8.0d0/3.0d0*PI*moi

   !       end do
   !    end do

   !    self%MOI_is = self%MOI_is/(self%M_surface_is*self%R_surface_is**2)

   ! END SUBROUTINE compute_MOI

!#######################################################################
   SUBROUTINE get_profiles(self)

      type(planet), intent(inout) :: self
      integer :: n, i, j
      real(8) :: grav, ener

      n = 0
      do i = 1, self%lay
         !Add one additional point because in each layer the first shell needs
         !to contribute at the bottom and the top to avoid gaps in the plots.
         !Don't do it for the inner core.
         n = n + self%layers(i)%shell_count
         if (i > 1) then
            n = n + 1
         end if
      end do

      if (.not. allocated(self%profiles)) then
         allocate (self%profiles(8, n))
      end if

      n = 1
      do i = 1, self%lay
         !Add a point at the bottom of each new layer containing the parameter
         !values at the bottom of the shell as otherwise there would be a gap
         !in the plots. Don't do it for the inner core.
         if (i > 1) then
            grav = G
            grav = grav*self%layers(i)%shells(1)%initials%real_vals(3)
            grav = grav/self%layers(i)%shells(1)%initials%real_vals(1)**2
            self%profiles(1, n) = self%layers(i)%shells(1)%initials%real_vals(1)
            self%profiles(2, n) = self%layers(i)%shells(1)%initials%real_vals(2)
            self%profiles(3, n) = self%layers(i)%shells(1)%initials%real_vals(4)
            self%profiles(4, n) = self%layers(i)%shells(1)%initials%real_vals(5)
            self%profiles(5, n) = self%layers(i)%shells(1)%initials%real_vals(3)
            self%profiles(6, n) = grav
            self%profiles(7, n) = self%layers(i)%shells(1)%initials%real_vals(16)
            self%profiles(8, n) = self%profiles(8, n - 1)
            n = n + 1
         end if

         do j = 1, self%layers(i)%shell_count
            self%profiles(1, n) = self%layers(i)%shells(j)%radius
            self%profiles(2, n) = self%layers(i)%shells(j)%temp
            self%profiles(3, n) = self%layers(i)%shells(j)%pres
            self%profiles(4, n) = self%layers(i)%shells(j)%dens
            self%profiles(5, n) = self%layers(i)%shells(j)%mass
            self%profiles(6, n) = self%layers(i)%shells(j)%gravity
            self%profiles(7, n) = self%layers(i)%shells(j)%MOI
            self%profiles(8, n) = self%layers(i)%shells(j)%dE_grav

            if (N > 1) then
               self%profiles(8, n) = self%profiles(8, n) + self%profiles(8, n - 1)
            end if

            n = n + 1
         end do
      end do
   END SUBROUTINE get_profiles

!#######################################################################
   SUBROUTINE compute_P_H2(self, average)
!Compute hydrogen partial pressure at the CMB according to Wu et al. 2018
!Note that in the paper it is stated, that this proceedure is only valid
!for P_H2 < 10MPa. However, Roskosz et al. 2013 have compared the Sievert's
!law to experimental data for nitrogen in Fe up to pressures of 20 GPa
!at temperatures of 2000-3000 K. They find very good agreement between
!the experiments and the predictions from Sivert's law. They find Nitrogen
!contents of up to ~14 wt% or ~40 mol%. At 3000 K and 20 GPa the density
!of iron is ~8.6 gcc. Using the ideal gas law to estimate the order of magnitude
!of the partial pressure P_N at these conditions yields  up to ~ 10³MPa.
!This implies that the Sievert's law extrapoaltes well to higher partial
!gas pressures. It is not clear if this is also true for hydrogen in iron.
!But given the simplicity of both system we assume that the hydrogen
!solubility in Fe extrapolates well to higher pressures just as for nitrogen.

      type(planet), intent(inout) :: self
      logical, intent(in), optional :: average
      logical :: average_dummy
      real(8) :: V_mantle, V_innermost_shell, &
                 N_H2O_innermost_shell, N_H2O_mantle, T_CMB, R_CMB, rho_H2, test
      integer :: i, j, lay
      print *, 'computing pH2 in layer ', self%lay
      if (.not. (present(average))) then
         average_dummy = .false.
      else
         average_dummy = average
      end if
!If mantle exists, compute H2 partial pressure in the lowermost shell
!of the mantle or the average H2 partial pressure in the entire mantle
      if (.not. size(self%layers(2)%shells) .eq. 0) then
         T_CMB = self%layers(2)%temp
         R_CMB = self%layers(2)%radius
         print *, 'Outer core exists'
      else
         print *, 'No outer core exists'
         T_CMB = self%layers(1)%temp
         R_CMB = self%layers(1)%radius
      end if

!Only core exists
      if (self%lay <= 2) then
         self%P_H2 = 0d0
         print *, 'No mantle exists'
         self%mantle_exists = .false.

      else
         ! No lower mantle
         if (self%layers(3)%shells(1)%volume .eq. 0d0) then
            print *, 'Setting lay to 4'
            lay = 4
         else
            lay = 3
         end if
         V_innermost_shell = 4d0/3d0*PI*(self%layers(lay)%shells(1)%radius**3 - &
                                         R_CMB**3)
         N_H2O_innermost_shell = self%layers(lay)%shells(1)%N_H2O
         N_H2O_innermost_shell = self%layers(lay)%shells(1)%indigenous_mass
         !N_H2O_innermost_shell = N_H2O_innermost_shell * 2d-1/(mH*2d0+mO)
         V_mantle = 4d0/3d0*PI*(self%layers(4)%radius**3 - R_CMB**3)
         N_H2O_mantle = 0d0

         !Comptue the average water content in the mantle
         if (average_dummy) then
            do i = 3, 4
               do j = 1, self%layers(i)%shell_count
                  N_H2O_mantle = N_H2O_mantle + self%layers(i)%shells(j)%N_H2O
               end do
            end do

            self%P_H2 = N_H2O_mantle*NA*kB*T_CMB/V_mantle
            rho_H2 = N_H2O_mantle*mH*2d0/V_mantle

            !Compute the partial pressure in the lowermost shell
         else
            self%P_H2 = N_H2O_innermost_shell*NA*kB*T_CMB/ &
                        V_innermost_shell
            rho_H2 = N_H2O_innermost_shell*mH*2d0/V_innermost_shell
         end if
!~         print *, 'rho_H2 =', rho_H2
         !Use van der Waal equation to estimate partial H2 pressure. Note that
         !using the IGL can yield deviations in the resulting hydrogen content
         !in the core of up to ~50% for 0-10 earth mass planets.
         self%P_H2 = rho_H2*Rgas*T_CMB/(2d0*mH - rho_H2*26.61d-6) - rho_H2**2*24.76d-3/(2d0*mH)**2
!~         print *, 'P_H2 =', self%P_H2
!~         print *, 'T_CMB =', T_CMB
         !This is a stupid way to prevent P_H2 < 0 in the case where eps_H2O = 0
         !and large X_impurity. Since impurities are currently only intruduced
         !for the anhydrous case this is kinda enough at this point. (But it's
         !still stupid!)
         self%P_H2 = self%P_H2*self%eps_H2O

      end if

      test = rho_H2*Rgas*T_CMB/(2d0*mH - rho_H2*26.61d-6) - rho_H2**2*24.76d-3/(2d0*mH)**2
!print *, 'N_H2O_innermost_shell =', N_H2O_innermost_shell
!print *, 'V_innermost_shell =', V_innermost_shell
!print *, 'P_H2 =', self%P_H2
!print *, 'test P_H2 =', test
!print *, 'deviation =', (sqrt(test)-sqrt(self%P_H2))/sqrt(self%P_H2)
   END SUBROUTINE compute_P_H2

!#######################################################################
   SUBROUTINE compute_xi_H_core(self)
!Compute molar hydrogen abundance in the core according to Wu et al. 2018

      type(planet), intent(inout) :: self
      real(8) :: T_CMB

      T_CMB = self%layers(2)%temp

      call compute_P_H2(self=self)

!Compute parameter x in Wu et al 2018: (x/2+1) Fe + x/2 H2O
!P0 is 1 bar
      self%xi_H_core = sqrt(self%P_H2/1d5)*exp((-31.8d3 - T_CMB*38.1d0)/(Rgas*T_CMB))

!self%xi_H_core = 2d-1

!Here the actual molar abundance xi_H of H in the core is computed
!With this the core composition is (1-xi_H) Fe + xi_H H
!It is defined as H/(H+Fe)
      self%xi_H_core = self%xi_H_core/(self%xi_H_core + 1d0)!*self%eps_H2O

      self%xi_H_core = 0d0
      print *, 'remember: xi_H_core is set to 0 for Venus!'
   END SUBROUTINE compute_xi_H_core

!#######################################################################
   SUBROUTINE remove_stuff(self)

      type(planet), intent(inout) :: self
      integer :: sh

!Note that when this routine is called, the shell has already been
!constructed which means the next shell is already initiated. The shell
!from which we whish to extract the content is the shell before the
!current shell (the current shell being the one that will be constructed
!in the subsequent integration step)
      sh = self%layers(self%lay)%shell_count - 1

      self%N_Mg = self%N_Mg - self%layers(self%lay)%shells(sh)%N_Mg
      self%N_Al = self%N_Al - self%layers(self%lay)%shells(sh)%N_Al
      self%N_Si = self%N_Si - self%layers(self%lay)%shells(sh)%N_Si
      self%N_Fe = self%N_Fe - self%layers(self%lay)%shells(sh)%N_Fe
      self%N_O = self%N_O - self%layers(self%lay)%shells(sh)%N_O
      self%N_S = self%N_S - self%layers(self%lay)%shells(sh)%N_S
      self%N_H2O = self%N_H2O - self%layers(self%lay)%shells(sh)%N_H2O
      self%N_H = self%N_H - self%layers(self%lay)%shells(sh)%N_H

      if (self%N_Fe .eq. 0.0d0 .and. self%N_Mg .eq. 0.0d0) then
         self%Mg_number_is = 0.0d0
      else
         self%Mg_number_is = self%N_Mg/(self%N_Mg + self%N_Fe)
      end if

      if (self%N_Si .eq. 0.0d0 .and. self%N_Mg .eq. 0.0d0) then
         self%Si_number_is = 0.0d0
      else
         self%Si_number_is = self%N_Si/(self%N_Mg + self%N_Si)
      end if

      self%M_H2O_is = self%N_H2O*mH2O
      self%M_H2O_is = self%M_H2O_is + self%N_H*(0.5d0*mO + mH)

   END SUBROUTINE remove_stuff

!#######################################################################
   SUBROUTINE add_stuff(self)

      type(planet), intent(inout) :: self
      integer :: sh

!Note that when this routine is called, the shell has already been
!constructed which means the shell count has been updated. The shell
!from which we whish to extract the content is the shell before the
!current shell (the current shell being the one that will be constructed
!in the subsequent integration step)
      sh = self%layers(self%lay)%shell_count - 1

      self%N_Mg = self%N_Mg + self%layers(self%lay)%shells(sh)%N_Mg
      self%N_Al = self%N_Al + self%layers(self%lay)%shells(sh)%N_Al
      self%N_Si = self%N_Si + self%layers(self%lay)%shells(sh)%N_Si
      self%N_Fe = self%N_Fe + self%layers(self%lay)%shells(sh)%N_Fe
      self%N_O = self%N_O + self%layers(self%lay)%shells(sh)%N_O
      self%N_S = self%N_S + self%layers(self%lay)%shells(sh)%N_S
      self%N_H2O = self%N_H2O + self%layers(self%lay)%shells(sh)%N_H2O
      self%N_H = self%N_H + self%layers(self%lay)%shells(sh)%N_H

      if (self%N_Fe .eq. 0.0d0 .and. self%N_Mg .eq. 0.0d0) then
         self%Mg_number_is = 0.0d0
      else
         self%Mg_number_is = self%N_Mg/(self%N_Mg + self%N_Fe)
      end if

      if (self%N_Si .eq. 0.0d0 .and. self%N_Mg .eq. 0.0d0) then
         self%Si_number_is = 0.0d0
      else
         self%Si_number_is = self%N_Si/(self%N_Mg + self%N_Si)
      end if

      self%M_H2O_is = self%N_H2O*mH2O
      self%M_H2O_is = self%M_H2O_is + self%N_H*(0.5d0*mO + mH)

   END SUBROUTINE add_stuff

!#######################################################################
   SUBROUTINE update_core_mass(self)
!Updates the core mass for a given bulk composition

      type(planet), intent(inout) :: self
      real(8) :: Q1, Q2, Q3, Q4, Q5, core_frac, FeMg, reldev, M
      real(8), dimension(5) :: molar_masses
      real(8), dimension(4) :: dummy, masses, dummy_1
      real(8), dimension(self%layer_dims(2)) :: wt_fractions_dummy
      integer :: i

      FeMg = (1e0 - self%Mg_number_should)/self%Mg_number_should
      masses = (/mFe, mS, mSi, mO/)
!Compute mole fraction of Mg in the mantle
      Q1 = 0d0
      do i = 1, size(self%fractions%axes(3)%real_array)
         Q1 = Q1 + self%fractions%axes(3)%real_array(i)* &
              (1d0 - self%Fe_number_layers(3))* &
              material_YMg(self%contents%axes(3)%int_array(i))
      end do

!Compute total normalized mass in the mantle
      Q2 = 0d0
      do i = 1, size(self%fractions%axes(3)%real_array)
         Q2 = Q2 + self%fractions%axes(3)%real_array(i)* &
              (1d0 - self%Fe_number_layers(3))* &
              material_YMg(self%contents%axes(3)%int_array(i))*mMg
      end do

      do i = 1, size(self%fractions%axes(3)%real_array)
         Q2 = Q2 + self%fractions%axes(3)%real_array(i)* &
              self%Fe_number_layers(3)* &
              material_YMg(self%contents%axes(3)%int_array(i))*mFe
      end do

      do i = 1, size(self%fractions%axes(3)%real_array)
         Q2 = Q2 + self%fractions%axes(3)%real_array(i)* &
              material_YSi(self%contents%axes(3)%int_array(i))*mSi
      end do

      do i = 1, size(self%fractions%axes(3)%real_array)
         Q2 = Q2 + self%fractions%axes(3)%real_array(i)* &
              material_YO(self%contents%axes(3)%int_array(i))*mO
      end do

      do i = 1, size(self%fractions%axes(3)%real_array)
         Q2 = Q2 + self%fractions%axes(3)%real_array(i)* &
              material_YH(self%contents%axes(3)%int_array(i))*mH
      end do

!Compute mole fraction of Fe in the mantle
      Q3 = 0d0
      do i = 1, size(self%fractions%axes(3)%real_array)
         Q3 = Q3 + self%fractions%axes(3)%real_array(i)* &
              self%Fe_number_layers(3)* &
              material_YMg(self%contents%axes(3)%int_array(i))
      end do

!Compute mole fraction of Fe in the core
      Q4 = self%xi_all_core(1)

!Compute total normalized mass in the core
      molar_masses = (/mFe, mH, mS, mSi, mO/)
      Q5 = 0d0

      do i = 1, size(self%xi_all_core)
         Q5 = Q5 + molar_masses(i)*self%xi_all_core(i)
      end do

!~ print *, 'Q1, Q2, Q3, Q4, Q5, Q6 =', Q1, Q2, Q3, Q4, Q5, 1d0 / mFe
!~ print *, 'Q4 / Q5 =', Q4 / Q5
!~ print *, 'xi all core =', self%xi_all_core
!~ print *, 'factor =', self%layers(1)%indigenous_mass / self%M_surface_should * (1d0 / mFe - Q4 / Q5)
!~ print *, 'M_ocean, M_tot =', self%M_ocean_should, self%M_surface_should
!~ print *, 'Mg mantle, Mg tot =', 1d0 - self%Fe_number_layers(3), self%Mg_number_should
!~ print *, 'fractions =', self%fractions%axes(3)%real_array

!Compute core mass fraction
      core_frac = (1d0 - self%M_ocean_should/self%M_surface_should)
      core_frac = core_frac*(Q3/Q2 - Q1/Q2*FeMg)
      core_frac = core_frac + self%layers(1)%indigenous_mass/self%M_surface_should* &
                  (1e0/mFe - Q4/Q5)
      core_frac = core_frac/(Q3/Q2 - Q4/Q5 - FeMg*Q1/Q2)
!~ print *, 'core_frac =', core_frac

      M = core_frac*self%M_surface_should/m_earth
!~ print *, 'Computed core mass in fortplanet:', M

!Decide which strategy is to be employed to probe the total core mass.
!By default the initial inner core mass in layer_masses is set to the
!total core mass assuming pure iron in the core. The outer core mass
!in layer_masses is set to the core mass if only outer core exists.
      reldev = abs(self%layers(1)%indigenous_mass/m_earth - self%layer_masses(1))/self%layer_masses(1)
!~ print *, 'IC mass =', self%layers(1)%indigenous_mass / m_earth
!~ print *, 'layermasses(1) =', self%layer_masses(1)
!Case 1: Total core mass already reached in inner core. Outer core mass
!must be set to inner core mass so that is skipped.
      if (reldev < eps_layer) then
!~ print *, 'Case 1'
         self%layer_masses(2) = self%layer_masses(1)

!Case 2: Both inner and outer core exist.
      else if (reldev >= eps_layer) then
!~ print *, 'Case 2'
!If redistribution of lighter elements is not accounted for:
!Total core mass must be updated.
!~ print *, 'iCS =', self%inner_core_segregation_model
         if (.not. self%inner_core_segregation_model) then
            self%layer_masses(2) = core_frac*self%M_surface_should/m_earth
!If redistribution of lighter elements is accoutned for:
!Outer core fractions need to be updated.
         else
            !THIS STUFF HERE DOES NOT WORK YET!
!~                 print *, 'layer_dims in update_core_mass =', self%layer_dims
!~                 print *, 'OC ambient mol fractions before =', self%ambient_mixtures(2)%fractions
!~                 print *, 'OC ambient wt fractions before =', self%ambient_mixtures(2)%weight_fractions
            !Convert material fractions to atomic mole fractions
            dummy = mat2at_core(xi=self%ambient_mixtures(2)%fractions, n_mats=4)
!~                 print *, 'dummy =', dummy
            !Convert atomic mole fractions to atomic wt fractions
            dummy = at2wt_core(at=dummy, m=masses, n_mats=4)
!~                 print *, 'dummy =', dummy
            !Rescale atomic wt fractions of S, Si, and O. Note that
            !the Fe content dummy(1) does not have an effect since the
            !fractions are normalized to unity anyways. Wt fractions are scaled
            !as wt_new = wt_old * M_Core / M_IC
            do i = 1, self%layer_dims(2) - 1
               dummy(i + 1) = dummy(i + 1)*self%layer_masses(2)/ &
                              (self%layers(1)%indigenous_mass/m_earth)
            end do
!~                 print *, 'dummy =', dummy
            !Compute new atomic mole fractions
            dummy = xi_general(eta=dummy, m=masses, n_mats=4)
!~                 print *, 'dummy =', dummy
            dummy_1 = at2mat_max(at=dummy, n_mats=4)
!~                 print *, 'dummy 1 =', dummy_1
!                self%ambient_mixtures(2)%fractions = dummy_1
!~                 print *, 'OC ambient wt fractions before fct =', self%ambient_mixtures(2)%weight_fractions

!~                 print *, 'OC ambient mol fractions after =', self%ambient_mixtures(2)%fractions
!~                 print *, 'OC ambient wt fractions after =', self%ambient_mixtures(2)%weight_fractions
            !call compute_mixture(self = self%ambient_mixtures(2))
            self%rho0_layers(2) = self%ambient_mixtures(2)%dens
!~                 print *, 'OC fractions before =', self%fractions%axes(2)%real_array
            !self%fractions%axes(2)%real_array = self%ambient_mixtures(2)%fractions
!~                 print *, 'OC fractions after =', self%fractions%axes(2)%real_array

         end if

      end if

!Case 3: No inner core exists. Nothing to be done.

   END SUBROUTINE update_core_mass

!#######################################################################
   SUBROUTINE get_layer_constraint(self, lay, skip)

      type(planet), intent(inout) :: self
      integer, intent(in) :: lay
      integer, intent(in) :: skip
!~ print *, 'get layerconstraints: skip / lay =', skip, lay
!Indigenous mass
      if (self%layer_constraints(lay + skip) == 0) then
         self%constraint_value_is = &
            self%layers(lay)%indigenous_mass

         !if the indigenous mass is probed in the next layer, it is currently
         !zero if it is checked from the bisection of the previous layer.
         !this has to be set manually as the next layer does not exist yet.
         if (skip .gt. 0) then

            self%constraint_value_is = 0.0d0
         end if

         self%constraint_value_should = &
            self%layer_masses(lay + skip)*m_earth

         self%direction = 1

!Enclosed mass
      else if (self%layer_constraints(lay + skip) == 1) then
         self%constraint_value_is = &
            self%layers(lay)%mass

!  self%constraint_value_should = &
!  (self%layer_masses(1)+self%layer_masses(2)+ &
!  self%layer_masses(4))*m_earth
         self%constraint_value_should = self%layer_masses(lay + skip)*m_earth

         self%direction = 1

!Radius
      else if (self%layer_constraints(lay + skip) == 2) then
         self%constraint_value_is = &
            self%layers(lay)%radius

         self%constraint_value_should = &
            self%layer_radii(lay + skip)*r_earth

         self%direction = 1

!Pressure
      else if (self%layer_constraints(lay + skip) == 3) then
         self%constraint_value_is = &
            self%layers(lay)%pres

         self%constraint_value_should = &
            self%layer_pres(lay + skip)

         self%direction = -1

!Temperature
      else if (self%layer_constraints(lay + skip) == 4) then
         self%constraint_value_is = &
            self%layers(lay)%temp

         self%constraint_value_should = &
            self%layer_temps(lay + skip)

         !For the inner core the transition is reached when the temperature
         !exceeds the melting temperature so the bisection direction is inverted
         if (self%lay == 1) then
            self%direction = 1
         else
            self%direction = -1
         end if

      end if

   END SUBROUTINE get_layer_constraint

!#######################################################################
   SUBROUTINE minor_bisection(self)

      type(planet), intent(inout) :: self
      integer :: sh, skip, i, skip_count, old_layer_constraint
      real(8) :: reldev, radius, mass, pres, temp, deltaT
      real(8) :: a, b, c, r

      call get_layer_constraint(self=self, lay=self%lay, skip=0)
! print *, 'Layer constraint should =', self%constraint_value_should
! print *, 'Layer constraint is =', self%constraint_value_is
! print *, 'layermasses =', self%layer_masses
      reldev = (self%constraint_value_should - self%constraint_value_is)/ &
               self%constraint_value_should
! print *, 'reldev, dir, eps =', reldev, self%direction, eps_layer
!Overshoot
      if (reldev*self%direction .lt. -eps_layer) then

   ! print *, 'Overshoot in minor bisection.'
   ! print *, 'Layer constraint type =', self%layer_constraints(self%lay)
   ! print *, 'Layer constraint should =', self%constraint_value_should
   ! print *, 'Layer constraint is =', self%constraint_value_is
         call remove_stuff(self=self)
         self%layers(self%lay)%overshoot = .true.
         sh = self%layers(self%lay)%shell_count

         !Two shells need to be reset here
         call reset_shell(self=self%layers(self%lay)%shells(sh))

         sh = sh - 1
!~   print *, 'Reconstruct shell =', sh
         self%layers(self%lay)%shell_count = sh

         call reset_shell(self=self%layers(self%lay)%shells(sh))
         call update_layer(self=self%layers(self%lay))

         !Compute the melting temperature of iron at current pressure
         if (self%lay == 1) then
            self%T_ICB = T_melt_iron(P=self%layers(self%lay)%pres, &
                                     X_Si=self%X_all_core(4), &
                                     X_O=self%X_all_core(5), &
                                     X_S=self%X_all_core(3))
            self%layer_temps(1) = self%T_ICB
!~                 print *, 'Update in minor bisection'
!~                 print *, 'Tm_Fe / P =', self%T_ICB, self%layers(self%lay)%pres*1d-9
!~                 print *, 'T =', self%layers(self%lay)%temp
         end if

         self%layers(self%lay)%dr = self%layers(self%lay)%dr*0.5d0
         self%layers(self%lay)%bisec = .true.
         self%layers(self%lay)%change_bisec = .false.
         self%subphase_refine_inhibit = .true.

         !If integration step size becomes too small, abort integration
         if (self%layers(self%lay)%dr .lt. 1.0d-10) then
            self%layer_iteration = .false.
            self%shell_iteration = .false.
         end if

!Layer transition reached
      else if (abs(reldev) .lt. eps_layer) then
   ! print *, 'Layer transition reached:', self%lay,'->', self%lay + 1
         if (self%lay == 1) then
!~   print *, 'layermasses(1) before update =', self%layer_masses(1)
            !call update_core_mass(self=self)
            !print *, 'layermasses after update core mass =', self%layer_masses
         end if
!~   ! updat coefficients for T-profile in mantle using core radius
!~   if(self%lay==1 .or. self%lay==2)then
!~         a = self%external_temp_profile(4,1)
!~         b = self%external_temp_profile(4,2)
!~         c = self%external_temp_profile(4,3)
!~         r = self%layers(self%lay)%radius * 1e-3
!~         print *, "update temp coeff =", self%external_temp_profile(4,:)
!~         self%external_temp_profile(4,1) = a - b * r + c * r**2
!~         print *, "update temp coeff =", self%external_temp_profile(4,:)
!~         self%external_temp_profile(4,2) = b - 2 * c * r
!~         print *, "update temp coeff =", self%external_temp_profile(4,:)
!~   end if
!~   print *, 'layermasses(1) after update =', self%layer_masses(1)
         self%layer_complete = .true.
         self%shell_iteration = .false.
         self%change_bisec = .true.

         self%subphase_refine_inhibit = .false.
         skip = 0

         !The construct method has already updated the shell count for the next
         !iteration. However, this iteration will never take place as the
         !shell that has just been constructed reached the minor constraint.
         !The shell count therefore has to be set back to this shell
         self%layers(self%lay)%shell_count = &
            self%layers(self%lay)%shell_count - 1

         !If this was the last layer, abort integration
         if (self%lay == size(self%contents%axes)) then
            self%layer_iteration = .false.
            ! print *, 'This was the last layer'

         else
            !Check if the condition for the next layer(s) is also met. In this case
            !skip the next layer and continue the integration with the
            !subsequent layer
            skip = 0
            skip_count = 0
            !Check layer transitions for all remaining layers
            do i = 1, size(self%contents%axes) - self%lay
!~     print *, 'i =', i
               call get_layer_constraint(self=self, lay=self%lay, skip=i)
               !   print *, 'Next layer constraint should =', self%constraint_value_should
               !   print *, 'Next layer constraint is =', self%constraint_value_is
               reldev = (self%constraint_value_should - self%constraint_value_is)/ &
                        self%constraint_value_should

               if (abs(reldev) .lt. eps_layer .or. reldev*self%direction .lt. -eps_layer) then
                  skip = skip + 1
               !   print *, 'Next layer transition reached aswell:', self%lay + i, '->', self%lay + i + 1
                  if (self%lay + skip == size(self%contents%axes)) then
                     self%layer_iteration = .false.
                     skip = skip - 1
                  end if

                  !Check again for layer 2 --> 3 because layer 3 has two layer
                  !constraints that could lead to 3 --> 4: pressure and mass
                  !but only relevant if layer 3 is not already skipped.
               else
                  if (self%lay <= 2) then
                     old_layer_constraint = self%layer_constraints(3)
                     self%layer_constraints(3) = 3
                     call get_layer_constraint(self=self, lay=self%lay, skip=i)
                     self%layer_constraints(3) = old_layer_constraint
                                 ! print *, 'Next layer constraint should =', self%constraint_value_should
                                 ! print *, 'Next layer constraint is =', self%constraint_value_is
                     reldev = (self%constraint_value_should - self%constraint_value_is)/ &
                              self%constraint_value_should

                     if (abs(reldev) .lt. eps_layer .or. reldev*self%direction .lt. -eps_layer) then
                        skip = skip + 1
                                 !   print *, 'Next layer transition reached aswell:', self%lay + i, '->', self%lay + i + 1
                        if (self%lay + skip == size(self%contents%axes)) then
                           self%layer_iteration = .false.
                           skip = skip - 1
                        end if
                     end if
                  end if
               end if
               !If one of the next layertransitions is not yet reached,
               !don't check for the subsequent ones as the iteration has
               !to continue in the next layer anyways. Also if it checks the
               !subsequent layers it's possible that a layer is not reached
               !but a later one is which leads to errors.
               if (skip < i) then
!~           print *, 'break'
                  exit
               end if
            end do

            do i = 1, skip + 1
!~           print *, 'Layer reached at M =', self%layers(self%lay)%mass/m_earth
               radius = self%layers(self%lay)%radius
               mass = self%layers(self%lay)%mass
               pres = self%layers(self%lay)%pres
               temp = self%layers(self%lay)%temp

               !If ICB is reached, update T_ICB
               if (self%lay == 1) then
                  self%T_ICB = temp
               end if

               !Compute temperature after temperature jump
               deltaT = self%temp_jumps(self%lay)
               temp = max(temp - deltaT, T_zero)
!~                 Update Si# and Fe# in the mantle from metal-silicate partitioning assuming
!~                 chemical equilibrium between core and mantle. If CS-model
!~                 is deactivated the mantle composition is computed from the inputs
!~                 Si_number_layers and Fe_number_layers
               if (self%core_segregation_model) then
                  if (self%lay == 9) then
                     print *, ''
                     print *, 'Doing the CS shit.... !!!!!!!'
                     print *, ''
                     print *, 'Si#, Fe# before:', self%Si_number_layers(3), self%Fe_number_layers(3)
                     !~                         print *, 'fractions =', self%fractions%axes(self%lay)%real_array
                     !~                         print *, 'fractions =', self%fractions%axes(3)%real_array
                     !~                         print *, 'xi_all core =', self%xi_all_core
                     self%Si_number_layers(3) = Si_number_mantle(self%T_CS, &
                                                                 self%P_CS, &
                                                                 self%xi_all_core)
                     self%Si_number_layers(4) = self%Si_number_layers(3)
                     print *, 'computed Si# in the mantle:', self%Si_number_layers(self%lay + 1)

                     self%Fe_number_layers(3) = Fe_number_mantle(self%T_CS, &
                                                                 self%P_CS, &
                                                                 self%xi_all_core)
                     !self%Fe_number_layers(3) = 1d-1
                     self%Fe_number_layers(4) = self%Fe_number_layers(3)
                     print *, 'computed Fe# in the mantle:', self%Fe_number_layers(self%lay + 1)
                  end if
               end if
               !Initiate next layer lay+1
               call init_layer(self=self%layers(self%lay + 1), &
                               n_mats=size(self%contents%axes(self%lay + 1)%int_array), &
                               contents=self%contents%axes(self%lay + 1)%int_array, &
                               fractions=self%fractions%axes(self%lay + 1)%real_array, &
                               r_in=radius, &
                               m=mass, &
                               T=temp, &
                               P=pres, &
                               tempType=self%tempType, &
                               rhoType=self%rhoType, &
                               adiabatType=self%adiabatType, &
                               q=self%q_layers(self%lay + 1), &
                               gammaG0=self%gammaG0_layers(self%lay + 1), &
                               eps_r=self%eps_r, &
                               rho0=self%rho0_layers(self%lay + 1), &
                               eps_T_zero=self%eps_T_zero, &
                               n_shells=n_shells_max_layer, &
                               eps_H2O=self%eps_H2O, &
                               eps_Al=self%eps_Al, &
                               lay=self%lay + 1, &
                               Si_number=self%Si_number_layers(self%lay + 1), &
                               Fe_number=self%Fe_number_layers(self%lay + 1), &
                               MOI=self%layers(self%lay)%MOI, & !MoI is currently the MoI from the last layer
                               omega=self%omega, &
                               xi_H=self%xi_H_layers(self%lay + 1), &
                               xi_Stv=self%xi_Stv_layers(self%lay + 1), &
                               X_impurity_0=self%X_impurity_0_layers(self%lay + 1), &
                               X_impurity_slope=self%X_impurity_slope_layers(self%lay + 1), &
                               external_temp_profile=self%external_temp_profile(self%lay + 1, :))

               self%n_shells_layers(self%lay) = self%layers(self%lay)%shell_count

               self%lay = self%lay + 1
            end do
         end if

!Nothing special happened
      else
         if (self%layers(self%lay)%change_bisec) then
            self%layers(self%lay)%bisec = .false.
         end if
      end if

   END SUBROUTINE minor_bisection

!########################################################################
   SUBROUTINE major_bisection(self)

      type(planet), intent(inout) :: self
      real(kind=8) :: param_is, param_should, eps
      integer, dimension(2) :: direction
      integer :: sh

      if (self%major_constraint == 'P_surface') then
         param_is = self%layers(self%lay)%pres
         param_should = self%P_surface_should
         !print *, "is, should = ", param_is, param_should
         eps = eps_P_surface
         direction = (/-1, -1/)

      elseif (self%major_constraint == 'T_surface') then
         param_is = self%layers(self%lay)%temp
         param_should = self%T_surface_should
         eps = eps_T_surface
         direction = (/-1, -1/)

      elseif (self%major_constraint == 'M_surface') then
         param_is = self%layers(self%lay)%mass
         param_should = self%M_surface_should
         eps = eps_M_surface
         direction = (/1, -1/)

      end if

!Overshoot
      if (direction(1)*(param_should - param_is)/param_should < direction(2)*eps &
          .or. param_is .lt. 0.0d0) then
! print *, ''
! print *, 'Overshoot in major bisection in layer', self%lay
! print *, 'Major constraint should =', param_should
! print *, 'Major constraint is =', param_is
         call remove_stuff(self=self)

         sh = self%layers(self%lay)%shell_count

         call reset_shell(self=self%layers(self%lay)%shells(sh))

         sh = sh - 1
         self%layers(self%lay)%shell_count = sh

         call reset_shell(self=self%layers(self%lay)%shells(sh))

         call update_layer(self=self%layers(self%lay))

         self%layers(self%lay)%overshoot = .true.
         self%layers(self%lay)%dr = self%layers(self%lay)%dr*0.5d0
         self%layers(self%lay)%bisec = .true.
         self%layers(self%lay)%change_bisec = .false.

!Surface reached
      elseif (abs(param_should - param_is)/param_should .lt. eps) then
         ! print *, 'Surface reached!'

         !The construct method has already updated the shell count for the next
         !iteration. However, this iteration will never take place as the
         !shell that has just been constructed reached the major constraint.
         !The shell count therefore has to be set back to this shell

         self%layers(self%lay)%shell_count = &
            self%layers(self%lay)%shell_count - 1
         self%layers(self%lay)%overshoot = .false.

         call update_layer(self=self%layers(self%lay))

         self%shell_iteration = .false.
         self%layer_iteration = .false.
         self%layers(self%lay)%temp = self%layers(self%lay)%temp - &
                                      self%temp_jumps(self%lay)

         self%n_shells_layers(self%lay) = self%layers(self%lay)%shell_count

      else
         if (self%layers(self%lay)%change_bisec) then
            self%layers(self%lay)%bisec = .false.
         end if
      end if

   END SUBROUTINE major_bisection

!########################################################################
   SUBROUTINE subphase_refine(self)

      type(planet), intent(inout) :: self
      type(shell) :: last_shell, current_shell
      logical :: trans, ug = .true., cm = .true.
      real(8) :: new_dr
      integer :: i, j, sh, layer_before, ph1, ph2

      trans = .false.
      sh = self%layers(self%lay)%shell_count

!Note that last and current shell here refer to the already constructed
!shells. As after construction the new shell is already initiated, the
!last two constructed shells are sh-1 and sh-2
      last_shell = self%layers(self%lay)%shells(sh - 2)
      current_shell = self%layers(self%lay)%shells(sh - 1)
!print *, 'trans before =', trans
!Check if between last shell and current shell a subphase transition
!has occured
!print *, 'contents =', last_shell%contents
      do i = 1, size(last_shell%contents)
         ph1 = last_shell%mixture%units(i)%phase
         ph2 = current_shell%mixture%units(i)%phase
         !print *,'ph1, ph2 =', ph1, ph2
         if (ph1 .ne. ph2) then
            trans = .true.
         end if
      end do
!print *, 'trans after =', trans
!If subphase transition is found, perform the shell refinement
      if (trans) then
         if (.not. self%subphase_refine_inhibit) then

            self%layers(self%lay)%overshoot = .true.
!~                 print *, 'Perform subphase refinement between shells', sh, sh-1

            call remove_stuff(self=self)
            !self%layers(self%lay)%overshoot = .true.

            call reset_shell(self=self%layers(self%lay)%shells(sh))
            !call reset_shell(self=self%layers(self%lay)%shells(sh-1))
!~                 print *, 'old sh =', sh
            sh = sh - 1
            self%layers(self%lay)%shell_count = sh
!~                 print *, 'new sh =', sh

            call reset_shell(self=self%layers(self%lay)%shells(sh))

            call update_layer(self=self%layers(self%lay))

            new_dr = self%layers(self%lay)%dr/real(self%subphase_res)

            do i = 1, self%subphase_res
!~                 print *, 'i =', i
               self%shell_iteration_count = self%shell_iteration_count + 1

               !Force dr to be the refined value
               self%layers(self%lay)%dr = new_dr
               sh = self%layers(self%lay)%shell_count
!~                         print *, 'Constructing shell:', sh
               call construct_layer(self=self%layers(self%lay))

               call update_shell(self=self%layers(self%lay)%shells(sh), update_grads=ug, &
                                 compute_mix=cm)
               !call update_layer(self=self%layers(self%lay))
               !print *, 'dens =', self%layers(self%lay)%shells(sh)%mixture%dens
               call merge_shells(self=self%layers(self%lay)%shells(sh + 1), &
                                 other=self%layers(self%lay)%shells(sh))

!                        print *, ''
!                        print *, 'Shell:', self%layers(self%lay)%shell_count

!                        print *, 'r=', self%layers(self%lay)%radius, &
!                        'm=', self%layers(self%lay)%mass, &
!                        'rho=', self%layers(self%lay)%dens, &
!                        'T=', self%layers(self%lay)%temp, &
!                        'P=', self%layers(self%lay)%pres, &
!                        'layer mass =', self%layers(self%lay)%indigenous_mass, &
!                        'dr =', new_dr
!                        print *, 'check rho =', self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)%dens

               call add_stuff(self=self)
               layer_before = self%lay

               !Only the major and minor bisection can set overshoot=.true.
               !If neither does, it must be false, so it is set to false. If the
               !integration that has been performed in this step did run into
               !an overshoot, it must be set to .true. for the next step, not the
               !current one.
               self%layers(self%lay)%overshoot = .false.
               call major_bisection(self=self)
               call minor_bisection(self=self)

               if (self%layers(self%lay)%overshoot) then
                  sh = self%layers(self%lay)%shell_count
                  self%layers(self%lay)%shell_count = sh
                  self%skip_bisec = .true.
                  exit
               end if

            end do
            !self%layers(self%lay)%bisec = .true.
            !self%layers(self%lay)%change_bisec = .false.
         end if
      end if

   END SUBROUTINE subphase_refine

!#######################################################################
   SUBROUTINE construct_planet(self, echo)

      type(planet), intent(inout) :: self
      logical, intent(inout), optional :: echo
      logical :: echo_dummy
      real(kind=8) :: olddens
      type(shell) :: lastshell, currentshell
      integer :: i, j, old_layer_constraint

      if (.not. present(echo)) then
         echo_dummy = .false.

      else
         echo_dummy = echo
      end if

      if (self%status == 'constructed') then
         print *, "This Planet has already been constructed. Use the &
         & reconstruct method to reconstruct the planet from seed with &
         &different specifications if you wish."

      else
         self%status = 'constructing'
         !self%lay = 1
         self%layer_iteration_count = 0
         self%shell_iteration_count = 0
         self%layer_iteration = .true.

         do while (self%layer_iteration)
   !   print *, '########################'
   !   print *, 'Layer:', self%lay
   !   print *, '########################'
   !   print *, 'n shell ocean =', self%layers(5)%shell_count
   !   print *, 'n shell outer core =', self%layers(2)%shell_count
   !   print *, 'actual layermasses:', self%layer_masses(1), self%layer_masses(2)
            
            self%layer_iteration_count = self%layer_iteration_count + 1
            self%shell_iteration = .true.
            self%change_bisec = .true.
            self%layer_complete = .false.
            self%subphase_refine_inhibit = .false.

            do while (self%shell_iteration)
               !   print *, 'Constructing shell:', self%layers(self%lay)%shell_count
               self%shell_iteration_count = self%shell_iteration_count + 1
               olddens = self%layers(self%lay)%dens

               self%skip_bisec = .false.
               call construct_layer(self=self%layers(self%lay))
      !  print *, ''
      !  print *, 'Shell:', self%layers(self%lay)%shell_count

      ! print *, 'r=', self%layers(self%lay)%radius, &
      !  'm=', self%layers(self%lay)%mass, &
      !  'rho=', self%layers(self%lay)%dens, &
      !  'T=', self%layers(self%lay)%temp, &
      !  'P (GPa)=', self%layers(self%lay)%pres*1d-9, &
      !  'layer mass =', self%layers(self%lay)%indigenous_mass, &
      !  'dr =', self%layers(self%lay)%dr
!      print *, 'check rho =', self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)%dens
               !print *, 'xi_H2O =', self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)%xi_H2O
               !print *, 'xi_H =', self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)%xi_H
!!!           print *, 'MOI factor =', self%layers(self%lay)%MOI/&
!!!           self%layers(self%lay)%mass/self%layers(self%lay)%radius**2
!      print *, 'grads =', self%layers(self%lay)%gradients
!      print *, 'gammaG0_layers =', self%gammaG0_layers
!~         print *, 'fracs =', self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)%fractions
!~         print *, 'conts =', self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)%contents
!~         print *, 'weight fracs =', self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)%weight_fractions

               call add_stuff(self=self)

               !Compute the melting temperature of iron at current pressure
               if (self%lay == 1) then
                  self%T_ICB = T_melt_iron(P=self%layers(self%lay)%pres, &
                                           X_Si=self%X_all_core(4), &
                                           X_O=self%X_all_core(5), &
                                           X_S=self%X_all_core(3))
                  self%layer_temps(1) = self%T_ICB
!~                         print *, 'Tm_Fe / P =', self%T_ICB, self%layers(self%lay)%pres*1d-9
!~                         print *, 'T =', self%layers(self%lay)%temp
               end if

               !Compute Olivine-Perovskite transition pressure
               self%P_MTZ = P_Ol_Pv(T=self%layers(self%lay)%temp)
               self%layer_pres(3) = self%P_MTZ

!~       !Only the major and minor bisection can set overshoot=.true.
               !If neither does, it must be false, so it is set to false. If the
               !integration that has been performed in this step did run into
               !an overshoot, it must be set to .true. for the next step, not the
               !current one.
               self%layers(self%lay)%overshoot = .false.

               !In the core the subphase refinement must be disabled in order
               !to prevent it to probe the Fe melting curve which is accounted
               !for as an actual layer transition.
               if (self%lay < 3) then
                  self%subphase_refine_inhibit = .true.
               end if

!~       if (self%layers(self%lay)%shell_count.gt.2)then
!~                 if(.not.self%subphase_refine_inhibit)then
!~                         call subphase_refine(self=self)
!~                 endif
!~       endif

               self%bisec = self%layers(self%lay)%bisec

               if (.not. self%skip_bisec) then
                  call major_bisection(self=self)
                  call minor_bisection(self=self)

                  !For the inner core check if Fe melting is reached
                  if (self%lay == 1) then
                     !If bisection has just been initiated for the the core mass,
                     !Do not check for melting. The bisection for the inner core
                     !has to be prioritiesed to the total core mass. If the total
                     !core mass is exceeded core integration needs to stop regardless
                     !of ICB is reached.
                     if (self%bisec .eqv. self%layers(self%lay)%bisec) then
                        ! print *, 'check melting'
                        old_layer_constraint = self%layer_constraints(1)
                        ! Temporarely set inner core constraint to temp and
                        ! check for layer transition
                        self%layer_constraints(1) = 4
                        call minor_bisection(self=self)

                        ! Reset inner core constraint to original value
                        self%layer_constraints(1) = old_layer_constraint
                        !self%layer_masses(1) = self%layer_masses(2)
                     end if
                  end if

                  !For the lower mantle check if pressure transition is reached
                  if (self%lay == 3) then
                     if (self%bisec .eqv. self%layers(self%lay)%bisec) then
                        old_layer_constraint = self%layer_constraints(1)
                        ! Temporarily set lower mantle constraint to pressure
                        ! and check for layer transition
                        self%layer_constraints(3) = 3
                        call minor_bisection(self=self)

                        ! Reset lower mantle constraint to to original value
                        self%layer_constraints(3) = old_layer_constraint
                        !Set LM mass to UM mass. This is important if an ocean
                        !is present. In this case the total mantle mass is used
                        !to define the mantle-ocean interface. It can happen either
                        !in the lower or upper mantle.
                        self%layer_masses(3) = self%layer_masses(4)
                     end if
                  end if

               end if

               if (self%layers(self%lay)%shell_count .gt. 2) then
                  if (.not. self%subphase_refine_inhibit) then
                     call subphase_refine(self=self)
                  end if
               end if

               if (self%layers(self%lay)%shell_count .gt. 1) then
                  lastshell = self%layers(self%lay)%shells(self%layers(self%lay)%shell_count - 1)

                  currentshell = self%layers(self%lay)%shells(self%layers(self%lay)%shell_count)
               end if

               if (self%shell_iteration_count >= n_shells_max_planet) then
                  self%layer_iteration = .false.
                  self%shell_iteration = .false.
                  print *, 'Total shell iteration limit exceeded.'
                  print *, 'r=', self%layers(self%lay)%radius, &
                     'm=', self%layers(self%lay)%mass, &
                     'rho=', self%layers(self%lay)%dens, &
                     'T=', self%layers(self%lay)%temp, &
                     'P=', self%layers(self%lay)%pres, &
                     'layer mass =', self%layers(self%lay)%indigenous_mass, &
                     'dr =', self%layers(self%lay)%dr

               end if

               if (self%layers(self%lay)%shell_count >= n_shells_max_layer) then
                  self%shell_iteration = .false.
                  self%layer_iteration = .false.
                  print *, 'Layer shell iteration limit exceeded in layer ', self%lay
                  print *, 'r=', self%layers(self%lay)%radius, &
                     'm=', self%layers(self%lay)%mass, &
                     'rho=', self%layers(self%lay)%dens, &
                     'T=', self%layers(self%lay)%temp, &
                     'P=', self%layers(self%lay)%pres, &
                     'layer mass =', self%layers(self%lay)%indigenous_mass, &
                     'dr =', self%layers(self%lay)%dr
               end if
            end do
         end do

         call update_layer(self=self%layers(self%lay))

         self%M_surface_is = self%layers(self%lay)%mass
         self%R_surface_is = self%layers(self%lay)%radius
         self%T_surface_is = self%layers(self%lay)%temp
         self%P_surface_is = self%layers(self%lay)%pres
         self%MOI_is = self%layers(self%lay)%MOI

         self%rho_mean = self%M_surface_is/(4.0d0/3.0d0*PI*self%R_surface_is**3)
         self%gravity = self%M_surface_is/self%R_surface_is**2*G

         !If planet has five layers, by convention, the last one is the surface ocean
         !Update the total ocean mass fraction in log
         if (self%lay == 5) then
            self%ocean_frac_is = self%layers(self%lay)%indigenous_mass/self%M_surface_is
            self%ocean_frac_is = log10(self%ocean_frac_is)
         end if

      end if

      do i = 1, self%lay
         do j = 1, self%layers(i)%shell_count
            self%n_shells = self%n_shells + 1
         end do
      end do

      ! print *, 'Integration terminated in layer ', self%lay
      ! print *, 'at shell ', self%n_shells

!Check if individual layers exist
      if (size(self%layers(1)%shells) .eq. 0) then
         self%inner_core_exists = .false.
      end if

      if (size(self%layers(2)%shells) .eq. 0) then
         self%outer_core_exists = .false.
      end if

   END SUBROUTINE construct_planet

END MODULE class_planet
