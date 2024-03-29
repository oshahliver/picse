
MODULE wrapper

   use functions
   use eosfort
   use run_params
   use class_unit
   use class_mixture
   use class_shell
   use class_layer
   use class_planet
   use class_weird_array
   use constants
   use phase

   implicit none

   private
   public :: load_eos_tables, create_planet

contains

!#######################################################################
!Load all EoS table files for subsequent table interpolations.
   SUBROUTINE load_eos_tables(table_dir)

      implicit none

      character(*), intent(in) :: table_dir
      integer :: n_tables = 10
      integer:: n_tables_hyd = 4

      ! if (verbose) then
      !    print *, "Loading EoS tables..."
      !    print *, "table directory:", table_dir
      ! end if

      call initiate(n_tables=n_tables, n_tables_hyd=n_tables_hyd, table_dir=table_dir)

   END SUBROUTINE load_eos_tables

!#######################################################################

   SUBROUTINE create_planet(T_center, P_center, fractions, &
                            contents, P_surface_should, T_surface_should, r_seed, layer_dims, &
                            tempType, rhoType, adiabatType, layerType, eps_r, layer_masses, &
                            layer_pres, layer_radii, temp_jumps, q_layers, gammaG0_layers, &
                            layer_constraints, Si_number_layers, Fe_number_layers, &
                            xi_all_core, X_all_core, R_surface_should, M_core_should, &
                            M_surface_should, M_ocean_should, subphase_res, spin, &
                            Mg_number_should, inner_core_segregation_type, core_segregation_type, &
                            pres_core_segregation, &
                            M_surface_is, &
                            R_surface_is, &
                            P_surface_is, &
                            T_surface_is, &
                            Mg_number_is, &
                            Si_number_is, &
                            Fe_count, &
                            Si_count, &
                            Mg_count, &
                            O_count, &
                            H2O_count, &
                            H_count, &
                            S_count,&
                            ocean_frac_is, &
                            moment_of_inertia, &
                            layer_properties, &
                            profiles, &
                            n_shells, &
                            n_shells_layers, &
                            out_frac, &
                            xi_Fe_mantle, &
                            Si_number_mantle, &
                            mantle_exists, &
                            inner_core_exists, &
                            outer_core_exists, &
                            E_grav_is,&
                            E_int_is, &
                            E_tot_is, &
                            L_int_is)

      implicit none

      ! Declare custom type variables
      type(weird_array) :: conts
      type(weird_array) :: add_conts
      type(weird_array) :: add_fracs
      type(weird_array) :: fracs
      type(planet) :: pl

      ! Declare input variables
      real(8), intent(in) :: P_center
      real(8), intent(in) :: T_center
      real(8), intent(in) :: fractions(:)
      real(8), intent(in) :: P_surface_should
      real(8), intent(in) :: T_surface_should
      real(8), intent(in) :: M_surface_should
      real(8), intent(in) :: M_ocean_should
      real(8), intent(in) :: Mg_number_should
      real(8), intent(in) :: R_surface_should
      real(8), intent(in) :: M_core_should
      real(8), intent(in) :: r_seed
      real(8), intent(in) :: eps_r

      real(8), intent(in) :: layer_masses(:)
      real(8), intent(in) :: layer_pres(:)
      real(8), intent(in) :: layer_radii(:)
      real(8), intent(in) :: temp_jumps(:)
      real(8), intent(in) :: q_layers(:)
      real(8), intent(in) :: gammaG0_layers(:)

      real(8), intent(in) :: Si_number_layers(:)
      real(8), intent(in) :: Fe_number_layers(:)

      integer, intent(in) :: tempType
      integer, intent(in) :: rhoType
      integer, intent(in) :: adiabatType
      integer, intent(in) :: layerType

      integer, intent(in) :: layer_dims(:)
      integer, intent(in) :: contents(:)

      ! Variable array size for core composition, this is probably going to cause some troubles!!!!!
      real(8), intent(in), dimension(5) :: xi_all_core
      real(8), intent(in), dimension(5) :: X_all_core

      integer, intent(in) :: layer_constraints(:)

      ! Declare optional input variables
      real(8), intent(in) :: spin
      integer, intent(in) :: subphase_res
      real(8), intent(in) :: pres_core_segregation
      integer, intent(in) :: inner_core_segregation_type
      integer, intent(in) :: core_segregation_type

      ! Declare output variables
      integer, intent(out) :: n_shells

      real(8), intent(out) :: H_count
      real(8), intent(out) :: H2O_count
      real(8), intent(out) :: Fe_count
      real(8), intent(out) :: Si_count
      real(8), intent(out) :: Mg_count
      real(8), intent(out) :: S_count
      real(8), intent(out) :: O_count

      real(8), intent(out) :: ocean_frac_is
      real(8), intent(out) :: M_surface_is
      real(8), intent(out) :: P_surface_is
      real(8), intent(out) :: T_surface_is
      real(8), intent(out) :: R_surface_is
      real(8), intent(out) :: Mg_number_is
      real(8), intent(out) :: Si_number_is
      real(8), intent(out) :: moment_of_inertia
      real(8), intent(out) :: Si_number_mantle
      real(8), intent(out) :: xi_Fe_mantle
      real(8), intent(out) :: E_tot_is, E_grav_is, E_int_is, L_int_is

      real(8), intent(out), dimension(10, 10) :: out_frac
      real(8), intent(out), dimension(size(layer_dims), 6) :: layer_properties
      integer, dimension(size(layer_dims)), intent(out) :: n_shells_layers
      real(8), dimension(8, 1000), intent(out) :: profiles

      logical, intent(out) :: mantle_exists
      logical, intent(out) :: outer_core_exists
      logical, intent(out) :: inner_core_exists

      integer :: n_layers, shell_idx
      integer :: i, j, c, ii
      real(8), dimension(n_params_integration) :: params
      real(8) :: major_vals_should(5)
      
      ! print *, "creating new planet..."
      ! TODO: Define these default values somewhere else!

      ! If no rotation is passed, set it to zero
      ! if (.not. present(spin)) then
      !    spin = 0e0
      ! end if

      ! if (.not. present(subphase_res)) then
      !    spin = 2
      ! end if

      ! if (.not. present(core_segregation_type)) then
      !    core_segregation_type = 0
      ! end if

      ! if (.not. present(core_segregation_type)) then
      !    inner_core_segregation_type = 0
      ! end if

      ! TODO: The default CS pressure should be computed self-consistently using
      ! a default model (the one used in Shah et al. 2022?)
      ! if (.not. present(pres_core_segregation) .and. core_segregation_type .gt. 0) then
      !    pres_core_segregation = 3e10
      ! end if

      ! Total number of planetary layers
      n_layers = size(layer_dims, 1)

      ! Initiate custom type arrays for material contents and mole fractions
      call init_weird_array(self=conts, ndim=n_layers, &
                            axes_lengths=layer_dims)

      call init_weird_array(self=fracs, ndim=n_layers, &
                            axes_lengths=layer_dims)

      ! Allocate values to custom type arrays
      c = 1
      do i = 1, n_layers
         do j = 1, layer_dims(i)
            fracs%axes(i)%real_array(j) = fractions(c)
            conts%axes(i)%int_array(j) = contents(c)
            c = c + 1
         end do
         !print *, '   conts(i) =', conts%axes(i)%int_array(:)

         !print *, 'fracs(i) =', fracs%axes(i)%real_array(:)
      end do

!Set out_frac to initial fractions because if integration is terminated
!before all layers were reached it would be set to zero in the outer core
!and in the next integration the desired material fractions would not
!be correctly adopted. Currently there are 5 possible layers:
!Inncer Core, Outer Core, Lower Mantle, Upper Mantle, Hydrosphere
      do i = 1, n_layers
         do j = 1, layer_dims(i)
            out_frac(i, j) = fracs%axes(i)%real_array(j)
         end do
      end do
      major_vals_should = (/M_surface_should, R_surface_should, P_surface_should, T_surface_should, M_surface_should/)
! Initiate custom type planet instance
      call init_planet(self=pl, T_center=T_center, P_center=P_center, &
                       R_seed=r_seed, contents=conts, fractions=fracs, tempType=tempType, rhoType=rhoType, &
                       adiabatType=adiabatType, layerType=layerType, eps_r=eps_r, layer_masses=layer_masses, &
                       layer_radii=layer_radii, layer_pres=layer_pres, temp_jumps=temp_jumps, &
                       q_layers=q_layers, gammaG0_layers=gammaG0_layers, n_layers=n_layers, &
                       eps_T_zero=0d0, layer_constraints=layer_constraints, Si_number_layers=Si_number_layers, &
                       Fe_number_layers=Fe_number_layers, &
                       T_surface_should=T_surface_should, P_surface_should=P_surface_should, &
                       omega=spin, subphase_res=subphase_res, xi_all_core=xi_all_core, &
                       X_all_core=X_all_core, P_CS=pres_core_segregation, core_segregation_model=core_segregation_type, &
                       M_surface_should=M_surface_should, M_ocean_should=M_ocean_should, &
                       Mg_number_should=Mg_number_should, M_core_should=M_core_should, &
                       inner_core_segregation_model=inner_core_segregation_type, R_surface_should=R_surface_should,&
                       major_vals_should = major_vals_should)
      
      call construct_planet(self=pl)
      call get_profiles(self=pl)
      call compute_L_int(self=pl)

      do i = 1, 8
         do j = 1, pl%n_shells + pl%lay - 1
            profiles(i, j) = pl%profiles(i, j)
         end do
      end do
      M_surface_is = pl%M_surface_is
      R_surface_is = pl%R_surface_is
      P_surface_is = pl%P_surface_is
      T_surface_is = pl%T_surface_is
      Mg_number_is = pl%Mg_number_is
      Si_number_is = pl%Si_number_is
      Si_count = pl%N_Si
      Mg_count = pl%N_Mg
      Fe_count = pl%N_Fe
      H2O_count = pl%N_H2O
      H_count = pl%N_H
      O_count = pl%N_O
      S_count = pl%N_S
      ocean_frac_is = pl%ocean_frac_is
      moment_of_inertia = pl%MOI_is / (pl%M_surface_is * pl%R_surface_is**2)
      n_shells = pl%n_shells
      n_shells_layers = pl%n_shells_layers
      Si_number_mantle = pl%Si_number_layers(3)
      xi_Fe_mantle = pl%Fe_number_layers(3)
      E_int_is = pl%E_int_is
      E_grav_is = pl%E_grav_is
      E_tot_is = E_int_is + E_grav_is
      L_int_is = pl%L_int_is

! Update layer fractions in case composition changed
      do i = 1, pl%lay
         if (size(pl%layers(i)%fractions) > 0) then
            do j = 1, layer_dims(i)
               out_frac(i, j) = pl%layers(i)%fractions(j)
            end do
         end if
      end do

! Update layer properties
      do i = 1, n_layers
         shell_idx = pl%layers(i)%shell_count
         if (shell_idx > 0) then
            params = pl%layers(i)%shells(shell_idx)%integration_parameters
            layer_properties(i, 1) = params(1)
            layer_properties(i, 2) = params(3)
            layer_properties(i, 3) = params(4)
            layer_properties(i, 4) = params(4)
            layer_properties(i, 5) = pl%layers(i)%indigenous_mass
            layer_properties(i, 6) = pl%layers(i)%radius
            !If layer exists, add bottom density. Else bottom and top density are
            !just the same. Distinction needs to be made here because if layer
            !does not exist a segmentation fault occurs if the bottom density
            !is attempted to be extracted.
            ! if (pl%layers(i)%shell_count .gt. 1) then ! .and. pl%layers(i)%shell_count .lt. 500) then
            !    layer_properties(i, 4) = pl%layers(i)%shells(1)%integration_parameters(4)
            ! else
            !    layer_properties(i, 4) = params(4)!pl%layers(i)%dens
            ! end if

         endif
      end do

      ! print *, 'NOTE: out_frac(5) is set to 1 manually in eosfort_wrapper! (No impurities in hydrosphere)'
      ! out_frac(5, 1) = 1d0

   END SUBROUTINE create_planet


END MODULE wrapper
