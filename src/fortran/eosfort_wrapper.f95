
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
public :: load_eos_tables, construct_new_planet, pure_curve


contains


!#######################################################################
!Load all EoS table files for subsequent table interpolations.
SUBROUTINE load_eos_tables()

implicit none

character(*), parameter :: table_dir = "./data/EoS_tables/"
integer :: n_tables=10, n_tables_hyd=4, o=1

print *, "table directory:", table_dir

call initiate(n_tables=n_tables, n_tables_hyd=n_tables_hyd, table_dir=table_dir)


END SUBROUTINE load_eos_tables

!#######################################################################
SUBROUTINE construct_new_planet(T_center, P_center, fractions, &
contents, P_surface_should, T_surface_should, r_seed, layer_dims, &
tempType, rhoType, adiabatType, layerType, eps_r, layer_masses, &
layer_pres, layer_radii, temp_jumps, q_layers, gammaG0_layers, &
layer_constraints, Si_number_layers, Fe_number_layers, eps_Al, &
eps_H2O, M_surface_is, R_surface_is, P_surface_is, T_surface_is, &
Mg_number_is, Si_number_is, Fe_count, Si_count, Al_count, Mg_count, &
O_count, H2O_count, H_count, S_count, ocean_frac_is, MOI, &
layer_properties, omega, xi_H_core, P_H2, xi_H_core_predicted, profiles,&
 n_shells, subphase_res, xi_Stv, n_shells_layers, out_frac, &
 X_impurity_0_layers, X_impurity_slope_layers, xi_all_core,&
 X_all_core, T_CS, P_CS, xi_Fe_mantle, Si_number_mantle, mantle_exists, &
 inner_core_exists, outer_core_exists, core_segregation_model,&
 M_surface_should, M_ocean_should, Mg_number_should, &
 inner_core_segregation_model, E_grav, E_int, external_temp_profile)

implicit none

real(8), intent(in) :: omega, xi_Stv
real(8), intent(inout) :: P_surface_should, T_surface_should, &
r_seed, eps_r, eps_Al, eps_H2O
real(8), dimension(2) :: Si_number_extremes
integer, intent(in) :: tempType, rhoType, adiabatType, &
layerType, layer_constraints(:), contents(:), subphase_res
real(8), intent(in) :: P_center, T_center, xi_H_core_predicted, P_CS
real(8), intent(in) :: layer_masses(:), layer_pres(:), &
layer_radii(:), temp_jumps(:), q_layers(:), gammaG0_layers(:), &
Si_number_layers(:), Fe_number_layers(:), X_impurity_0_layers(:),&
X_impurity_slope_layers(:)
real(8), intent(in) :: fractions(:)
integer, intent(in) :: layer_dims(:)
integer, intent(out) :: n_shells
type(weird_array) :: conts, add_conts, add_fracs, fracs
type(planet) :: pl
real(8), intent(in) :: M_surface_should, M_ocean_should, &
Mg_number_should
real(8), dimension(5), intent(in) :: xi_all_core, X_all_core
real(8), intent(out), dimension(10, 10) :: out_frac
integer :: n_layers, n_additionals, i, j, c
real(8), intent(out) :: M_surface_is, R_surface_is, T_surface_is, &
P_surface_is, Si_number_is, Mg_number_is, Fe_count, Si_count, Al_count,&
Mg_count, O_count, H2O_count, ocean_frac_is, MOI, xi_H_core, P_H2, &
H_count, S_count, T_CS, xi_Fe_mantle, Si_number_mantle, E_grav, E_int
real(8), intent(out), dimension(size(layer_dims), 6) :: layer_properties
real(8), dimension(8,1000), intent(out) :: profiles
integer, dimension(size(layer_dims)), intent(out) :: n_shells_layers
logical, intent(out) :: mantle_exists, inner_core_exists, outer_core_exists
logical, intent(in) :: core_segregation_model, inner_core_segregation_model
real(8), intent(in), dimension(5,3) :: external_temp_profile(:,:)

n_layers = size(layer_dims, 1)


!~ if(.not.present(external_temp_profile))then
!~ 	external_temp_profile = &
!~ 	reshape((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 , 0d0, 0d0, &
!~ 	0d0 , 0d0, 0d0, 0d0/), shape(external_temp_profile))
!~ endif

!n_additionals = size(additional_dims, 1)

call init_weird_array(self=conts, ndim=n_layers, &
axes_lengths=layer_dims)

call init_weird_array(self=fracs, ndim=n_layers, &
axes_lengths=layer_dims)

!call init_weird_array(self=add_conts, ndim=n_additionals, &
!axes_lengths=additional_dims)

!call init_weird_array(self=add_fracs, ndim=n_additionals, &
!axes_lengths=additional_dims)
c = 1

!The iron number in the core is defined via the H content
!In order to avoid numerical issues in the table interpolation
!the Fe number must be strictly smaller than 1. For pure iron
!a value very close to one is set here.
!Fe_number_layers(1) = min(1d0-xi_H_core_predicted, 1d0-1d-10)
!Fe_number_layers(2) = min(1d0-xi_H_core_predicted, 1d0-1d-10)
do i=1, n_layers
	do j=1, layer_dims(i)
		fracs%axes(i)%real_array(j) = fractions(c)
		conts%axes(i)%int_array(j) = contents(c)
		c = c + 1
	enddo
	!print *, 'conts(i) =', conts%axes(i)%int_array(:)
	!print *, 'fracs(i) =', fracs%axes(i)%real_array(:)
enddo

!~ print *, 'fractions in eosfort_wrapper:'
!~ do i=1, n_layers
!~ 	print *, 'layer ', i, ':', fracs%axes(i)%real_array
!~ enddo

do i=1, n_layers
  Si_number_extremes =  Si_number_max(Mg_number = 1.0d0 - Fe_number_layers(i), &
  contents=conts%axes(i)%int_array)
enddo

!c = 1
!do i=1, n_additionals
!	do j=1, additional_dims(i)
!		add_fracs%axes(i)%real_array(j) = additional_fractions(c)
!		add_conts%axes(i)%int_array(j) = additional_contents(c)
!		c = c + 1
!	enddo
!enddo

!Set out_frac to initial fractions because if integration is terminated
!before all layers were reached it would be set to zero in the outer core
!and in the next integration the desired material fractions would not
!be correctly adopted. Currently there are 5 possible layers
do i=1, n_layers
	do j=1, layer_dims(i)
		out_frac(i,j) = fracs%axes(i)%real_array(j)
	enddo
enddo

call init_planet(self=pl, T_center=T_center, P_center=P_center, &
R_seed=r_seed, contents=conts, fractions=fracs, tempType=tempType, rhoType=rhoType, &
adiabatType=adiabatType, layerType=layerType, eps_r=eps_r, layer_masses=layer_masses, &
layer_radii = layer_radii, layer_pres = layer_pres, temp_jumps=temp_jumps, &
q_layers=q_layers, gammaG0_layers=gammaG0_layers, n_layers=n_layers, &
eps_T_zero=0d0, layer_constraints=layer_constraints, Si_number_layers=Si_number_layers,&
Fe_number_layers=Fe_number_layers, eps_Al=eps_Al, eps_H2O=eps_H2O, &
T_surface_should=T_surface_should, P_surface_should=P_surface_should, &
omega=omega, xi_H_core_predicted=xi_H_core_predicted, subphase_res=subphase_res,&
xi_Stv=xi_Stv, X_impurity_0_layers=X_impurity_0_layers, &
X_impurity_slope_layers=X_impurity_slope_layers, xi_all_core=xi_all_core,&
X_all_core=X_all_core, P_CS=P_CS, core_segregation_model = core_segregation_model, &
M_surface_should = M_surface_should, M_ocean_should = M_ocean_should, &
Mg_number_should = Mg_number_should, &
inner_core_segregation_model = inner_core_segregation_model,&
external_temp_profile=external_temp_profile)

print *, "check2"
!~ print *, 'fractions before construct in eosfort_wrapper:', pl%fractions%axes(2)%real_array
call construct_planet(self=pl)
!~ print *, 'fractions after construct in eosfort_wrapper:', pl%fractions%axes(2)%real_array

!~ call compute_xi_H_core(self=pl)
call get_profiles(self=pl)
call compute_E_grav(self=pl)
call compute_E_int(self=pl)

do i=1, 8
	do j=1, pl%n_shells + pl%lay - 1
		profiles(i,j) = pl%profiles(i,j)
	enddo
enddo
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
Al_count = pl%N_Al
O_count = pl%N_O
S_count = pl%N_S
ocean_frac_is = pl%ocean_frac_is
MOI = pl%MOI_is
P_H2 = pl%P_H2
xi_H_core = pl%xi_H_core
n_shells = pl%n_shells
n_shells_layers = pl%n_shells_layers
T_CS = pl%T_CS
Si_number_mantle = pl%Si_number_layers(3)
xi_Fe_mantle = pl%Fe_number_layers(3)
mantle_exists = pl%mantle_exists
inner_core_exists = pl%inner_core_exists
outer_core_exists = pl%outer_core_exists
E_grav = pl%E_grav
E_int = pl%E_int


!~ print *, 'eosfort_wrapper: integrated up to layer ', pl%lay
!~ print *, 'out_frac before update in eosfort_wrapper =', out_frac(2,:)
do i=1, pl%lay
	if (size(pl%layers(i)%fractions) > 0)then
		do j=1, layer_dims(i)
			out_frac(i,j) = pl%fractions%axes(i)%real_array(j)
		enddo
	endif
enddo

!~ print*, 'out fracs in eosfort_wrapper =', out_frac(3,:)
 
do i=1, n_layers
	layer_properties(i,1) = pl%layers(i)%pres
	layer_properties(i,2) = pl%layers(i)%temp
	layer_properties(i,3) = pl%layers(i)%dens
	!If layer exists add bottom density else bottom and top density are
	!just the same. Distinction needs to be made here because if layer
	!does not exist a segmentation fault occurs if the bottom density
	!is attempted to be extracted.
	if(pl%layers(i)%shell_count .gt. 4 .and. pl%layers(i)%shell_count .lt. 500)then
		layer_properties(i,4) = pl%layers(i)%shells(1)%dens
	
	else
		layer_properties(i,4) = pl%layers(i)%dens
	endif
	layer_properties(i,5) = pl%layers(i)%indigenous_mass
	layer_properties(i,6) = pl%layers(i)%radius
enddo

print *, 'NOTE: out_frac(5) is set to 1 manually in eosfort_wrapper!'
out_frac(5,1) = 1d0
END SUBROUTINE construct_new_planet

!#######################################################################
SUBROUTINE pure_curve(N, ll, T_surface, P_surface, P_min, P_max)

integer, intent(in) :: ll, N
real(8), intent(in) :: T_surface, P_surface
real(8), intent(in) :: P_min, P_max
real(8), dimension(N) :: P_center
type(planet), dimension(N) :: planets
integer, dimension(1,1) :: contents
real(8), dimension(1,1) :: fractions
real(8) :: delta
integer :: i

contents(1,1) = ll
fractions (1,1) = 1d0

!Create log pressure difference between two adjacent points
delta = (log10(P_max)-log10(P_min))/(N-1)

!Compute central pressure for each point
do i=1, N
	P_center(i) = 10**(log10(P_min) + (i-1)*delta)
enddo

END SUBROUTINE pure_curve

END MODULE wrapper
