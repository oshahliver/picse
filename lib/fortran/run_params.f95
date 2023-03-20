MODULE run_params

   implicit none

!Relative uncertainties allowed for probing various parameters
   real(8), parameter :: eps_layer = 1d-6
   real(8), parameter :: eps_P_surface = 1d-4
   real(8), parameter :: eps_M_surface = 1d-4
   real(8), parameter :: eps_T_surface = 1d-2

!Maximum number of shells per layer and per planet
   integer, parameter :: n_shells_max_layer = 150
   integer, parameter :: n_shells_max_planet = 5000

!Number of integration parameters
   integer, parameter :: n_params_integration = 8

!Zero temperature for layertransition
   real(8) :: T_zero = -1d4

!Number of eos tables
   integer, parameter :: n_tables = 10

!Zero isothermal bulk modulus.
!This is to avoid spurious behaviours especially at high temperatures
!where KT can become negative in some cases. The zero value for KT has
!to be chosen large enough to avoid dTdr to get too large. If it does,
!the zero temperature can be reached even in one integration step. This
!is supposed to be avoided by the adaptive integration step size, but
!if the gradient during one integration changes too much, this estimation
!is not sufficient in some cases to avoid T to drop drastically in one step
   real(8) :: KT_zero = 5d10

!Default number of output parameters for table interpolation
!Parameters are: rho, dTdP_S, dPdrho, alpha, xi_H2O, X_impurity
   integer :: n_out = 6

!max FeO content in the mantle
   real(8), parameter :: xi_Fe_mantle_max = 1d0

!Butcher-table for the 4th order RK solver
   real(8), dimension(4, 4), parameter :: a_list = &
                                          transpose(reshape((/0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, &
                                                              0.0, 0.0, 0.0, 0.0, 1.0, 0.0/), shape(a_list)))
   real(8), dimension(4), parameter :: b_list = (/1./6., 1./3., 1./3., 1./6./)
   real(8), dimension(4), parameter :: c_list = (/0., 1./2., 1./2., 1./)

!Specify for which tables variable iron content exists (1) and for
!which tables a hydrated state exists (1)
   integer, dimension(n_tables) :: eos_tables_variable_Fe = (/0, 1, 0, 1, 1, 1, 1, 0, 0, 0/)
   integer, dimension(n_tables) :: eos_tables_hydrated = (/0, 0, 1, 1, 1, 1, 1, 0, 0, 0/)

!List containing the table numbers that has to be taken for interpolation
!in the hydrated case. E.g. if ll=1, table 1 is taken anyways. If ll=4
!then table 11 is taken which contains the hydrated (Mg,Fe)O data.
   integer, dimension(n_tables) :: which_eos_tables_hydrated = (/1, 2, 11, 12, 13, 14, 15, 8, 9, 10/)

!if variable Fe content, each component is translated to use the right
!eos table
!integer, dimension(8) :: eos_translation = (/1,2,3,9,10,11,12,8/)

END MODULE run_params
