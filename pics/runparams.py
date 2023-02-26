# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 15:07:45 2018

@author: os18o068
"""
from pics.physicalparams import m_earth, r_earth

dur = 0.0
time_round = 3
mass_round = 4
pres_round = 4
radius_round = 3
temp_round = 2
dens_round = 3
layerbisection_limit = 50
shell_iteration_limit = 1000  # max number of shell iterations before interupt
iteration_limit = 0  # max number of iteration steps per planet to match R and M
coreiteration_limit = 0
pressureiteration_limit = 0
bisectionstep_limit = 15
bisection_counter_list = []
bisectionstep_counter_list = []

supported_base_types = ["telluric", "aqua"]

initial_predictor_keys = [
    "M_surface_should",
    "Mg_number_should",
    "ocean_fraction_should",
    "Si_number_mantle",
    "T_surface_should",
    "Fe_number_mantle",
    "fractions",
]


fortplanet_keys_translator = [
    "t_center",
    "p_center",
    "contents",
    "fractions",
    "layer_dims",
    "layer_masses",
    "layer_pres",
    "layer_radii",
    "layer_constraints",
    "temp_jumps",
    "q_layers",
    "gammag0_layers",
    "si_number_layers",
    "fe_number_layers",
    "r_seed",
    "adiabattype",
    "rhotype",
    "layertype",
    "temptype",
    "eps_r",
    "p_surface_should",
    "t_surface_should",
    "spin",
    "subphase_res",
    "xi_all_core",
    "x_all_core",
    "pres_core_segregation",
    "core_segregation_type",
    "m_ocean_should",
    "m_surface_should",
    "mg_number_should",
    "inner_core_segregation_type",
    "r_surface_should",
    "m_core_should",
]

fortplanet_output_keys = [
    "M_surface_is",
    "R_surface_is",
    "P_surface_is",
    "T_surface_is",
    "Mg_number_is",
    "Si_number_is",
    "Fe_count",
    "Si_count",
    "Mg_count",
    "O_count",
    "H2O_count",
    "H_count",
    "S_count",
    "ocean_fraction_is",
    "moment_of_inertia_is",
    "layer_properties_dummy",
    "profiles",
    "shell_count",
    "shell_count_layers",
    "fractions",
    "x_Fe_mantle",
    "Si_number_mantle",
    "mantle_exists",
    "inner_core_exists",
    "outer_core_exists",
    "gravitational_energy",
    "internal_energy",
]


fortplanet_input_keys = [
    "T_center",
    "P_center",
    "contents",
    "fractions",
    "layer_dims",
    "layer_masses",
    "layer_pressures",
    "layer_radii",
    "layer_constraints",
    "temperature_jumps",
    "debye_exponents_layers",
    "grueneisen_gammas_layers",
    "Si_number_layers",
    "Fe_number_layers",
    "seed_radius",
    "adiabat_type",
    "density_type",
    "layer_type",
    "temperature_type",
    "eps_r",
    "P_surface_should",
    "T_surface_should",
    "spin",
    "subphase_resolution",
    "x_all_core",
    "eta_all_core",
    "core_segregation_pressure",
    "core_segregation_type",
    "M_ocean_should",
    "M_surface_should",
    "Mg_number_should",
    "inner_core_segregation_type",
    "R_surface_should",
    "M_core_should",
]

parameter_mask = [
    "Mg_number",
    "radius",
    "moment_of_inertia",
    "Si_number",
    "Si_number_mantle",
    "x_Fe_mantle",
    "pres_core_segregation",
    "temp_core_segregation",
    "w_S_outer_core",
    "w_Si_outer_core",
    "w_O_outer_core",
    "x_FeS_outer_core",
    "x_FeSi_outer_core",
    "x_FeO_outer_core",
    "T0_CMB",
    "radius_core",
    "radius_inner_core",
    "radius_HMI",
    "temperature_TBL",
    "x_SiO2_mantle",
    "x_FeO_mantle",
    "temp_center",
    "pres_center",
    "log_oxygen_fug",
    "delta_temp_MTZ",
    "mass_water",
    "mass_water_mantle",
    "Fe_number_outer_core",
    "core_mass_frac",
    "temp_CMB",
    "pres_CMB",
    "dens_CMB",
    "temp_ICB",
    "pres_ICB",
    "dens_ICB",
    "temp_MTZ",
    "pres_MTZ",
    "dens_MTZ",
    "x_MgO_mantle",
    "radius_MTZ",
    "x_S_outer_core",
    "x_S_mantle",
    "S_to_Fe",
    "Si_to_Fe",
    "Mg_to_Fe",
    "S_to_Si",
    "mass",
    "inner_core_mass_frac",
    "water_mass_frac",
    "pres_surface",
    "w_liquid",
    "w_supercrit",
    "w_solid",
    "h_liquid",
    "h_supercrit",
    "h_solid",
    "w_liquid_to_w_H2O",
    "w_supercrit_to_w_H2O",
    "w_solid_to_w_H2O",
    "vol_liquid",
    "vol_supercrit",
    "vol_solid",
    "gravity",
    "hydro_structure",
    "temp_HMI",
    "pres_HMI",
    "dens_HMI",
    "bulk_structure",
]

# water phase regions
water_phases = ["ice Ih", "high pressure Mazevet", "high pressure French", "IAPWS"]

# Coefficients for external temperature profiles from Cammarano et al. 2006
external_temp_profiles = [  # Cold
    [
        [1200, 0, 0],
        [1200, 0, 0],
        [1200, -0.5, -0.0004655],
        [1200, -0.5, -0.0004655],
        [273, 0, 0],
    ],
    # Hot1
    [
        [1476.0, -0.01106, -0.00006433],
        [1476.0, -0.01106, -0.00006433],
        [1437.0, 0.02399, -0.00003742],
        [1437.0, 0.02399, -0.00003742],
        [273.0, 0, 0],
    ],
    # Hot 2
    [
        [1888, -0.01253, -0.00008785],
        [1888, -0.01253, -0.00008785],
        [1453, -0.00301, -0.0000274],
        [1453, -0.00301, -0.0000274],
        [273.0, 0, 0],
    ],
]

# Different types of hydrosphere structures

hydro_structures = {
    "a": 21020,
    "b": 2120,
    "c": 1210,
    "d": 1201,
    "e": 210,
    "f": 1020,
    "g": 120,
    "h": 202,
    "i": 20,
    "j": 2,
    "k": 0,
}

param_labels = [
    r"$\T \ [K]$",
    r"$\P \ [Pa]$",
    r"$\rho \ [kg/m^3]$",
    r"$\ (dT/dP)_S \ [K \ Pa^{-1}]$",
    r"$\ (dP/d \rho )_T \ [m^2 / s^2]$",
    r"$\alpha_{th} \ [K^{-1}]$",
    r"$c_P \ [\rm J \ kg^{-1} K^{-1}]$",
    r"$\rho \ [kg/m^3]$",
    r"$\rho \ [kg/m^3]$",
]


plot_units = {
    "Pa": 1.0,
    "kPa": 1.0e3,
    "MPa": 1.0e6,
    "GPa": 1.0e9,
    "TPa": 1.0e12,
    "bar": 1.0e5,
    "kbar": 1.0e8,
    "Mbar": 1.0e11,
    "g": 1.0e-3,
    "kg": 1.0,
    "t": 1.0e3,
    "kt": 1.0e6,
    "Mt": 1.0e9,
    "m_earth": m_earth,
    "r_earth": r_earth,
    "mm": 1.0e-3,
    "cm": 1.0e-2,
    "m": 1.0,
    "km": 1.0e3,
    "Mm": 1.0e6,
    "K": 1.0,
    "kg/m3": 1.0,
    "gcc": 1.0e-3,
}

dens_max = 2.0e4
dens_min = 5.0e2
eps_P = 1.0e-4
eps_h = 1.0e-2
eps_r = 0.25
eps_layerfix = 1.0e-6
eps_layer = eps_layerfix
# eps_majorconstrain = 1.0e-4
eps_FeSi = 1.0e-6
eps_Mtot = 1.0e-6
eps_Rtot = 1.0e-6
eps_Psurf = 1.0e-3
eps_Tsurf = 1.0e-3
eps_Pspace = 1.0e-5
eps_Mg_number = 1.0e-6
acc_Mg_number = 1.0e-3
acc_M_surface = 1.0e-3
acc_T_surface = 1.0e-3

eos_pres_trans = 1.0e9

# pressure switch for H2O table eos
pressure_switch_H2O = 1.0e10

# max allowed iron content of the mantles
xi_Fe_mantle_max = 5e-1

# some plot parameters
plotrad_min = 0.0
plotdens_min = 0.0
plotpres_min = 0.0
plotmass_min = 0.0
plottemp_min = 0.0
plotrad_max = 1.4e0
plotdens_max = 1.5e4
plotpres_max = 4.0e2
plotmass_max = 1.5
plottemp_max = 6.0e3


material_plot_list_fort = [
    r"$\rm H_2O$",
    r"$\rm Fe$",
    "SiO2",
    r"$\rm MgO$",
    "MgSiO3",
    "Mg2SiO4",
    "Mg2Si2O6",
    "FeS",
    "Mg(OH)2",
]

# define colors for parameters
pres_color = (0.4, 0.5, 1.0)
dens_color = (0.8, 0.75, 0.0)
mass_color = (0.2, 0.8, 0.2)
temp_color = (1.0, 0.2, 0.2)
vesc_color = (0.5, 0.5, 0.5)
grav_color = (0.5, 0.5, 0.0)
moi_color = (0.0, 0.0, 0.0)

param_colors = [
    temp_color,
    dens_color,
    grav_color,
    pres_color,
    mass_color,
    moi_color,
    vesc_color,
]

pure_curves_color_list = [
    (0.5, 0.5, 1),  # water
    (0.5, 0.5, 0.5),  # Fe
    (0, 0.5, 0.5),  # SiO2
    (1, 0.5, 0.5),  # MgO
]

layerColors = {
    "liquid": (0.25, 0.5, 1),
    "supercritical": (0.5, 0.25, 1),
    "solid": (0.75, 0.8, 1),
    "vapor": (1, 0, 0),
    "gas": (0, 1, 0),
    "inner core": (0.75, 0.75, 0.75),
    "outer core": (1.0, 0.5, 0.25),
    "inner mantle": (0.4, 0.25, 0),
    "outer mantle": (0.5, 0.3, 0.0),
}

layerCodes = {
    1: "supercritical",
    2: "solid",
    0: "liquid",
    3: "inner core",
    4: "outer core",
    5: "inner mantle",
    6: "outer mantle",
}

# define color scheme for fancy plots
color_list = [
    (0.7, 0.2, 0.6),
    (0.6, 0.4, 1.0),
    (0.2, 0.4, 1.0),
    (0.2, 0.6, 0.9),
    (0.2, 0.8, 0.8),
    (0.2, 0.8, 0.4),
    (0.6, 0.8, 0.4),
    (0.6, 0.6, 0.2),
    (0.8, 0.4, 0.2),
    (1.0, 0.2, 0.2),
    (1.0, 0.5, 0.5),
    (0.7, 0.7, 0.7),
    (0.5, 0.5, 0.5),
    (0.2, 0.2, 0.2),
    (0.0, 0.0, 0.0),
]

# define parameter limits for plots
param_lims = [
    [0, 10000],  # T
    [0, 1.0e11],  # P
    [0, 10000],  # rho
    [0, 1.0e-7],  # dPdT_S
    [0, 1.0e8],  # dPdrho_T
    [0, 1.0e-4],  # alpha
    [0, 5000],  # cp
    [0.0, 0],  # s
    [0, 0],  # u
    [0, 5],  # phase
]

# define plot color for materials
# Conventions are:
# BM (0), MGD (1),
# liquid water (0), Mg2SiO4 (1), Mg2Si2O6 (2), Fe2SiO4 (3), Fe2Si2O6 (4), MgO (5),
# MgSiO3 (6), FeO (7), FeSiO3 (8), Fe (9), FeS (10), water ice VII (11)
material_colors = []

# suffix for different output files
suffix = {
    "planet": ".pl",
    "layer": ".lay",
    "shell": ".shl",
}  # .shl to avoid confusion with .sh shell scripts

lwdht = 1.0
laytranslwdht = 1.2
tlfntsize = 8
lfntsize = 10
gridalpha = 0.25
legendalpha = 0.75

# color for axis grid on plots
grid_color = (0.6, 0.6, 1.0)

# background color for figures
background_color = (1.0, 1.0, 1.0)

plot_params = {
    "lwdt": lwdht,
    "laytranslwdht": laytranslwdht,
    "tlfntsize": tlfntsize,
    "lfntsize": lfntsize,
    "gridcolor": grid_color,
    "gridalpha": gridalpha,
    "legendalpha": legendalpha,
    "legendedgecolor": "None",
    "backgroundcol": background_color,
    "pad": 5,
    "majorticklen": 10,
    "minorticklen": 5,
    "lwdthgridminor": 1,
    "lwdthgridmajor": 2,
}
