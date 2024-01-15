# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 09:02:20 2018

@author: os18o068
"""
# physical constants in SI
Rgas = 8.3144598  # universal gas constant
G = 6.67408e-11  # gravity
c_light = 2.99792458e8  # speed of light in m/s
kB = 1.38065e-23  # Boltzmann constant
kBcgs = kB * 1.0e7
sigmaSB = 5.670374e-8  # Stefan-Boltzmann constant
hPl = 6.62607015e-34  # Planck constant in J s
m_p = 1.6726219e-27  # proton mass in kg
NA = 6.0221409e23  # Avogadro8's constant
au = 1.49598e11  # astronomical unit in m
yr = 365.26  # year in days
day = 24.0 * 3600  # day in seconds
E_grav_unit = 2.241835488e32 # grav. bind. energy of uniform density earth

r_sun = 695510.0e3  # solar radius in m
T_sun = 5778  # solar surface temperature in K
m_sun = 1.989e30  # solar mass in kg

# Mg# and Si# from Sotin 2007, table 2, solar b
# for these values Al, Ca and Ni are represented by Si, Mg and Fe
Mg_number_solar = 0.534
Si_number_solar = 0.504

Si_earth = [0.47, 0.56]  # Allowed range for bulk composition from Sotin et al. 2007
Mg_earth = [0.47, 0.53]  # Allowed range for bulk composition from Sotin et al. 2007
SFe_earth = [0.081, 0.12]
SFe_solar = [0.2, 0.3]
SFe_poor = [0.0, 0.0]

# Bulk DMM composition according to Workmann 2005 or Allègre 1995
# Composition is given as wt% of relevant oxides
# Convention is: MgO, SiO2, FeO, CaO, Al2O3, NiO
bulk_DMM_lower_Workmann = [37.71, 44.7, 8.03, 3.17, 3.98, 0.24]
bulk_DMM_upper_Workmann = [38.73, 44.9, 8.18, 3.54, 4.44, 0.25]
bulk_DMM_Allegre = [37.77, 46.117, 7.485, 3.232, 4.09, 0.25]

y_oxides = [1, 1, 1, 1, 2, 1]

# molar masses in kg mol-1 for different elements:
mFe = 55.845e-3
mO = 15.999e-3
mH = 1.00794e-3
mHe = 4.0026e-3
mMg = 24.305e-3
mS = 32.065e-3
mSi = 28.0855e-3
mN = 14.0067e-3
mAr = 39.948e-3
mAl = 26.9815e-3
mCl = 35.453e-3
mCa = 40.078e-3
mNi = 58.6934e-3
mC = 12.0107e-3

# molar mass of air accounting for nitrogen, oxygen and argon
mAir = 2 * mN * 0.78 + 2 * mO * 0.21 + mAr * 0.01

# molar mass for some minerals and other substances
mH2O = 2 * mH + mO  # Water
mOl = 2 * mMg + mSi + 4 * mO  # iron free Olivine
mPer = mMg + mO  # Periclase
mBr = mMg + 2 * (mO + mH)  # Brucite
mEn = 2 * mMg + 2 * mSi + 6 * mO  # iron free Enstatite
mPerov = mMg + mSi + 3 * mO  # iron free Perovskite
mWs = mMg + mO  # iron free Wuestite
mStv = mSi + 2 * mO  # stishovite

# upper <-> lower transition parameters (Sotin 2007)
Ptrans0 = 2.5e10
Ttrans0 = 800.0
alphatrans = -0.0017e9

T_triple = 273.16  # triple point temperature
P_triple = 611.657  # triple point pressure in Pa
T_critical = 647.096
P_critical = 22.064e6

# high pressure ice transition parameters (Sotin 2007)
Ttrans0_hpi = 355.0
Ptrans0_hpi = 2.216e9
alphatrans_hpi = 0.534e9
beta_hpi = 5.22

# surface parameters for the earth (standardized values)
P_earth = 1.013e5
r_earth = 6.371e6
m_earth = 5.9722e24
V_ocean_earth = 1.3e18  # Total volume of the oceans on earth in m3

r_jupiter = 71.492e6
m_jupiter = 1.898e27

# temperature and pressure lower bound values for planetary interior to avoid
# spurious behaviour in the integration. These values do not apply for the
# atmopshere as both T and P must be allowed to become smaller
T_zero = 50.0
P_zero = 10.0
# density of atomic hydrogen at T=2.7 K and P=1.0e-11 Pa
rho_zero = 1.0e-17

# Temperature jumps at CMB, MTZ and Crust
temperature_jumps = [0.0, 800.0, 300.0, 1200.0]

# MR data for solar system planets
# convention is: mercury, venus, earth, moon, mars, jupiter, europa, Io,
# Ganymede, saturn, Titan, uranus, neptune, pluto
names_solar = [
    "Mercury",
    "Venus",
    "Earth",
    "Moon",
    "Mars",
    "Jupiter",
    "Europa",
    "Io",
    "Ganymede",
    "Saturn",
    "Enceladus",
    "Titan",
    "Uranus",
    "Neptune",
    "Pluto",
]

abbrevations_solar = [
    "M",
    "V",
    "E",
    "m",
    "M",
    "J",
    "e",
    "i",
    "g",
    "S",
    "e",
    "t",
    "U",
    "N",
    "P",
]

M_solar = [
    0.0553,
    0.815,
    1.0,
    0.0123,
    0.107,
    317.83,
    0.008035,
    0.015,
    0.0248,
    95.159,
    0.000018,
    0.0225,
    14.536,
    17.147,
    0.002,
]

R_solar = [
    0.3829,
    0.9499,
    1.0,
    0.2727,
    0.532,
    10.97,
    0.245,
    0.2859,
    0.4135,
    9.14,
    0.0395,
    0.4037,
    3.981,
    3.865,
    0.186,
]


MoI_solar = [
    0.346,
    None,
    0.3307,
    0.3929,
    0.3662,
    0.2756,
    0.346,
    0.37824,
    0.3115,
    0.22,
    0.3305,
    0.341,
    0.23,
    0.23,
]

delta_MoI_solar = [
    0.014,
    None,
    0.0,
    0.0009,
    0.0017,
    0.0006,
    0.005,
    0.00022,
    0.0028,
    None,
    0.0025,
    None,
    None,
    None,
]

T_solar = {
    "sun": 5778.0,
    "mercury": 340,
    "venus": 737,
    "earth": 287,
    "moon": 250,
    "mars": 210,
    "jupiter": 165.0,
    "europa": 102.0,
    "io": 110.0,
    "callisto": 134.0,
    "ganymede": 110.0,
    "saturn": 134.0,
    "enceladus": 75.0,
    "titan": 93.7,
    "uranus": 76.0,
    "neptune": 72.0,
    "pluto": 44.0,
}

M_solar = {
    "sun": 333.0e3,
    "mercury": 0.0553,
    "venus": 0.815,
    "earth": 1.0,
    "moon": 0.0123,
    "mars": 0.107,
    "jupiter": 317.83,
    "europa": 0.008035,
    "io": 0.015,
    "callisto": 0.018,
    "ganymede": 0.0248,
    "saturn": 95.159,
    "enceladus": 0.000018,
    "titan": 0.0225,
    "uranus": 14.536,
    "neptune": 17.147,
    "pluto": 0.002,
}

hosts_solar = {
    "sun": "",
    "mercury": "sun",
    "venus": "sun",
    "earth": "sun",
    "moon": "earth",
    "mars": "sun",
    "jupiter": "sun",
    "europa": "jupiter",
    "io": "jupiter",
    "callisto": "jupiter",
    "ganymede": "jupiter",
    "saturn": "sun",
    "enceladus": "saturn",
    "titan": "saturn",
    "uranus": "sun",
    "neptune": "sun",
    "pluto": "sun",
}


R_solar = {
    "sun": 109.0,
    "mercury": 0.3829,
    "venus": 0.9499,
    "earth": 1.0,
    "moon": 0.2727,
    "mars": 0.532,
    "jupiter": 10.97,
    "europa": 0.245,
    "io": 0.2859,
    "callisto": 0.378,
    "ganymede": 0.4135,
    "saturn": 9.14,
    "enceladus": 0.0395,
    "titan": 0.4037,
    "uranus": 3.981,
    "neptune": 3.865,
    "pluto": 0.186,
}

ecc_solar = {
    "sun": None,
    "mercury": 0.2056,
    "venus": 0.0068,
    "earth": 0.0167,
    "moon": 0.0549,
    "mars": 0.0934,
    "jupiter": 0.0484,
    "europa": 0.009,
    "io": 0.0041,
    "callisto": 0.0074,
    "ganymede": 0.0013,
    "saturn": 0.0541,
    "enceladus": 0.0047,
    "titan": 0.0288,
    "uranus": 0.0472,
    "neptune": 0.0086,
    "pluto": 0.2488,
}

# Distance to primary in 1.0e6 km
axes_solar = {
    "sun": 0.0,
    "mercury": 57.9,
    "venus": 108.2,
    "earth": 149.6,
    "moon": 0.384,
    "mars": 227.9,
    "jupiter": 778.6,
    "europa": 0.67,
    "io": 0.422,
    "callisto": 1.88,
    "ganymede": 1.07,
    "saturn": 1433.0,
    "enceladus": 0.0257,
    "titan": 1.222,
    "uranus": 2872.0,
    "neptune": 4495.0,
    "pluto": 5906.0,
}


MoI_solar = {
    "sun": 0.07,
    "mercury": 0.346,
    "venus": 0.337,
    "earth": 0.3307,
    "moon": 0.3929,
    "mars": 0.3662,
    "jupiter": 0.2756,
    "europa": 0.346,
    "io": 0.37824,
    "callisto": 0.3549,
    "ganymede": 0.3185,
    "saturn": 0.22,
    "enceladus": 0.3305,
    "titan": 0.3414,
    "uranus": 0.23,
    "neptune": 0.23,
    "pluto": 0.31,
}

# orbital period in days
orbit_solar = {
    "sun": 0.0,
    "mercury": 87.97,
    "venus": 225.0,
    "earth": yr,
    "moon": 27.32,
    "mars": 1.881 * yr,
    "jupiter": 11.86 * yr,
    "europa": 3.551,
    "io": 42 / 24,
    "callisto": 16.69,
    "ganymede": 7.155,
    "saturn": 29.46 * yr,
    "enceladus": 1.37,
    "titan": 15.945,
    "uranus": 84.01 * yr,
    "neptune": 164.8 * yr,
    "pluto": 248.59,
}


# rotation period in days
rotation_solar = {
    "sun": 0.0,
    "mercury": 87.97,
    "venus": 243,
    "earth": day,
    "moon": 27.32,
    "mars": 24.6 / 24,
    "jupiter": 0.41,
    "europa": 3.551,
    "io": 42 / 24,
    "callisto": 16.69,
    "ganymede": 7.155,
    "saturn": 0.44,
    "enceladus": 1.37,
    "titan": 15.945,
    "uranus": -0.71833,
    "neptune": 0.67,
    "pluto": 6.39,
}


# in units 1.0e-6
J2_solar = {
    "sun": None,
    "mercury": None,
    "venus": None,
    "earth": 0.001083,  # Yoder 1995
    "moon": 203.4,
    "mars": None,
    "jupiter": None,
    "europa": 435.5,
    "io": 1859.5,
    "callisto": 32.7,
    "ganymede": 127.53,
    "saturn": None,
    "enceladus": 5435.2,
    "titan": 31.8,
    "uranus": None,
    "neptune": None,
    "pluto": 0.282,
}  # McKinnon 1989

# in units 1.0e-6
C22_solar = {
    "sun": None,
    "mercury": None,
    "venus": None,
    "earth": 1.5744,  # Yoder 1995
    "moon": None,
    "mars": None,
    "jupiter": None,
    "europa": 131.5,
    "io": 558.8,
    "callisto": 10.2,
    "ganymede": 38.26,
    "saturn": None,
    "enceladus": None,
    "titan": 9.983,
    "uranus": None,
    "neptune": None,
    "pluto": None,
}

delta_MoI_solar = {
    "sun": None,
    "mercury": 0.014,
    "venus": 0.024,
    "earth": 0.0,
    "moon": 0.0009,
    "mars": 0.0005,
    "jupiter": 0.0006,
    "europa": 0.005,
    "io": 0.00022,
    "callisto": None,
    "ganymede": 0.0028,
    "saturn": None,
    "enceladus": 0.0025,
    "titan": 0.0005,
    "uranus": None,
    "neptune": None,
    "pluto": None,
}

# in units 1.0e-6
delta_J2_solar = {
    "sun": None,
    "mercury": None,
    "venus": None,
    "earth": None,
    "moon": None,
    "mars": None,
    "jupiter": None,
    "europa": 8.2,
    "io": 2.7,
    "callisto": 0.8,
    "ganymede": 2.9,
    "saturn": None,
    "enceladus": None,
    "titan": 0.4,
    "uranus": None,
    "neptune": None,
    "pluto": None,
}


# values in terrestrial units
M_trappist1 = [1.02, 1.16, 0.3, 0.77, 0.93, 1.15, 0.33]
R_trappist1 = [1.12, 1.1, 0.78, 0.92, 1.05, 1.15, 0.77]
M_trappist1_error = [0.15, 0.14, 0.04, 0.08, 0.08, 0.1, 0.06]
R_trappist1_error = [0.03, 0.03, 0.02, 0.03, 0.03, 0.03, 0.03]
names_trappist1 = ["b", "c", "d", "e", "f", "g", "h"]

# values in jovian units
names_others = [
    "K-011 b",
    "K-011 c",
    "K-029 c",
    "K-048 b",
    "K-060 c",
    "K-097 b",
    "K-098 b",
    "K-099 b",
    "K-102 d",
    "K-114 c",
    "K-138 c",
    "K-138 d",
    "K-128 b",
]

M_others = [
    0.006,
    0.009,
    0.01259,
    0.0124,
    0.01211,
    0.011,
    0.011,
    0.019,
    0.012,
    0.009,
    0.0062,
    0.002014,
    0.0023,
]
R_others = [
    0.161,
    0.256,
    0.28,
    0.168,
    0.17,
    0.132,
    0.178,
    0.132,
    0.105,
    0.143,
    0.107,
    0.1081,
    0.128,
]
M_others_error = [
    (0.004, 0.003),
    (0.009, 0.005),
    (0.00387, 0.00406),
    (0.00661, 0.00661),
    (0.00255, 0.00255),
    (0.006, 0.006),
    (0.005, 0.005),
    (0.004, 0.004),
    (0.006, 0.006),
    (0.002, 0.002),
    (0.00602, 0.00352),
    (0.002121, 0.001218),
    (0.00151, 0.00079),
]
R_others_error = [
    (0.003, 0.004),
    (0.004, 0.005),
    (0.018, 0.018),
    (0.009, 0.009),
    (0.013, 0.013),
    (0.012, 0.012),
    (0.02, 0.02),
    (0.007, 0.007),
    (0.004, 0.004),
    (0.016, 0.016),
    (0.006, 0.006),
    (0.0067, 0.0067),
    (0.008, 0.008),
]

solar_system = [M_solar, R_solar, names_solar]


# fit parameters for EOS
# Conventions are:
# BM (0), MGD (1),
# liquid water (0), Mg2SiO4 (1), Mg2Si2O6 (2), Fe2SiO4 (3), Fe2Si2O6 (4), MgO (5),
# MgSiO3 (6), FeO (7), FeSiO3 (8), Fe (9), FeS (10), water ice VII (11)
# additional materials can easily be added by adding the corresponding values to
# the appropriate arrays and specifying which EOS should be used for this
# material. If a new EOS has to be used for a material it can be implemented
# in PIMPeos.py
material_list_fort = [
    "H2O",
    "Fe",
    "SiO2",
    "MgO",
    "MgSiO3",
    "Mg2SiO4",
    "Mg2Si2O6",
    "FeS",
    "FeO",
    "FeSi",
]

material_list = [
    "H20",  # 0
    "Mg2SiO4",  # 1
    "Mg2Si2O6",  # 2
    "Fe2SiO4",  # 3
    "Fe2Si2O6",  # 4
    "MgO",  # 5
    "MgSiO3",  # 6
    "FeO",  # 7
    "FeSiO3",  # 8
    "Fe",  # 9
    "FeS",  # 10
    "Mg(OH)2",  # 11
    "alpha-(Mg,Fe)2SiO4",  # 12
    "beta-(Mg,Fe)2SiO4",  # 13
    "gamma-(Mg,Fe)2SiO4",  # 14
    "post-(Mg,Fe)SiO3",  # 15
    "SiO2 Stv",  # 16
    "SiO2 CaCl2-type",  # 17
    "FeSi",  # 18
    "Fe(l)",  # 19
    "Fe3C",  # 20
]

gas_material_list = ["H", "H2", "He", "O2", "air", "H2O", "N"]

material_plot_list = [
    r"\rm H_20",
    r"\rm Mg_2SiO_4",
    r"\rm Mg_2Si_2O_6",
    r"\rm Fe_2SiO_4",
    r"\rm Fe_2Si_2O_6",
    r"\rm MgO",
    r"\rm MgSiO_3",
    r"\rm FeO",
    r"\rm FeSiO_3",
    r"\rm Fe",
    r"\rm FeS",
    r"\rm Mg(OH)_2",
    r"\rm \alpha-(Mg,Fe)_2SiO_4",
    r"\rm \beta-(Mg,Fe)_2SiO_4",
    r"\rm \gamma-(Mg,Fe)_2SiO_4",
    r"\rm post-(Mg,Fe)SiO_3",
    r"\rm SiO_2 \ Stv",
    r"\rm CaCl_2-type SiO_2",
    r"\rm FeSi",
    r"\rm Fe(l)",
    r"\rm Fe_3C",
]

gas_material_plot_list = [
    r"$\rm H$",
    r"$\rm He$",
    r"$\rm O_2$",
    r"$ \rm H_20$",
    r"$\rm N$",
]

# molar mass in kg for each substance
molar_mass_list = [
    mH * 2 + mO,
    mMg * 2 + mSi + mO * 4,
    mMg * 2 + mSi * 2 + mO * 6,
    mFe * 2 + mSi + mO * 4,
    mFe * 2 + mSi * 2 + mO * 6,
    mMg + mO,
    mMg + mSi + mO * 3,
    mFe + mO,
    mFe + mSi + mO * 3,
    mFe,
    mFe + mS,
    mMg + 2 * mO + 2 * mH,
    None,
    None,
    None,
    mMg + mSi + mO * 3,
    mSi + 2 * mO,
    mSi + 2 * mO,
    mFe + mSi,
    mFe,
    mFe * 3 + mC,
]

# adhesion coefficient for van der waals eos in J m3 mol-2
a_VDW_list = [None, 24.76e-3, 3.45e-3, 137.8e-3, 135.8e-3, 557.29e-3, 137.0e-3]

# covolume coefficient for van der waals eos in m3 mol-1
b_VDW_list = [None, 26.61e-6, 23.7e-6, 31.8e-6, 36.4e-6, 31.0e-6, 38.7e-6]


# for hydrogen, the VDW eos yields reasonable results for P < 10MPa at
# 300 K (error of a few %) and very good results for higher temperatures or
# lower pressure
gas_molar_mass_list = [mH, 2 * mH, mHe, 2 * mO, mAir, 2 * mH + mO, mN]
gas_gamma_list = [1.36, 1.36, 1.67, 1.3, 1.4, 1.39]

# ambient bulk modulus in Pa for each material
K0_list = [
    2.2e9,  # H2O
    128.0e9,  # Mg2SiO4
    111.0e9,  # Mg2Si2O6
    128.0e9,  # Fe2SiO4
    111.0e9,  # Fe2Si2O6
    160.2e9,  # MgO Dorogokupets 2010
    256.7e9,  # MgSiO3, Wolf 2015
    137.8e9,  # B8-FeO, Fischer et al. 2011
    263.0e9,  # FeSiO3
    174.0e9,  # epsilon-Fe (hcp), Belonoshko 2010
    135.0e9,  # FeS Alibert 2014
    43.4e9,  # Mg(OH)2
    128.0e9,
    170.0e9,
    198.0e9,
    205.4e9,  # post-MgSiO3, Sun 2018
    302.0e9,  # SiO2 Stv, Fischer 2018
    341.0e9,  # CaCl2-type Silica, Fischer 2018
    230.6e9,  # FeSi Fischer 2014
    174.0e9,  # epsilon-Fe (hcp), Belonoshko 2010
    311.1e9,  # Fe3C (fcc), Takahashi 2019
]


K0prime_list = [
    4.0,  # H2O
    4.3,  # MgSi2O4
    7.0,  # Mg2Si2O6
    4.3,  # Fe2SiO4
    7.0,  # Fe2Si2O6
    3.99,  # MgO Dorogokupets 2010
    4.044,  # MgSiO3, Wolf 2015
    4.0,  # B8-FeO, Fischer et al. 2011
    3.9,  # FeSiO3
    5.3,  # epsilon-Fe (hcp), Belonoshko 2010
    6.0,  # FeS Alibert 2014
    5.4,  # Mg(OH)2
    4.3,
    4.6,
    5.3,
    5.069,  # post-MgSiO3, Sun 2018
    5.24,  # SiO2 Stv, Fischer 2018
    3.2,  # CaCl2-type Silica, Fischer 2018
    4.17,  # FeSi Fischer 2014
    5.3,  # epsilon-Fe (hcp), Belonoshko 2010
    3.4,  # Fe2C (fcc), Takahashi 2019
]

T0_list = [
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
    300.0,
]
aT_list = [
    0.0,
    2.832e-5,
    2.86e-5,
    2.832e-5,
    2.86e-5,
    None,
    2.302e-5,
    None,
    None,
    None,
    None,
    7.3e-5,
    2.832e-5,
    2.832e-5,
    2.832e-5,
    2.406e-5,
    None,
    None,
    None,
    None,
    None,
]
bT_list = [
    0.0,
    7.58e-9,
    7.2e-9,
    7.58e-9,
    7.2e-9,
    None,
    9.441e-9,
    None,
    None,
    None,
    None,
    3.6e-8,
    7.58e-9,
    7.58e-9,
    7.58e-9,
    1.034e-8,
    None,
    None,
    None,
    None,
    None,
]
cT_list = [
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    None,
    0.0,
    None,
    None,
    None,
    None,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    None,
    None,
    None,
    None,
    None,
]
aP_list = [
    None,
    -0.016e9,
    0.0,
    -0.016e9,
    0.0,
    None,
    -0.04e9,
    None,
    None,
    None,
    None,
    -0.011e9,
    -0.016e9,
    -0.016e9,
    -0.016e9,
    -0.033e9,
    None,
    None,
    None,
    None,
    None,
]

# Transition coeffs pv -> ppv from Townsend 2016 fig. 3
# The fit is P(T) = P0(GPa)+ a(GPa/K)*T(K)
# coeffs = [P0, a]
phase_transition_coeffs_MgSiO3 = [1.01e2, 8.93e-3]

# for Mg2SiO4 different phases are taken into account
# convention is: alpha, beta, gamma
# only K0, K0' and rho0 are assumed to change for the different phases. The
# temperature dependence on K0T and rho0T is assumed to be the same for all
# phases
phase_list_Mg2SiO4 = ["alpha", "beta", "gamma"]

K0prime_list_Mg2SiO4 = [4.3, 4.6, 5.3]

aP_list_Mg2SiO4 = [-0.016e9, -0.018e9, -0.021e9]

# From Ye et al 2009
aT_list_Mg2SiO4 = [2.19e-5, 1.6e-5, 2.54e-5]

# From Ye et al 2009
bT_list_Mg2SiO4 = [2.4e-8, 2.2e-8, 1.22e-8]

# From Ye et al 2009
cT_list_Mg2SiO4 = [0.0, 0.0, 0.0]

# From Ye eta al 2009
# Change of aP and bT as function of water content assuming linear change
daTdH2O_list_Mg2SiO4 = [-4.6e-6, -6.7e-6, -9.7e-6]

dbTdH2O_list_Mg2SiO4 = [0.9e-8 * 0.5, 2.5e-8 * 0.5, 2.6e-8 * 0.5]

# fitted coefficients for linear phase transition pressures (Pa) as function of
# temperature (°C) for alpha-beta and beta-gamma transition in the pure and
# anhydrous Mg2SiO4 system
phase_transition_coeffs_Mg2SiO4 = [[1.159e1, 1.999e-3], [1.089e1, 5.679e-3]]

# fit parameters for water content in saturated Mg2SiO4 as function of
# pressure (GPa), temperature (°C), iron number (mol%)

saturation_coeffs_Mg2SiO4 = [
    [
        [
            [-0.003059202294896886, 0.03497819014827574],
            [0.0012983169164147625, 0.0012182638756681875],
        ],
        [
            [1.952572561001825e-06, -2.3875892100643497e-05],
            [-7.538595758202829e-07, 0.0],
        ],
    ],
    [
        [[0.11918296124562064, -0.0017502866503096384], [-0.0040531056078208505, 0.0]],
        [[-5.392791752515525e-05, 0.0], [1.892085047859246e-06, 0.0]],
    ],
    [
        [[-0.05671829751163112, 0.023029510635385924], [0.00793727381682779, 0.0]],
        [[2.740249961349819e-05, 0.0], [-4.069716435073367e-06, 0.0]],
    ],
]

# ambient density and bulk modulus depend linearly on iron content and
# water content. Here are the coefficients for the linear fits by Mao 2016 eq.8
# pressure derivative of bulk modulus is assumed to be constant due to a lack
# of consitent experimental data
# for the density for alpha, beta, gamma phase
# conventions are: rho0, C_Fe, C_H20
# rho0(X_Fe, X_H2O) = coeff[0] + coeff[1]*X_Fe + coeff[2]*X_H2O
rho0_coeffs_Mg2SiO4 = [
    [3.222e3, 0.012e3, -0.049e3],
    [3.468e3, 0.013e3, -0.047e3],
    [3.562e3, 0.015e3, -0.048e3],
]

# for the bulk modulus of alpha, beta gamma phase from Mao 2016 (eq. 10,
# eq. 11 and eq. 12 respectively)
# conventions are: K0, C_Fe, C_H20, C_Fe,H2O
K0_coeffs_Mg2SiO4 = [
    [128.0e9, 0.1e9, -3.8e9, -0.4e9],
    [170.0e9, 0.0e9, -12.4e9, 0.5e9],
    [185.0e9, 0.36e9, -11.8e9, 0.0e9],
]


# for the shear modulus of alpha, beta gamma phase from Mao 2016
# conventions are: G0, C_Fe, C_H2O, C_Fe,H2O
G0_coeffs_Mg2SiO4 = [
    [81.6, -0.3, -2.2, -0.52],
    [114.0, -0.8, -9.7, 0.5],
    [120.4, -0.1, -5.6, -0.9],
]

# coefficients to fit relative density change of Mg(OH)2 with respect to
# MgO. The form of the fit is:
#                   Delta(T, P) = c1*log(P) + c2*log(P)**2 + c3*log(P)**3
#                               + c4*log(P)**4 + c5*T + c6*T*log(P)
#                               + c7*T*log(P)**2 + c8*T*log(P)**3
#                               + c9

delta_rho_MgO_coeffs = [
    1.3283,
    -2.8464e-1,
    2.66913e-2,
    -9.2556e-4,
    4.8976e-4,
    -1.9528e-4,
    2.9641e-05,
    -1.4748e-06,
    -1.9716,
]

# for Perovskite and post-perovskite the effect of iron is computed by scaling
# the results in table1 from Sun 2018 (from Wolf 2015 for 13 mol% Fe) to
# the iron content in the material
# compute deltas for all MGD parameters for 1 mol% Fe
dq_pv = (0.56 - 1.5) / 13
dgamma0_pv = (1.4 - 1.54) / 13
dthetaD0_pv = (1000 - 950) / 13
dK0prime_pv = 0.0
dd0_pv = 31.0  # use linear scaling from Sinmyo 2014 table 2 (at 26 GPa, 2073 K)
dK0_pv = (238.4e9 - 256.7e9) / 13

# the same is done for post-perovskite from Shieh 2005 for 9 mol% Fe in
# comparison with 0 mol% Fe from Sun 2018
# only values for V0, K0 and K0prime are given, the thermal pressure
# parameters are taken from perovskite
dq_ppv = dq_pv
dgamma0_ppv = dgamma0_pv
dthetaD0_ppv = dthetaD0_pv
dK0_ppv = (255.0e9 - 205.4e9) / 9.0
dK0prime_ppv = (5.069 - 3.7) / 9.0
dd0_ppv = 33.6  # use linear scaling from Mao 2014 fig. 6

# Hydration of post-Perovskite from Townsend 2015.
# Only effect on rho0, K0 and K0prime have bee

# Note: in Sotin 2007 there is a wrong value for FeS, i.e. 4.9. Comparing with
# the reference for FeS a value of 5.9 seems to represent the measurements
# to better agreement
rho0_list = [
    1.0e3,  # H2O
    3.222e3,  # Mg2SiO4
    3.215e3,  # Mg2Si2O6
    4.404e3,  # Fe2SiO4
    4.014e3,  # Fe2Si2O6
    3.583e3,  # MgO Dorogokupets 2010
    4.107e3,  # MgSiO3, Sun 2018
    5.989e3,  # B8-FeO, Fischer et al. 2011
    5.178e3,  # FeSiO3
    8.334e3,  # epsilon-Fe (hcp), Belonoshko 2010
    4.9e3,  # FeS Alibert 2014
    2.323e3,  # Mg(OH)2
    3.222e3,  # alpha-Ol
    3.468e3,  # beta-Ol
    3.562e3,  # gamma-Ol
    4.020e3,  # post-MgSiO3, adjusted to match density contrast from Dorfman 2014
    4.287e3,  # SiO2 Stv, Fischer et al. 2018
    4.287e3,  # CaCl2-type Silica, Fischer et al. 2018
    6.543e3,  # FeSi Fischer 2014
    8.334e3,  # epsilon-Fe (hcp), Belonoshko 2010
    8.01e3,  # Fe3C (fcc), Takahashi 2019
]

gamma0_list = [
    None,  # H2O
    None,  # Mg2SiO4
    None,  # Mg2Si2O6
    None,  # Fe2SiO4
    None,  # Fe2Si2O6
    1.524,  # MgO Dorogokupets 2010
    1.54,  # MgSiO3, Wolf 2015
    1.45,  # FeO Alibert 2014
    1.96,  # FeSiO3
    2.434,  # epsilon-Fe (hcp), Belonoshko 2010
    1.36,  # FeS Alibert 2014
    0.43,  # Mg(OH)2
    None,  # alpha
    None,  # beta
    None,  # gamma
    1.495,  # post-MgSiO3, Sun 2018
    1.71,  # SiO2 Stv, Fischer 2018
    2.14,  # CaCl2-type Silica, Fischer 2018
    1.3,  # FeSi Fischer 2014
    2.434,  # epsilon-Fe (hcp), Belonoshko 2010
    1.06,  # Fe3C (fcc), Takahashi 2019
]


thetaD0_list = [
    None,  # H2O
    None,  # Mg2SiO4
    None,  # Mg2Si2O6
    None,  # Fe2SiO4
    None,  # Fe2Si2O6
    599 / 0.775,  # MgO Dorogokupets 2010
    950,  # MgSiO3, Wolf 2015
    430.0,  # FeO Alibert 2014
    1017.0,  # FeSiO3
    227.0,  # Fe, Belonoshko 2010
    998.0,  # FeS Alibert 2014
    1470.0,  # Mg(OH)2
    None,
    None,
    None,
    995.0,  # post-MgSiO3, Sun 2018
    1109.0,  # SiO2 Stv, Fischer 2018
    1109.0,  # CaCl2-type Silica, Fischer 2018
    417.0,  # FeSi Fischer 2014
    227.0,  # epsilon-Fe (hcp), Belonoshko 2010
    314.0,  # Fe3C (fcc), Takahashi 2019
]

q_list = [
    None,
    None,
    None,
    None,
    None,
    1.65,  # MgO Dorogokupets 2010
    1.5,  # MgSiO3, Wolf 2015
    3.0,  # FeO, Alibert 2014
    2.5,
    0.489,  # Fe, Belonoshko 2010
    0.91,  # FeS Alibert 2014
    1.0,  # Mg(OH)2
    None,
    None,
    None,
    1.97,  # post-MgSiO3, Sun 2018
    1.0,  # SiO2 Stv, Fischer 2018
    1.0,  # CaCl2-type Silica, Fischer 2018
    1.7,  # FeSi Fischer 2014
    0.489,  # epsilon-Fe (hcp), Belonoshko 2010
    1.92,  # Fe2C (fcc), Takahashi 2019
]


# for Mg(0H)2 different phases are taken into account. Values are taken from
# the ab initio data of Hermann 2016
# convention is: alpha (=Br-P3), beta (=Br-P41212), gamma (=Per)
# where we assume in the case of the dissociated phase gamma that all water
# will be transportet to the top and only pure MgO remains
phase_list_Brucite = ["alpha", "beta", "gamma"]

K0_list_Brucite = [43.4e9, 67.3e9, K0_list[5]]

K0prime_list_Brucite = [5.4, 4.9, K0prime_list[5]]

rho0_list_Brucite = [2323.0, 2544.0, rho0_list[5]]

n_list = [3, 7, 10, 7, 10, 2, 5, 2, 5, 1, 2, 3, 5, 7, 7, 7, 3, 3, 2, 1, 4]

# Si, O (Fischer et al. 2015)
a_KD = [1.3, 0.6]
b_KD = [-1.35e4, -3.8e3]
c_KD = [0.0, 22.0]

# Si, O (Fischer et al. 2015)
a_KDi = [0.6, 0.1]
b_KDi = [-11.7e3, -2.2e3]
c_KDi = [0.0, 5.0]
eki = [[0.0, -0.11], [-0.06, -0.12]]

# Transition pressure from non-ideal to ideal FeO-Fe solution
# Schaefer et al. 2017
P_FeO_Fe_trans = 50e9

# Margules parameters for Fe-FeO solutions from Frost et al. 2010
Margules_FeO_Fe = [[83307.0, 8.978, 0.09], [135943.0, 31.122, 0.059]]

# Note: conventions here are consistent with the fortran code!
# H2O, FeH, SiO2, MgO, MgSiO3, Mg2SiO4, Mg2Si2O6, FeS, FeO
material_YMg = [0, 1, 0, 1, 1, 2, 2, 0, 0, 0, 0]
material_YSi = [0, 0, 1, 0, 1, 1, 2, 0, 0, 1, 0]
material_YO = [1, 0, 2, 1, 3, 4, 6, 0, 1, 0, 0]
material_YH = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
material_YS = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
material_YFe = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]

EOSstring_list = ["BM3", "MGD", "Vinet", "Bel", "water EOS", "IGL", "const", "VR"]
EOS_type = [4, 0, 0, 0, 0, 1, 1, 1, 1, 3, 1, 0, 0, 0, 0, 1, 1, 1, 1, 3, 1]
gas_EOS_type = [5, 5, 5, 5, 5]
