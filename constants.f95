MODULE constants

use class_weird_array

implicit none

!Generic universal constants
real(kind=8), parameter :: PI = 3.141592653589793
real(kind=8), parameter :: G = 6.67408E-11 !Gravity
real(kind=8), parameter :: Rgas = 8.3144598 !Gas constant
real(kind=8), parameter :: c_light = 2.99792458d8 !Speed of light
real(kind=8), parameter :: kB = 1.38065d-23 !Boltzmann constant
real(kind=8), parameter :: NA = 6.0221409d23 !Avogadro number

!Molar masses in kg
!Atomic species
real(kind=8), parameter :: mFe = 55.845d-3
real(kind=8), parameter :: mO = 15.999d-3
real(kind=8), parameter :: mH = 1.00794d-3
real(kind=8), parameter :: mHe = 4.0026d-3
real(kind=8), parameter :: mMg = 24.305d-3
real(kind=8), parameter :: mS = 32.065d-3
real(kind=8), parameter :: mSi = 28.0855d-3
real(kind=8), parameter :: mN = 14.0067d-3
real(kind=8), parameter :: mAr = 39.948d-3
real(kind=8), parameter :: mAl = 26.9815d-3
real(kind=8), parameter :: mCl = 40.078d-3

!Molecular substances
real(kind=8), parameter :: mH2O = 2.0d0*mH+mO
real(kind=8), parameter :: mOl = 2*mMg+mSi+4*mO
real(kind=8), parameter :: mPer = mMg + mO
real(kind=8), parameter :: mBr = mPer + mH2O
real(kind=8), parameter :: mEn = 2*mMg+2*mSi+6*mO
real(kind=8), parameter :: mPv = mMg+mSi+3*mO
real(kind=8), parameter :: mStv = mSi + mO*2

!Parameter for metal-silicate partitioning
!Si, O
real(8), dimension(2) :: a_KD = (/1.3d0, 0.6d0/) 
real(8), dimension(2) :: b_KD = (/-1.35d4, -3.8d3/)
real(8), dimension(2) :: c_KD = (/0d0, 22d0/)

real(8), dimension(2) :: a_KD_i = (/0.6d0, 0.1d0/) 
real(8), dimension(2) :: b_KD_i = (/-11.7d3, -2.2d3/)
real(8), dimension(2) :: c_KD_i = (/0d0, 5d0/)
real(8), dimension(2,2) :: eki = reshape((/0d0, -0.06d0, -0.11d0, -0.12d0/), shape(eki))

!~ eki(1,1) = 0d0
!~ eki(1,2) = -0.11d0
!~ eki(2,1) = -0.06d0
!~ eki(2,2) = -0.12d0
!Note that masses are in gram pre mole in Fischer et al. 2015
real(8), dimension(2) :: m_partitioning = (/mO * 1e-3, mSi * 1e-3/)

!Planetary paramters
!Masses & Radii
real(kind=8), parameter :: m_earth = 5.9722d24 !Earth mass in kg
real(kind=8), parameter :: r_earth = 6.371d6 !Earth radius in m

integer, parameter :: n_eos_tables = 10

!Primary composition parameters for different materials
!Note that the parameters are given for refractories and metals.
!In the current version of the model Fe is incorporated at the Mg site
!and Al at the Si site. Volatiles are treated via the secondary
!composition parameters For this reason in this treatment, YO and YH
!are 0 for pure water.
!
!Convention is:
!H2O, FeH, SiO2, MgO, MgSiO3, Mg2SiO4, Mg2Si2O6, FeS, FeO
integer, dimension(n_eos_tables) :: material_YMg = (/0,0,0,1,1,2,2,0,0,0/)
integer, dimension(n_eos_tables) :: material_YSi = (/0,0,1,0,1,1,2,0,0,1/)
integer, dimension(n_eos_tables) :: material_YO =  (/1,0,2,1,3,4,6,0,1,0/)
integer, dimension(n_eos_tables) :: material_YH =  (/2,0,0,0,0,0,0,0,0,0/)
integer, dimension(n_eos_tables) :: material_YS =  (/0,0,0,0,0,0,0,1,0,0/)
integer, dimension(n_eos_tables) :: material_YFe =  (/0,1,0,0,0,0,0,1,1,1/)

!Average density correction as function of H content x.
!rh(P,T,x) = rho(P,T, x=0) + drho*x
!This is only used for FeHx so drho = 0 for all other materials
real(8), dimension(n_eos_tables) :: drho_H_content = (/0d0, -1.26d3, 0d0, 0d0, 0d0, 0d0,&
										0d0, 0d0, 0d0, 0d0/)

!some EoS parameters for the individual materials
!0 means no value given
real, dimension(n_eos_tables) :: material_gammaG0 = (/0d0, 1.36d0, 1.71d0, 1.45d0, &
										 1.675d0, 0d0, 0d0, 1.2d0, 1.45d0, 1.3d0/)
real, dimension(n_eos_tables) :: material_q = (/0d0, 0.489d0, 1d0, 3d0, 1.39d0, &
									0d0, 0d0, 0.91d0,3d0,1.7d0/)


!Phase transition parameters for H2O
real(8) :: T_triple_H2O = 273.16d0
real(8) :: P_triple_H2O = 611.657d0
real(8) :: T_critical_H2O = 647.096d0
real(8) :: P_critical_H2O = 22.064d6
real(8) :: rho_critical_H2O = 322d0
real(8) :: P0_H2O = 1.01325d5
real(8) :: T0_H2O = 273.153d0
real(8) :: T0_melt = 273.153d0
real(8) :: s0 = 189.13d0 !j kg-1 K-1
real(8), dimension(3) :: LtoVII_coeffs_Dunaeva2010 = (/-2.29294290d-01, &
 4.11453047d+01,  5.46863218d+02/)
real(8), dimension(5) :: coeffs_French2010 = (/-2.26536372e-07,  &
2.64778698e-04, -1.19895845e-01,  2.37743503e+01, 3.96728823e+02/)

!Phase transition parameters for MgO
real(8) :: T_triple1_br = 794.5d0
real(8) :: T_triple2_br = 1087.0d0
real(8) :: P_triple1_br = 29.67d9
real(8) :: P_triple2_br = 33.65d9

real(8), dimension(5) :: bruce1 = (/-7.1999d-2, 8.20957d0, -3.46793d2, &
6.48875d3, -4.507d4/)
real(8), dimension(5) :: bruce2 = (/-1.01963d1, 1.26708d3, -5.90348d4, &
1.22208d6, -9.48251d6/)
real(8), dimension(2) :: bruce3 = (/-14.9321d0, 1589.76d0/)
real(8), dimension(5) :: bruce4 = (/-8.7169d-4, 1.1958d-1, -5.57672d0, &
8.55575d1, 1.084d3/)
real(8), dimension(3) :: bruce5 = (/-31.1131d0, 252.469d0, 863.148d0/)


END MODULE constants
