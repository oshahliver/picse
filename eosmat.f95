MODULE class_unit

use functions
use eosfort
use constants
use phase

implicit none

type unit

integer :: ll, phase
real(8) :: temp, pres, molar_mass
real(8) :: dens, Fe_number, FeMg, X_H2O, eps_H2O, eps_Al, X_Al
real(8) :: K_isoth, dPdrho, alpha_th, xi_H2O, xi_Al, xi_Fe, &
xi_AlSi, xi_AlMg, xi_H
logical :: saturation
logical :: constant_Fe = .true.

end type unit

contains


SUBROUTINE init_unit(self, P, T, eps_H2O, eps_Al, xi_Fe, ll, xi_H)

type(unit), intent(inout) :: self
integer, intent(in) :: ll
real(8), intent(in) :: P, T, eps_H2O, eps_Al, xi_Fe
real(8), intent(inout), optional :: xi_H

self%pres=P
self%temp=T
self%eps_H2O = eps_H2O
self%eps_Al = eps_Al
self%xi_Fe = xi_Fe
self%ll = ll
self%phase = 1

if(present(xi_H))then
	self%xi_H = xi_H
else
	self%xi_H = 0d0
endif

END SUBROUTINE init_unit

!#######################################################################
SUBROUTINE compute_unit(self, T, P, order, alloc)

type(unit), intent(inout) :: self
integer, intent(inout), optional :: order
integer :: order_dummy
real(8), intent(in) :: T, P
real(8) :: T_dummy, P_dummy
integer, dimension(n_out) :: which
real(8), dimension(n_out) :: res
real(8), dimension(3) :: vals
logical, intent(inout), optional :: alloc
logical :: alloc_dummy
integer :: i, ll
do i=1, n_out
	which(i) = i
enddo
vals(1) = T
vals(2) = P
P_dummy = P
T_dummy = T
!For FeH_x, the table needs xi_Fe = 1-xi_H in order to compute the eos
!parameters consistently for a given hydrogen content. However, in the
!shell xi_Fe is used to define the ratio Fe/Mg, which is not the same thing.
!So for FeH_x, the iron content for the table interpolation must be 1-xi_H
!without changing xi_Fe = 1 in order to avoid counting Mg in the core.
if (self%ll == 2)then
	vals(3) = 1d0 - self%xi_H
else
	vals(3) = self%xi_Fe
endif
!vals(4) = self%eps_Al
!vals(5) = self%eps_H2O

if(.not.present(alloc))then
  alloc_dummy = .false.
else
  alloc_dummy = alloc
endif
if(present(order))then
  order_dummy=order
else

  order_dummy = 1
endif

self%temp=T
self%pres=P
call get_phase(ll=self%ll, P=P_dummy, T=T_dummy, ph=self%phase, xiFe=vals(3))
call compute(which=which, n_out=n_out, order=order_dummy, &
alloc=alloc_dummy, vals=vals, res=res, ll=self%ll, eps_H2O = self%eps_H2O)
!Here some parameters are updated
!The density is corrected for the H content. This is only relevant
!For FeHx in the core and only valid for x << 1
self%dens = res(1)
self%dPdrho = res(3)
self%alpha_th = res(4)

!For pure water, the water content must be set to 1
if(self%ll==1)then
	self%X_H2O = 1d0
	self%xi_Al = 0d0

elseif(self%ll==2)then
	!In pure iron, the water equivalent content is measured via xi_H
	!and xi_H2O is set to zero to make sure that no H2O is counted in
	!the core but only H which can then be converted into H2O equivalent
	self%xi_H2O = 0d0
	
	self%X_H2O = compute_eta(xi=self%xi_H2O, m1=mFe, m2=mH2O)
	self%xi_Al = 0d0

else
	!Here the water and aluminum contents are extracted from the table.
	!In general these parameters depend on T,P,xi_Fe,eps_H2O,eps_Al
	self%X_H2O = res(5)
	self%xi_Al = res(6)
	!print *, 'll, ph, P, T, X_H2O =', self%ll, self%phase, self%pres*1d-9, self%temp,  self%X_H2O
	
endif

!print *, 'dens/dPdrho/alpha in unit =', self%dens, self%dPdrho, self%alpha_th

END SUBROUTINE compute_unit

END MODULE class_unit

!#######################################################################
!#######################################################################
MODULE class_mixture

use class_unit

implicit none

type mixture

integer :: lay, n_mats
integer, dimension(:), allocatable :: contents, YO, YSi, YMg, YS, YFe
real(8), dimension(:), allocatable :: fractions, densities, weight_fractions
real(8), dimension(:), allocatable :: dPdrhos, alpha_ths, xi_Fe, &
xi_AlSi, xi_AlMg, xi_H2O, xi_Al, X_H2O, molar_masses
type(unit), dimension(:), allocatable :: units
integer, dimension(:), allocatable :: additionals
real(8) :: temp, pres, dens, eps_H2O, dPdrho
real(8) :: K_isoth, alpha_th, SiMg, FeMg, eps_Al, Fe_number, &
Si_number, xi_H, xi_Stv
logical :: saturation, force_bisection
logical :: constant_Fe = .true.

end type mixture

contains

!#######################################################################
SUBROUTINE init_mixture(self, n_mats, fractions, contents, P, T, &
eps_H2O, eps_Al, alloc, Fe_number, lay, &
Si_number, xi_H, xi_Stv)

implicit none

type(mixture), intent(inout) :: self
integer, intent(in) :: n_mats
integer, intent(inout), optional :: lay
integer, dimension(n_mats), intent(in) :: contents
real(8), dimension(n_mats), intent(in) :: fractions
real(8), intent(in) :: P, T, eps_H2O, eps_Al
real(8), intent(inout), optional :: Fe_number, Si_number, xi_H, xi_Stv
integer :: i, order
logical, intent(inout), optional :: alloc
logical :: alloc_dummy, negative=.false.

order=1

if(.not.present(alloc))then
  alloc_dummy = .false.
else
  alloc_dummy = alloc
endif

!~ print *, '		init mixture'
!To avoid some shit check if stuff is already allocated. If yes don't 
!try to do it again because more shit will happen. Don't touch my allocated
!stuff man!
if(.not.allocated(self%weight_fractions))then
	allocate(self%weight_fractions(n_mats))
	allocate(self%contents(n_mats))
	allocate(self%fractions(n_mats))
	allocate(self%units(n_mats))
	allocate(self%densities(n_mats))
	allocate(self%alpha_ths(n_mats))
	allocate(self%dPdrhos(n_mats))
	allocate(self%YO(n_mats))
	allocate(self%YSi(n_mats))
	allocate(self%YMg(n_mats))
	allocate(self%YS(n_mats))
	allocate(self%YFe(n_mats))
	allocate(self%X_H2O(n_mats))
	allocate(self%xi_Al(n_mats))
	allocate(self%xi_AlMg(n_mats))
	allocate(self%xi_AlSi(n_mats))
	allocate(self%xi_H2O(n_mats))
	allocate(self%xi_Fe(n_mats))
	allocate(self%molar_masses(n_mats))
endif

self%contents = contents
self%fractions = fractions
self%n_mats = n_mats
!~ print *, 'fracs in init mixture =', self%fractions
do i=1, n_mats
  self%YO(i) = material_YO(self%contents(i))
  self%YMg(i) = material_YMg(self%contents(i))
  self%YSi(i) = material_YSi(self%contents(i))
  self%YS(i) = material_YS(self%contents(i))
  self%YFe(i) = material_YFe(self%contents(i))
enddo

self%force_bisection = .false.
self%eps_H2O = eps_H2O
self%eps_Al = eps_Al
self%pres = P
self%temp = T
self%dens = 0.0d-10
self%Fe_number = Fe_number
self%Si_number = Si_number
self%xi_Stv = xi_Stv

if(present(lay))then
  self%lay = lay
else
  self%lay = 1
endif
if(present(xi_H))then
	self%xi_H = xi_H
else
	self%xi_H = 0.0d0
endif
if(present(xi_Stv))then
	self%xi_Stv = xi_Stv
else
	self%xi_Stv = 0.0d0
endif
if(.not.Si_number==1.0d0)then
  self%SiMg = self%Si_number/(1.0d0 - self%Si_number)
else
  self%SiMg = 1.0d10
endif
if(.not.Fe_number==1.0d0)then
  self%FeMg = self%Fe_number/(1.0d0 - self%Fe_number)
else
  self%FeMg = 1.0d10
endif
!Here the iron partitioning can be taken into account but is not currently
!For the time being the iron content is just equally distributed over
!the different materials

!self%xi_Fe = compute_xiFei(n_mats=size(self%contents), &
!P=self%pres, &
!T=self%temp, &
!Fe_number=self%Fe_number, &
!contents=self%contents)
self%xi_Fe = self%Fe_number

!print *, 'vals in init mixture =', self%temp, self%pres
do i=1, n_mats
  call init_unit(self=self%units(i), T=T, P=P, eps_H2O=self%eps_H2O, &
  eps_Al=self%eps_Al, xi_Fe=self%xi_Fe(i), ll=self%contents(i), &
  xi_H = self%xi_H)
  call compute_unit(self=self%units(i), T=self%temp, P=self%pres, &
  order=order, alloc=alloc_dummy)
  self%densities(i) = self%units(i)%dens
  self%alpha_ths(i) = self%units(i)%alpha_th
  self%dPdrhos(i) = self%units(i)%dPdrho
  self%X_H2O(i) = self%units(i)%X_H2O
  self%xi_Al(i) = self%units(i)%xi_Al
  !Compute dry molar masses
  !For water impurity, do not account for iron
	if (self%contents(i) == 1)then
		self%molar_masses(i) = mO + 2 * mH
	else
		self%molar_masses(i) = self%YMg(i) * ((1.0d0 - self%xi_Fe(i)) * mMg + &
		self%xi_Fe(i) * mFe) + self%YSi(i) * ((1.0d0 - self%xi_Al(i)) * mSi + &
		self%xi_Al(i) * mAl) + self%YO(i) * mO + &
		self%YS(i) * mS + &
		self%YFe(i) * mFe
	endif
  self%xi_AlSi(i) = 0.0d0
  self%xi_AlMg(i) = 0.0d0
  
  !Compute the molar water content in material i given it's mass fraction
  self%xi_H2O(i) = compute_xi(eta=self%X_H2O(i), m1=self%molar_masses(i), m2=mH2O)
  !print *, 'X_H2O, xi_H2O =', self%X_H2O(i), self%xi_H2O(i)

enddo
!~ print *, '		fractions in init mix 1 =', self%fractions!, self%weight_fractions
call update_mixture_fractions(self=self)
!~ print *, '		fractions in init mix 2 =', self%fractions!, self%weight_fractions
END SUBROUTINE init_mixture

!#######################################################################
SUBROUTINE update_mixture_fractions(self)

type(mixture), intent(inout) :: self
real(8), dimension(self%n_mats-2) :: additional
integer :: i
!print *, 'Fe# =', self%Fe_number

!Here additional contains only one component. In principle it could
!contain as many additional materials as desired.
if (size(self%fractions)>2)then
	do i=1, self%n_mats-2
		additional(i) = self%fractions(2+i)
	enddo
else
	additional(:) = 0d0
endif

if (self%lay>2)then
	call compute_abundance_vector(SiMg=self%SiMg, &
	FeMg=self%FeMg, &
	n_mats=size(self%contents), &
	YSii=self%YSi, &
	YMgi=self%YMg, &
	xiFei=self%xi_Fe, &
	xiH2Oi=self%xi_H2O, &
	xiAlMgi=self%xi_AlMg, &
	xiAlSii=self%xi_AlSi, &
	abundances=self%fractions, &
	contents=self%contents, &
	additional = additional)
endif

self%weight_fractions = eta_general(self%fractions, self%molar_masses, &
size(self%fractions))

!~ print *, '---'
!~ print *, 'contents in mix =', self%contents(:)
!~ print *, 'fractions in mix =', self%fractions(:)
!~ print *, 'weight fractions in mix =', self%weight_fractions(:)
!~ print *, 'molar masses =', self%molar_masses

END SUBROUTINE update_mixture_fractions

!#######################################################################
SUBROUTINE update_mean_mixtures(self)

type(mixture), intent(inout) :: self
integer :: i

!Update properties
self%dens = 0.0d0
self%dPdrho = 0.0d0
self%alpha_th = 0.0d0

!~ print *, 'contents =', self%contents
!~ print *, 'fractions =', self%fractions
!~ print *, ' weight fractions =', self%weight_fractions 

!Use linear mixing law to compute density, dPdrho and thermal expansion
do i=1, size(self%contents)
  self%dens = self%dens + self%weight_fractions(i)/self%densities(i)
  self%dPdrho = self%dPdrho + self%weight_fractions(i)/self%dPdrhos(i)/ &
				self%densities(i)**2
  self%alpha_th = self%alpha_th + self%weight_fractions(i)/self%densities(i)* &
					self%alpha_ths(i)
enddo

self%dens = 1.0d0/self%dens
self%dPdrho = self%dPdrho*self%dens**2
self%dPdrho = 1.0d0/self%dPdrho

self%K_isoth = self%dPdrho * self%dens

self%alpha_th = self%alpha_th* self%dens

!print *, 'stuff in update mean:', self%densities(:), self%dPdrhos(:), self%alpha_ths(:)
!print *, 'results in update mean:', self%dens, self%dPdrho, self%alpha_th

END SUBROUTINE update_mean_mixtures

!#######################################################################
SUBROUTINE compute_mixture(self, T, P)

implicit none

type(mixture), intent(inout) :: self
real(kind=8), intent(in), optional :: T, P
integer :: i, order

if (present(P))then
  self%pres = P
endif

if (present(T))then
  self%temp = T
endif

order=1
!Compute all parameters for individual materials
do i=1, size(self%contents)
  call compute_unit(self=self%units(i), T=self%temp, P=self%pres, order=order)
  self%densities(i) = self%units(i)%dens
  self%alpha_ths(i) = self%units(i)%alpha_th
  self%dPdrhos(i) = self%units(i)%dPdrho
  self%X_H2O(i) = self%units(i)%X_H2O
  self%xi_Al(i) = self%units(i)%xi_Al
enddo

call update_mean_mixtures(self=self)

END SUBROUTINE compute_mixture


!SUBROUTINE update_mixture_fractions(self, newfractions)

!type(mixture), intent(inout) :: self
!real(kind=8), dimension(size(self%contents)), intent(in) :: newfractions

!self%fractions = newfractions

!call update_mean_mixtures(self=self)

!END SUBROUTINE update_mixture_fractions

!########################################################################
SUBROUTINE update_mixture(self, T, P)

type(mixture), intent(inout) :: self
real(kind=8), intent(in) :: T, P
integer :: i

self%temp = T
self%pres = P

!Update individual material units
do i=1, size(self%contents)
  self%units(i)%temp=T
  self%units(i)%pres=P
enddo

END SUBROUTINE update_mixture

!#######################################################################
SUBROUTINE print_mixture(self, indv)

implicit none

type(mixture), intent(inout) :: self
logical, intent(in), optional :: indv
logical :: individual
integer :: i

if(present(indv))then
	individual=indv
else
	individual = .false.
endif

print *, ''
print *, 'Mixture properties:'
print *, 'P/T:', self%pres, self%temp
print *, 'rho:', self%dens
print *, 'dPdrho:', self%dPdrho
print *, 'K_isoth:', self%K_isoth
print *, 'Contents:', self%contents(:)
print *, 'Fractions:', self%fractions(:)

if(individual)then
	do i=1, size(self%contents)
	  print *, '----------'
	  print *, 'material', i,':', self%units(i)%dens, self%units(i)%K_isoth, self%units(i)%alpha_th
	enddo
endif
END SUBROUTINE print_mixture


END MODULE class_mixture
