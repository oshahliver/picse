MODULE functionspy

   ! use constants
   ! use run_params
   ! use eosfort
   use linalg
   use constants

   implicit none
   private
   public :: construct_abundance_matrix, compute_abundance_vector, compute_core_mass, get_core_mass_q_vector

contains

!########################################################################
   ! subroutine progress_bar(j, n)

   !    integer :: j, k, perc, dec1, dec2
   !    real(8) :: jj
   !    integer, intent(in) :: n
   !    character(len=72) :: bar = "???.??% |                                                   |"

   !    jj = real(j)/real(n)
   !    perc = int(jj*100)
   !    dec1 = int(((jj*1d2) - real(perc))*1d1)
   !    dec2 = int(((jj*1d2) - real(perc))*1d2)
   !    dec2 = dec2 - 10*dec1

   !    write (unit=bar(1:3), fmt="(i3)") perc
   !    write (unit=bar(5:5), fmt="(i1)") dec1
   !    write (unit=bar(6:6), fmt="(i1)") dec2

   !    do k = 1, perc/2
   !       bar(9 + k:9 + k + 1) = "=>"
   !    end do

   !    write (unit=6, fmt="(a1, a72)", advance="no") char(13), bar

   !    if (jj /= 100d0) then
   !       flush (unit=6)

   !    else

   !       write (unit=6, fmt=*)

   !    end if
   !    return

   ! end subroutine progress_bar

!    function fnc(x) result(y)
!       implicit none
!       real(8), intent(in) :: X
!       real(8) :: y

!       y = x

!    end function fnc

! ! ########################################################################
!    function MOI_integrand_linear(r1, r2, rho1, rho2) result(moi)

!       implicit none

!       real(8), intent(in) :: r1, r2, rho1, rho2
!       real(8) :: moi
!       real(8) :: dr, a, moi1, moi2

!       dr = r2 - r1
!       a = (rho2 - rho1)/dr

!       moi1 = 1d0/5d0*rho1*r1**5 + 1d0/6d0*a*r1**6 - 1d0/5d0*a*r1*r1**5
!       moi2 = 1d0/5d0*rho1*r2**5 + 1d0/6d0*a*r2**6 - 1d0/5d0*a*r1*r2**5

!       moi = moi2 - moi1

!    end function MOI_integrand_linear

! !########################################################################
!    function eta_general(xi, m, n_mats) result(y)
! !Compute weight fractions from mole fractions of different components.

!       implicit none

!       integer :: i, n_mats
!       real(8), dimension(n_mats) :: y
!       real(8) :: m_tilde
!       real(8), intent(in), dimension(n_mats) :: xi, m

!       y = 0d0
!       m_tilde = 0d0

!       do i = 1, n_mats
!          m_tilde = m_tilde + m(i)*xi(i)
!       end do

!       do i = 1, n_mats
!          y(i) = xi(i)*m(i)/m_tilde
!       end do

!    end function eta_general

! !########################################################################
!    function xi_general(eta, m, n_mats) result(y)
! !Compute mole fractions from weight fractions of different components.

!       implicit none

!       integer :: i, n_mats
!       real(8), dimension(n_mats) :: y
!       real(8) :: dummy
!       real(8), intent(in), dimension(n_mats) :: eta, m

!       print *, 'eta =', eta

!       dummy = 0d0
!       do i = 1, n_mats
!          dummy = dummy + eta(i)/m(i)
!       end do

!       do i = 1, n_mats
!          y(i) = eta(i)/m(i)/dummy
!       end do

!    end function xi_general

! !#######################################################################
!    function at2mat(at, n_mats) result(y)
! !Convert atomic mole fractions to molecular mole fractions. The atomic
! !composition of the molecular species is given by the matrix N and is
! !fixed to represent the outer core composition at this point. In future
! !versions it could be passed as an argument to allow for a more diverse
! !compositions.

!       integer, intent(in) :: n_mats
!       real(8), dimension(n_mats) :: y, x
!       real(8), dimension(n_mats, n_mats) :: matrix, N
!       integer :: i, j, k
!       real(8), dimension(n_mats), intent(in) :: at

! !Initiate matrix
!       do i = 1, n_mats
!          do j = 1, n_mats
!             matrix(i, j) = 0d0
!             if (i == j) then
!                N(i, j) = 1d0
!             else
!                N(i, j) = 0d0
!             end if

!             if (i == 1) then
!                N(i, j) = 1d0
!             end if

!          end do
!       end do

! !Compute matrix elements
!       do i = 1, n_mats
!          x(i) = 0d0
!          do j = 1, n_mats
!             do k = 1, n_mats
!                matrix(i, j) = matrix(i, j) + N(k, j)*at(i)
!             end do
!             matrix(i, j) = matrix(i, j) - N(i, j)
!          end do
!       end do

!       x(1) = 1d0
!       call solve_linear_system(n=n_mats, l=n_mats, matrix=matrix, b=x, x=y)

!    end function at2mat

! !#######################################################################
!    function at2mat_max(at, n_mats) result(y)
! !Convert atomic mole fractions to molecular mole fractions. The atomic
! !composition of the molecular species is given by the matrix N and is
! !fixed to represent the outer core composition at this point. In future
! !versions it could be passed as an argument to allow for a more diverse
! !compositions. For a composition of the form Fe, FeM1, FeM2, ... the
! !maximum concentration of lighter elements is acheived for Fe = 0. If
! !Fe < 0. the Fe content is set to zero and the other fractions rescaled
! !to satisfy the normalization criteria.

!       integer, intent(in) :: n_mats
!       real(8) :: norm
!       real(8), dimension(n_mats) :: y, x
!       real(8), dimension(n_mats, n_mats) :: matrix, N
!       integer :: i, j, k
!       real(8), dimension(n_mats), intent(in) :: at

! !Initiate matrix
!       do i = 1, n_mats
!          do j = 1, n_mats
!             matrix(i, j) = 0d0
!             if (i == j) then
!                N(i, j) = 1d0
!             else
!                N(i, j) = 0d0
!             end if

!             if (i == 1) then
!                N(i, j) = 1d0
!             end if
!          end do
!       end do

! !Compute matrix elements
!       do i = 1, n_mats
!          x(i) = 0d0
!          do j = 1, n_mats
!             do k = 1, n_mats
!                matrix(i, j) = matrix(i, j) + N(k, j)*at(i)
!             end do
!             matrix(i, j) = matrix(i, j) - N(i, j)
!          end do
!       end do

!       do i = 1, n_mats
!          matrix(1, i) = 1d0
!       end do

!       x(1) = 1d0
!       call solve_linear_system(n=n_mats, l=1, matrix=matrix, b=x, x=y)

! !Check if too much lighter elements are present and rescale molecular
! !mole fractions to satisfy normalization criteria.
!       if (y(1) .lt. 0d0) then
!          y(1) = 0d0
!          norm = sum(y)
!          do i = 1, n_mats
!             y(i) = y(i)/norm
!          end do
!       end if

!    end function at2mat_max

! !#######################################################################
!    function mat2at_core(xi, n_mats) result(y)
! !Takes the material fractions (Fe, FeS, FeSi, FeO) and converts it to
! !atomic fractions xiFe, xiS, xiSi, xiO.

!       integer, intent(in) :: n_mats
!       integer :: i
!       real(8), dimension(n_mats), intent(in) :: xi
!       real(8), dimension(n_mats) :: y, dummy

!       dummy(1) = 1d0 - (sum(xi) - xi(1))
!       do i = 1, n_mats - 1
!          dummy(i + 1) = xi(i + 1)
!       end do

!       do i = 1, n_mats
!          y(i) = dummy(i)/sum(dummy)
!       end do

!    end function mat2at_core

! !#######################################################################
!    function at2wt_core(at, m, n_mats) result(y)
! !Takes atomic mole fractions xiFe, xiS, xiSi, xiO and converts
! !them to weight fractions (Fe, S, Si, O)

!       integer, intent(in) :: n_mats
!       integer :: i
!       real(8), dimension(n_mats), intent(in) :: at
!       real(8), dimension(4) :: m
!       real(8), dimension(n_mats) :: y
!       real(8) :: m_tilde

!       m = (/mFe, mS, mSi, mO/)
!       m_tilde = 0d0

!       do i = 1, n_mats
!          m_tilde = m_tilde + at(i)*m(i)
!       end do
!       print *, 'm_tilde =', m_tilde
!       do i = 1, n_mats
!          y(i) = at(i)*m(i)/m_tilde
!       end do

!    end function at2wt_core

! !#######################################################################
!    function wt2mol_core(wt, xiH) result(y)
! !Give X_S, X_Si, X_O as wt frac and xiH as mole frac to compute the mole
! !fracs of FeHx, FeS, FeSi and FeO in the core.

!       implicit none

!       real(8), intent(in) :: xiH
!       real(8), dimension(3), intent(in) :: wt
!       real(8), dimension(4) :: y
!       real(8), dimension(4, 4) :: matrix
!       integer :: i, j

! !Constuct matrix

!    end function wt2mol_core

! !########################################################################
!    function compute_xi(eta, m1, m2) result(y)
! !Computes mole fractions from mass fraction of species 2 where species 1
! !has molar abundance of (1-xi)

!       implicit none

!       real(kind=8) :: y
!       real(kind=8), intent(in) :: m1, m2, eta

!       y = m1*eta/(eta*m1 - eta*m2 + m2)

!    end function compute_xi

! !########################################################################
!    function compute_eta(xi, m1, m2) result(y)
! !Compute mass fraction from mole fraction of species 2 where species 1
! !has a weight abundance of (1-eta)

!       implicit none

!       real(kind=8) :: y
!       real(kind=8), intent(in) :: xi, m1, m2

!       y = xi*m2/((1.0d0 - xi)*m1 + xi*m2)

!    end function compute_eta

! !#######################################################################
!    function epsilon_ki(T, k, i) result(eps)
! !Compute epsilon parameter for metal-silicate partitioning according to
! !Fischer et al. 2015.
!       real(8), intent(in) :: T
!       integer, intent(in) :: k, i
!       real(8) :: eps

!       eps = eki(k, i)*m_partitioning(i)/0.242d0*1873d0/T
!       eps = eps - m_partitioning(i)/55.85d0 + 1d0

!    end function epsilon_ki

! !#######################################################################
!    function partition_coefficient(T, P, ll) result(KD)
! !Compute metal-silicate partition coefficient of O or Si at T and P from
! !Fischer et al. 2015.

!       implicit none

!       real(8), intent(in) :: T, P
!       integer, intent(in) :: ll
!       real(8) :: KD

! !Note that in Fischer et al. 2015 the pressure is in GPa
!       KD = 10**(a_KD(ll) + b_KD(ll)/T + c_KD(ll)*P*1d-9/T)

!    end function partition_coefficient

! !#######################################################################
!    function partition_coefficient_i(T, P, i, xi) result(KD)
! !Metal-silicate partition coefficient of Si and O according to Fischer
! !et al. 2015. Non-ideal effects for multicomponent fluids are accounted
! !for using the epsilon-formalism.
!       real(8), intent(in) :: T, P
!       integer, intent(in) :: i
!       integer :: k
!       real(8), dimension(2) :: xi_part
!       real(8), dimension(5), intent(in) :: xi
!       real(8) :: KD, delta, eps

!       xi_part(1) = xi(4) !Silicate
!       xi_part(2) = xi(5) !Oxygen
!       KD = a_KD_i(i) + b_KD_i(i)/T + c_KD_i(i)*P*1d-9/T
! !~ print *, 'KD =', KD
!       eps = epsilon_ki(T, i, i)

!       KD = KD + eps*log(1d0 - xi_part(i))/2.303d0
! !~ print *, 'KD =', KD
! !~ print *, 'xi =', xi
! !~ print *, 'eps =', eps
!       do k = 1, 2
!          if (k == i) then
!             cycle
!          else
!             eps = epsilon_ki(T, k, i)
! !~                 print *, 'eps =', eps
!             delta = 1d0/2.303d0*eps*xi_part(k)
! !~                 print *, 'delta =', delta
!             delta = delta*(1d0 + log(1d0 - xi_part(k))/xi_part(k) - &
!                            1d0/(1d0 - xi_part(i)))
! !~                 print *, 'delta =', delta
!             KD = KD + delta
!          end if
!       end do
! !~ print *, 'KD =', KD
! !~ print *, '----'
!       do k = 1, 2
!          if (k == i) then
!             cycle
!          else
!             eps = epsilon_ki(T, k, i)
! !~                 print *, 'eps =', eps
!             delta = 1d0/2.303d0*eps*xi_part(k)**2*xi_part(i)
! !~                 print *, 'delta =', delta
!             delta = delta*(1d0/(1d0 - xi_part(i)) + 1d0/(1d0 - xi_part(k)) + &
!                            xi_part(i)/(2d0*(1d0 - xi_part(i))**2) - 1d0)
! !~                 print *, 'delta =', delta
!             KD = KD - delta
!          end if
!       end do
! !~ print *, 'KD =', KD
!       KD = 10d0**KD

!    end function partition_coefficient_i

! !########################################################################
!    function xi_FeO(T, P, xi) result(xiFeO)
! !Compute mole fraction of FeO in the mantle at given T, P and Fe and O
! !content in the core assuming chemical equ. between core and mantle from
! !Fischer et al. 2015.

!       implicit none

!       real(8), intent(in) :: T, P
!       real(8), dimension(5), intent(in) :: xi
!       real(8) :: xiFeO, KD

!       KD = partition_coefficient_i(T, P, 2, xi)
! !~ print *, 'KD_O =', KD
!       xiFeO = xi(1)*xi(5)/KD

!    end function xi_FeO

! !#######################################################################
!    function xi_SiO2(T, P, xi) result(xiSiO2)
! !Compute mole fraction of siO2 in the mantle at given T, P and Fe and O
! !and Si content in the core assuming chemical equ. between core and mantle
! !from Fischer et al. 2015.

!       implicit none

!       real(8), intent(in) :: T, P
!       real(8), dimension(5), intent(in) :: xi
!       real(8) :: xiSiO2, KD, xiFeO

!       xiFeO = xi_FeO(T, P, xi)
!       KD = partition_coefficient_i(T, P, 1, xi)
! !~ print *, 'xiFeO =', xiFeO
! !~ print *, 'KD_Si =', KD
!       xiSiO2 = xi(4)*xiFeO/(xi(1)*KD)

!       if (isnan(xiSiO2)) then
!          xiSiO2 = 0d0
!       end if

!    end function xi_SiO2

! !#######################################################################
!    function compute_XH2O_melt(T) result(XH2O)
! !Computes H2O content in silicate melt as function of temperature
! !according to Fei 2020 (Fig. 3). The result is in wt fraction.

!       implicit none

!       real(8), intent(in) :: T
!       real(8) :: XH2O

!       XH2O = (10d0 + 250d0*(1d3/T - 0.408d0))/100d0

!    end function compute_XH2O_melt

! !########################################################################
!    function Si_number_mantle(T, P, xi) result(Si_number)
! !Compute Si# in the mantle from Fe, O and Si content in the core and the
! !P and T at the CMB assuming chemical equ. between core and mantle from
! !Fischer et al. 2015.

!       implicit none

!       real(8), intent(in) :: T, P
!       real(8), dimension(5), intent(in) :: xi
!       real(8) :: Si_number, SiMg, xiSiO2, xiFeO, Fe_num

!       xiSiO2 = xi_SiO2(T, P, xi)
!       xiFeO = xi_FeO(T, P, xi)
!       Fe_num = xiFeO/(1d0 - xiSiO2)
!       Fe_num = min(Fe_num, xi_Fe_mantle_max)
!       Fe_num = max(Fe_num, 0d0)
!       SiMg = xiSiO2/(1d0 - xiFeO - xiSiO2)

! !~ print *, 'P, T, xi =', P, T, xi
! !~ print *, 'xiSiO2 =', xiSiO2
! !~ print *, 'xiFeO =', xiFeO
!       Si_number = SiMg/(1d0 + SiMg)
!       Si_number = min(Si_number, 1d0/(2d0 - Fe_num)*0.9999d0)
!       Si_number = max(Si_number, 1d0/(3d0 - 2d0*Fe_num)*1.0001d0)
!    end function Si_number_mantle

! !#######################################################################
!    function Fe_number_mantle(T, P, xi) result(Fe_number)
! !Compute Si# in the mantle from Fe, O and Si content in the core and the
! !P and T at the CMB assuming chemical equ. between core and mantle from
! !Fischer et al. 2015.

!       implicit none

!       real(8), intent(in) :: T, P
!       real(8), dimension(5), intent(in) :: xi
!       real(8) :: Fe_number, SiMg, xiSiO2, xiFeO

!       xiSiO2 = xi_SiO2(T, P, xi)
!       xiFeO = xi_FeO(T, P, xi)
!       Fe_number = xiFeO/(1d0 - xiSiO2)
! !~ print *, 'P, T, xi =', P, T, xi
! !~ print *, 'xiSiO2 =', xiSiO2
! !~ print *, 'xiFeO =', xiFeO
! !~ Fe_number = min(Fe_number, xi_Fe_mantle_max)
! !~ Fe_number = max(Fe_number, 0d0)

!       if (Fe_number < 0d0 .or. Fe_number > xi_Fe_mantle_max) then
!          Fe_number = xi_Fe_mantle_max
!       end if

!    end function Fe_number_mantle

! !########################################################################
!    function Si_number_max(Mg_number, contents) result(y)

!       implicit none

!       real(8), dimension(2) :: y
!       real(8) :: mn, mx
!       integer, intent(in) :: contents(:)
!       real(8), intent(in) :: Mg_number
!       real(8), allocatable :: Si(:), Mg(:)
!       integer :: i, n

!       n = size(contents)

!       allocate (Si(n), Mg(n))

!       do i = 1, n
!          Si(i) = material_YSi(contents(i))
!          Mg(i) = material_YMg(contents(i))*Mg_number
!       end do

! !Compute Si# for given material
!       do i = 1, n
!          Si(i) = Si(i)/Mg(i)/(1d0 + Si(i)/Mg(i))
!       end do

!       y(1) = minval(Si)
!       y(2) = maxval(Si)

! !print *, 'Si# can take values between:', y(:)

!    end function Si_number_max

! !########################################################################
!    function xi_stv_lim(YMg, YSi, xiFe, xiH2O, SiMg, xiAlSi, xiAlStv, &
!                        xiH2OStv, xiAlMg) result(y)
! !Computes the mole fractions xi_stv in the limiting case where only
! !one other phase given by YMg, YSi is present

!       implicit none

!       real(kind=8) :: y, denom
!       real(kind=8), intent(in) :: xiFe, xiH2O, SiMg, xiAlSi, &
!                                   xiAlStv, xiH2OStv, xiAlMg
!       integer, intent(in) :: YSi, YMg

!       y = YSi*(1.0d0 - xiH2O)*(1.0d0 - xiAlSi)
!       y = y - YMg*(1.0d0 - xiH2O)*(1.0d0 - xiAlMg)*(1.0d0 - xiFe)*SiMg

!       denom = YSi*(1.0d0 - xiH2O)*(1.0d0 - xiAlSi)
!       denom = denom - (1.0d0 - xiH2OStv)*(1.0d0 - xiAlStv)
!       denom = denom - YMg*(1.0d0 - xiH2O)*(1.0d0 - xiAlMg)*(1.0d0 - xiFe)*SiMg

!       y = y/denom
!       y = max(y, 0.0d0)

!    end function xi_stv_lim

! !########################################################################
!    function compute_xiFei(Fe_number, P, T, n_mats, contents) result(y)

!       implicit none

!       integer, intent(in) :: n_mats
!       real(8), intent(in) :: Fe_number, P, T
!       real(8), dimension(n_mats) :: y
!       real(8), dimension(n_mats - 1) :: DFe
!       integer, dimension(n_mats), intent(in) :: contents
!       integer :: i

! !Compute iron partitioning coefficient
! !Here the fit results to experimental data will be called

! !Compute partitioning coefficient of iron between substance i and i+1
! !Convention is: DFe = xi_Fe_(i+1)/xi_Fe_i
!       if (n_mats > 1) then
!          do i = 1, n_mats - 1
!             DFe(i) = 1.0d0

!             !If 3 materials are present, it is assumed that the third one is Stv
!             !which contains no iron. In general it can also be another phase
!             !for which a different value of DFe could be set.
!             if (i == 3) then
!                DFe(i) = 0.0d0
!             end if

!          end do

!          y(1) = n_mats*Fe_number/(DFe(1) + 1.0d0)

!          do i = 1, n_mats - 1
!             y(i + 1) = n_mats*Fe_number*DFe(i)/(DFe(i) + 1.0d0)
!          end do

!       else

!          y(1) = Fe_number

!       end if

!    end function compute_xiFei

! !########################################################################
!    function xiAli(xiFe, P, T, n_mats) result(y)

!       implicit none

!       integer, intent(in) :: n_mats
!       real(kind=8), intent(in) :: P, T
!       real(kind=8), dimension(n_mats), intent(in) :: xiFe
!       real(kind=8), dimension(n_mats) :: y

! !Compute Al content in each phase at given P, T and for the Fe content
! !in the phase that has been computed at P, T via the function xiFei
!       y(1) = 0.0d0
!       y(2) = 0.0d0
!       y(3) = 0.0d0

!    end function xiAli

! !########################################################################
!    function Al_partitioning_SiMg(xiAl, P, T, n_mats) result(y)

!       implicit none

!       real(kind=8), intent(in) :: P, T
!       integer, intent(in) :: n_mats
!       real(kind=8), dimension(n_mats) :: xiAl
!       real(kind=8), dimension(n_mats, 2) :: y
!       real(kind=8), dimension(n_mats) :: DAl_SiMg
!       integer :: i

!       do i = 1, n_mats
!          if (i < 3) then
!             DAl_SiMg(i) = 0.5d0

!          else
!             DAl_SiMg(i) = 1.0d10
!          end if
!       end do

!       do i = 1, n_mats
!          y(i, 1) = xiAl(i)/(1.0d0 + DAl_SiMg(i))
!          y(i, 2) = xiAl(i)*(1.0d0 - 1.0d0/(1.0d0 + DAl_SiMg(i)))
!       end do

!    end function Al_partitioning_SiMg

! !#######################################################################
!    function compute_xiH2O(xiFei, xiAli, P, T, n_mats, contents) result(y)

!       implicit none

!       integer, intent(in) :: n_mats
!       integer, dimension(n_mats), intent(in) :: contents
!       real(kind=8), dimension(n_mats), intent(in) :: xiFei, xiAli
!       real(kind=8), intent(in) :: P, T
!       real(kind=8), dimension(n_mats) :: y
!       integer :: i

!       do i = 1, n_mats

!          !Here, call a function to compute the water content as function of
!          !P, T, Fe, Al for each phase

!          y(i) = 0.0d0
!       end do

!    end function compute_xiH2O

! !#######################################################################
!    SUBROUTINE update_all_fractions(X_new, weight_fractions, fractions, &
!                                    molar_masses, n_mats)
! !Uses the newly computed weight fraction of the impurity and updates all
! !other weight fractions and all mole fractions. new_last_frac is the
! !weight fraction of the impurity (last material)

!       real(8), intent(in) :: X_new
!       real(8) :: xi_new
!       integer, intent(in) :: n_mats
!       real(8), intent(inout), dimension(n_mats) :: weight_fractions, fractions
!       real(8), intent(in), dimension(n_mats) :: molar_masses
!       real(8), dimension(n_mats) :: olds
!       real(8) :: sum_others, m_others, m_impurity, delta
!       integer :: i

! !First update the weight fractions
!       olds = weight_fractions
!       weight_fractions(n_mats) = X_new

! !Compute normalization factor
!       sum_others = 0d0
!       do i = 1, n_mats - 1
!          sum_others = sum_others + olds(i)
!       end do

! !Compute weight fractions
!       do i = 1, n_mats - 1
!          weight_fractions(i) = olds(i)*(1d0 - X_new)/sum_others
!       end do

! !Compute the normalized mass of the non-impurities for the conversion
! !to mole fraction
!       m_others = 0d0

!       do i = 1, 2
! !The layer itself does not contain the right fractions so they have to
! !be taken from the shells. Since the first two material fractions will
! !be constant we can just take the first shell.
!          delta = molar_masses(i)
!          m_others = m_others + delta*fractions(i)
!       end do

! !Normalized mass of the non-impurities
!       m_others = m_others/(fractions(1) + fractions(2))

! !Mass of the impurity (last material)
!       m_impurity = molar_masses(n_mats)

! !Compute the mole fraction of the impurity
!       xi_new = compute_xi(eta=X_new, m1=m_others, m2=m_impurity)

! !Use this to compute all the other mole fractions
!       olds = fractions
!       fractions(n_mats) = xi_new

! !Compute normalization factors
!       sum_others = 0d0
!       do i = 1, n_mats - 1
!          sum_others = sum_others + olds(i)
!       end do

! !Compute mole fractions
!       do i = 1, n_mats - 1
!          fractions(i) = olds(i)*(1d0 - xi_new)/sum_others
!       end do

!    END SUBROUTINE update_all_fractions


!#######################################################################
!    SUBROUTINE compute_P_H2(self, average)
! !Compute hydrogen partial pressure at the CMB according to Wu et al. 2018
! !Note that in the paper it is stated, that this proceedure is only valid
! !for P_H2 < 10MPa. However, Roskosz et al. 2013 have compared the Sievert's
! !law to experimental data for nitrogen in Fe up to pressures of 20 GPa
! !at temperatures of 2000-3000 K. They find very good agreement between
! !the experiments and the predictions from Sivert's law. They find Nitrogen
! !contents of up to ~14 wt% or ~40 mol%. At 3000 K and 20 GPa the density
! !of iron is ~8.6 gcc. Using the ideal gas law to estimate the order of magnitude
! !of the partial pressure P_N at these conditions yields  up to ~ 10³MPa.
! !This implies that the Sievert's law extrapoaltes well to higher partial
! !gas pressures. It is not clear if this is also true for hydrogen in iron.
! !But given the simplicity of both system we assume that the hydrogen
! !solubility in Fe extrapolates well to higher pressures just as for nitrogen.

!       type(planet), intent(inout) :: self
!       logical, intent(in), optional :: average
!       logical :: average_dummy
!       real(8) :: V_mantle, V_innermost_shell, &
!                  N_H2O_innermost_shell, N_H2O_mantle, T_CMB, R_CMB, rho_H2, test
!       integer :: i, j, lay
!       print *, 'computing pH2 in layer ', self%lay
!       if (.not. (present(average))) then
!          average_dummy = .false.
!       else
!          average_dummy = average
!       end if
! !If mantle exists, compute H2 partial pressure in the lowermost shell
! !of the mantle or the average H2 partial pressure in the entire mantle
!       if (.not. size(self%layers(2)%shells) .eq. 0) then
!          T_CMB = self%layers(2)%temp
!          R_CMB = self%layers(2)%radius
!          print *, 'Outer core exists'
!       else
!          print *, 'No outer core exists'
!          T_CMB = self%layers(1)%temp
!          R_CMB = self%layers(1)%radius
!       end if

! !Only core exists
!       if (self%lay <= 2) then
!          self%P_H2 = 0d0
!          print *, 'No mantle exists'
!          self%mantle_exists = .false.

!       else
!          ! No lower mantle
!          if (self%layers(3)%shells(1)%volume .eq. 0d0) then
!             print *, 'Setting lay to 4'
!             lay = 4
!          else
!             lay = 3
!          end if
!          V_innermost_shell = 4d0/3d0*PI*(self%layers(lay)%shells(1)%radius**3 - &
!                                          R_CMB**3)
!          N_H2O_innermost_shell = self%layers(lay)%shells(1)%N_H2O
!          N_H2O_innermost_shell = self%layers(lay)%shells(1)%indigenous_mass
!          !N_H2O_innermost_shell = N_H2O_innermost_shell * 2d-1/(mH*2d0+mO)
!          V_mantle = 4d0/3d0*PI*(self%layers(4)%radius**3 - R_CMB**3)
!          N_H2O_mantle = 0d0

!          !Comptue the average water content in the mantle
!          if (average_dummy) then
!             do i = 3, 4
!                do j = 1, self%layers(i)%shell_count
!                   N_H2O_mantle = N_H2O_mantle + self%layers(i)%shells(j)%N_H2O
!                end do
!             end do

!             self%P_H2 = N_H2O_mantle*NA*kB*T_CMB/V_mantle
!             rho_H2 = N_H2O_mantle*mH*2d0/V_mantle

!             !Compute the partial pressure in the lowermost shell
!          else
!             self%P_H2 = N_H2O_innermost_shell*NA*kB*T_CMB/ &
!                         V_innermost_shell
!             rho_H2 = N_H2O_innermost_shell*mH*2d0/V_innermost_shell
!          end if
! !~         print *, 'rho_H2 =', rho_H2
!          !Use van der Waal equation to estimate partial H2 pressure. Note that
!          !using the IGL can yield deviations in the resulting hydrogen content
!          !in the core of up to ~50% for 0-10 earth mass planets.
!          self%P_H2 = rho_H2*Rgas*T_CMB/(2d0*mH - rho_H2*26.61d-6) - rho_H2**2*24.76d-3/(2d0*mH)**2
! !~         print *, 'P_H2 =', self%P_H2
! !~         print *, 'T_CMB =', T_CMB
!          !This is a stupid way to prevent P_H2 < 0 in the case where eps_H2O = 0
!          !and large X_impurity. Since impurities are currently only intruduced
!          !for the anhydrous case this is kinda enough at this point. (But it's
!          !still stupid!)
!          self%P_H2 = self%P_H2*self%eps_H2O

!       end if

!       test = rho_H2*Rgas*T_CMB/(2d0*mH - rho_H2*26.61d-6) - rho_H2**2*24.76d-3/(2d0*mH)**2
! !print *, 'N_H2O_innermost_shell =', N_H2O_innermost_shell
! !print *, 'V_innermost_shell =', V_innermost_shell
! !print *, 'P_H2 =', self%P_H2
! !print *, 'test P_H2 =', test
! !print *, 'deviation =', (sqrt(test)-sqrt(self%P_H2))/sqrt(self%P_H2)
!    END SUBROUTINE compute_P_H2

!#######################################################################
!    SUBROUTINE compute_xi_H_core(self)
! !Compute molar hydrogen abundance in the core according to Wu et al. 2018

!       type(planet), intent(inout) :: self
!       real(8) :: T_CMB

!       T_CMB = self%layers(2)%temp

!       call compute_P_H2(self=self)

! !Compute parameter x in Wu et al 2018: (x/2+1) Fe + x/2 H2O
! !P0 is 1 bar
!       self%xi_H_core = sqrt(self%P_H2/1d5)*exp((-31.8d3 - T_CMB*38.1d0)/(Rgas*T_CMB))

! !self%xi_H_core = 2d-1

! !Here the actual molar abundance xi_H of H in the core is computed
! !With this the core composition is (1-xi_H) Fe + xi_H H
! !It is defined as H/(H+Fe)
!       self%xi_H_core = self%xi_H_core/(self%xi_H_core + 1d0)!*self%eps_H2O

!       self%xi_H_core = 0d0
!       print *, 'remember: xi_H_core is set to 0 for Venus!'
!    END SUBROUTINE compute_xi_H_core


!########################################################################
   subroutine construct_abundance_matrix(SiMg, n_mats, &
                                         YMgi, YSii, matrix, b, xiFei, additional)

!Constructs linear system of equations to solve for the individual
!molar abundances for fixed atomic abundances and materials in a
!mixture.

      implicit none

      integer, intent(in) :: n_mats
      integer, dimension(n_mats), intent(in) :: YMgi, YSii
      real(kind=8), dimension(n_mats), intent(in) :: xiFei
      real(8), intent(in) :: additional(:)
      real(kind=8), intent(in) :: SiMg
      real(kind=8), dimension(n_mats, n_mats), intent(out) :: matrix
      real(8), intent(out) :: b(:, :)
      integer :: i, j

      do j = 1, n_mats
         ! matrix(1, j) = (1.0d0 - xiH2Oi(j))*(SiMg*(1.0d0 - xiAlMgi(j))* &
         !                                     (1.0d0 - xiFei(j))*YMgi(j) - (1.0d0 - xiAlSii(j))*YSii(j))
         matrix(1, j) = SiMg * (1.0d0 - xiFei(j))*YMgi(j) - YSii(j)
         matrix(2, j) = 1.0d0

      end do

!With the default compositional parameters, i.e. Fe/Mg and Si/Mg only
!two components can be uniquely constrained. If more than two components,
!i.e. n_mats>2, are given, then the molar abundances of the additional
!components must be fixed beforehand. They will not influence the ratio
!of the first two abundances but the absolut values as the sum of all
!abundances must be equal to one.

      if (n_mats > 2) then
         do i = 3, n_mats
            do j = 1, n_mats
               matrix(i, j) = 0d0

            end do
            b(i, 1) = additional(i - 2)
            matrix(i, i) = 1d0
         end do
      end if

      b(1, 1) = 0d0
      b(2, 1) = 1d0

   end subroutine construct_abundance_matrix

!########################################################################
   subroutine compute_abundance_vector(SiMg, n_mats, &
                                       YMgi, YSii, xiFei, abundances, &
                                       contents, additional)

!Here the molar abundances of the individual components of a mixture
!are computed using atomic abundances and the different materials as
!inputs. The 'additional' parameter contains additional components of
!the mixture which are not constrained by the main inputs (i.e. Fe and
!Si contents). It is an optional argument and can contain as many
!additional components as desired. These components will not be taken
!into account in the computation of the abundances even if they contain
!Fe, Si or Mg!

      implicit none

      integer, intent(in) :: n_mats
      integer, dimension(n_mats), intent(in) :: YMgi, YSii, contents
      real(kind=8), dimension(n_mats), intent(in) :: xiFei
      real(kind=8), intent(in) :: SiMg
      real(kind=8), dimension(n_mats, n_mats) :: matrix
      real(kind=8), dimension(n_mats), intent(out) :: abundances
      real(kind=8), dimension(n_mats, 1) :: abundances_dummy
      real(8), dimension(n_mats, 1) :: b_vec
      real(8), intent(in), optional :: additional(:)
      real(8), allocatable :: additional_dummy(:)
      integer :: i
     
      if (n_mats > 2) then
         allocate (additional_dummy(n_mats - 2))
         if (present(additional)) then
            additional_dummy(:) = additional

         else
            additional_dummy(:) = 0d0
         end if

      else
         allocate (additional_dummy(n_mats))
         additional_dummy(1:) = 0d0

      end if

      do i = 1, n_mats
         abundances_dummy(i, 1) = 0.0d0
      end do

!In the core the first material is pure iron. There no compositional
!gradients are currently allowed and the fractions need not be computed
      if (contents(1) == 2 .or. contents(1) == 1) then
!pass

      elseif (size(contents) == 1) then
         abundances(1) = 1d0

      else
         call construct_abundance_matrix(SiMg=SiMg, n_mats=n_mats, &
                                         YMgi=YMgi, YSii=YSii, xiFei=xiFei, &
                                         matrix=matrix, b=b_vec, &
                                         additional=additional_dummy)

         call gauss_elimination(a=matrix, sol=abundances_dummy, b=b_vec)

         do i = 1, n_mats
            abundances(i) = abundances_dummy(i, 1)
         end do
      end if

   end subroutine compute_abundance_vector


!###################################################################################
   function compute_core_mass(M_tot, M_ocean, FeMg, SiMg, Fe_numbers, xi_all_core, contents,&
       inner_core_mass_fraction, inner_core_mass, additional, mode) result(core_mass)

      real(8), intent(in), optional :: inner_core_mass_fraction, inner_core_mass
      real(8), intent(in) :: M_tot, M_ocean, FeMg, SiMg
      real(8) :: core_mass
      integer :: n_mats, lay
      real(8), dimension(5) :: Q
      real(8), dimension(5), intent(in) :: xi_all_core
      real(8), intent(in) :: Fe_numbers(:)
      real(8) :: Fe_number_mantle
      integer, intent(in) :: contents(:)
      real(8), intent(in), optional :: additional(:)
      real(8) :: core_fraction
      integer, intent(in) :: mode
      real(8), allocatable :: fractions(:)
      integer, allocatable :: ymg(:), ysi(:)
      integer :: i

      ! By default, the lower mantle will be used
      ! TODO. Change to input argument
      lay = 3
      n_mats = size(contents)

      ! Allocate the arrays
      allocate(fractions(n_mats))
      allocate(ymg(n_mats))
      allocate(ysi(n_mats))
      
      ! Get the stoichometric for the relevant materials
      do i=1, n_mats
         ymg(i) = material_YMg(contents(i))
         ysi(i) = material_YSi(contents(i))
      enddo

      ! Compute the fractions in the mantle
      call compute_abundance_vector(SiMg=SiMg,&
       n_mats=n_mats, &
       YMgi=ymg, &
       YSii=ysi, &
       xiFei=Fe_numbers, &
       abundances=fractions, &
       contents=contents,&
       additional=additional)

      ! Fe number is assumed homogeneous in mantle so it does not matter which one we take here
      Fe_number_mantle = Fe_numbers(1)

      ! Compute the coefficients for the core mass computation based on the composition
      Q = get_core_mass_q_vector(Fe_number_mantle, xi_all_core, fractions, contents)

      ! Compute core mass from the absolute inner core mass
      if (mode == 1 .and. present(inner_core_mass)) then
         core_fraction = (1d0 - M_ocean / M_tot)
         core_fraction = core_fraction * (Q(3) / Q(2) - Q(1) / Q(2) * FeMg)
         core_fraction = core_fraction + inner_core_mass / M_tot * (1e0 / mFe - Q(4) / Q(5))
         core_fraction = core_fraction / (Q(3) / Q(2) - Q(4) / Q(5) - FeMg * Q(1) / Q(2))
         
      ! Compute core mass from the inner core mass fraction
      elseif (mode == 2 .and. present(inner_core_mass_fraction)) then
         core_fraction = (1d0 - M_ocean / M_tot)
         core_fraction = core_fraction * (Q(3) / Q(2) - Q(1) / Q(2) * FeMg)
         core_fraction = core_fraction / ((Q(3) / Q(2) - Q(4) / Q(5) - FeMg * Q(1) / Q(2)) &
          + inner_core_mass_fraction * (Q(4) / Q(5) - 1e0 / mFe))
      endif
      core_mass = M_tot * core_fraction
   end function compute_core_mass

! ##################################################################################
   function get_core_mass_q_vector(Fe_number, xi_all_core, fractions, contents) result(Q)
      real(8), dimension(5) :: Q
      real(8), dimension(5) :: masses = (/mFe, mH, mS, mSi, mO/)
      integer :: i, mat
      real(8) :: frac
      real(8), intent(in) :: Fe_number
      real(8), intent(in) :: fractions(:), xi_all_core(:)
      integer, intent(in) :: contents(:)
      ! Compute mole fraction of Mg in the mantle
      Q(1) = 0d0
      do i=1, size(fractions)
         mat = contents(i)
         frac = fractions(i)
         Q(1) = Q(1) + frac * (1d0 - Fe_number) * material_YMg(mat)
      enddo
      ! Compute total normalized mass in the mantle
      Q(2) = 0d0

      ! Contribution of Mg
      do i=1, size(fractions)
         mat = contents(i)
         frac = fractions(i)
         Q(2) = Q(2) + frac * (1d0 - Fe_number) * material_YMg(mat) * mMg
      enddo
      ! Contribution of Fe
      do i=1, size(fractions)
         mat = contents(i)
         frac = fractions(i)
         Q(2) = Q(2) + frac * Fe_number * material_YMg(mat) * mFe
      enddo
      ! Contribution of Si
      do i=1, size(fractions)
         mat = contents(i)
         frac = fractions(i)
         Q(2) = Q(2) + frac * material_YSi(mat) * mSi
      enddo
      ! Contribution of O
      do i=1, size(fractions)
         mat = contents(i)
         frac = fractions(i)
         Q(2) = Q(2) + frac * material_YO(mat) * mO
      enddo
      ! Contribution of H
      do i=1, size(fractions)
         mat = contents(i)
         frac = fractions(i)
         Q(2) = Q(2) + frac * material_YH(mat) * mH
      enddo
      ! Compute mole fraction of Fe in the mantle
      Q(3) = 0d0
      do i=1, size(fractions)
         mat = contents(i)
         frac = fractions(i)
         Q(3) = Q(3) + frac * Fe_number * material_YMg(mat)
      enddo
      ! Compute total fraction of Fe in the core
      Q(4) = xi_all_core(1)

      ! Compute total normalized mass in the core
      Q(5) = 0d0
      do i=1, size(xi_all_core)
         Q(5) = Q(5) + masses(i) * xi_all_core(i)
      enddo

   end function get_core_mass_q_vector


END MODULE functionsPy
