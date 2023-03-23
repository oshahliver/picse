MODULE functions

   use constants
   use run_params
   use eosfort
   use linalg

   implicit none

contains

!########################################################################
   subroutine progress_bar(j, n)

      integer :: j, k, perc, dec1, dec2
      real(8) :: jj
      integer, intent(in) :: n
      character(len=72) :: bar = "???.??% |                                                   |"

      jj = real(j)/real(n)
      perc = int(jj*100)
      dec1 = int(((jj*1d2) - real(perc))*1d1)
      dec2 = int(((jj*1d2) - real(perc))*1d2)
      dec2 = dec2 - 10*dec1

      write (unit=bar(1:3), fmt="(i3)") perc
      write (unit=bar(5:5), fmt="(i1)") dec1
      write (unit=bar(6:6), fmt="(i1)") dec2

      do k = 1, perc/2
         bar(9 + k:9 + k + 1) = "=>"
      end do

      write (unit=6, fmt="(a1, a72)", advance="no") char(13), bar

      if (jj /= 100d0) then
         flush (unit=6)

      else

         write (unit=6, fmt=*)

      end if
      return

   end subroutine progress_bar

!########################################################################
   function MOI_integrand_linear(r1, r2, rho1, rho2) result(moi)

      implicit none

      real(8), intent(in) :: r1, r2, rho1, rho2
      real(8) :: moi
      real(8) :: dr, a, moi1, moi2

      dr = r2 - r1
      a = (rho2 - rho1)/dr

      moi1 = 1d0/5d0*rho1*r1**5 + 1d0/6d0*a*r1**6 - 1d0/5d0*a*r1*r1**5
      moi2 = 1d0/5d0*rho1*r2**5 + 1d0/6d0*a*r2**6 - 1d0/5d0*a*r1*r2**5

      moi = moi2 - moi1

   end function MOI_integrand_linear

!########################################################################
   function eta_general(xi, m, n_mats) result(y)
!Compute weight fractions from mole fractions of different components.

      implicit none

      integer :: i, n_mats
      real(8), dimension(n_mats) :: y
      real(8) :: m_tilde
      real(8), intent(in), dimension(n_mats) :: xi, m

      y = 0d0
      m_tilde = 0d0

      do i = 1, n_mats
         m_tilde = m_tilde + m(i)*xi(i)
      end do

      do i = 1, n_mats
         y(i) = xi(i)*m(i)/m_tilde
      end do

   end function eta_general

!########################################################################
   function xi_general(eta, m, n_mats) result(y)
!Compute mole fractions from weight fractions of different components.

      implicit none

      integer :: i, n_mats
      real(8), dimension(n_mats) :: y
      real(8) :: dummy
      real(8), intent(in), dimension(n_mats) :: eta, m

      print *, 'eta =', eta

      dummy = 0d0
      do i = 1, n_mats
         dummy = dummy + eta(i)/m(i)
      end do

      do i = 1, n_mats
         y(i) = eta(i)/m(i)/dummy
      end do

   end function xi_general

!#######################################################################
   function at2mat(at, n_mats) result(y)
!Convert atomic mole fractions to molecular mole fractions. The atomic
!composition of the molecular species is given by the matrix N and is
!fixed to represent the outer core composition at this point. In future
!versions it could be passed as an argument to allow for a more diverse
!compositions.

      integer, intent(in) :: n_mats
      real(8), dimension(n_mats) :: y, x
      real(8), dimension(n_mats, n_mats) :: matrix, N
      integer :: i, j, k
      real(8), dimension(n_mats), intent(in) :: at

!Initiate matrix
      do i = 1, n_mats
         do j = 1, n_mats
            matrix(i, j) = 0d0
            if (i == j) then
               N(i, j) = 1d0
            else
               N(i, j) = 0d0
            end if

            if (i == 1) then
               N(i, j) = 1d0
            end if

         end do
      end do

!Compute matrix elements
      do i = 1, n_mats
         x(i) = 0d0
         do j = 1, n_mats
            do k = 1, n_mats
               matrix(i, j) = matrix(i, j) + N(k, j)*at(i)
            end do
            matrix(i, j) = matrix(i, j) - N(i, j)
         end do
      end do

      x(1) = 1d0
      call solve_linear_system(n=n_mats, l=n_mats, matrix=matrix, b=x, x=y)

   end function at2mat

!#######################################################################
   function at2mat_max(at, n_mats) result(y)
!Convert atomic mole fractions to molecular mole fractions. The atomic
!composition of the molecular species is given by the matrix N and is
!fixed to represent the outer core composition at this point. In future
!versions it could be passed as an argument to allow for a more diverse
!compositions. For a composition of the form Fe, FeM1, FeM2, ... the
!maximum concentration of lighter elements is acheived for Fe = 0. If
!Fe < 0. the Fe content is set to zero and the other fractions rescaled
!to satisfy the normalization criteria.

      integer, intent(in) :: n_mats
      real(8) :: norm
      real(8), dimension(n_mats) :: y, x
      real(8), dimension(n_mats, n_mats) :: matrix, N
      integer :: i, j, k
      real(8), dimension(n_mats), intent(in) :: at

!Initiate matrix
      do i = 1, n_mats
         do j = 1, n_mats
            matrix(i, j) = 0d0
            if (i == j) then
               N(i, j) = 1d0
            else
               N(i, j) = 0d0
            end if

            if (i == 1) then
               N(i, j) = 1d0
            end if
         end do
      end do

!Compute matrix elements
      do i = 1, n_mats
         x(i) = 0d0
         do j = 1, n_mats
            do k = 1, n_mats
               matrix(i, j) = matrix(i, j) + N(k, j)*at(i)
            end do
            matrix(i, j) = matrix(i, j) - N(i, j)
         end do
      end do

      do i = 1, n_mats
         matrix(1, i) = 1d0
      end do

      x(1) = 1d0
      call solve_linear_system(n=n_mats, l=1, matrix=matrix, b=x, x=y)

!Check if too much lighter elements are present and rescale molecular
!mole fractions to satisfy normalization criteria.
      if (y(1) .lt. 0d0) then
         y(1) = 0d0
         norm = sum(y)
         do i = 1, n_mats
            y(i) = y(i)/norm
         end do
      end if

   end function at2mat_max

!#######################################################################
   function mat2at_core(xi, n_mats) result(y)
!Takes the material fractions (Fe, FeS, FeSi, FeO) and converts it to
!atomic fractions xiFe, xiS, xiSi, xiO.

      integer, intent(in) :: n_mats
      integer :: i
      real(8), dimension(n_mats), intent(in) :: xi
      real(8), dimension(n_mats) :: y, dummy

      dummy(1) = 1d0 - (sum(xi) - xi(1))
      do i = 1, n_mats - 1
         dummy(i + 1) = xi(i + 1)
      end do

      do i = 1, n_mats
         y(i) = dummy(i)/sum(dummy)
      end do

   end function mat2at_core

!#######################################################################
   function at2wt_core(at, m, n_mats) result(y)
!Takes atomic mole fractions xiFe, xiS, xiSi, xiO and converts
!them to weight fractions (Fe, S, Si, O)

      integer, intent(in) :: n_mats
      integer :: i
      real(8), dimension(n_mats), intent(in) :: at
      real(8), dimension(4) :: m
      real(8), dimension(n_mats) :: y
      real(8) :: m_tilde

      m = (/mFe, mS, mSi, mO/)
      m_tilde = 0d0

      do i = 1, n_mats
         m_tilde = m_tilde + at(i)*m(i)
      end do
      print *, 'm_tilde =', m_tilde
      do i = 1, n_mats
         y(i) = at(i)*m(i)/m_tilde
      end do

   end function at2wt_core

!#######################################################################
   function wt2mol_core(wt, xiH) result(y)
!Give X_S, X_Si, X_O as wt frac and xiH as mole frac to compute the mole
!fracs of FeHx, FeS, FeSi and FeO in the core.

      implicit none

      real(8), intent(in) :: xiH
      real(8), dimension(3), intent(in) :: wt
      real(8), dimension(4) :: y
      real(8), dimension(4, 4) :: matrix
      integer :: i, j

!Constuct matrix

   end function wt2mol_core

!########################################################################
   function compute_xi(eta, m1, m2) result(y)
!Computes mole fractions from mass fraction of species 2 where species 1
!has molar abundance of (1-xi)

      implicit none

      real(kind=8) :: y
      real(kind=8), intent(in) :: m1, m2, eta

      y = m1*eta/(eta*m1 - eta*m2 + m2)

   end function compute_xi

!########################################################################
   function compute_eta(xi, m1, m2) result(y)
!Compute mass fraction from mole fraction of species 2 where species 1
!has a weight abundance of (1-eta)

      implicit none

      real(kind=8) :: y
      real(kind=8), intent(in) :: xi, m1, m2

      y = xi*m2/((1.0d0 - xi)*m1 + xi*m2)

   end function compute_eta

!#######################################################################
   function epsilon_ki(T, k, i) result(eps)
!Compute epsilon parameter for metal-silicate partitioning according to
!Fischer et al. 2015.
      real(8), intent(in) :: T
      integer, intent(in) :: k, i
      real(8) :: eps

      eps = eki(k, i)*m_partitioning(i)/0.242d0*1873d0/T
      eps = eps - m_partitioning(i)/55.85d0 + 1d0

   end function epsilon_ki

!#######################################################################
   function partition_coefficient(T, P, ll) result(KD)
!Compute metal-silicate partition coefficient of O or Si at T and P from
!Fischer et al. 2015.

      implicit none

      real(8), intent(in) :: T, P
      integer, intent(in) :: ll
      real(8) :: KD

!Note that in Fischer et al. 2015 the pressure is in GPa
      KD = 10**(a_KD(ll) + b_KD(ll)/T + c_KD(ll)*P*1d-9/T)

   end function partition_coefficient

!#######################################################################
   function partition_coefficient_i(T, P, i, xi) result(KD)
!Metal-silicate partition coefficient of Si and O according to Fischer
!et al. 2015. Non-ideal effects for multicomponent fluids are accounted
!for using the epsilon-formalism.
      real(8), intent(in) :: T, P
      integer, intent(in) :: i
      integer :: k
      real(8), dimension(2) :: xi_part
      real(8), dimension(5), intent(in) :: xi
      real(8) :: KD, delta, eps

      xi_part(1) = xi(4) !Silicate
      xi_part(2) = xi(5) !Oxygen
      KD = a_KD_i(i) + b_KD_i(i)/T + c_KD_i(i)*P*1d-9/T
!~ print *, 'KD =', KD
      eps = epsilon_ki(T, i, i)

      KD = KD + eps*log(1d0 - xi_part(i))/2.303d0
!~ print *, 'KD =', KD
!~ print *, 'xi =', xi
!~ print *, 'eps =', eps
      do k = 1, 2
         if (k == i) then
            cycle
         else
            eps = epsilon_ki(T, k, i)
!~                 print *, 'eps =', eps
            delta = 1d0/2.303d0*eps*xi_part(k)
!~                 print *, 'delta =', delta
            delta = delta*(1d0 + log(1d0 - xi_part(k))/xi_part(k) - &
                           1d0/(1d0 - xi_part(i)))
!~                 print *, 'delta =', delta
            KD = KD + delta
         end if
      end do
!~ print *, 'KD =', KD
!~ print *, '----'
      do k = 1, 2
         if (k == i) then
            cycle
         else
            eps = epsilon_ki(T, k, i)
!~                 print *, 'eps =', eps
            delta = 1d0/2.303d0*eps*xi_part(k)**2*xi_part(i)
!~                 print *, 'delta =', delta
            delta = delta*(1d0/(1d0 - xi_part(i)) + 1d0/(1d0 - xi_part(k)) + &
                           xi_part(i)/(2d0*(1d0 - xi_part(i))**2) - 1d0)
!~                 print *, 'delta =', delta
            KD = KD - delta
         end if
      end do
!~ print *, 'KD =', KD
      KD = 10d0**KD

   end function partition_coefficient_i

!########################################################################
   function xi_FeO(T, P, xi) result(xiFeO)
!Compute mole fraction of FeO in the mantle at given T, P and Fe and O
!content in the core assuming chemical equ. between core and mantle from
!Fischer et al. 2015.

      implicit none

      real(8), intent(in) :: T, P
      real(8), dimension(5), intent(in) :: xi
      real(8) :: xiFeO, KD

      KD = partition_coefficient_i(T, P, 2, xi)
!~ print *, 'KD_O =', KD
      xiFeO = xi(1)*xi(5)/KD

   end function xi_FeO

!#######################################################################
   function xi_SiO2(T, P, xi) result(xiSiO2)
!Compute mole fraction of siO2 in the mantle at given T, P and Fe and O
!and Si content in the core assuming chemical equ. between core and mantle
!from Fischer et al. 2015.

      implicit none

      real(8), intent(in) :: T, P
      real(8), dimension(5), intent(in) :: xi
      real(8) :: xiSiO2, KD, xiFeO

      xiFeO = xi_FeO(T, P, xi)
      KD = partition_coefficient_i(T, P, 1, xi)
!~ print *, 'xiFeO =', xiFeO
!~ print *, 'KD_Si =', KD
      xiSiO2 = xi(4)*xiFeO/(xi(1)*KD)

      if (isnan(xiSiO2)) then
         xiSiO2 = 0d0
      end if

   end function xi_SiO2

!#######################################################################
   function compute_XH2O_melt(T) result(XH2O)
!Computes H2O content in silicate melt as function of temperature
!according to Fei 2020 (Fig. 3). The result is in wt fraction.

      implicit none

      real(8), intent(in) :: T
      real(8) :: XH2O

      XH2O = (10d0 + 250d0*(1d3/T - 0.408d0))/100d0

   end function compute_XH2O_melt

!########################################################################
   function Si_number_mantle(T, P, xi) result(Si_number)
!Compute Si# in the mantle from Fe, O and Si content in the core and the
!P and T at the CMB assuming chemical equ. between core and mantle from
!Fischer et al. 2015.

      implicit none

      real(8), intent(in) :: T, P
      real(8), dimension(5), intent(in) :: xi
      real(8) :: Si_number, SiMg, xiSiO2, xiFeO, Fe_num

      xiSiO2 = xi_SiO2(T, P, xi)
      xiFeO = xi_FeO(T, P, xi)
      Fe_num = xiFeO/(1d0 - xiSiO2)
      Fe_num = min(Fe_num, xi_Fe_mantle_max)
      Fe_num = max(Fe_num, 0d0)
      SiMg = xiSiO2/(1d0 - xiFeO - xiSiO2)

!~ print *, 'P, T, xi =', P, T, xi
!~ print *, 'xiSiO2 =', xiSiO2
!~ print *, 'xiFeO =', xiFeO
      Si_number = SiMg/(1d0 + SiMg)
      Si_number = min(Si_number, 1d0/(2d0 - Fe_num)*0.9999d0)
      Si_number = max(Si_number, 1d0/(3d0 - 2d0*Fe_num)*1.0001d0)
   end function Si_number_mantle

!#######################################################################
   function Fe_number_mantle(T, P, xi) result(Fe_number)
!Compute Si# in the mantle from Fe, O and Si content in the core and the
!P and T at the CMB assuming chemical equ. between core and mantle from
!Fischer et al. 2015.

      implicit none

      real(8), intent(in) :: T, P
      real(8), dimension(5), intent(in) :: xi
      real(8) :: Fe_number, SiMg, xiSiO2, xiFeO

      xiSiO2 = xi_SiO2(T, P, xi)
      xiFeO = xi_FeO(T, P, xi)
      Fe_number = xiFeO/(1d0 - xiSiO2)
!~ print *, 'P, T, xi =', P, T, xi
!~ print *, 'xiSiO2 =', xiSiO2
!~ print *, 'xiFeO =', xiFeO
!~ Fe_number = min(Fe_number, xi_Fe_mantle_max)
!~ Fe_number = max(Fe_number, 0d0)

      if (Fe_number < 0d0 .or. Fe_number > xi_Fe_mantle_max) then
         Fe_number = xi_Fe_mantle_max
      end if

   end function Fe_number_mantle

!########################################################################
   function Si_number_max(Mg_number, contents) result(y)

      implicit none

      real(8), dimension(2) :: y
      real(8) :: mn, mx
      integer, intent(in) :: contents(:)
      real(8), intent(in) :: Mg_number
      real(8), allocatable :: Si(:), Mg(:)
      integer :: i, n

      n = size(contents)

      allocate (Si(n), Mg(n))

      do i = 1, n
         Si(i) = material_YSi(contents(i))
         Mg(i) = material_YMg(contents(i))*Mg_number
      end do

!Compute Si# for given material
      do i = 1, n
         Si(i) = Si(i)/Mg(i)/(1d0 + Si(i)/Mg(i))
      end do

      y(1) = minval(Si)
      y(2) = maxval(Si)

!print *, 'Si# can take values between:', y(:)

   end function Si_number_max

!########################################################################
   function xi_stv_lim(YMg, YSi, xiFe, xiH2O, SiMg, xiAlSi, xiAlStv, &
                       xiH2OStv, xiAlMg) result(y)
!Computes the mole fractions xi_stv in the limiting case where only
!one other phase given by YMg, YSi is present

      implicit none

      real(kind=8) :: y, denom
      real(kind=8), intent(in) :: xiFe, xiH2O, SiMg, xiAlSi, &
                                  xiAlStv, xiH2OStv, xiAlMg
      integer, intent(in) :: YSi, YMg

      y = YSi*(1.0d0 - xiH2O)*(1.0d0 - xiAlSi)
      y = y - YMg*(1.0d0 - xiH2O)*(1.0d0 - xiAlMg)*(1.0d0 - xiFe)*SiMg

      denom = YSi*(1.0d0 - xiH2O)*(1.0d0 - xiAlSi)
      denom = denom - (1.0d0 - xiH2OStv)*(1.0d0 - xiAlStv)
      denom = denom - YMg*(1.0d0 - xiH2O)*(1.0d0 - xiAlMg)*(1.0d0 - xiFe)*SiMg

      y = y/denom
      y = max(y, 0.0d0)

   end function xi_stv_lim

!########################################################################
   function compute_xiFei(Fe_number, P, T, n_mats, contents) result(y)

      implicit none

      integer, intent(in) :: n_mats
      real(8), intent(in) :: Fe_number, P, T
      real(8), dimension(n_mats) :: y
      real(8), dimension(n_mats - 1) :: DFe
      integer, dimension(n_mats), intent(in) :: contents
      integer :: i

!Compute iron partitioning coefficient
!Here the fit results to experimental data will be called

!Compute partitioning coefficient of iron between substance i and i+1
!Convention is: DFe = xi_Fe_(i+1)/xi_Fe_i
      if (n_mats > 1) then
         do i = 1, n_mats - 1
            DFe(i) = 1.0d0

            !If 3 materials are present, it is assumed that the third one is Stv
            !which contains no iron. In general it can also be another phase
            !for which a different value of DFe could be set.
            if (i == 3) then
               DFe(i) = 0.0d0
            end if

         end do

         y(1) = n_mats*Fe_number/(DFe(1) + 1.0d0)

         do i = 1, n_mats - 1
            y(i + 1) = n_mats*Fe_number*DFe(i)/(DFe(i) + 1.0d0)
         end do

      else

         y(1) = Fe_number

      end if

   end function compute_xiFei

!########################################################################
   function xiAli(xiFe, P, T, n_mats) result(y)

      implicit none

      integer, intent(in) :: n_mats
      real(kind=8), intent(in) :: P, T
      real(kind=8), dimension(n_mats), intent(in) :: xiFe
      real(kind=8), dimension(n_mats) :: y

!Compute Al content in each phase at given P, T and for the Fe content
!in the phase that has been computed at P, T via the function xiFei
      y(1) = 0.0d0
      y(2) = 0.0d0
      y(3) = 0.0d0

   end function xiAli

!########################################################################
   function Al_partitioning_SiMg(xiAl, P, T, n_mats) result(y)

      implicit none

      real(kind=8), intent(in) :: P, T
      integer, intent(in) :: n_mats
      real(kind=8), dimension(n_mats) :: xiAl
      real(kind=8), dimension(n_mats, 2) :: y
      real(kind=8), dimension(n_mats) :: DAl_SiMg
      integer :: i

      do i = 1, n_mats
         if (i < 3) then
            DAl_SiMg(i) = 0.5d0

         else
            DAl_SiMg(i) = 1.0d10
         end if
      end do

      do i = 1, n_mats
         y(i, 1) = xiAl(i)/(1.0d0 + DAl_SiMg(i))
         y(i, 2) = xiAl(i)*(1.0d0 - 1.0d0/(1.0d0 + DAl_SiMg(i)))
      end do

   end function Al_partitioning_SiMg

!#######################################################################
   function compute_xiH2O(xiFei, xiAli, P, T, n_mats, contents) result(y)

      implicit none

      integer, intent(in) :: n_mats
      integer, dimension(n_mats), intent(in) :: contents
      real(kind=8), dimension(n_mats), intent(in) :: xiFei, xiAli
      real(kind=8), intent(in) :: P, T
      real(kind=8), dimension(n_mats) :: y
      integer :: i

      do i = 1, n_mats

         !Here, call a function to compute the water content as function of
         !P, T, Fe, Al for each phase

         y(i) = 0.0d0
      end do

   end function compute_xiH2O

!#######################################################################
   SUBROUTINE update_all_fractions(X_new, weight_fractions, fractions, &
                                   molar_masses, n_mats)
!Uses the newly computed weight fraction of the impurity and updates all
!other weight fractions and all mole fractions. new_last_frac is the
!weight fraction of the impurity (last material)

      real(8), intent(in) :: X_new
      real(8) :: xi_new
      integer, intent(in) :: n_mats
      real(8), intent(inout), dimension(n_mats) :: weight_fractions, fractions
      real(8), intent(in), dimension(n_mats) :: molar_masses
      real(8), dimension(n_mats) :: olds
      real(8) :: sum_others, m_others, m_impurity, delta
      integer :: i

!First update the weight fractions
      olds = weight_fractions
      weight_fractions(n_mats) = X_new

!Compute normalization factor
      sum_others = 0d0
      do i = 1, n_mats - 1
         sum_others = sum_others + olds(i)
      end do

!Compute weight fractions
      do i = 1, n_mats - 1
         weight_fractions(i) = olds(i)*(1d0 - X_new)/sum_others
      end do

!Compute the normalized mass of the non-impurities for the conversion
!to mole fraction
      m_others = 0d0

      do i = 1, 2
!The layer itself does not contain the right fractions so they have to
!be taken from the shells. Since the first two material fractions will
!be constant we can just take the first shell.
         delta = molar_masses(i)
         m_others = m_others + delta*fractions(i)
      end do

!Normalized mass of the non-impurities
      m_others = m_others/(fractions(1) + fractions(2))

!Mass of the impurity (last material)
      m_impurity = molar_masses(n_mats)

!Compute the mole fraction of the impurity
      xi_new = compute_xi(eta=X_new, m1=m_others, m2=m_impurity)

!Use this to compute all the other mole fractions
      olds = fractions
      fractions(n_mats) = xi_new

!Compute normalization factors
      sum_others = 0d0
      do i = 1, n_mats - 1
         sum_others = sum_others + olds(i)
      end do

!Compute mole fractions
      do i = 1, n_mats - 1
         fractions(i) = olds(i)*(1d0 - xi_new)/sum_others
      end do

   END SUBROUTINE update_all_fractions

!########################################################################
   subroutine construct_abundance_matrix(SiMg, FeMg, n_mats, &
                                         YMgi, YSii, xiH2Oi, matrix, b, xiFei, xiAlSii, xiAlMgi, additional)

!Constructs linear system of equations to solve for the individual
!molar abundances for fixed atomic abundances and materials in a
!mixture.

      implicit none

      integer, intent(in) :: n_mats
      integer, dimension(n_mats), intent(in) :: YMgi, YSii
      real(kind=8), dimension(n_mats), intent(in) :: xiH2Oi, xiAlSii, xiAlMgi, &
                                                     xiFei
      real(8), intent(in) :: additional(:)
      real(kind=8), intent(in) :: SiMg, FeMg
      real(kind=8), dimension(n_mats, n_mats), intent(out) :: matrix
      real(8), intent(out) :: b(:, :)
      integer :: i, j, k

      do j = 1, n_mats
         matrix(1, j) = (1.0d0 - xiH2Oi(j))*(SiMg*(1.0d0 - xiAlMgi(j))* &
                                             (1.0d0 - xiFei(j))*YMgi(j) - (1.0d0 - xiAlSii(j))*YSii(j))

         !matrix(2,j) = (1.0d0-xiH2Oi(j))*(FeMg*(1.0d0-xiAlMgi(j))*&
         ! (1.0d0-xiFei(j))*YMgi(j) - (1.0d0-xiAlMgi(j))*xiFei(j)*YMgi(j))

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

!b(1,1)=-xiStv*(1.0d0-xiH2Oi(n_mats))*(SiMg*&
!(1.0d0-xiAlMgi(n_mats))*(1.0d0-xiFei(n_mats))*YMgi(n_mats) - &
!(1.0d0-xiAlSii(n_mats))*YSii(n_mats))

!-xiStv

!print *, 'matrix'
!do i=1, n_mats - 1
      ! print*, matrix(i,:)
!enddo

   end subroutine construct_abundance_matrix

!########################################################################
   subroutine compute_abundance_vector(SiMg, FeMg, n_mats, &
                                       YMgi, YSii, xiH2Oi, xiFei, xiAlSii, xiAlMgi, abundances, &
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
      real(kind=8), dimension(n_mats), intent(in) :: xiH2Oi, xiAlSii, xiAlMgi, &
                                                     xiFei
      real(kind=8), intent(in) :: SiMg, FeMg
      real(kind=8), dimension(n_mats, n_mats) :: matrix
      real(kind=8), dimension(n_mats), intent(out) :: abundances
      real(kind=8), dimension(n_mats, 1) :: abundances_dummy
      real(8), dimension(n_mats, 1) :: b_vec
      real(8), intent(in), optional :: additional(:)
      real(8), allocatable :: additional_dummy(:)
      integer :: i, j

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

!~ print *, 'contents =', contents(:)
!~ print *, 'xiAlSii =', xiAlSii(:)
!~ print *, 'xiAlMgi =', xiAlMgi(:)
!~ print *, 'FeMg =', FeMg
!~ print *, 'SiMg =', SiMg
!~ print *, 'xiH2Oi =', xiH2Oi
!~ print *, 'xiFei =', xiFei
!~ print *, 'additional =', additional_dummy

!In the core the first material is pure iron. There no compositional
!gradients are currently allowed and the fractions need not be computed
      if (contents(1) == 2 .or. contents(1) == 1) then
!pass

      elseif (size(contents) == 1) then
         abundances(1) = 1d0

      else

         call construct_abundance_matrix(SiMg=SiMg, FeMg=FeMg, n_mats=n_mats, &
                                         YMgi=YMgi, YSii=YSii, xiH2Oi=xiH2Oi, xiFei=xiFei, xiAlSii=xiAlSii, &
                                         xiAlMgi=xiAlMgi, matrix=matrix, b=b_vec, &
                                         additional=additional_dummy)

!call mat_print('matrix', matrix)
!print *, 'b_vec =', b_vec

         call gauss_elimination(a=matrix, sol=abundances_dummy, b=b_vec)

!print *, 'abundances =', abundances_dummy(:,1)

         do i = 1, n_mats
            abundances(i) = abundances_dummy(i, 1)
         end do
      end if

   end subroutine compute_abundance_vector

!########################################################################
   subroutine gradients(r, y, nmat, ll, gammaG0, tempType, grads, &
                        q, d0, adiabatType, eps_T_zero, xi_Fe, eps_H2O, &
                        eps_Al, omega, phi, fractions, weight_fractions, &
                        composition_gradients, molar_masses, lay, xi_H, &
                        external_temp_profile)

      implicit none

      integer, intent(in) :: nmat, tempType, adiabatType, lay
      real(8), intent(in) :: r, gammaG0, d0, eps_T_zero, q, xi_Fe, &
                             eps_H2O, eps_Al, omega, phi, xi_H
      real(8), intent(in), dimension(n_params_integration) :: y
      real(8), intent(inout), dimension(nmat) :: fractions, weight_fractions
      real(8), intent(in), dimension(nmat) :: composition_gradients, &
                                              molar_masses
      integer, intent(in), dimension(nmat) :: ll
      real(8), intent(in), dimension(3), optional :: external_temp_profile
      real(8) :: dPdr, dTdr, drhodr, dmdr, gammaG, dMOIdr, dXdr, dXdm
      real(8), intent(out), dimension(n_params_integration) :: grads
      real(8) :: P, m, T, d, dPdrho_mean, alpha_mean, KS_mean, KT_mean
      real(8) :: rho_mean, delta_T_zero, T_zero, xi_H2O, xi_Al, X
      real(8), dimension(n_out) :: params_individual
      integer, dimension(n_out) :: params
      real(8), dimension(3) :: vals
      integer :: i, order


      order = 1
!ommit y(5) which is the MoI and only passively integrated
      P = y(1)
      m = y(2)
      T = y(3)
      d = y(4)
      X = y(6)

! If an external temperature profile is followed extract the temperature
! from the data instead of the integration step. Just replace the temp
! value from the integration step by the one from the data.
!~ if(present(external_temp_profile))then
!~         T = extract_temperature(r, external_temp_profile, lay, y)
!~ endif

!X here is only the impurity mass fraction (case nmat > 2). The other
!fractions depend on X and must be updated here in order to correctly
!compute the gradients with the correct fractions. This is only done
!for the mantle.
      if (nmat > 2 .and. lay > 2 .and. lay < 5) then
         call update_all_fractions(X_new=X, weight_fractions=weight_fractions, &
                                   fractions=fractions, molar_masses=molar_masses, n_mats=size(fractions))
      end if

!At low T and high P the EoS have trouble and can lead to spurious
!behaviour in the code. Therefore I truncate the computaion of the
!thermal profile if T decreases below T_zero. T_zero is defined such
!that gamma-Olivine yields still meaningful results in the EoS. At
!lower T an isothermal profile is envoked
      delta_T_zero = 45.0e0*(P*1.0e-9 - 20.0e0)*eps_T_zero
      T_zero = 100.0e0 + delta_T_zero

      do i = 1, n_out
         params(i) = i
      end do

      vals(1) = T
      vals(2) = P
!vals(4) = eps_H2O
!vals(5) = eps_Al

!Compute mean quantities
      dPdrho_mean = 0.0d0
      alpha_mean = 0.0d0
      rho_mean = 0.0d0
      KT_mean = 0.0d0

      do i = 1, nmat
         if (ll(i) == 2) then
            vals(3) = 1d0 - xi_H
         else
            vals(3) = xi_Fe
         end if

         call compute(vals=vals, which=params, n_out=n_out, ll=ll(i), res=params_individual, &
                      order=order, alloc=.false., eps_H2O=eps_H2O)

         !To avoid spurious behaviour if water reaches vapor state
         if (ll(i) == 1) then
            params_individual(1) = max(500d0, params_individual(1))

         end if

         rho_mean = rho_mean + weight_fractions(i)/params_individual(1)
         dPdrho_mean = dPdrho_mean + weight_fractions(i)/params_individual(3)/params_individual(1)**2
         alpha_mean = alpha_mean + params_individual(4)*weight_fractions(i)/params_individual(1)
      end do

      rho_mean = 1.0e0/rho_mean
      dPdrho_mean = dPdrho_mean*rho_mean**2
      dPdrho_mean = 1.0e0/dPdrho_mean

      KT_mean = max(rho_mean*dPdrho_mean, KT_zero)

      alpha_mean = alpha_mean*rho_mean

!Compute the radial gradients of all structure parameters:

!pressure
      dPdr = rho_mean*(omega/r*COS(phi) - G*m/r**2)

!enclosed mass
      dmdr = 4.0d0*PI*r**2*rho_mean

!If nmat > 2 the last material will change (impurity)
!composition gradient of last material.
      dXdm = composition_gradients(nmat)
      dXdr = dmdr*dXdm

!temperature
      if (tempType == 0) then
         dTdr = 0.0d0

      elseif (tempType == 1) then

         !Assumes that water is not mixed with other materials
         if (ll(1) == 1) then

            !For low temperatures no reliable EoS for water is implemented
            !Switch to isothermal profile there as we are not interested
            !in modelling this region for the planets anyways
            if (T .lt. 0e0) then
               dTdr = -1.0e-6

            else
               !Compute adiabatic gradient as dTdr = dTdP*dPdr
               dTdr = params_individual(2)*dPdr
            end if

         else
            if (T .lt. T_zero .and. T .lt. 1000.0e0) then
               dTdr = -1.0e-6

            else
               !Compute Grueneisenparameter at current density
               gammaG = gammaG0*(rho_mean/d0)**(-q)
               if (adiabatType == 0) then
                  dTdr = dPdr*gammaG*T/KT_mean

               elseif (adiabatType == 1) then
                  !Compute adiabatic bulk modulus
                  KS_mean = KT_mean*(1.0d0 + gammaG*alpha_mean*T)

                  !Check if alpha is nan and take KT instead of KS
                  if (alpha_mean /= alpha_mean) then
                     KS_mean = KT_mean
                  end if

                  if (KS_mean .lt. 1.0e9 .or. KS_mean .gt. 1.0e12) then
                     if (KT_mean .lt. KT_zero) then
                        KT_mean = KT_zero
                     end if

                     KS_mean = KT_mean

                  end if

                  !Compute radial adiabatic gradient
                  dTdr = dPdr*gammaG*T/KS_mean
               end if
            end if
         end if
      end if

!density
      drhodr = dPdr/dPdrho_mean

!moment of inertia
      dMOIdr = 8d0/3d0*PI*rho_mean*r**4

      grads(1) = dPdr
      grads(2) = dmdr
      grads(3) = dTdr
      grads(4) = drhodr
      grads(5) = dMOIdr
      grads(6) = dXdr
      grads(7) = - G * m / r * dmdr ! gravitational energy

      ! TODO. implement self-consistent prescription for heat capacity
      grads(8) = 10000 * T * dmdr ! internal energy

!~ print *, 'fractions =', weight_fractions(:)
!~ print *, 'vals =', vals(:)
!~ print *, 'means =', rho_mean, dPdrho_mean, alpha_mean, KT_mean
!~ print *, 'grads in grads =', grads(:)

   end subroutine gradients
!########################################################################

!########################################################################
   subroutine integrateRK(r, y_in, y_out, r_start, h, nmat, &
                          fractions, tempType, ll, gammaG0, q, d0, adiabatType, grads, &
                          eps_T_zero, xi_Fe, eps_H2O, eps_Al, omega, phi, weight_fractions, &
                          composition_gradients, molar_masses, lay, xi_H, external_temp_profile)

      implicit none

      integer :: i, j, k, s
      integer, intent(in) ::  nmat, tempType, adiabatType, lay
      integer, intent(in), dimension(nmat) :: ll
      real(8), intent(in), dimension(3), optional :: external_temp_profile
      real(8), intent(inout), dimension(nmat) :: fractions, weight_fractions
      real(8), intent(in), dimension(nmat) :: composition_gradients, &
                                              molar_masses
      real(8), intent(in) :: r_start, h, gammaG0, q, d0, eps_T_zero, &
                             xi_Fe, eps_H2O, eps_Al, omega, phi, xi_H
      real(8), intent(in), dimension(n_params_integration) :: y_in
      real(8), intent(out), dimension(n_params_integration) :: y_out
      real(8), intent(out), dimension(n_params_integration) :: grads
      real(8), dimension(4, n_params_integration) :: k_list
      real(8), dimension(n_params_integration) :: ki_list
      real(8), intent(out) :: r
!~ print *, 'y_in =', y_in
      r = r_start
      do i = 1, 4
!        print *, "i=", i
         do k = 1, size(y_out)
            ki_list(k) = y_in(k)
         end do
!        print *, "ki_list =", ki_list

         do j = 1, i - 1
            do s = 1, size(y_out)
               ki_list(s) = ki_list(s) + h*k_list(j, s)*a_list(i, j)
            end do
!                print *, "ki_list =", ki_list
         end do
         r = r_start + h*c_list(i)

         do k = 1, size(y_out)
            y_out(k) = ki_list(k)
         end do
!        print *, "grads before f=", grads

         !call f(grads, y_out)
         call gradients(r, y_out, nmat, ll, gammaG0, tempType, &
                        grads, q, d0, adiabatType, eps_T_zero, xi_Fe, eps_H2O, eps_Al, &
                        omega, phi, fractions, weight_fractions, &
                        composition_gradients=composition_gradients, &
                        molar_masses=molar_masses, lay=lay, xi_H=xi_H, &
                        external_temp_profile=external_temp_profile)
         !call f(grads, y_out)

         do k = 1, size(y_out)
            k_list(i, k) = grads(k)
         end do
      end do

      do i = 1, size(y_out)
         y_out(i) = y_in(i)

      end do

      do i = 1, size(y_out)
         do j = 1, 4
            y_out(i) = y_out(i) + h*k_list(j, i)*b_list(j)
         end do
      end do

! If an external temperature profile is followed extract the temperature
! from the data instead of the integration step. Just replace the temp
! value from the integration step by the one from the data.
!~ if(present(external_temp_profile))then
!~         y_out(3) = extract_temperature(r, external_temp_profile, lay, y_out)
!~ endif

!~ print *, "y_out =", y_out

!~ if (y_out(3)<10.0)then
!~         y_out(3)=10.0
!~ endif

   end subroutine integrateRK
!########################################################################

!
   FUNCTION extract_temperature(radius, coeffs, lay, y) result(temp)

      implicit none

      real(8), intent(in) :: radius
      integer, intent(in) :: lay
      real(8), intent(in), dimension(n_params_integration) :: y
      real(8), dimension(3), intent(in) :: coeffs
      real(8) :: temp

      temp = coeffs(1) + coeffs(2)*(radius*1e-3) + coeffs(3)*(radius*1e-3)**2

! Temperature jump between ocean and ice layer
      if (radius .gt. 1.53d6) then
         temp = temp - 1.5d1
      end if

! Convective --> conductive transition in the mantle
      if (lay == 4 .or. lay == 3) then
         if (radius .gt. 1350d3) then
            temp = 1.769d4 - 12.05*radius*1d-3
            !print *, "temp, radius =", temp, radius * 1d-3
         end if

         temp = max(temp, 273d0)
      end if

      temp = max(temp, 5d1)

   END FUNCTION extract_temperature

END MODULE functions

