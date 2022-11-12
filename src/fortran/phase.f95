MODULE phase

   use constants
   use run_params

   implicit none

contains

!#######################################################################
   FUNCTION T_solidus_H2O_high(P) result(T)

      real(8), intent(in) :: P
      real(8) :: T
      integer :: i

      T = 0d0
      if (P < 1e9) then
         do i = 1, size(LtoVII_coeffs_Dunaeva2010)
            T = T + LtoVII_coeffs_Dunaeva2010(i)* &
                (P*1d-9)**(size(LtoVII_coeffs_Dunaeva2010) - i)
         end do

      else
         do i = 1, size(coeffs_French2010)
            T = T + coeffs_French2010(i)* &
                (P*1d-9)**(size(coeffs_French2010) - i)
         end do
      end if

   END FUNCTION T_solidus_H2O_high

!#######################################################################
   FUNCTION T_solidus_H2O(P) result(T)

      real(8), intent(in) :: P
      real(8) :: T

      if (P >= 5d9) then
         T = T_solidus_H2O_high(P)

      end if

   END FUNCTION T_solidus_H2O

!#######################################################################
   FUNCTION T_melt_MgSiO3(P) result(T)
!Compute melting temperature of silicate-perovskite from Belonoshko 2005
      real(8), intent(in) :: P
      real(8) :: T, a, b, T0, P0

      a = 4.6d9
      b = 0.33d0
      T0 = 1831d0
      P0 = 0d0

      T = T0*(1d0 + (P - P0)/a)**b

   END FUNCTION T_melt_MgSiO3

!#######################################################################
   FUNCTION T_melt_CaSiO3(P) result(T)
!Compute melting temperature of calcium-perovskite from Braithwaite &
!Stixrude 2019
      real(8), intent(in) :: P
      real(8) :: T, a, b, c, T0, P0

      a = 48.8d9
      b = 0.413d0
      c = 860d9
      T0 = 4020d0
      P0 = 44.2d9

      T = T0*(1d0 + (P - P0)/a)**b*exp(-(P - P0)/c)

   END FUNCTION T_melt_CaSiO3

!#######################################################################
   FUNCTION T_liquidus_pyrolite(P) result(T)

      real(8), intent(in) :: P
      real(8) :: T, T0, a, c

      T0 = 1940d0
      a = 29d9
      c = 1.9d0

      T = T0*(1d0 + P/a)**(1d0/c)

   END FUNCTION T_liquidus_pyrolite

!#######################################################################
   FUNCTION T_melt_iron(P, X_Si, X_O, X_S) result(T)
!Compute melting temperature of iron as a function of pressure according
!to Li et al. 2020. The effect of impurities is modelled for S, Si and O
!according to Andrault et al. 2016.
      real(8), intent(in) :: P, X_Si, X_O, X_S
      real(8) :: T
      real(8) :: a = 23d9, b = 2.26d0, T0 = 1811d0, P0 = 1d5

!Melting temperature of pure Fe from Li et al. 2020
      T = T0*((P - P0)/a + 1d0)**(1d0/b)

!Effect of impurities from Andrault 2016 (there as K/wt%)
      T = T - 3d3*X_Si - 5d3*X_O - 1d4*X_S

   END FUNCTION T_melt_iron

!#######################################################################
   SUBROUTINE get_phase(ll, T, P, ph, xiFe)

      implicit none

      integer, intent(in) :: ll
      real(8), intent(in) :: xiFe
      real(8), intent(inout) :: T, P
      integer, intent(out) :: ph
      real(8) :: T_trans, P_trans, P_trans1, P_trans2, pres
      integer :: i, j

!~ print *, 'll / T / P in get_phase =', ll, T, P
!Water

!Magnesiowüstit
      if (ll == 4 .or. ll == which_eos_tables_hydrated(4)) then
         pres = P*1d-9
         if (P .lt. P_triple2_br) then
            T_trans = 0d0
            do i = 1, size(bruce4)
               T_trans = T_trans + bruce4(i)*pres**(size(bruce4) - i)
            end do
            if (T_trans .lt. T) then
               ph = 3
            else
               T_trans = 0d0
               do i = 1, size(bruce1)
                  T_trans = T_trans + bruce1(i)*pres**(size(bruce1) - i)
               end do
               if (T_trans .lt. T) then
                  ph = 1

               else
                  if (P < P_triple1_br) then
                     ph = 2

                  else
                     T_trans = 0d0
                     do i = 1, size(bruce2)
                        T_trans = T_trans + bruce2(i)*pres**(size(bruce2) - i)
                     end do
                     if (T_trans < T) then
                        ph = 3

                     else
                        ph = 2
                     end if

                  end if
               end if

            end if
         else
            ph = 3
         end if

!Perovksite
      elseif (ll == 5 .or. ll == which_eos_tables_hydrated(5)) then
         P_trans = 1.01d11 + 8.93e6*T
         if (P .lt. P_trans) then
            ph = 1
         else
            ph = 2
         end if

!Olivine
      elseif (ll == 6 .or. ll == which_eos_tables_hydrated(6)) then
         P_trans1 = P_alpha_beta(T, xiFe)
         P_trans2 = P_beta_gamma(T, xiFe)
         if (P < P_trans1) then
            ph = 1

         elseif (P >= P_trans1) then

            if (P_trans1 < P_trans2) then
               if (P <= P_trans2) then
                  ph = 2
               else
                  ph = 3
               end if
            elseif (P_trans1 >= P_trans2) then
               ph = 3
            end if
         end if

!Iron
      elseif (ll == 2) then
         T_trans = T_melt_iron(P, 0d0, 0d0, 0d0)
         if (T > T_trans) then
            ph = 2
         else
            ph = 1
         end if

      else
         ph = 1
      end if

   END SUBROUTINE get_phase

!#######################################################################
   function P_alpha_beta(T, xiFe) result(P)

      implicit none

      real(8), intent(in) :: T
      real(8) :: P, xiFe

      P = 1.159d1 + 1.999d-3*(T - 273.15d0) - 13d0*xiFe
      P = P*1d9 - 13d9*xiFe

   end function P_alpha_beta

!#######################################################################
   function P_beta_gamma(T, xiFe) result(P)

      implicit none

      real(8), intent(in) :: T
      real(8) :: P, xiFe

      P = 1.089d1 + 5.679d-3*(T - 273.15d0) - 20d0*xiFe
      P = P*1d9

   end function P_beta_gamma

!#######################################################################
   function P_Ol_Pv(T) result(P)

      implicit none

      real(8), intent(in) :: T
      real(8) :: P, P0 = 25d0, T0 = 800d0, alpha = -0.0017d0

      P = P0 !+ alpha*(T-T0)
      P = P*1e9

   end function P_Ol_Pv

END MODULE phase
