MODULE integrator

   use constants
   use fortfunctions
   use run_params
   use eosfort
   use linalg

   implicit none

contains

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


END MODULE integrator

