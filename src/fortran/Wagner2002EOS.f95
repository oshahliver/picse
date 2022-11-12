subroutine delta(the, bb, del, a, result)

   implicit None

   real(kind=8), intent(in) :: the, bb, del, a
   real(kind=8), intent(out) :: result

   result = the**2 + bb*((del - 1.0)**2)**a

end subroutine delta

subroutine theta(tau, a, del, beta, result)

   implicit None

   real(kind=8), intent(in) :: tau, a, del, beta
   real(kind=8), intent(out):: result

   result = (1.-tau) + a*((del - 1.)**2)**(1.0/(2.0*beta))

end subroutine theta

subroutine psi(tau, c, d, del, result)

   implicit None

   real(kind=8), intent(in) :: tau, c, d, del
   real(kind=8), intent(out) :: result

   result = EXP(-c*(del - 1.0)**2 - d*(tau - 1.0)**2)

end subroutine psi

subroutine ddeltadtau(the, bb, del, a, result)

   implicit None

   real(kind=8), intent(in) :: the, bb, del, a
   real(kind=8), intent(out) :: result

   result = -2.*the

end subroutine ddeltadtau

subroutine ddeltaddelta(the, aa, bb, beta, del, a, result)

   implicit None

   real(kind=8), intent(in) :: the, a, bb, beta, del, aa
   real(kind=8), intent(out) :: result

   result = 2*the*aa/beta*(del - 1.0)**(1./beta - 1.0) + 2*a*bb*(del - 1.0)**(2.*a - 1.)

end subroutine ddeltaddelta

subroutine ddeltaddelta_2(the, aa, bb, beta, del, a, result)

   implicit None

   real(kind=8), intent(in) :: the, a, bb, beta, del, aa
   real(kind=8), intent(out) :: result

   result = 2.*aa**2./beta**2.*(del - 1.)**(2./beta - 2.) + 2.*the*aa/beta* &
            (1./beta - 1.)*(del - 1.)**(1./beta - 2.) + 2.*a*bb*(2.*a - 1.)*(del - 1.)**(2.*a - 2.)

end subroutine ddeltaddelta_2

subroutine dpsidtau(tau, cc, dd, del, result)

   implicit None

   real(kind=8), intent(in) :: tau, cc, dd, del
   real(kind=8), intent(out) :: result
   real(kind=8) :: ps

   call psi(tau, cc, dd, del, ps)

   result = -2.*dd*(tau - 1.)*ps

end subroutine dpsidtau

subroutine dpsiddeltadtau(tau, cc, dd, del, result)

   implicit None

   real(kind=8), intent(in) :: tau, cc, dd, del
   real(kind=8), intent(out) :: result
   real(kind=8) :: ps

   call psi(tau, cc, dd, del, ps)

   result = 4.*cc*dd*(del - 1.)*(tau - 1.)*ps

end subroutine dpsiddeltadtau

subroutine ddeltaddeltadtau(the, aa, bb, beta, del, result)

   implicit None

   real(kind=8), intent(in) :: the, aa, bb, beta, del
   real(kind=8), intent(out) :: result

   result = -2.*aa/beta*(del - 1.)**(1./beta - 1.)

end subroutine ddeltaddeltadtau

subroutine dpsiddelta(tau, cc, dd, del, result)

   implicit None

   real(kind=8), intent(in) :: tau, cc, dd, del
   real(kind=8), intent(out) :: result
   real(kind=8) :: ps

   call psi(tau, cc, dd, del, ps)

   result = -2.*cc*(del - 1.)*ps

end subroutine dpsiddelta

subroutine dpsiddelta_2(tau, cc, dd, del, result)

   implicit None

   real(kind=8), intent(in) :: tau, cc, dd, del
   real(kind=8), intent(out) :: result
   real(kind=8) :: ps

   call psi(tau, cc, dd, del, ps)

   result = -2.*cc*ps*(1.-(del - 1.)**2*2.*cc)

end subroutine dpsiddelta_2

subroutine ddeltadtau_2(result)

   implicit None

   real(kind=8), intent(out) :: result

   result = 2.0

end subroutine ddeltadtau_2

subroutine dpsidtau_2(tau, cc, dd, del, result)

   implicit None

   real(kind=8), intent(in) :: tau, cc, dd, del
   real(kind=8), intent(out) :: result
   real(kind=8) :: ps

   call psi(tau, cc, dd, del, ps)

   result = 2.*dd*ps*((tau - 1.)**2*2*dd - 1.)

end subroutine dpsidtau_2
