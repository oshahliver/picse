# this is a simple python program that uses a fortran subroutine
import numpy as np
from matplotlib import pyplot as plt
from picse.materials import mazevet
import time
from picse.utils.function_tools import functionTools as ftool

# -----------------------------------------------------------------------------
# functions


def pressure(d=None, T=None, **kwargs):
    # NOTE: here d is in gcc and output is in Mbar
    P = mazevet.compute(d / 1000.0, T)[0]
    # convert to Pa
    return P * 1.0e6 * 1.0e5


def density(P=None, T=None, **kwargs):
    """Compute density in kg/m3"""
    return ftool.bisec(
        whicharg="d",
        f=pressure,
        T=T,
        y=P,
        a=1.0e-3,
        b=1.0e6,
        identity="density",
        acc=1.0e-10,
        **kwargs
    )


def F(T=None, P=None, d=None):
    """Compute Helmholty free energy"""
    if d == None:
        d = density(P=P, T=T)

    return mazevet.compute(d / 1000.0, T)[5]


def U(T=None, P=None, d=None):
    """Compute internal energy"""
    if d == None:
        d = density(P=P, T=T)

    return mazevet.compute(d / 1000.0, T)[6] * 1000.0


def s_spec(T=None, P=None, d=None):
    pass


def u_spec(P=None, T=None, d=None, **kwargs):
    """compute specific internal energy"""
    if d == None:
        d = density(P=P, T=T)

    return mazevet.compute(d / 1000.0, T)[3] * 1000.0 * 1.0e-7


def cp_spec(T=None, P=None, d=None, acc=1.0e-4):
    """Compute isobaric specific heat capacity in J/K"""
    return ftool.deriv(f=u_spec, whicharg="T", x0=T, P=P, acc=acc)


def alpha_th_p(T=None, P=None, d=None):
    """Compute isobaric thermal expansion coefficient"""
    if d == None:
        d = density(P=P, T=T)

    drhodT_P = ftool.deriv(f=density, whicharg="T", x0=T, P=P, acc=1.0e-4)
    return -1.0 / d * drhodT_P


def dlnTdlnP(d=None, t=None, **kwargs):
    # derivative of log(T) with respect ot log(P) at constant entropy
    return 1.0 / mazevet.compute(d / 1000.0, t)[1]


def dlnPdlnrho(t=None, d=None, **kwargs):
    # derivative of log(P) with respect to log(rho) at constant entropy
    return mazevet.compute(d / 1000.0, t)[2]


def dPdrho_T(T=None, P=None, d=None):
    if d == None:
        d = density(P=P, T=T)

    dPdrho = ftool.deriv(f=pressure, whicharg="d", x0=d, T=T, acc=1.0e-6)
    return dPdrho


def dPdrho_S(T=None, P=None, d=None, **kwargs):
    """computes pressure derivative with respect to rho at constant entropy
    from the logarithmic derivative
    """
    if d == None:
        d = density(P=P, T=T)

    # compute logarithmic derivative
    # remember: pres is in Mbar and dens in gcc in Mazevet
    derlog = mazevet.compute(d / 1000.0, T)[2]

    # convert Mbar to Pa in the denuminator
    derlog = derlog * 1.0e6 * 1.0e5

    # convert gcc to kg m-3 in the denominator
    derlog = derlog / 1000.0

    # convert logarithmic derivative to linear derivative
    derlog = d / P * derlog

    return derlog


def dTdP_S(d=None, T=None, P=None, type=1, **kwargs):
    """Computes adiabatic gradient"""
    if type == 0:
        # if no density is given, compute it for the PT at hand
        if d == None:
            d = density(P=P, T=T)

        # compute isobaric derivative of specific internal energy with
        # respect to density at given P-T-rho
        dudrho = ftool.deriv(f=u_spec, whicharg="d", x0=d, t=T, acc=1.0e-4)

        # compute isobaric derivative of specific entropy with respect to
        # density at given pressure
        dsdrho = 1.0 / (T * d ** 2) * (d ** 2 * dudrho - P)

        # compute dT/dP_s
        der = -1.0 / d ** 2 / dsdrho
        return der

    elif type == 1:
        if d == None:
            d = density(P=P, T=T)

        alpha = alpha_th_p(T=T, P=P, d=d)
        cp = cp_spec(T=T, P=P, d=d)

        return alpha * T / (cp * d)


def together(T=None, P=None):
    T = max(T, 300)

    d = density(T=T, P=P)
    pres, dlnpdlrho, dlnPdlnrho, u_spec, c_v, E, U = mazevet.compute(d / 1000.0, T)

    # convert Mbar to Pa in the denuminator
    dlnPdlnrho = dlnPdlnrho * 1.0e6 * 1.0e5

    # convert gcc to kg m-3 in the denominator
    dlnPdlnrho = dlnPdlnrho / 1000.0

    # convert logarithmic derivative to linear derivative
    dPdrho = dPdrho_T(T=T, P=P, d=d)

    # compute specific isobaric heat capacity
    cp = cp_spec(T=T, P=P, d=d)

    # compute specific internal energy
    u_spec = u_spec * 1000.0 * 1.0e-7

    # compute specific entropy
    s_spec = 0.0

    # compute isobaric thermal expansion
    alpha = alpha_th_p(T=T, P=P, d=d)

    # compute dTdP_S
    dTdP_S = alpha * T / (cp * d)

    return d, dTdP_S, dPdrho, alpha, cp, s_spec, u_spec


def test(T=None, P=None):
    t0 = time.time()
    a, b, c, d, e, f, g = together(P=P, T=T)
    t = time.time()
    print("elapsed time for simultaneous evaluation:", t - t0)

    t0 = time.time()
    dens = density(T=T, P=P)
    dPdrho = dPdrho_T(T=T, P=P, d=dens)
    dTdP = dTdP_S(T=T, P=P, d=dens)
    cp = cp_spec(T=T, P=P, d=dens)
    alpha = alpha_th_p(T=T, P=P, d=dens)
    u = u_spec(T=T, P=P, d=dens)
    s = s_spec(T=T, P=P, d=dens)
    t = time.time()
    print("elapsed time for seprate evaluation:", t - t0)

    print("d:", a, dens)
    print("dTdP_S:", b, dTdP)
    print("dPdhro_T:", c, dPdrho)
    print("alpha:", d, alpha)
    print("cp:", e, cp)
    print("s:", f, s)
    print("u:", g, u)


def compute(what="dens", **kwargs):
    if what == "dens":
        return density(**kwargs)

    elif what == "pres":
        return pressure(**kwargs)

    # other parameter outputs can be added here after they have been specified
    # above...
    # elif which == 'other parameter':
    #    return other parameter stuff


# -----------------------------------------------------------------------------
# testing
# this part is normally commented out
"""
d0 = 1500.
testpres = (compute('pres', d = d0, t = 278.5))
testdens = (compute('dens', p = testpres, t = 278.5, acsc = 1.0e-5))
print ('test pressure:', testpres)
print ('test density: ',testdens)
print ('original density: ', d0)
print ('deviation:',round((d0 - testdens)/d0 *100,6),'%')
"""
