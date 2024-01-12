"""This file contains the routines to set up and extract the core properties.

The core is devided into two sublayers:

-- Inner core: Fe
-- Outer core: Fe + FeS + FeO + FeSi

The transition from the inner to the outer core is calculated from the liquidus of
iron-alloys using the results from 
"""

from picse.utils import fortfunctions
from picse.physicalparams import (
    material_list_fort,
    material_YMg,
    material_YSi,
    material_YO,
    material_YH,
    material_YS,
    material_YFe,
    mSi,
    mS,
    mO,
    mH,
    mFe,
    mMg,
    r_earth,
    m_earth,
    sigmaSB,
    G,
    mH2O,
    material_list,
)


def compute_core_mass(props, n=3, M_IC=None, inner_core_mass_fraction = None, mode = 1):
    """Computes the core mass of a planet at given total mass, composition and
    value for Mg#
    """

    # print (contents, Mg_number, M_surface, Mg_number_mantle, SiMg, M_ocean, xi_H_core)
    Mg_number_mantle = min(1.0 - props["Fe_number_mantle"], 0.9999999999)
    FeMg_mantle = (1.0 - Mg_number_mantle) / Mg_number_mantle
    FeMg = (1.0 - props["Mg_number_should"]) / props["Mg_number_should"]
    SiMg = props["Si_number_mantle"] / (1.0 - props["Si_number_mantle"])

    core_mass = fortfunctions.functionspy.compute_core_mass(m_tot = props["M_surface_should"],
                                                            m_ocean = props["M_surface_should"]*10**props["ocean_fraction_should"],
                                                            femg = FeMg,
                                                            simg = SiMg,
                                                            fe_numbers = [1.0 - Mg_number_mantle for i in props["contents"][n]],
                                                            xi_all_core = props["x_all_core"],
                                                            contents = props["contents"][n],
                                                            inner_core_mass_fraction = inner_core_mass_fraction,
                                                            inner_core_mass = M_IC,
                                                            mode = mode,
                                                            additional = [])
    return core_mass

    # Compute the fractions in the mantle
    fractions = fortfunctions.functionspy.compute_abundance_vector(
        simg=SiMg,
        femg=FeMg_mantle,
        n_mats=len(props["contents"][n]),
        ymgi=[material_YMg[i - 1] for i in props["contents"][n]],
        ysii=[material_YSi[i - 1] for i in props["contents"][n]],
        # xih2oi=[0.0 for i in props["contents"][n]],
        xifei=[1.0 - Mg_number_mantle for i in props["contents"][n]],
        # xialsii=[0.0 for i in props["contents"][n]],
        # xialmgi=[0.0 for i in props["contents"][n]],
        contents=props["contents"][n],
        additional=[],
    )

    print ("core mass test =", core_mass)
    # Count
    # Note that the indices are shifted by one because the original definition
    # of the arrays comes from the fortran code.
    Q1 = sum(
        [
            fractions[i] * Mg_number_mantle * material_YMg[props["contents"][n][i] - 1]
            for i in range(len(props["contents"][n]))
        ]
    )

    # Compute total normalized mass in the mantle
    Q2 = (
        sum(
            [
                fractions[i]
                * Mg_number_mantle
                * material_YMg[props["contents"][n][i] - 1]
                for i in range(len(props["contents"][n]))
            ]
        )
        * mMg
        + sum(
            [
                fractions[i]
                * (1.0 - Mg_number_mantle)
                * material_YMg[props["contents"][n][i] - 1]
                for i in range(len(props["contents"][n]))
            ]
        )
        * mFe
        + sum(
            [
                fractions[i] * material_YSi[props["contents"][n][i] - 1]
                for i in range(len(props["contents"][n]))
            ]
        )
        * mSi
        + sum(
            [
                fractions[i] * material_YO[props["contents"][n][i] - 1]
                for i in range(len(props["contents"][n]))
            ]
        )
        * mO
        + sum(
            [
                fractions[i] * material_YH[props["contents"][n][i] - 1]
                for i in range(len(props["contents"][n]))
            ]
        )
        * mH
    )

    Q3 = sum(
        [
            fractions[i]
            * (1.0 - Mg_number_mantle)
            * material_YMg[props["contents"][n][i] - 1]
            for i in range(len(props["contents"][n]))
        ]
    )

    Q4 = props["x_all_core"][0]
    m_core = [mFe, mH, mS, mSi, mO]
    Q5 = sum(
        [m_core[i] * props["x_all_core"][i] for i in range(len(props["x_all_core"]))]
    )

    if not M_IC == None and inner_core_mass_fraction == None:
        # Compute core mass fraction with total inner core mass
        core_frac = 1.0 - 10 ** props["ocean_fraction_should"]
        core_frac *= Q3 / Q2 - Q1 / Q2 * FeMg
        core_frac += M_IC / props["M_surface_should"] * (1.0 / mFe - Q4 / Q5)
        core_frac /= Q3 / Q2 - Q4 / Q5 - FeMg * Q1 / Q2

    elif not inner_core_mass_fraction == None and M_IC == None:
        # Compute core mass fraction with inner core mass fraction
        core_frac = 1.0 - 10 ** props["ocean_fraction_should"]
        core_frac *= Q3 / Q2 - Q1 / Q2 * FeMg
        # core_frac += M_IC / props["M_surface_should"] * (1.0 / mFe - Q4 / Q5)
        core_frac /= (Q3 / Q2 - Q4 / Q5 - FeMg * Q1 / Q2) + inner_core_mass_fraction * (Q4 / Q5 - 1. / mFe)
    else:
        raise ValueError("Both inner_core_mass_fraction and M_IC are passed to function.")
    return core_frac * props["M_surface_should"]
