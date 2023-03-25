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


def compute_core_mass(props, n=2, M_IC=0.0):
    """Computes the core mass of a planet at given total mass, composition and
    value for Mg#
    """

    # print (contents, Mg_number, M_surface, Mg_number_mantle, SiMg, M_ocean, xi_H_core)
    Mg_number_mantle = min(1.0 - props["Fe_number_mantle"], 0.9999999999)
    FeMg_mantle = (1.0 - Mg_number_mantle) / Mg_number_mantle
    FeMg = (1.0 - props["Mg_number_should"]) / props["Mg_number_should"]
    SiMg = props["Si_number_mantle"] / (1.0 - props["Si_number_mantle"])

    # Compute the fractions in the mantle
    fractions = fortfunctions.functionspy.compute_abundance_vector(
        simg=SiMg,
        femg=FeMg_mantle,
        n_mats=len(props["contents"][n]),
        ymgi=[material_YMg[i - 1] for i in props["contents"][n]],
        ysii=[material_YSi[i - 1] for i in props["contents"][n]],
        xih2oi=[0.0 for i in props["contents"][n]],
        xifei=[1.0 - Mg_number_mantle for i in props["contents"][n]],
        xialsii=[0.0 for i in props["contents"][n]],
        xialmgi=[0.0 for i in props["contents"][n]],
        contents=props["contents"][n],
        additional=[],
    )
    

    # Compute mole fraction of Mg in the mantle
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

# Compute mole fraction of Fe in the mantle
    Q3 = sum(
        [
            fractions[i]
            * (1.0 - Mg_number_mantle)
            * material_YMg[props["contents"][n][i] - 1]
            for i in range(len(props["contents"][n]))
        ]
    )

# Compute mole fraction of Fe in the outer core
    Q4 = props["x_all_core"][0]
    m_core = [mFe, mH, mS, mSi, mO]
    
# Compute total normilized mass in the outer core
    Q5 = sum(
        [m_core[i] * props["x_all_core"][i] for i in range(len(props["x_all_core"]))]
    )
    # print ("xi_all_core in compute core mass:", props["x_all_core"])
    # print ("fractions in mantle =", fractions)
    # print ("Q1, Q2, Q3, Q4, Q5 =", Q1, Q2, Q3, Q4, Q5)
    # print ("Q3/Q2, Q4 / Q5 =", Q3/Q2, Q4/Q5)
    # print ("FeMG =", FeMg)
    # print("inner core mass =", M_IC)
    core_frac = 1.0 - 10 ** props["ocean_fraction_should"]
    # print ("core_frac =", core_frac)
    core_frac *= Q3 / Q2 - Q1 / Q2 * FeMg
    # print ("core_frac =", core_frac)
    core_frac += M_IC / props["M_surface_should"] * (1.0 / mFe - Q4 / Q5)
    # print ("core_frac =", core_frac)
    core_frac /= Q3 / Q2 - Q4 / Q5 - FeMg * Q1 / Q2
    # print ("core_frac =", core_frac)
    return core_frac * props["M_surface_should"]
