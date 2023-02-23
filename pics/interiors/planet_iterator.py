#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:03:10 2019

@author: oshah
"""
# from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# from pics.utils import functionTools as ftool
from matplotlib import pyplot as plt

# import matplotlib.patches as patches
from sklearn.linear_model import LinearRegression
from pics.utils import matrix_generator as mgen
import numpy as np

# from pics.utils import readPREM, logTrans, plotTools
import random
from pics.physicalparams import (
    m_earth,
    r_earth,
    Mg_number_solar,
    M_solar,
    R_solar,
    M_trappist1,
    R_trappist1,
    M_trappist1_error,
    R_trappist1_error,
    names_trappist1,
    names_solar,
    M_others,
    R_others,
    names_others,
    M_others_error,
    R_others_error,
    mH,
    mFe,
    mMg,
    mO,
    mSi,
    NA,
    m_jupiter,
    r_jupiter,
    mH2O,
    mOl,
    G,
    sigmaSB,
)

from pics.runparams import suffix, plot_params
from pics.runparams import grid_color, background_color, color_list

image_loc = "/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/"

plt.rcParams["axes.axisbelow"] = False


class Toolkit:
    def __init__(self, a=0, b=0):
        self.a = a
        self.b = b
        self.bisection = False
        self.delta = 0.0
        self.oldvals = []
        self.iteration = False
        self.number = 0
        self.iterator_specs_keys = [
            "what",
            "how",
            "val_should",
            "predictor",
            "all_val_should_weights",
            "all_howval_weights",
            "acc",
            "iterationLimit",
            "deltaType",
            "unpredictable",
        ]

    def get_specs(self, planet, omit=[]):
        iterator_specs = {
            "what": ["M_surface", "T_surface"],  # --> target properties
            "how": ["P_center", "T_center"],  # --> adjustable properties
            "val_should": [
                planet.M_surface_should * m_earth,
                planet.T_surface_should,
            ],  # --> target values
            "predictor": ["linear", "linear"],  # --> no effect at this point
            "all_val_should_weights": [
                "log",
                "log",
            ],  # --> log or lin extrapolation for targets
            "all_howval_weights": [
                "exp",
                "exp",
            ],  # --> exp or lin prediction for adjustables
            "acc": [1e-3, 1e-2],  # --> desired relative accuracies
            "iterationLimit": 20,  # --> max. number of iterations
            "deltaType": 0,  # --> mode for initial adjustment of adjustables
            "unpredictable": False,  # --> no effect at this point
        }

        if planet.label == "telluric":
            pass

        else:
            pass
            # raise NotImplementedError(
            #     "Iterator does not support base type {} yet".format(planet.label)
            # )

        planet.iterator_specs = iterator_specs
        return iterator_specs

    def all_params(self, planet):
        try:
            all_what = {
                "Mg_number_is": planet.Mg_number_is,
                "T_surface_is": planet.T_surface_is,
                "P_surface_is": planet.P_surface_is,
                "M_surface_is": planet.M_surface_is,
                "M_ocean_is": planet.M_ocean_is,
                "ocean_fraction_is": planet.ocean_fraction_is,
                # "xi_H_core": planet.xi_H_core,
                # "xi_FeO_mantle": planet.xi_FeS_core,
                # "E_tot_is": planet.E_tot_is,
                # "L_int_is": planet.L_int_is,
            }

        # If layer 4 does not exist, no ocean is there
        except IndexError:
            all_what = {
                "Mg_number_is": planet.Mg_number_is,
                "T_surface_is": planet.T_surface_is,
                "P_surface_is": planet.P_surface_is,
                "M_surface_is": planet.M_surface_is,
                "M_ocean_is": 0.0,
                "ocean_fraction_is": -10,
                # "x_H_core": planet.xi_H_core,
                # "xi_FeO_mantle": planet.xi_S_core,
                # "E_tot_is": planet.E_tot_is,
                # "L_int_is": planet.L_int_is,
            }

        except ValueError:
            all_what = {
                "Mg_number_is": planet.Mg_number_is,
                "T_surface_is": planet.T_surface_is,
                "P_surface_is": planet.P_surface_is,
                "M_surface_is": planet.M_surface_is,
                "M_ocean_is": planet.M_ocean_is,
                "ocean_fraction_is": -10,
                # "xi_H_core": planet.xi_H_core,
                # "xi_FeO_mantle": planet.xi_S_core,
                # "E_tot_is": planet.E_tot_is,
                # "L_int_is": planet.L_int_is,
            }
        try:
            all_how = {
                "M_core": planet.layer_masses[1] + planet.layer_masses[0],
                "T_center": planet.T_center,
                "P_center": planet.P_center,
                "M_outer_mantle": planet.layer_properties[3]["indigenous_mass"]
                / m_earth,
            }

        except IndexError:
            try:
                all_how = {
                    "M_core": planet.layer_masses[1] + planet.layer_masses[0],
                    "T_center": planet.T_center,
                    "P_center": planet.P_center,
                    "M_outer_mantle": 0.0,
                }

            except IndexError:
                all_how = {
                    "M_core": planet.layer_masses[0],
                    "T_center": planet.T_center,
                    "P_center": planet.P_center,
                    "M_outer_mantle": 0.0,
                }

        return all_what, all_how

    # def bisect(self, val_is, val_should, acc, direction, predictor=True, log=False):
    #     reldev = (val_should - val_is) / val_should
    #     # if the original criterion is no longer met, that means that the
    #     # probed parameter has been overshoot. Then parameter delta must be
    #     # reduced and the previous step undone

    #     if direction[0] * direction[1] * reldev < -acc:
    #         print("overshoot")
    #         self.iteration = True
    #         if predictor == "linear":
    #             x1 = self.oldhowvals[-2]
    #             x2 = self.oldhowvals[-1]

    #             y1 = self.oldwhatvals[-2]
    #             y2 = self.oldwhatvals[-1]

    #             # Extract values for passive parameter
    #             p1 = [self.oldpassivevals[i][-2] for i in range(len(self.passives))]
    #             p2 = [self.oldpassivevals[i][-1] for i in range(len(self.passives))]

    #             # compute predictor slope for passive parameter
    #             self.passive_slope = [
    #                 (p2[i] - p1[i]) / (x2 - x1) for i in range(len(self.passives))
    #             ]

    #             if log:
    #                 x1 = np.log10(x1)
    #                 x2 = np.log10(x2)
    #                 y1 = np.log10(y1)
    #                 y2 = np.log10(y2)
    #                 val_should = np.log10(val_should)

    #             # compute predictor slope and intercept
    #             slope = (y2 - y1) / (x2 - x1)
    #             self.delta = (val_should - y2) / slope

    #         else:
    #             self.oldhowvals.pop(-1)
    #             self.oldwhatvals.pop(-1)

    #             # perform bisection step
    #             self.delta = self.delta * 0.5
    #             self.bisection = True
    #             print("starting bisection")

    #         # print ('overshoot detected')

    #     elif abs(reldev) <= acc:
    #         self.iteration = False
    #         # print ('desired precission reached')

    #     else:
    #         self.iteration = True

    def iterate(
        self,
        planet=None,
        **kwargs
        #
        # what=["M_surface"],
        # how=["P_center"],
        # val_should=[1.0 * m_earth],
        # acc=[1.0e-3],
        # predictor=["linear"],
        # iterationLimit=25,
        # update_val_should=False,
        # echo=False,
        # updateType=0,
        # unpredictable=True,
        # deltaType=0,
        # write_learning_set=False,
        # start_number=0,
        # log=False,
        # test_data_dir="test_set_00",
        # passives=[],
        # test=False,
        # passive_predictors=[],
        # all_val_should_weights=["lin"],
        # all_howval_weights=["exp"],
        # deltas=None,
    ):
        """Iteratively re-constructs a given planetary object changing a
        specified parameter (how) in each iteration to match a given other
        parameter (what). By default the total mass is probed by iteratively
        adjusting the central pressure of the planet. The desired value
        for the probed parameter is set by the val_should argument and the
        desired precission at which the val_should and val_is match is set by
        the argument acc. This method can be iteratively used to fix several
        surface parameters at once. For instance surface temperature and
        surface pressure can be fixed by first running the iterate methode for
        the temperature and then for the pressure and repeat hthis until the
        precission for both variables is reached. This does, however, sometimes
        require large numbers of iteration steps and of course not all composition
        specifications an total masses can be expected to give a solution for
        the desired surface pressure and temperature.
        Note, if what = P_surface_is or what = T_surface_is, then the
        respective parameter should not be used as majorconstraint in the
        planet.construct() methode to avoid conflicts between the two.

        Note: in order to constrain total mass and surface pressure at the same
        time it is recommended to proceed as follows:

        Use P_surface_is as majorconstraint for the planet integration and
        what = M_surface_is, how = P_center. This methode is very stable and
        efficient.

        In principle it is possible, to use M_surface_is as majorconstraint and
        what = P_surface_is, how = p_center, but as the surface pressure is
        very sensitive to small changes in P_center, the iteration is very
        inefficient and can only  find approximate solutions (i.e. reldev >~
        10%)
        """

        if not planet.status == "very much alive":
            raise TypeError(
                "You have passed an unconstructed planet object to the iterator"
            )

        else:
            pass

        # Set default iterator specifications
        specs = self.get_specs(planet)

        # Update all specifications that are passed manually by the user
        if "iterator_specs" in kwargs:
            for key, val in kwargs["iterator_specs"].items():
                if key in self.iterator_specs_keys:

                    specs.update({key: val})

                else:
                    raise KeyError("Invalid iterator specification given")

        passives = []
        passive_predictors = []
        self.iteration = True
        self.bisection = False
        self.passives = passives
        n_points = 2 ** len(specs["how"])

        for i in range(len(specs["what"])):
            w = specs["what"][i]
            if w == "M_surface" or w == "M_ocean":
                planet.initials[w + "_should"] = specs["val_should"][i] / m_earth

            elif w == "E_tot":
                planet.initials[w + "_should"] = specs["val_should"][i] / (
                    G * m_earth**2 / r_earth * 3 / 5
                )

            elif w == "L_int":
                planet.initials[w + "_should"] = specs["val_should"][i] / (
                    4 * np.pi * r_earth**2 * sigmaSB * 300**4
                )

            else:
                planet.initials[w + "_should"] = specs["val_should"][i]

        # print ('planet fractions in iterate =', planet.fractions)
        # perform first construction to determine iteration direction for
        # bisection algorithm
        self.oldhowvals = []
        self.oldwhatvals = []
        self.oldshouldvals = []
        self.oldpassivevals = [[] for i in passives]
        self.oldpassivehowvals = [[] for i in passives]
        self.delta = [None for i in specs["what"]]
        newval = [None for i in specs["what"]]

        # generate dict containing all possible what and how parameters
        all_what, all_how = self.all_params(planet)

        sanity_borders = {
            "M_core": [1.0e-10, 10.0],  # in earth masses
            "T_center": [100.0, 25000.0],  # in K
            "P_center": [1.0e7, 5.0e13],  # in Pa
            "M_outer_mantle": [1.0e-6, 10.0],  # in earth masses
        }

        # these values are fixed based on runtime experience and have not been
        # optimized by any means
        if specs["deltaType"] == 0:
            initial_deltas = {
                "M_core": 0.1,
                "T_center": 0.25,
                "P_center": 0.25,
                "M_outer_mantle": 0.1,
            }

        elif specs["deltaType"] == 1:
            initial_deltas = {
                "M_core": 0.01,
                "T_center": 0.05,
                "P_center": 0.05,
                "M_outer_mantle": 0.01,
            }

        self.oldhowvals.append([all_how[h] for h in specs["how"]])
        self.oldwhatvals.append([all_what[w + "_is"] for w in specs["what"]])
        self.oldshouldvals.append([v for v in specs["val_should"]])

        for i in range(len(passives)):
            if passives[i] == "xi_H_core" or passives[i] == "xi_FeO_mantle":
                self.oldpassivevals[i].append(all_what[passives[i]])
                self.oldpassivehowvals[i].append(None)

            else:
                self.oldpassivevals[i].append(all_what[passives[i] + "_is"])
                self.oldpassivehowvals[i].append(all_how["P_center"])

        self.passive_slope = [0.0 for i in passives]

        val_is = [all_what[w + "_is"] for w in specs["what"]]

        reldev = [
            (specs["val_should"][i] - val_is[i]) / specs["val_should"][i]
            for i in range(len(specs["what"]))
        ]
        direction = [[0, 0] for w in specs["what"]]

        self.already_met = [0 for i in range(len(specs["what"]))]
        print("val_should =", specs["val_should"])
        print("initial howval =", [all_how[h] for h in specs["how"]])
        print("initial whatval =", [all_what[w + "_is"] for w in specs["what"]])

        # print("initial passiveval is=", self.oldpassivevals)

        # Construct initial grid points for prediction
        # First grid point is input planet.
        # Second grid point is both howvals predicted
        # Third and forth grid point is only one howval predicted
        for p in range(2 ** (len(specs["how"])) - 1):
            # Predict both howvals
            for i in range(len(specs["what"])):
                print("\nProbing ", specs["what"][i])
                if specs["what"][i] == "Mg_number":
                    direction[i] = [-1, None]
                    if specs["how"][i] == "M_core":

                        # initial core mass guess too large -> iteration direction down
                        if direction[i][0] * reldev[i] < -specs["acc"][i]:
                            direction[i][1] = -1
                            self.delta[i] = -initial_deltas[specs["how"][i]] * (
                                planet.layer_masses[1] + planet.layer_masses[0]
                            )

                        # initial core mass guess too small -> iteration direction up
                        elif direction[i][0] * reldev[i] > specs["acc"][i]:
                            direction[i][1] = 1
                            self.delta[i] = initial_deltas[specs["how"][i]] * (
                                planet.layer_masses[1] + planet.layer_masses[0]
                            )

                        # accuracy already met
                        else:
                            print("Condition already met")
                            print("reldev =", reldev[i])
                            self.delta[i] = (
                                initial_deltas[specs["how"][i]]
                                * all_how[specs["how"][i]]
                                / 2.0
                            )
                            self.already_met[i] = 1

                    elif specs["how"][i] == "M_outer_mantle":
                        # initial outer mantle mass guess too small -> iteration direction up
                        if direction[i][0] * reldev[i] < -specs["acc"][i]:
                            direction[i][1] = -1
                            self.delta[i] = (
                                initial_deltas[specs["how"][i]] * planet.layer_masses[1]
                            )

                        # initial outer mantle mass guess too large -> iteration direction down
                        elif direction[i][0] * reldev[i] > specs["acc"][i]:
                            direction[i][1] = 1
                            self.delta[i] = (
                                -initial_deltas[specs["how"][i]]
                                * planet.layer_masses[1]
                            )

                        # accuracy already met
                        else:
                            # print ('condition already met')
                            self.delta[i] = 0.0
                            self.already_met[i] = 1

                    elif specs["how"][i] == "P_center":
                        # Mg# too small
                        # initial central pressure too low -> iteration direction up
                        if direction[i][0] * reldev[i] < -specs["acc"][i]:
                            direction[i][1] = -1
                            self.delta[i] = (
                                initial_deltas[specs["how"][i]] * planet.P_center
                            )

                        elif direction[i][0] * reldev[i] > specs["acc"][i]:
                            direction[i][1] = 1
                            self.delta[i] = (
                                -initial_deltas[specs["how"][i]] * planet.P_center
                            )

                        else:
                            self.delta[i] = 0.0
                            self.already_met[i] = 1

                elif (
                    specs["what"][i] == "T_surface"
                    or specs["what"][i] == "P_surface"
                    or specs["what"][i] == "M_surface"
                    or specs["what"][i] == "M_ocean"
                    or specs["what"][i] == "ocean_frac"
                    or specs["what"][i] == "E_tot"
                    or specs["what"][i] == "L_int"
                ):
                    # State dependency of target on probed parameter
                    # negative means target increases with increasing parameter
                    if (
                        specs["what"][i] == "M_ocean"
                        or specs["what"][i] == "ocean_frac"
                    ):
                        direction[i] = [-1, None]

                    else:
                        direction[i] = [-1, None]

                    # initial guess for central value too low
                    if direction[i][0] * reldev[i] < -specs["acc"][i]:
                        direction[i][1] = -1
                        self.delta[i] = (
                            initial_deltas[specs["how"][i]] * all_how[specs["how"][i]]
                        )

                    # initial guess for central value too high
                    elif direction[i][0] * reldev[i] > specs["acc"][i]:
                        direction[i][1] = 1
                        self.delta[i] = (
                            -initial_deltas[specs["how"][i]] * all_how[specs["how"][i]]
                        )

                    # accuracy already met
                    else:
                        print("Condition already met")
                        print("reldev =", reldev[i])
                        self.delta[i] = (
                            initial_deltas[specs["how"][i]]
                            * all_how[specs["how"][i]]
                            / 2.0
                        )
                        self.already_met[i] = 1

                    print("initial delta =", self.delta[i])

                # print("all_how = ", all_how)
                # print("how =", specs["how"])
                all_how[specs["how"][i]] += self.delta[i]

                # Add random fluctuation to avoid singular matrix if the initial
                # values of two data points are identical
                rand = random.random()
                # print("random fluc =", rand)
                all_how[specs["how"][i]] *= 1.0 + rand * 1e-2

                # force newval to stay within the defined value ranges for the
                # currently iterated parameter
                newval[i] = min(
                    all_how[specs["how"][i]], sanity_borders[specs["how"][i]][1]
                )
                newval[i] = max(newval[i], sanity_borders[specs["how"][i]][0])

                if abs(reldev[i]) <= specs["acc"][i]:
                    print(
                        "\n --> Desired precission for {} reached after initial iteration.".format(
                            specs["what"][i]
                        )
                    )

                # M_core has to be treated seperately here as it is not an attribute
                # of the Planet class but rather the first entry of layer_masses of a
                # Planet.Planet object. Also the planets properties only need to be
                # updated if the condition is NOT already met
                if specs["how"][i] == "M_core":
                    planet.initials["layer_masses"][0] = newval[i]
                    planet.initials["layer_masses"][1] = newval[
                        i
                    ]  # * (1. - planet.initials['inner_core_frac'])

                elif specs["how"][i] == "M_outer_mantle":
                    planet.initials["layer_masses"][3] = newval[i]

                else:
                    planet.initials[specs["how"][i]] = newval[i]
            """
            print ('vals should =', val_should)
            print ('vals is =', val_is)
            print ('initial howval =', [all_how[h] for h in how])
            print ('initial whatval =', [all_what[w+'_is'] for w in what])
            print ('initial reldev =', reldev)
            print ('initial deltas =', self.delta)
            """

            if sum(self.delta) == 0.0:
                print("All parameters already satisfied")
                self.iteration = False

            # For the last grid point don't perform integration but only update
            # nitial values. The integration will be performed in the first
            # iteration below.
            if p + 1 < 2 ** (len(specs["how"])) - 1:
                print("\n creating additional grid point {}".format(p + 1))
                planet.reset()
                planet.update_initials()
                planet.construct()

                all_what, all_how = self.all_params(planet)
                val_is = [all_what[w + "_is"] for w in specs["what"]]
                reldev = [
                    (specs["val_should"][i] - val_is[i]) / specs["val_should"][i]
                    for i in range(len(specs["what"]))
                ]

                self.oldhowvals.append([n for n in newval])
                self.oldwhatvals.append([v for v in val_is])

                # print("initial howvals =", self.oldhowvals)
                # print("initial whatvals =", self.oldwhatvals)
                # print("initial newvals =", newval)

        count = 0
        exitcode = 0
        while self.iteration:
            count += 1
            print("\n iteration", count)
            # print ('layer_masses =', planet.layer_masses)
            planet.reset()
            # print ('layer_masses =', planet.layer_masses)
            planet.update_initials()
            # print ('layer_masses =', planet.layer_masses)
            planet.construct()
            # print ('layer_masses =', planet.layer_masses)
            all_what, all_how = self.all_params(planet)
            val_is = [all_what[w + "_is"] for w in specs["what"]]
            reldev = [
                (specs["val_should"][i] - val_is[i]) / specs["val_should"][i]
                for i in range(len(specs["how"]))
            ]

            if len(self.passives) > 0:
                self.passive_slope = []

                for i in range(len(passives)):
                    if passives[i] == "xi_H_core" or passives[i] == "xi_FeO_mantle":
                        self.oldpassivevals[i].append(all_what[passives[i]])
                        self.oldpassivehowvals[i].append(None)

                    else:
                        self.oldpassivevals[i].append(all_what[passives[i] + "_is"])
                        self.oldpassivehowvals[i].append(all_how["P_center"])

                    # Extract dependant parameter for passive parameter
                    x1 = self.oldhowvals[-2][passive_predictors[i]]
                    x2 = self.oldhowvals[-1][passive_predictors[i]]

                    # Extract values for passive parameter
                    p1 = self.oldpassivevals[i][-2]
                    p2 = self.oldpassivevals[i][-1]

                    if np.isnan(p2):
                        p2 = p1

                    # compute predictor slope for passive parameter
                    try:
                        self.passive_slope.append((p2 - p1) / (x2 - x1))

                    except ZeroDivisionError:
                        self.passive_slope.append(0.0)

            # print ('x1, x2 =', x1, x2)
            # print ('p1, p2 =', p1, p2)
            # print ('passiveslope =', self.passive_slope)

            self.oldhowvals.append([n for n in newval])
            self.oldwhatvals.append([v for v in val_is])

            # print ('howvals =', self.oldhowvals)
            # print ('whatvals =', self.oldwhatvals)

            # Size of coefficient vector for multi-linear predictor
            n_coeffs = 2 ** len(specs["how"])

            # print("oldwhatvals =", self.oldwhatvals)
            linsys = mgen.LinearSystem(
                self.oldwhatvals[-(2 ** len(specs["how"])) :],
                self.oldhowvals[-(2 ** len(specs["how"])) :],
                target=specs["val_should"],
                param_weights=specs["all_val_should_weights"],
                data_weights=specs["all_howval_weights"],
            )

            # print("body =", linsys.body)
            # print("vec =", linsys.vec)

            try:
                # Solve linear system to get coefficients a_ij
                # M_ki a_ij = x_ki
                # k: number of whatvals used for predictor
                # i: number of terms in predictor step (depends on model)
                # j: number of howvals which have to be predicted
                # a = np.linalg.solve(matrix, x)
                linsys.create_model()
                linsys.predict()
                newval = linsys.pred

            except np.linalg.LinAlgError as err:
                if "Singular matrix" in str(err):
                    print("Singular matrix")
                    """The reason for a singular matrix can be the fact that
                    one of the parameters is already met in which case the
                    delta is zero and hence y1 = y2. In this case adopt 
                    separate linear extrapolation for the other parameters
                    keeping the ones that are already met constant.
                    """
                    ind = []

                    # gather indices for parameters that are not already met
                    for i in range(len(self.delta)):
                        if self.delta[i] != 0.0:
                            ind.append(i)
                            print("i =", i, "not met")
                            # perform sperate linear extrapolation for those parameters
                            # while keeping the others constant
                            x1 = self.oldhowvals[-2][i]
                            x2 = self.oldhowvals[-1][i]

                            y1 = self.oldwhatvals[-2][i]
                            y2 = self.oldwhatvals[-1][i]

                            try:
                                slope = (y2 - y1) / (x2 - x1)
                                self.delta[i] = (specs["val_should"][i] - y2) / slope

                            # In case the iteration has reached a bottleneck where
                            # it does not proceed anymore, the delta will be zero zero and
                            # hence the new slope aswell. In this case just finish the
                            # iteration and try to match the current parameter in the
                            # next sweep
                            except ZeroDivisionError:
                                pass

                            newval[i] = max(
                                self.oldhowvals[-1][i] + self.delta[i],
                                sanity_borders[specs["how"][i]][0],
                            )

                        else:
                            newval[i] = x2

                else:
                    pass

            # print("newvals =", newval)
            # print ('pred =', linsys.pred)
            for i in range(len(specs["how"])):
                # newval[i] = min(newval[i], sanity_borders[how[i]][1])
                # newval[i] = max(newval[i], sanity_borders[how[i]][0])
                if (
                    newval[i] >= sanity_borders[specs["how"][i]][1]
                    or newval[i] <= sanity_borders[specs["how"][i]][0]
                ):
                    self.iteration = False
                    exitcode = 1
                    print(
                        "WARNING: Ceiling for parameter "
                        + specs["how"][i]
                        + " reached."
                    )
                    print("vapor =", planet.vapor_reached)
                    print("T / P =", planet.T_surface_is, planet.P_surface_is * 1e-5)
            # print ('a =', a)
            if len(self.passives) > 0:
                self.passive_delta = [
                    max(
                        (
                            newval[passive_predictors[i]]
                            - self.oldhowvals[-1][passive_predictors[i]]
                        )
                        * self.passive_slope[i],
                        1.0e-10,
                    )
                    for i in range(len(passives))
                ]

                for p in range(len(passives)):
                    if not self.passive_slope[p] == 0.0 and not np.isnan(
                        self.passive_slope[p]
                    ):
                        newpassiveval = [
                            self.oldpassivevals[i][-1] + self.passive_delta[i]
                            for i in range(len(passives))
                        ]

                    else:
                        newpassiveval = [0.0 for i in range(len(passives))]

            # Sometimes iteration enters a bottle neck. In this case abort!
            for i in range(len(specs["how"])):
                if (
                    newval[i] == self.oldhowvals[-1][i]
                    and self.oldhowvals[-1][i] == self.oldhowvals[-2][i]
                ):
                    self.iteration = False
                    exitcode = 2
                    print(
                        "WARNING: Bottleneck for"
                        + specs["how"][i]
                        + " reached. Aborting iteration..."
                    )

                if np.isnan(newval[i]):
                    self.iteration = False
                    print(
                        "WARNING: " + specs["how"][i] + " is NaN. Aborting iteration..."
                    )
                    exitcode = 3
            """          
            try:
                print ('newpassiveval =', newpassiveval)
                print ('passivedelta =', self.passive_delta)
                
            except UnboundLocalError:
                pass
            """

            print("new how =", newval)
            print("new what =", val_is)
            print("val_should =", specs["val_should"])
            print("reldev =", reldev)

            for i in range(len(specs["how"])):
                if specs["how"][i] == "M_core":
                    planet.initials["layer_masses"][0] = newval[i]
                    planet.initials["layer_masses"][1] = newval[
                        i
                    ]  # * (1. - planet.initials['inner_core_frac'])
                    print("setting IC mass to:", newval[i])
                elif specs["how"][i] == "M_outer_mantle":
                    planet.initials["layer_masses"][3] = newval[i]

                else:
                    planet.initials[specs["how"][i]] = newval[i]

            for i in range(len(passives)):
                if passives[i] == "xi_H_core" or passives[i] == "xi_FeO_mantle":
                    planet.initials[passives[i] + "_predicted"] = newpassiveval[i]

                else:
                    pass
                    # planet.initials[passives[i]] = newpassiveval[i]

            self.already_met = [0 for i in range(len(specs["how"]))]
            acc_reached_count = 0
            for h in range(len(specs["how"])):
                if abs(reldev[h]) < specs["acc"][h]:
                    acc_reached_count += 1
                    self.already_met[h] = 1

                else:
                    pass

            if acc_reached_count == len(specs["how"]):
                print(
                    "\n --> Desired precission for all parameters reached after {} iterations!".format(
                        count
                    )
                )
                print("---------------------")
                print("relative deviations =", reldev)
                print("desired accuracies = {}".format(specs["acc"]))
                print("---------------------")
                self.iteration = False

            if count >= specs["iterationLimit"]:
                self.iteration = False
                print("WARNING: iteration limit reached after", count, "iterations !")

        return exitcode
