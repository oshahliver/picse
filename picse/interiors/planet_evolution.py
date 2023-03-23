import numpy as np
from picse.physicalparams import sigmaSB, G, m_earth, r_earth, T_zero
from picse.utils.function_tools import functionTools as ftool
import sys, os


class Toolkit:
    def __init__(self):
        self.iterator = None

    def get_specs(self, planet):
        pass

    def evolve(self, planet, iterator=None, iterator_specs={}):
        if not planet.status == "very much alive":
            raise TypeError(
                "You have passed an unconstructed planet object to the iterator"
            )

        else:
            pass

        # The thermal evolution model requires the total mass and the total
        # energy of the planet to be adjusted. in case the planet was
        # created using different boundary conditions, the iterator specs
        # must be updated here for the thermal evolution model.
        iterator_specs.update(
            {
                "what": ["M_surface", "ener_tot"],  # --> target properties
                "how": ["P_center", "T_center"],  # --> adjustable properties
            }
        )

        eps = 2e-1
        t = 0.0
        E = planet.ener_tot_is

        albedo = 0.3
        L_in = 4 * np.pi * r_earth ** 2 * 1366 * (1 - albedo) / 4

        # Define the integrator for the energy balance equation
        # Note. Different models could be implemented here.
        def gradient(y=[], source=0.0, **kwargs):
            """
            Computes rate of change of total energy based on the luminosity of the planet.
            """
            # Update the boundary conditions for the current time step
            iterator_specs.update(
                {"val_should": [planet.M_surface_should * m_earth, y[0]]}
            )
            iterator.iterate(planet=planet, iterator_specs=iterator_specs)
            newgrad = (
                -4
                * np.pi
                * planet.R_surface_is ** 2
                * sigmaSB
                * planet.T_surface_is ** 4
            )
            return [(newgrad + source) / (3 / 5 * G * m_earth ** 2 / r_earth)]

        # Integration loop
        for i in range(25):
            print("\n########################")
            print(f"Time step {i+1}")
            print("Surface temperature =", planet.T_surface_is)

            dt = (
                min(
                    abs(
                        planet.ener_grav
                        * 3
                        / 5
                        * G
                        * m_earth ** 2
                        / r_earth
                        / planet.luminosity
                    ),
                    abs(
                        planet.ener_int
                        * 3
                        / 5
                        * G
                        * m_earth ** 2
                        / r_earth
                        / planet.luminosity
                    ),
                )
                * eps
            )
            print("dt = {:e}".format(dt / (365 * 24 * 3600)))

            sys.stdout = open(os.devnull, "w")
            t, E, = ftool.integrate(
                dydx=gradient,
                y=[E],
                N=1,
                start=t,
                end=t + dt,
                source=L_in,
                order=2,
                noisy=True,
                whicharg="t",
            )
            sys.stdout = sys.__stdout__

            if planet.T_surface_is <= T_zero:
                print(f"WARNING: Lower bound for surface temperature reached at {round(planet.T_surface_is, 3)} K.")
                break

class TimeLine:
    def __init__(self, specs={}):
        self.objects = {"planet": None, "atmosphere": None}

    def add(self, object):
        """Adds an evolvable object to the timeline
        """
        self.objects.update({"planet": object})

    def create(self, planet, specs={}):
        """ Creates the thermal evolution time line for the 
        planet according to the user defined specifications.
        """
        pass
