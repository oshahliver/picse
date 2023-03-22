import numpy as np
from picse.physicalparams import sigmaSB, G, m_earth, r_earth
import sys, os


class Toolkit:
    def __init__(self):
        self.iterator = None

    def get_specs(self, planet):
        return {}

    def evolve(self, planet, iterator=None, iterator_specs={}):
        if not planet.status == "very much alive":
            raise TypeError(
                "You have passed an unconstructed planet object to the iterator"
            )

        else:
            pass

        # Compute time step
        eps = 1e-1

        # Initial conditions
        E = planet.ener_tot_is  # Total energy
        R = planet.R_surface_is  # Radius of the sphere
        T = planet.T_surface_is  # Surface temperature of the sphere
        L_in = 0.0

        # Integration loop
        for i in range(15):
            print("\n########################")
            print(f"Time step {i+1}")
            print("Surface temperature =", planet.T_surface_is)
            print ("Luminostiy =", planet.luminosity)
            # Compute outgoing radiation
            L_out = planet.luminosity

            # Compute RK2 step
            E_0 = E
            L_out_0 = L_out
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
            E_mid = E_0 + (L_in - L_out_0) * dt / 2 / (
                3 / 5 * G * m_earth ** 2 / r_earth
            )

            # Compute radius and temperature at midpoint
            # Here the mass of the planet can change too if mass loss
            # or accretion processes are considered.
            # temporarely supress all prints
            sys.stdout = open(os.devnull, "w")

            iterator_specs.update(
                {"val_should": [planet.M_surface_should * m_earth, E_mid]}
            )
            iterator.iterate(planet=planet, iterator_specs=iterator_specs)
            R_mid = planet.R_surface_is
            T_mid = planet.T_surface_is

            # Compute final energy
            L_out_mid = 4 * np.pi * R_mid ** 2 * sigmaSB * T_mid ** 4
            E += (L_in - L_out_mid) * dt / (3 / 5 * G * m_earth ** 2 / r_earth)

            # Compute new radius and temperature
            iterator_specs.update(
                {"val_should": [planet.M_surface_should * m_earth, E]}
            )
            iterator.iterate(planet=planet, iterator_specs=iterator_specs)
            sys.stdout = sys.__stdout__

            R = planet.R_surface_is
            T = planet.T_surface_is

            if planet.T_surface_is <= 200:
                print("WARNING: Surface temperature too low")
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
