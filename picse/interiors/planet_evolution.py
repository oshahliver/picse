import numpy as np
from picse.physicalparams import sigmaSB, G, m_earth, r_earth, T_zero, year, day
from picse.utils.function_tools import functionTools as ftool
import sys, os
import copy
from picse.utils.file_tools import internal_data
from alive_progress import alive_bar


class Toolkit:
    def __init__(self):
        self.iterator = None

    def get_specs(self, planet):
        self.evolver_specs = {
            "order": 2,
            "start": 0,
            "end": 1e7,
            "acc": 1e-5,
            "tag": "",
            "eps": 0.25,
            "write":False,
            "source":0,
            "out_path":os.getcwd()
        }

    def evolve(
        self,
        planet,
        iterator=None,
        iterator_specs={},
        evolver_specs = {},
    
    ):
        if not planet.status == "very much alive":
            raise TypeError(
                "You have passed an unconstructed or spurious planet object to the iterator"
            )

        else:
            pass

        # Set default iterator specifications
        self.get_specs(planet)

        # Update all specifications that are passed manually by the user
        for key, val in evolver_specs.items():
            if key in self.evolver_specs.keys():
                self.evolver_specs.update({key: val})

            else:
                raise KeyError(f"Invalid evolver specification <{key}> given")

        # The thermal evolution model requires the total mass and the total
        # energy of the planet to be adjusted. in case the planet was
        # created using different boundary conditions, the iterator specs
        # must be updated here for the thermal evolution model.
        new_iterator_specs = copy.deepcopy(iterator_specs)
        new_iterator_specs.update(
            {
                "what": ["M_surface", "ener_tot"],  # --> target properties
                "how": ["P_center", "T_center"],  # --> adjustable properties
                "all_val_should_weights":["log", "lin"], # --> energy must be lin!
                "all_howval_weights":["exp", "lin"], # --> energy must be lin!
            }
        )

        if not new_iterator_specs["what"][0] == "M_surface":
            raise ValueError ("For thermal evolution models the first target property must be the total mass.")

        t = self.evolver_specs["start"]
        y = [planet.ener_tot_is]

        iterate = True
        iteration_count = 0

        # print(f"\nInitial state:")
        # print("Surface temperature (K) =", planet.T_surface_is)
        # print("Radius (r_earth) = {:e}".format(round(planet.R_surface_is / r_earth, 3)))
        # print ("time (yr) = {:e}".format(t / (year*day)))

        if self.evolver_specs["write"]:
            data = {key:[] for key, val in internal_data.labels.items()}
            internal_data.add_to_timeline(planet, data)
            data.update({"time":[self.evolver_specs["start"]]})
        
        with alive_bar(
            manual=True,
            title=f"Creating time line {self.evolver_specs['tag']}",
            bar="bubbles",
            spinner="pulse",
            dual_line=True,
        ) as bar:
            # Define the integrator for the energy balance equation
            # Note. Different models could be implemented here.
            def gradient(y=[], source=0.0, **kwargs):
                """
                Computes rate of change of total energy based on the luminosity of the planet.
                """
                # Update the boundary conditions for the current time step
                new_iterator_specs.update(
                    {"val_should": [planet.M_surface_should * m_earth, y[0]]}
                )
                try:
                    old_text = kwargs["same_text"]
                except KeyError:
                    old_text = ""

                iterator.iterate(
                    planet=planet,
                    iterator_specs=new_iterator_specs,
                    progress_bar=bar,
                    old_text=old_text,
                )
                newgrad = (
                    -4
                    * np.pi
                    * planet.R_surface_is ** 2
                    * sigmaSB
                    * planet.T_surface_is ** 4 * .2
                )

                return [(newgrad + source) / (3 / 5 * G * m_earth ** 2 / r_earth)]

            # Integration loop
            while iterate:
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
                    * self.evolver_specs["eps"]
                )

                if t < self.evolver_specs["end"]:
                    dt = min(dt, self.evolver_specs["end"] - t)

                else:
                    pass

                sys.stdout = open(os.devnull, "w")
                text = f"Time step: {iteration_count + 1} (Iteration {0})"

                # Perform the next time integration step
                t, y = ftool.integrate(
                    dydx=gradient,
                    y=y,
                    start=t,
                    end=t + dt,
                    source=self.evolver_specs["source"],
                    order=self.evolver_specs["order"],
                    noisy=False,
                    whicharg="t",
                    same_text=text,
                )

                bar.text = text
                if self.evolver_specs["write"]:
                    internal_data.add_to_timeline(planet, data)
                    data["time"].append(t)

                if planet.T_surface_is <= T_zero:
                    print(
                        f"\nWARNING: Lower bound for surface temperature reached at {round(planet.T_surface_is, 3)} K."
                    )
                    break

                # Check if end time is reached
                if not self.evolver_specs["end"] is None:
                    reldev = (t - self.evolver_specs["end"]) / self.evolver_specs["end"]
                    # End time reached
                    if abs(reldev) < self.evolver_specs["acc"]:
                        iterate = False

                iteration_count += 1

                sys.stdout = sys.__stdout__
                bar(((t - self.evolver_specs["start"]) / (self.evolver_specs["end"] - self.evolver_specs["start"])))
        
        if self.evolver_specs["write"]:
            return data
        else:
            return 

class TimeLine:
    def __init__(self, specs={}):
        self.objects = {"planet": None}

    def get_specs(self):
        pass

    def create(self, planet, planetary_specs={}, iterator_specs={}):
        """ Creates the thermal evolution time line for the 
        planet according to the user defined specifications.
        """
        
        pass
