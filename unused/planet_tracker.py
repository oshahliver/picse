from picse.interiors import planet_creator
from picse.interiors import planet_iterator
import numpy as np
from picse.physicalparams import sigmaSB


def energy_loss(temperature, radius):
    """Compute the energy lost by the planet per unit time
    due to black body radiation.

    Parameters:
    -----------
    luminosity : float
        The surface luminosity of the planet.
    radius : float
        The radius of the planet.

    Returns:
    --------
    float
        The energy lost by the planet per unit time.
    """
    return 4 * np.pi * radius ** 2 * sigmaSB * temperature ** 4


class Engine:
    def __init__(self):
        pass


class BaseEngine(Engine):
    def __init__(self):
        Engine.__init__(self)

    def run(self):
        pass

    def stop(self):
        pass


class GravityEngine(Engine):
    def __init__(self):
        Engine.__init__(self)

    def run(self):
        pass

    def stop(self):
        pass


class ThermoEngine(Engine):
    def __init__(self):
        Engine.__init__(self)

    def driver(
        self,
        planet,
        iterator,
        total_mass,
        bulk_composition,
        surface_temperature,
        surface_pressure,
        timestep,
        final_time,
    ):
        """Evolve a planet in time using the energy balance model.

        Parameters:
        -----------
        total_mass : float
            The total mass of the planet in kg.
        bulk_composition : dict
            A dictionary with the mass fractions of each component
            (e.g. {"Fe": 0.3, "Mg": 0.4, "Si": 0.3}).
        surface_temperature : float
            The initial surface temperature of the planet in K.
        surface_pressure : float
            The initial surface pressure of the planet in Pa.
        surface_luminosity : float
            The surface luminosity of the planet in W.
        timestep : float
            The timestep for the simulation in years.
        final_time : float
            The final time for the simulation in years.

        Returns:
        --------
        numpy.ndarray
            An array with the time values.
        numpy.ndarray
            An array with the surface temperature values over time.
        """
        # Compute the initial planetary structure
        boundary_conditions = {
            "total_mass": total_mass,
            "bulk_composition": bulk_composition,
            "surface_temperature": surface_temperature,
            "surface_pressure": surface_pressure,
        }

        radius, temperature = planet.R_surface_is, planet.T_surface_is
        # Initialize the time and temperature arrays
        time_array = np.arange(0, final_time + timestep, timestep)
        temperature_array = np.zeros_like(time_array)
        temperature_array[0] = surface_temperature

        # Compute the energy lost by the planet per unit time
        energy_lost = energy_loss(surface_temperature, radius)

        # Time loop
        for i in range(1, len(time_array)):
            # Compute the new surface temperature using the energy balance equation
            new_surface_temperature = temperature[-1] - energy_lost * timestep / (
                4 * np.pi * radius ** 2 * sigmaSB
            )

            # Update the surface temperature and store it in the temperature array
            temperature_array[i] = new_surface_temperature

            # Compute the new planetary structure
            planet = planet_creator.TelluricPlanet()
            planet.construct()
            iterator.iterate(planet=planet)

            radius = planet.R_surface_is

        # Use vectorization to compute the energy loss for all time steps at once
        energy_lost_array = np.full_like(temperature_array, energy_lost)
        energy_loss_time = (
            energy_lost_array * timestep / (4 * np.pi * radius ** 2 * sigmaSB)
        )
        temperature_array[1:] -= energy_loss_time[:-1]

        return time_array, temperature_array

    def run(self, planet, time_series):
        t = time_series.time_start

        # call the driver on the planet
        self.driver(planet)

    def stop(self):
        pass


class TimeSeries:
    def __init__(self, specs={}):
        self.time_start = 0.0

    def set_up(self, t1, eps_t=0.5):
        self.time_end = t1
        self.eps_t = eps_t


class PlanetaryHistory:
    def __init__(self):
        # List to store all thermal trajectories for a planet
        self.trajectories = []

    def set_up(self, planetary_specs={}, evolution_specs=[], iterator_specs={}):
        """Set up the planetary properties and the evolution track specifications"""
        # Create an instance of the iterator toolkit
        self.iterator = planet_iterator.Toolkit()

        self.trajectories = [TimeSeries(specs=evs) for evs in evolution_specs]
        self.ancestor = planet_creator.TelluricPlanet(**planetary_specs)
        self.ancestor.construct()
        self.iterator.iterate(planet=self.ancestor, **iterator_specs)
        self.ancestor.check_convergence()

        if self.ancestor.converged:
            pass

        else:
            print(r"Warning: there was an issue with the creation of the ancestor.")

    def create(self, engine):

        for ts in self.trajectories:
            engine.run(self.ancestor, ts)

    def add_trajectory(self, ts):
        self.trajectories.append(ts)

    def remove_trajectory(self, tag):
        pass
