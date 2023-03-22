import numpy as np

class Toolkit():
    def __init__(self):
        self.iterator = None

    def get_specs(self, planet):
        return {}

    def evolve(self, planet, **kwargs):
        if not planet.status == "very much alive":
            raise TypeError(
                "You have passed an unconstructed planet object to the iterator"
            )

        else:
            pass

        # Set default iterator specifications
        specs = self.get_specs(planet)

        if "iterator_specs" in kwargs:
            for key, val in kwargs["iterator_specs"].items():
                if key in self.iterator_specs_keys:

                    specs.update({key: val})

                else:
                    raise KeyError("Invalid iterator specification given")
        
        

class TimeLine():
    def __init__(self, specs = {}):
        self.objects = {"planet":None, "atmosphere":None}
    
    def add(self, object):
        """Adds an evolvable object to the timeline
        """
        self.objects.update({"planet":object})

    def create(self, planet, specs = {}):
        """ Creates the thermal evolution time line for the 
        planet according to the user defined specifications.
        """
        pass

