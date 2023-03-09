from pics.utils import test_interface as ti
import numpy as np
import inspect

# List of input, output and other variable names specified elsewhere in the code
input_keys = ["a", "b", "c"]
output_keys = ["d", "e", "f", "x", "y", "z"]
other_keys = ["v1", "v2", "v3"]


class MySim:
    def __init__(self, a=1.0, b=2.0):
        # Define input parameters
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)

        for key in args:
            setattr(self, key, values[key])

        # Do some calculations to define some more input parameters
        self.c = self.a + self.b

        # Define output parameters from the sim
        for key in output_keys:
            setattr(self, key, 0.0)

        # Define parameters calculated from the output parameters of the sim
        for key in other_keys:
            setattr(self, key, 0.0)

    def run(self):
        # Gather input variables in dict
        kwargs = dict((name, getattr(self, name)) for name in input_keys)

        # Run the actual simulation
        output = ti.interface.do_some_science_stuff(**kwargs)

        # Assign the output variables from the simulation
        self.d, self.e, self.f, self.x, self.y, self.z = output  # --> this is ugly!!!

    def update(self):
        self.v1 = 2 * self.d
        self.v2 = self.e + self.d**2
        self.v3 = np.sqrt(self.f) + self.v1


sim = MySim(a=5.5, b=2.3)
sim.run()
sim.update()

print("The simulation results are:")
print([sim.__dict__[o] for o in output_keys])
print([sim.__dict__[o] for o in other_keys])
