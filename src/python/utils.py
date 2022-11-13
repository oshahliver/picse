from pathlib import Path

def get_project_root() -> Path:
	return Path(__file__).parent.parent.parent


class Parameters():
    def __init__(self, default_values):
        self.default_values = default_values
        self.allowed_keys = self.default_values.keys()

    def set_values(self, **kwargs):
        """ Set planetary parameters to custom values.
        """
        for key, value in kwargs.items():
            if key in self.allowed_keys:
                setattr(self, key, value)

            else:
                raise KeyError("Invalid keyargument <{}> passed".format(key))


    def set_default_values(self):
        """ Set planetary parameters to default values.
        """
        for key, value in self.default_values.items():
            setattr(self, key, value)


class PlanetaryParams(Parameters):
    def __init__(self):
        
        self.default_values = {"Si_number_mantle":0.4,
        "Fe_number_mantle":0.0, 
        "M_surface_should": 1.0,
        "Mg_number_should":0.5,
        "T_surface_should":300.,
        "P_surface_should":1e5,
        "ocean_frac_should":0.,
        "contents":[[2], [2,9,9,9,9], [4,5], [6,7]],
        "fractions":[[1.], [1., 0., 0., 0., 0.], [.5, .5], [.5, .5]],
        "layermasses":[.5, .5, 0., 100.]
        }

        Parameters.__init__(self, self.default_values)            


class RunParams(Parameters):
    def __init__(self, **kwargs):
        self.default_values ={"adiabat_type":0,
        "layer_constraints":[1,1,3,1]}

        Parameters.__init__(self, self.default_values) 