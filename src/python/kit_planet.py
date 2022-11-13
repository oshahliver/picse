# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import copy


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
    def __init__(self, type = "telluric", **kwargs):
        
        if type == "telluric":
            self.default_values = {"Si_number_mantle":0.4,
            "Fe_number_mantle":0.0, 
            "M_surface_should": 1.0,
            "Mg_number_should":0.5,
            "T_surface_should":300.,
            "P_surface_should":1e5,
            "ocean_fraction_should":0.,
            "contents":[[2], [2,9,9,9,9], [4,5], [6,7]],
            "fractions":[[1.], [1., 0., 0., 0., 0.], [.5, .5], [.5, .5]],
            "layermasses":[.5, .5, 0., 100.]
            }
        
        
        elif type == "aqua":
            self.default_values = {"Si_number_mantle":0.4,
            "Fe_number_mantle":0.0, 
            "M_surface_should": 1.0,
            "Mg_number_should":0.5,
            "T_surface_should":300.,
            "P_surface_should":1e5,
            "ocean_fraction_should":0.1,
            "contents":[[2], [2,9,9,9,9], [4,5], [6,7], [1]],
            "fractions":[[1.], [1., 0., 0., 0., 0.], [.5, .5], [.5, .5], [1.]],
            "layermasses":[.5, .5, 0., 0., 100.]
            }


        elif type == "inferno":
            self.default_values = {"Si_number_mantle":0.4,
            "Fe_number_mantle":0.0, 
            "M_surface_should": 1.0,
            "Mg_number_should":0.5,
            "T_surface_should":1000.,
            "P_surface_should":1e5,
            "ocean_fraction_should":0.,
            "contents":[[2], [2,9,9,9,9], [4,5], [6,7]],
            "fractions":[[1.], [1., 0., 0., 0., 0.], [.5, .5], [.5, .5]],
            "layermasses":[.5, .5, 0., 100.]
            }

        elif type == "ice":
            self.default_values = {"Si_number_mantle":0.4,
            "Fe_number_mantle":0.0, 
            "M_surface_should": 1.0,
            "Mg_number_should":0.5,
            "T_surface_should":200.,
            "P_surface_should":1e5,
            "ocean_fraction_should":0.1,
            "contents":[[2], [2,9,9,9,9], [4,5], [6,7]],
            "fractions":[[1.], [1., 0., 0., 0., 0.], [.5, .5], [.5, .5]],
            "layermasses":[.5, .5, 0., 100.]
            }
        
        else:
            
            if not "default_values" in kwargs:
                raise KeyError("If you don't specify the planet type you must provide all parameters!")
            
            else:
                self.default_values = {"Si_number_mantle":0.4,
                "Fe_number_mantle":0.0, 
                "M_surface_should": 1.0,
                "Mg_number_should":0.5,
                "T_surface_should":300.,
                "P_surface_should":1e5,
                "ocean_fraction_should":0.,
                "contents":[[2], [2,9,9,9,9], [4,5], [6,7]],
                "fractions":[[1.], [1., 0., 0., 0., 0.], [.5, .5], [.5, .5]],
                "layermasses":[.5, .5, 0., 100.]
                }
                
                self.default_values.update(kwargs["default_values"])

        Parameters.__init__(self, self.default_values)


class RunParams(Parameters):
    def __init__(self, **kwargs):
        self.default_values ={"adiabat_type":0,
        "layer_constraints":[1,1,3,1]}

        Parameters.__init__(self, self.default_values) 


class AllInputParams(RunParams, PlanetaryParams):
    def __init__(self):
        pass
    

class Planet():
    def __init__(self, label = "A random planet", **kwargs):
        #self.default_values = {}
        #self.default_values.update(planetary_params.default_values)
        #self.default_values.update(run_params.default_values)
        self.default_values = None
        self.label = label
       
    def set_values(self, default = False, **kwargs): 
        omit_keys = ["default_values", "allowed_keys", "label"]
        
        if "planetary_params" in kwargs:
            for key, value in kwargs["planetary_params"].__dict__.items():
                if not key in omit_keys:
                    setattr(self, key, value)
                    
        if "run_params" in kwargs:
            for key, value in kwargs["run_params"].__dict__.items():
                if not key in omit_keys:
                    setattr(self, key, value)
        
        # Set the given values as the default for this planetary object
        if default:
            self.default_values = copy.deepcopy(self.__dict__)
    
    
    def update_values(self):
        self.default_values = copy.deepcopy(self.__dict__)
    
    def update_core(self):
        pass
    
    def update_composition(self):
        pass
    
    def update_bulk(self):
        pass
    
    def update_mantle(self):
        pass

    def update_hydrosphere(self):
        pass

    def update_atmosphere(self):
        pass

    def update_orbit(self):
        pass

    def update(self, default = False):
        """ Computes all dependant planetary parameters
        """
        self.update_composition()
        self.update_bulk()
        self.update_core()
        self.update_mantle()
        self.update_hydrosphere()

        # Set the given values as the default for this planetary object
        if default:
            self.default_values = copy.deepcopy(self.__dict__)
    
    
    def reset(self):
        """ Resets all planetary parameters to the values stored in default_values
        """
        for key, value in self.__dict__.items():
            if not key == "default_values":
                setattr(self, key, self.default_values[key])
                
                
    def show(self, style = 0):
        """ Prints out a simple overview of all relevant planetary parameters.
        """
        show_parameters = {"Mg_number_should", "ocean_fraction_should", "T_surface_should", "Fe_number_mantle"}
    
        print("\n####################\nPlanet Label: {}\n".format(self.label))
        
        if style == 0:
            for sp in show_parameters:
                print ("{a} = {b}".format(a = sp, b = self.__dict__[sp]))
                
    
    def dump(self, traget):
        pass
    
    
    def write(self, target, format = "csv"):
        pass
    
    
    def load(self, file):
        pass
    
    
    def construct(self):
        pass
     

class TelluricPlanet(Planet):
    def __init__(self):
        Planet.__init__(self, label = "telluric")
        pp = PlanetaryParams(type = self.label)
        rp = RunParams(type = self.label)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params = pp, run_params = rp, default = True)        


class AquaPlanet(Planet):
    def __init__(self):
        Planet.__init__(self, label = "aqua")
        pp = PlanetaryParams(type = self.label)
        rp = RunParams(type = self.label)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params = pp, run_params = rp, default = True)       


class InfernoPlanet(Planet):
    def __init__(self):
        Planet.__init__(self, label = "inferno")
        pp = PlanetaryParams(type = self.label)
        rp = RunParams(type = self.label)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params = pp, run_params = rp, default = True)   


class IcePlanet(Planet):
    def __init__(self):
        Planet.__init__(self, label = "ice")
        pp = PlanetaryParams(type = self.label)
        rp = RunParams(type = self.label)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params = pp, run_params = rp, default = True)   


class CustomPlanet(Planet):
    def __init__(self, planetary_parameters = {}, run_parameters = {}):
    
        default_values1 = {"ocean_fraction_should":0.25,
                                     "Mg_number_should":0.6,
                                     "T_surface_should":320,
                                     "Fe_number_mantle":0.15}


        default_values2 = {"adiabat_type":0}
        
        default_values1.update(planetary_parameters)
        Planet.__init__(self, label = "custom")
        pp = PlanetaryParams(type = self.label, default_values = default_values1)
        rp = RunParams(type = self.label, default_values = default_values2)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params = pp, run_params = rp, default = True) 

# pp = PlanetaryParams()
# rp = RunParams()

# pp.set_values(**pp.default_values)
# rp.set_values(**rp.default_values)

# pl = Planet()
# pl.set_values(planetary_params = pp, run_params = rp, default = True)


# pl.Si_number_mantle = 6

# #print (pl.default_values)
# pl.reset()

# pl.update(default = True)
# pl.show()

# pl1 = AquaPlanet()
# pl2 = TelluricPlanet()
# pl3 = InfernoPlanet()
# pl4 = IcePlanet()
# pl5 = CustomPlanet(planetary_parameters={"Fe_number_mantle":0.2})

# pl1.show()
# pl2.show()
# pl3.show()
# pl4.show()
# pl5.show()

