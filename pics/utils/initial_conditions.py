import pickle
import numpy as np
import sys
from pics.utils.context_manager import temp_path
from pics.utils.internal_data import get_predictor_model


def predict_initials(
    M_surface_should=1.0,
    Mg_number_should=0.5,
    ocean_fraction_should=0.0,
    Si_number_mantle=0.4,
    T_surface_should=300.0,
    Fe_number_mantle=0.0,
    fractions=[],
    model_name="models_with_ocean.pkl",
):

    xi_Fe_core = fractions[1][0]
    
    models = get_predictor_model(model_name)

    # Loop over parameters to be predicted
    for i in range(len(models)):

        # Central temperature
        if i == 0:
            params = [
                [np.log10(M_surface_should)],
                [Mg_number_should],
                # [xi_Fe_core],
                [T_surface_should],
                [ocean_fraction_should],
            ]

        # Central pressure
        elif i == 1:
            params = [
                [np.log10(M_surface_should)],
                [Mg_number_should],
                # [xi_Fe_core],
                [Fe_number_mantle],
                # [of]
            ]

        # Core mass
        elif i == 2:
            pc = models[1].result[0]
            params = [
                [M_surface_should],
                [Mg_number_should],
                [xi_Fe_core],
                [Fe_number_mantle],
                [Si_number_mantle],
                # [pc],
                # [of]
            ]

        models[i].evaluate(np.array(params))

    return [
        10 ** models[0].result[0],
        10 ** models[1].result[0],
        models[2].result[0],
    ]
