from importlib import resources
from pics.utils import file_manager as fman
import pkg_resources
import pickle
import pandas as pd
import csv
import io
import numpy as np

import os
import pandas as pd


def create_data(planets):
    labels = {
        "mass": "M_surface_is",
        "radius": "R_surface_is",
        "mg_number": "Mg_number_is",
        "si_number": "Si_number_is",
        "si_number_mantle": "Si_number_mantle",
        "fe_number_mantle": "Fe_number_mantle",
        "temp_surface": "T_surface_is",
        "pres_surface": "P_surface_is",
        "temp_center": "T_center",
        "pres_center": "P_center",
        "core_mass_fraction": "core_mass_fraction_is",
        "ocean_mass_fraction": "ocean_fraction_is",
    }

    data = {}
    for key, val in labels.items():

        data.update({key: [getattr(pl, val, None) for pl in planets]})

    return pd.DataFrame(data)


def check_type(s):
    try:
        float(s)
        if "." in s:
            return float(s)
        else:
            return int(s)

    except ValueError:
        return s


def read_sample(path, specs={}, comment="#"):
    # Read data
    with open(path) as fd:
        df = pd.read_csv(fd, delimiter=",", comment=comment)

    # Read in the comments and store them as a dictionary
    with open(path, "r") as f:
        comments = {}
        for line in f:
            if line.startswith(comment):
                key, value = line.strip(comment).strip().split(":")
                comments[key.strip()] = value.strip()

    metadata = {}
    current_key = None

    with open(path, "r") as f:
        for line in f:
            if line.startswith(comment):
                if line[1] != "!":
                    current_key = line.strip()[1:].strip()
                    metadata[current_key] = {}
                else:
                    if current_key is None:
                        continue
                    key_value = line.strip()[2:].strip().split(":")
                    if "--" in key_value[1]:
                        val = [check_type(x.strip()) for x in key_value[1].split("--")]
                    else:
                        val = check_type(key_value[1].strip())

                    metadata[current_key][key_value[0].strip()] = val

    return df, metadata


def write_sample(data, path, meta={}, specs={}, comment="#"):
    if len(list(meta.keys())) != 0:
        metadata = meta
        add_meta = True

    else:
        add_meta = False

    allowed_specs = {
        "conflict": ["ommit", "add", "overwrite"],
        "format": ["csv", "bin"],
        "create_log": ["False", "True", "true", "false", True, False],
    }
    allowed_keys = list(allowed_specs.keys())
    allowed_vals = list(allowed_specs.values())

    # Check validity of passed specifications via the specs key argument
    invalid_keys = set(specs.keys()) - set(allowed_keys)
    if invalid_keys:
        raise ValueError(f"Invalid keys found: {invalid_keys}")

    # Check if all values provided for allowed keys are valid
    invalid_vals = {
        k: v for k, v in specs.items() if k in allowed_keys and v not in allowed_vals
    }
    if invalid_vals:
        raise ValueError(f"Invalid values found: {invalid_vals}")

    # Set default values for all spec keys that are not passed
    # By convention the first element of the possibility list is the default
    for key, val in allowed_vals.items():
        if not key in specs.keys():
            specs.update({key: val[0]})

    if add_meta:
        # Open a file for writing and write metadata as comments
        with open(path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile, delimiter=",")
            for key, value in metadata.items():
                if isinstance(value, dict):
                    string = f"{comment} {key}:"
                    writer.writerow([string])
                    for k, v in value.items():
                        # check type of value
                        if isinstance(v, list):
                            var_string = "--".join(str(x) for x in v)
                        else:
                            var_string = str(v)
                        string = f"{comment}!{k}: {var_string}"
                        writer.writerow([string])
                else:

                    if isinstance(value, list):
                        var_string = "--".join(str(x) for x in value)
                    else:
                        var_string = str(value)
                    string = f"{comment} {key}: {var_string}"
                    writer.writerow([string])

            # Write DataFrame to the file without adding an empty row
            data.to_csv(csvfile, index=False)

    else:
        data.to_csv(path)

    # Create a log file for the data set from which the simulation can
    # directly be reproduced
    if specs["create_log"] in ["True", "true", True]:
        pass


def get_eos_dir():
    """Get the absoulte path of the package directory containing the
    EoS tables.
    """
    with resources.path("pics.data", "EoS_tables") as tabdir:
        return tabdir


def get_predictor_dir():
    """Get the pkl models to predict initial conditions for given set
    of boundary conditions
    """
    with resources.path("pics.data", "initial_conditions") as moddir:
        return moddir


def get_training_dir():
    """Get the directory containing the training data"""
    with resources.path("pics.data", "training_data") as moddir:
        return moddir


def get_predictor_model(filename):
    file_path = "{}/{}".format(get_predictor_dir(), filename)

    with open(file_path, "rb") as fp:
        model = renamed_load(fp)

    return model


def get_training_data(filename):
    file_path = "{}/{}".format(get_training_dir(), filename)

    with open(file_path, "rb") as fp:
        df = pd.read_csv(file_path)

    return df


class RenameUnpickler(pickle.Unpickler):
    """Implements solution to change pickle module suggested by bossylobster:
    https://stackoverflow.com/questions/2121874/python-pickling-after-changing-a-modules-directory/2121918#2121918

    This can be used to unpickle a pickled object if the location of the module it was
    created with changed and allows reuse of the same pickle files in distributed packages.
    """

    def find_class(self, module, name):
        renamed_module = module
        if module == "regression":
            renamed_module = "pics.utils.regression"

        return super(RenameUnpickler, self).find_class(renamed_module, name)


def renamed_load(file_obj):
    return RenameUnpickler(file_obj).load()
