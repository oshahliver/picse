from importlib import resources
import pickle
import io


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


def get_predictor_model(filename):
    file_path = "{}/{}".format(get_predictor_dir(), filename)

    with open(file_path, "rb") as fp:
        model = renamed_load(fp)

    return model


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
