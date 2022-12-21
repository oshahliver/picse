from setuptools import find_packages
from numpy.distutils.core import setup, Extension
import site
import os

PACKAGE_NAME = "pics"
LIB_DIR = "lib"
SRC_DIR = "lib"

# Generate object files for static linking:
# TODO: Execute command to compile .f90 files

# Specify the object files for static linking
# Must be pre-compiled at this point
extra_objects = [
    "{}/fortran/my_types.o".format(LIB_DIR),
    "{}/fortran/run_params.o".format(LIB_DIR),
    "{}/fortran/constants.o".format(LIB_DIR),
    "{}/fortran/LinAlg.o".format(LIB_DIR),
    "{}/fortran/class_table_new.o".format(LIB_DIR),
    "{}/fortran/phase.o".format(LIB_DIR),
    "{}/fortran/eosfort.o".format(LIB_DIR),
    "{}/fortran/functions.o".format(LIB_DIR),
    "{}/fortran/eosmat.o".format(LIB_DIR),
    "{}/fortran/fortshell.o".format(LIB_DIR),
    "{}/fortran/fortlayer.o".format(LIB_DIR),
    "{}/fortran/fortplanet.o".format(LIB_DIR),
]


ext = Extension(
    name="{}.utils.wrapper".format(PACKAGE_NAME),
    sources=["{}/fortran/eosfort_wrapper.f95".format(SRC_DIR)],
    extra_objects=extra_objects,
    # libraries=['three'], # --> does not work!
    f2py_options=["--quiet"],
)


setup(
    name=PACKAGE_NAME,
    version="0.0.1",
    # package_dir={"": "foo"}, # --> does not include .py files!
    install_requires=[
        "astropy",
        "plotly",
        "tabulate",
        "scipy",
        "scikit-learn",
        "numpy",
        "pandas",
        "seafreeze",
        "iapws",
        "matplotlib",
    ],
    packages=find_packages(),
    ext_modules=[ext],
    #   # py_modules = ['main', 'bar.myclass'], # --> does not work!
    #   optional=os.environ.get('CIBUILDWHEEL', '0') != '1',
)
