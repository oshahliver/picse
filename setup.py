from setuptools import find_packages
from numpy.distutils.core import setup, Extension
import subprocess

PACKAGE_NAME = "picse"
LIB_DIR = "lib"
SRC_DIR = "lib"

# Generate object files for static linking:
subprocess.run(["make", "static"])

# Specify the object files for static linking
# Must be pre-compiled at this point
fort_objects = [
    "my_types",
    "run_params",
    "constants",
    "LinAlg",
    "class_table_new",
    "phase",
    "eosfort",
    "integrator",
    "utils",
    "eosmat",
    "fortshell",
    "fortlayer",
    "fortplanet",
]


ext1 = Extension(
    name="{}.utils.fortplanet".format(PACKAGE_NAME),
    sources=["{}/fortran/eosfort_wrapper.f95".format(SRC_DIR)],
    extra_objects=["{}/fortran/{}.o".format(LIB_DIR, fo) for fo in fort_objects],
    f2py_options=["--quiet"],
)

ext2 = Extension(
    name="{}.utils.fortutils".format(PACKAGE_NAME),
    sources=["{}/fortran/utils.f95".format(SRC_DIR)],
    extra_objects=["{}/fortran/{}.o".format(LIB_DIR, fo) for fo in ["constants", "LinAlg"]],
    f2py_options=["--quiet"],
)

setup(
    name=PACKAGE_NAME,
    version="0.0.1",
    install_requires=[
        "astropy",
        "plotly",
        "kaleido",
        "tabulate",
        "scipy",
        "scikit-learn",
        "numpy",
        "pandas",
        "seafreeze",
        "iapws",
        "keras",
        "tensorflow",
        "matplotlib",
        "alive_progress",
    ],
    packages=find_packages(),
    ext_modules=[ext1, ext2],
    # include_package_data=True,
    package_dir={"picse": "picse"},
    package_data={"": ["*.tab", "*.pkl", "*.csv"]},
    #   optional=os.environ.get('CIBUILDWHEEL', '0') != '1',
)

# remove static library files
subprocess.run(["make", "clean"])
