from setuptools import find_packages
from numpy.distutils.core import setup, Extension
import subprocess
import shutil

PACKAGE_NAME = "picse"
LIB_DIR = "lib"
SRC_DIR = "lib"

# Remove previous build directory if it exists
shutil.rmtree('build', ignore_errors=True)

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
    "functions",
    "functionspy",
    "eosmat",
    "fortshell",
    "fortlayer",
    "fortplanet",
]

# Create extension for the fortplanet module
ext1 = Extension(
    name="{}.utils.fortplanet".format(PACKAGE_NAME),
    sources=["{}/fortran/eosfort_wrapper.f95".format(SRC_DIR)],
    extra_objects=["{}/fortran/{}.o".format(LIB_DIR, fo) for fo in fort_objects],
    # libraries=['three'], # --> does not work!
    f2py_options=["--quiet"],
)

# Create extension for the fortfunctions module
ext2 = Extension(
    name="{}.utils.fortfunctions".format(PACKAGE_NAME),
    sources=["{}/fortran/functionspy.f95".format(SRC_DIR)],
    extra_objects=["{}/fortran/{}.o".format(LIB_DIR, fo) for fo in ["LinAlg", "constants"]],
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
    #   # py_modules = ['main', 'bar.myclass'], # --> does not work!
    #   optional=os.environ.get('CIBUILDWHEEL', '0') != '1',
)


subprocess.run(["make", "clean"])