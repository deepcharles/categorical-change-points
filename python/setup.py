from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
from pathlib import Path

def get_eigen_include()->str:
    """Helper function to determine the Eigen include path.

    The purpose of this function is to postpone importing includeigen
    until it is actually installed.
    """
    import includeigen
    return includeigen.get_include()

def get_numpy_include()->str:
    """Helper function to determine the numpy include path.

    The purpose of this function is to postpone importing includeigen
    until it is actually installed.
    """
    import numpy as np
    return np.get_include()

# Collect paths to cpp sources
cpp_sources = list(Path("src/categoricalfpop").resolve().glob("*.cpp"))
cpp_sources += list(Path("../cpp").resolve().glob("*.cpp"))
cpp_sources = [str(path) for path in cpp_sources]  # convert to strings

ext_modules = [
    Pybind11Extension(
        "categoricalfpop",
        sorted(cpp_sources),
        cxx_std=14,
    ),
]

if __name__ == "__main__":
    setup(
        ext_modules=ext_modules,
        include_dirs=[get_numpy_include(), get_eigen_include()],
    )
