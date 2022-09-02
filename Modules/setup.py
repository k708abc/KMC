from setuptools import setup
from Cython.Build import cythonize

# setup(ext_modules=cythonize("kmc_functions.pyx"))
setup(ext_modules=cythonize("atoms_recalculate.pyx"))
setup(ext_modules=cythonize("lattice_form.pyx"))
