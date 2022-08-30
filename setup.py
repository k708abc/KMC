from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("Rejection_free_kmc.pyx"))
