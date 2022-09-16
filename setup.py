from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import Cython.Compiler.Options

Cython.Compiler.Options.annotate = True

ext_modules = [
    Extension("Modules.kmc_functions", ["Modules/kmc_functions.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))
