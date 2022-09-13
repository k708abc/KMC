from setuptools import setup
from Cython.Build import cythonize
import Cython.Compiler.Options

Cython.Compiler.Options.annotate = True

setup(ext_modules=cythonize("atoms_recalculate.pyx"), annotate=True)
setup(ext_modules=cythonize("cal_rates.pyx"), annotate=True)
setup(ext_modules=cythonize("Calc_grid_index.pyx"), annotate=True)
setup(ext_modules=cythonize("deposition.pyx"), annotate=True)
setup(ext_modules=cythonize("event_collection.pyx"), annotate=True)
setup(ext_modules=cythonize("find_candidates.pyx"), annotate=True)
setup(ext_modules=cythonize("growth_mode_determination.pyx"), annotate=True)
setup(ext_modules=cythonize("kmc_functions.pyx"), annotate=True)
setup(ext_modules=cythonize("lattice_form.pyx"), annotate=True)
setup(ext_modules=cythonize("recording.pyx"), annotate=True)
setup(ext_modules=cythonize("rejection_free_choose.pyx"), annotate=True)
