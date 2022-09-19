from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import Cython.Compiler.Options

Cython.Compiler.Options.annotate = True

ext_modules = [
    Extension("kmc_functions", ["kmc_functions.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))

"""
ext_modules = [
    Extension("kmc_functions", ["kmc_functions.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))

ext_modules = [
    Extension("lattice_form", ["lattice_form.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))


ext_modules = [
    Extension("rejection_free_choose", ["rejection_free_choose.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))

ext_modules = [
    Extension("InputParameter", ["InputParameter.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))


ext_modules = [
    Extension("find_candidates", ["find_candidates.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))


ext_modules = [
    Extension("vector_operation", ["vector_operation.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))


ext_modules = [
    Extension("event_collection", ["event_collection.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))

ext_modules = [
    Extension("deposition", ["deposition.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))



ext_modules = [
    Extension("cal_rates", ["cal_rates.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))



ext_modules = [
    Extension("atoms_recalculate", ["atoms_recalculate.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))

ext_modules = [
    Extension("Calc_grid_index", ["Calc_grid_index.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))


ext_modules = [
    Extension("search_element", ["search_element.pyx"]),
]
setup(ext_modules=cythonize(ext_modules, annotate=True))


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
"""
