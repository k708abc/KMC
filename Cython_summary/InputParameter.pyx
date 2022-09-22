# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False
# cython: cdivision=True

import yaml
from libc.math cimport round

cdef class Params:
    def __init__(self, str filename):
        if filename:
            with open(filename) as f:
                input_yaml = yaml.safe_load(f)
            for k, v in input_yaml.items():
                setattr(self, k, v)

    cdef double temperature_eV(self):
        return self.temperature * 8.617e-5

    cdef int atoms_in_BL(self):
        return int(2 * (self.cell_size_xy) ** 2)

    cdef double dep_rate_atoms_persec(self):
        return self.dep_rate / 60 * self.atoms_in_BL()

    cdef double total_time(self):
        return self.dep_time * 60

    cdef double interval(self):
        return self.total_time() * self.img_per / 100

    cdef int total_atoms(self):
        return int(round(self.dep_rate_atoms_persec() * self.total_time()))

    cdef int rec_num_atom_interval(self):
        return int(round(self.total_atoms() * self.img_per / 100))

    cdef int num_one_layer(self):
        return self.cell_size_xy**2

    cdef int num_grids(self):
        return self.num_one_layer() * self.z_max()

    cdef double unit_height(self):
        return 3 * (self.distance_intra + self.distance_inter)

    cdef int z_max(self):
        return self.cell_size_z * 6 -1


