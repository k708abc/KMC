# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

cdef class Params:
    cdef double temperature_eV(self)
    cdef int atoms_in_BL(self)
    cdef double dep_rate_atoms_persec(self)
    cdef double total_time(self)
    cdef double interval(self)
    cdef int total_atoms(self)
    cdef int rec_num_atom_interval(self)
    cdef int num_one_layer(self)
    cdef int num_grids(self)
    cdef double unit_height(self)
    cdef int z_max(self)
    cdef public int cell_size_xy, cell_size_z, first_put_num, max_workers_val, repeat_val, trans_num
    cdef public double dep_rate, dep_time, img_par, distance_inter, distance_intra, img_per, prefactor
    cdef public double repeat_E1_diff, repeat_E1_end, repeat_E1_start, repeat_E2_diff, repeat_E2_end, repeat_E2_start, temperature
    cdef public str comments, record_name
    cdef public bint first_put_check, start_from_middle, trans_check
    cdef public dict energies_binding
    cdef public dict energies_diffusion
    cdef public list repeat_combos