# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False
from libcpp.vector cimport vector

cdef class Params:
    cdef int cell_size_xy, cell_size_z, first_put_num, max_workers_val, repeat_val, trans_num
    cdef double dep_rate, dep_time, img_par, distance_inter, distance_intra, img_per, prefactor
    cdef double repeat_E1_diff, repeat_E1_end, repeat_E1_start, repeat_E2_diff, repeat_E2_end, repeat_E2_start, temperature
    cdef str comments, record_name
    cdef bint first_put_check, start_from_middle, trans_check
    cdef dict energies_binding
    cdef dict energies_diffusion
    cdef list repeat_combos

cdef class common_functions:
    cdef double start_time, total_event_time, elapsed_time, prog_time, z_intra, z_inter
    cdef int n_atoms, n_events, rec_num_atoms, n_events_perdep, event_number, setting_value, unit_length, total_time
    cdef int z_units, unit_height, z_max, num_one_layer, num_grids, index
    cdef int move_atom, target
    cdef bint height_change_rem, height_change_add
    cdef vector[vector[int]] pos_rec
    cdef vector[double] time_rec
    cdef vector[double] cov_rec
    cdef vector[int] n_events_rec
    cdef vector[int] num_atoms_rec
    cdef vector[vector[double]] lattice
    cdef vector[vector[int]] bonds
    cdef vector[int] atom_set
    cdef vector[vector[int]] event
    cdef vector[vector[double]] event_time
    cdef vector[double] event_time_tot
    cdef vector[vector[int]] index_list
    cdef vector[vector[int]] diffuse_candidates
    cdef vector[int] highest_atom
    cdef vector[double] energy_bonding
    cdef vector[double] energy_diffuse
    cdef vector[int] related_atoms
    cdef vector[int] related_atoms_move
    cdef Params init_value
    cpdef loop(self)
    cpdef start_setting(self)
    cdef update_progress(self)
    cpdef start_rejection_free(self)
    cdef energy_summarize(self)
    cdef bint height_check_add(self, int pos)
    cdef bint height_check_remove(self, int pos)
    cdef update_events(self)
    cdef recalculate_total_rate(self)
    cdef double compare_rate(self)
    cdef update_after_deposition(self, int dep_pos)
    cdef int deposition(self)
    cdef rejection_free_deposition(self)
    cpdef put_first_atoms_rf(self)
    cdef rejection_free_event(self)
    cdef event_progress(self)
    cpdef pickle_dump(self, obj, path)
    cpdef parameter_record(self)
    cdef rejection_free_loop(self)
    cpdef record_position(self)
    cpdef end_of_loop(self)




