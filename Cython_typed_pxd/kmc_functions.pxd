# cython: language_level=3, boundscheck=False, wraparound=False

from InputParameter cimport Params


cdef class common_functions:
    cdef Params init_value
    cdef double start_time, total_event_time, elapsed_time, prog_time
    cdef int n_atoms, n_events, rec_num_atoms, n_events_perdep, event_number, setting_value, unit_length, total_time
    cdef bint height_change_rem, height_change_add
    cdef int move_atom, target
    cdef list time_rec
    cdef list cov_rec
    cdef list n_events_rec
    cdef list num_atoms_rec
    cdef list event_time_tot
    cdef list energy_bonding
    cdef list energy_diffuse
    cdef list related_atoms
    cdef list lattice
    cdef list bonds
    cdef list atom_set
    cdef list event
    cdef list event_time
    cdef list diffuse_candidates
    cdef list highest_atom
    cdef list index_list
    cdef list pos_rec
    #
    """
    cdef list[double] time_rec
    cdef list[double] cov_rec
    cdef list[int] n_events_rec
    cdef list[int] num_atoms_rec
    cdef list[double] event_time_tot
    cdef list[double] energy_bonding
    cdef list[double] energy_diffuse
    cdef list[int] related_atoms
    cdef list[list[double]] lattice
    cdef list[list[int]] bonds
    cdef list[int] atom_set
    cdef list[list[int]] event
    cdef list[list[double]] event_time
    cdef list[list[int]] diffuse_candidates
    cdef list[int] highest_atom
    cdef list[list[int]] index_list
    cdef list[list[int]] pos_rec
    """
    #
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
    cdef pickle_dump(self, obj, path)
    cdef parameter_record(self)
    cdef rejection_free_loop(self)
    cdef record_position(self)
    cpdef end_of_loop(self)