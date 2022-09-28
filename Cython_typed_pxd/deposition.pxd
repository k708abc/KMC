# cython: language_level=3, boundscheck=False, wraparound=False

cdef list find_candidates(list atom_set, list bonds, list index_list)

cdef int dep_position(list candidate)

cdef int deposit_an_atom(list atom_set, list bonds, list index_list)