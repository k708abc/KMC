# cython: language_level=3, boundscheck=False, wraparound=False

cdef list find_aboves(list indexes, list index_list, int unit_length)

cdef list find_lower_sites(list indexes, list index_list, int unit_length)

cdef list find_shares(list sites_list, list bonds)

cdef list find_candidates(list bonds, list index_list, int atom_x, int atom_y, int atom_z, int unit_length, int z_max)

