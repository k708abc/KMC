# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False


from libcpp.vector cimport vector

cdef vector[int] find_aboves(vector[int] indexes, vector[vector[int]] index_list, int unit_length)
cdef vector[int] find_lower_sites(vector[int] indexes, vector[vector[int]] index_list, int unit_length)
cdef vector[int] find_shares(vector[int] sites_list, vector[vector[int]] bonds)
cdef vector[int] find_candidates(vector[vector[int]] bonds, vector[vector[int]] index_list, int atom_x, int atom_y, int atom_z, int unit_length, int z_max)
