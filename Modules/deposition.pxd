# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from libcpp.vector cimport vector

cdef vector[int] find_candidates(vector[int] atom_set, vector[vector[int]] bonds, vector[vector[int]] index_list)
cdef int dep_position(vector[int] candidate)
cdef int deposit_an_atom(vector[int] atom_set, vector[vector[int]] bonds, vector[vector[int]] index_list)