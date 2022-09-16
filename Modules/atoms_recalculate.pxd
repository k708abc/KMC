# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from libcpp.vector cimport vector

cdef vector[int] recalculate(
    int target,
    vector[int] atom_set,
    vector[vector[int]] bonds,
    vector[vector[int]] diffuse_candidates,
    bint height_change,
    int trans_val,
    vector[vector[int]] index_list,
    int unit_length
    )