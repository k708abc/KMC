# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from libc.stdlib cimport rand, RAND_MAX
from libcpp.vector cimport vector
from vector_operation cimport remove_duplicate


cdef vector[int] find_candidates(vector[int] atom_set, vector[vector[int]] bonds, vector[vector[int]] index_list):
    cdef vector[int] candidate
    cdef vector[int] atom_index
    cdef int condition, bond, index
    for index in range(atom_set.size()):
        condition = atom_set[index]
        atom_index = index_list[index]
        if (atom_index[2] == 0) and (condition == 0):
            candidate.push_back(index)
        elif condition == 1:
            for bond in bonds[index]:
                if atom_set[bond] == 0:
                    candidate.push_back(bond)
    return remove_duplicate(candidate)


cdef int dep_position(vector[int] candidate):
    cdef int N = candidate.size()
    return candidate[int(rand()/RAND_MAX * N)]


cdef int deposit_an_atom(vector[int] atom_set, vector[vector[int]] bonds, vector[vector[int]] index_list):
    cdef vector[int] caindidate 
    cdef int dep_pos
    candidate = find_candidates(atom_set, bonds, index_list)
    dep_pos = dep_position(candidate)
    return dep_pos
