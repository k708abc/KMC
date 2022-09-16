# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

import random
import math
from libcpp.vector cimport vector

cdef vector[int] find_candidates(vector[int] atom_set, vector[vector[int]] bonds, vector[vector[int]] index_list):
    cdef vector[int] candidate = []
    cdef vector[int] atom_index
    cdef int condition, bond, index
    for index in range(len(atom_set)):
        condition = atom_set[index]
        atom_index = index_list[index]
        if condition == 1:
            pass
        elif atom_index[2] == 0:
            candidate.push_back(index)
        else:
            for bond in bonds[index]:
                if atom_set[bond] == 1:
                    candidate.push_back(index)
                    break
    return list(set(candidate))


cdef int dep_position(vector[int] candidate):
    return candidate[math.floor(random.random() * len(candidate))]


cdef int deposit_an_atom(vector[int] atom_set, vector[vector[int]] bonds, vector[vector[int]] index_list):
    cdef vector[int] caindidate 
    cdef int dep_pos
    candidate = find_candidates(atom_set, bonds, index_list)
    dep_pos = dep_position(candidate)
    return dep_pos
