import random
import math

cdef list find_candidates(list atom_set, list bonds, list index_list):
    cdef list candidate = []
    cdef tuple atom_index
    cdef int condition, bond
    for index in range(len(atom_set)):
        condition = atom_set[index]
        atom_index = index_list[index]
        if condition == 1:
            pass
        elif atom_index[2] == 0:
            candidate.append(atom_index)
        else:
            for bond in bonds[index]:
                if atom_set[bond] == 1:
                    candidate.append(index)
                    break
    return list(set(candidate))


cdef int dep_position(list candidate):
    return candidate[math.floor(random.random() * len(candidate))]


cpdef int deposit_an_atom(list atom_set, list bonds, list index_list):
    cdef list caindidate 
    candidate = find_candidates(atom_set, bonds, index_list)
    return dep_position(candidate)
