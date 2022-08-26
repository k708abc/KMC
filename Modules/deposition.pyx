import random
import math

cdef list find_candidates(dict atom_set, bonds):
    cdef list candidate = []
    cdef tuple atom_index, bond
    cdef int condition
    for atom_index, condition in atom_set.items():
        if condition == 1:
            pass
        elif atom_index[2] == 0:
            candidate.append(atom_index)
        else:
            for bond in bonds[atom_index]:
                if atom_set[bond] == 1:
                    candidate.append(atom_index)
                    break
    return list(set(candidate))


cdef tuple dep_position(list candidate):
    return candidate[math.floor(random.random() * len(candidate))]


cpdef tuple deposit_an_atom(dict atom_set, bonds):
    cdef list caindidate 
    candidate = find_candidates(atom_set, bonds)
    return dep_position(candidate)
