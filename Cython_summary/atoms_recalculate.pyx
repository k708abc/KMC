# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from Calc_grid_index cimport grid_num
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
):
    cdef vector[int] candidate, t_list
    cdef int i, fill, bond, t_x, t_y, t_z
    for fill in diffuse_candidates[target]:
        if atom_set[fill] == 1:
            candidate.push_back(fill)
    candidate.push_back(target)
    for bond in bonds[target]:
        if atom_set[bond] == 0:
            for fill in diffuse_candidates[bond]:
                if atom_set[fill] == 1:
                    candidate.push_back(fill)
    if height_change:
        t_list = index_list[target]
        t_x = t_list[0]
        t_y = t_list[1]
        t_z = t_list[2]
        for i in range(trans_val, -1, -1):
            if atom_set[grid_num(t_x, t_y, i, unit_length)] == 1:
                candidate.push_back(grid_num(t_x, t_y, i, unit_length))
    return candidate
