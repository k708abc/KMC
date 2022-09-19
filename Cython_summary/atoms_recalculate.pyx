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
    cdef vector[int] candidate
    cdef int i, fill, bond, t_x, t_y, t_z
    candidate = [fill for fill in diffuse_candidates[target] if atom_set[fill] == 1] + [
        target
    ]
    for bond in bonds[target]:
        if atom_set[bond] == 0:
            candidate = candidate + [
                fill for fill in diffuse_candidates[bond] if atom_set[fill] == 1
            ]
    if height_change:
        t_x, t_y, t_z = index_list[target]
        candidate = candidate + [
            grid_num(t_x, t_y, i, unit_length)
            for i in range(trans_val, -1, -1)
            if atom_set[grid_num(t_x, t_y, i, unit_length)] == 1
        ]
    return candidate
