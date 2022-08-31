from Modules.Calc_grid_index import grid_num

cpdef list recalculate(
    int target,
    list atom_set,
    list bonds,
    list diffuse_candidates,
    bint height_change,
    int trans_val,
    list index_list,
    int unit_length
):
    cdef list candidate, 
    cdef int i, fill, bond, t_x, t_y, t_z
    candidate = [fill for fill in diffuse_candidates[target] if atom_set[fill] == 1] + [
        target
    ]
    for bond in bonds[target]:
        if atom_set[bond] == 0:
            candidate += [
                fill for fill in diffuse_candidates[bond] if atom_set[fill] == 1
            ]
    if height_change:
        tx, t_y, t_z = index_list[target]
        candidate += [
            grid_index(t_x, t_y, i, unit_length)
            for i in range(trans_val, -1, -1)
            if atom_set[grid_index(t_x, t_y, i, unit_length)] == 1
        ]
    return candidate
