cpdef list recalculate(
    tuple target,
    dict atom_set,
    dict bonds,
    dict diffuse_candidates,
    bint height_change,
    int trans_val,
):
    cdef list candidate, 
    cdef tuple fill, bond
    cdef int i
    candidate = [fill for fill in diffuse_candidates[target] if atom_set[fill] == 1] + [
        target
    ]
    for bond in bonds[target]:
        if atom_set[bond] == 0:
            candidate += [
                fill for fill in diffuse_candidates[bond] if atom_set[fill] == 1
            ]

    if height_change:
        candidate += [
            (target[0], target[1], i)
            for i in range(trans_val, -1, -1)
            if atom_set[(target[0], target[1], i)] == 1
        ]
    return candidate
