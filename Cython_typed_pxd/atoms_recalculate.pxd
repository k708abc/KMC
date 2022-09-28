# cython: language_level=3, boundscheck=False, wraparound=False

cdef list recalculate(
    int target,
    list atom_set,
    list bonds,
    list diffuse_candidates,
    bint height_change,
    int trans_val,
    list index_list,
    int unit_length
)