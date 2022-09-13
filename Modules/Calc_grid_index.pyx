# cython: language_level=3, boundscheck=False, wraparound=False


cpdef int grid_num(int x, int y, int z, int unit_length):
    return unit_length**2 * z + unit_length * y + x
