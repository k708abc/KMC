# cython: language_level=3, boundscheck=False, wraparound=False


cdef list form_first_3BL(double z_intra, double z_inter, int unit_length)
cdef list lattice_full_layers(double unit_height, int unit_length, int z_max, list lattice_first, int num_grids)
cdef list neighbor_points(
    int x, int y, int z, int unit_length, int z_max
)
cdef list search_bond(int unit_length, int z_max, int num_grids)
cdef tuple lattice_form(input_params)
