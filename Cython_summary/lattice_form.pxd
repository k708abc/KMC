# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False
from libcpp.vector cimport vector


cdef vector[vector[vector[double]]] form_first_3BL(double z_intra, double z_inter, int unit_length)
cdef vector[vector[double]] lattice_full_layers(double unit_height, int unit_length, int z_max, vector[vector[vector[double]]] lattice_first, int num_grids)
cdef vector[int] neighbor_points(
    int x, int y, int z, int unit_length, int z_max
)
cdef vector[vector[int]] search_bond(int unit_length, int z_max, int num_grids)
#
cdef vector[vector[double]] lattice_form_lattice(
    int unit_length, 
    int z_units, 
    double z_intra, 
    double z_inter, 
    double unit_height, 
    int z_max,
    int num_one_layer,
    int num_grids
    )
cdef vector[vector[int]] lattice_form_bonds(int unit_length, int num_grids, int z_max)
cdef vector[int] lattice_form_atom_set(int num_grids)
cdef vector[vector[int]] lattice_form_event(int num_grids)
cdef vector[vector[double]] lattice_form_event_time(int num_grids)
cdef vector[double] lattice_form_event_time_tot(int num_grids)
cdef vector[vector[int]] lattice_form_index_list(int unit_length, int num_grids, int z_max)
cdef vector[vector[int]] lattice_form_diffuse_candidates(
    int unit_length,
    int num_grids,
    int z_max, 
    vector[vector[int]] bonds, 
    vector[vector[int]] index_list
    )
cdef vector[int] lattice_form_highest_atom(int num_one_layer)