# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False
from libcpp.vector cimport vector
from libcpp.pair cimport pair

cdef double total_energy(
    vector[int] atom_set, vector[vector[int]] bonds, int target, vector[double] energy_bonding, vector[double] energy_diffuse, vector[int] highest_atom, vector[vector[int]] index_list, int unit_length
)

cdef vector[int] find_filled_sites(vector[int] atom_set, vector[int] indexes)

cdef vector[int] find_empty_sites(vector[int] atom_set, vector[int] indexes)

cdef pair[vector[int], bint] check_cluster(int c_target, vector[int] atom_set, vector[vector[int]] bonds, vector[int] cluster_start, vector[vector[int]] index_list)

cdef vector[int] judge_isolation(
    vector[int] atom_set,
    vector[vector[int]] bonds,
    int target,
    vector[int] nn_atoms,
    vector[int] events,
    vector[vector[int]] index_list
)

cdef pair[vector[int], vector[double]] possible_events(
    vector[int] atom_set,
    vector[vector[int]] bonds,
    int target,
    double pre,
    double kbt,
    double energy,
    vector[vector[int]] diffuse_candidates,
    vector[vector[int]] index_list
)

cdef pair[vector[int], vector[int]] site_events(
    vector[int] atom_set,
    vector[vector[int]] bonds,
    int target,
    int unit_length,
    double pre,
    double kbt,
    vector[double] energy_bonding,
    vector[double] energy_diffuse,
    vector[vector[int]] diffuse_candidates,
    vector[int] highest_atom,
    vector[vector[int]] index_list
)