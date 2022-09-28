# cython: language_level=3, boundscheck=False, wraparound=False
from InputParameter cimport Params

cdef double total_energy(
    list atom_set, list bonds, int target, list energy_bonding, list energy_diffuse, list highest_atom, list index_list, int unit_length
)

cdef list find_filled_sites(list atom_set, list indexes)

cdef list find_empty_sites(list atom_set, list indexes)

cdef tuple check_cluster(int c_target, list atom_set, list bonds, list cluster_start, list index_list)

cdef list judge_isolation(
    list atom_set,
    list bonds,
    int target,
    list nn_atoms,
    list events,
    list index_list
)

cdef tuple possible_events(
    list atom_set,
    list bonds,
    int target,
    Params params,
    double energy,
    list diffuse_candidates,
    list index_list
)

cdef tuple site_events(
    list atom_set,
    list bonds,
    int target,
    Params params,
    list energy_bonding,
    list energy_diffuse,
    list diffuse_candidates,
    list highest_atom,
    list index_list
)