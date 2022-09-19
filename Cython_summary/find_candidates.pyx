# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False
import collections
from Calc_grid_index cimport grid_num
from libcpp.vector cimport vector
from vector_operation cimport remove_duplicate

cdef vector[int] find_aboves(vector[int] indexes, vector[vector[int]] index_list, int unit_length):
    cdef vector[int] above_list = []
    cdef int x, y, z, new_index, index
    for index in indexes:
        x, y, z = index_list[index]
        new_index = grid_num(x, y, z + 1, unit_length)
        above_list.push_back(new_index)
    return above_list


cdef vector[int] find_lower_sites(vector[int] indexes, vector[vector[int]] index_list, int unit_length):
    cdef vector[int] below_list = []
    cdef int x, y, z, new_index, index
    for index in indexes:
        x, y, z = index_list[index]
        new_index = grid_num(x, y, z - 1, unit_length)
        below_list.push_back(new_index)
    return below_list


cdef vector[int] find_shares(vector[int] sites_list, vector[vector[int]] bonds):
    cdef vector[int] all_bonds = [], sfared_sites
    cdef int index, count, site
    for index in sites_list:
        all_bonds.insert(all_bonds.end(), bonds[index].begin(), bonds[index].end())
    counts = collections.Counter(all_bonds)
    shared_sites = [site for site, count in counts.items() if count >= 2]
    return shared_sites


cdef vector[int] find_candidates(vector[vector[int]] bonds, vector[vector[int]] index_list, int atom_x, int atom_y, int atom_z, int unit_length, int z_max):
    cdef vector[int] nnn_sites, nn_sites, above_layer, lower_layer, nn_vert, shared, shared_nnn_nighbor, candidates
    cdef int index
    index = grid_num(atom_x, atom_y, atom_z, unit_length)
    nnn_sites = [
        grid_num((atom_x - 1) % unit_length, atom_y, atom_z, unit_length),
        grid_num(atom_x, (atom_y - 1) % unit_length, atom_z, unit_length),
        grid_num((atom_x + 1) % unit_length, atom_y, atom_z, unit_length),
        grid_num(atom_x, (atom_y + 1) % unit_length, atom_z, unit_length),
        grid_num((atom_x + 1) % unit_length, (atom_y - 1) % unit_length, atom_z, unit_length),
        grid_num((atom_x - 1) % unit_length, (atom_y + 1) % unit_length, atom_z, unit_length),
    ]
    nn_sites = bonds[index]
    if atom_z == 0:
        above_layer = find_aboves(nn_sites, index_list, unit_length)
        shared = find_shares(above_layer, bonds)
        lower_layer = []
        nn_vert = []
    elif atom_z == 1:
        above_layer = find_aboves(nnn_sites, index_list, unit_length)
        lower_layer = []
        shared = []
        nn_vert = bonds[grid_num(atom_x, atom_y, atom_z + 1, unit_length)]
    elif atom_z == z_max:
        above_layer = []
        lower_layer = find_lower_sites(nn_sites, index_list, unit_length)
        nn_vert = []
        shared = find_shares(lower_layer, bonds)
    elif atom_z == z_max - 1:
        above_layer = []
        lower_layer = find_lower_sites(nnn_sites, index_list, unit_length)
        nn_vert = bonds[grid_num(atom_x, atom_y, atom_z - 1, unit_length)]
        shared = []
    elif atom_z % 2 == 1:
        above_layer = find_aboves(nnn_sites, index_list, unit_length)
        lower_layer = find_lower_sites(nn_sites, index_list, unit_length)
        nn_vert = bonds[grid_num(atom_x, atom_y, atom_z + 1, unit_length)]
        shared = find_shares(lower_layer, bonds)
    elif atom_z % 2 == 0:
        above_layer = find_aboves(nn_sites, index_list, unit_length)
        lower_layer = find_lower_sites(nnn_sites, index_list, unit_length)
        nn_vert = bonds[grid_num(atom_x, atom_y, atom_z - 1, unit_length)]
        shared = find_shares(above_layer, bonds)
    #
    shared_nnn_nighbor = find_shares(nnn_sites, bonds)
    #
    candidates.insert(candidates.end(), nnn_sites.begin(), nnn_sites.end())
    candidates.insert(candidates.end(), nn_sites.begin(), nn_sites.end())
    candidates.insert(candidates.end(), above_layer.begin(), above_layer.end())
    candidates.insert(candidates.end(), lower_layer.begin(), lower_layer.end())
    candidates.insert(candidates.end(), nn_vert.begin(), nn_vert.end())
    candidates.insert(candidates.end(), shared.begin(), shared.end())
    candidates.insert(candidates.end(), shared_nnn_nighbor.begin(), shared_nnn_nighbor.end())
    candidates = remove_duplicate(candidates)
    return candidates
