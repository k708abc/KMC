# cython: language_level=3, boundscheck=False, wraparound=False

import collections
from Calc_grid_index cimport grid_num

cdef list find_aboves(list indexes, list index_list, int unit_length):
    cdef list above_list
    cdef int x, y, z, new_index, index
    above_list = []
    for index in indexes:
        x, y, z = index_list[index]
        new_index = grid_num(x, y, z + 1, unit_length)
        above_list.append(new_index)
    return above_list


cdef list find_lower_sites(list indexes, list index_list, int unit_length):
    cdef list below_list = []
    cdef int x, y, z, new_index, index
    for index in indexes:
        x, y, z = index_list[index]
        new_index = grid_num(x, y, z - 1, unit_length)
        below_list.append(new_index)
    return below_list


cdef list find_shares(list sites_list, list bonds):
    cdef list all_bonds = [], sfared_sites
    cdef int index, count, site
    for index in sites_list:
        all_bonds += bonds[index]
    counts = collections.Counter(all_bonds)
    shared_sites = [site for site, count in counts.items() if count >= 2]
    return shared_sites


cdef list find_candidates(list bonds, list index_list, int atom_x, int atom_y, int atom_z, int unit_length, int z_max):
    cdef list nnn_sites, nn_sites, above_layer, lower_layer, nn_vert, shared, shared_nnn_nighbor, candidates
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
    candidates = (
        nnn_sites
        + nn_sites
        + above_layer
        + lower_layer
        + nn_vert
        + shared
        + shared_nnn_nighbor
    )

    return list(set(candidates))
