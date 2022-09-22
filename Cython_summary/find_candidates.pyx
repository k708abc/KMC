# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False


from Calc_grid_index cimport grid_num
from libcpp.vector cimport vector
from vector_operation cimport remove_duplicate

cdef vector[int] find_aboves(vector[int] indexes, vector[vector[int]] index_list, int unit_length):
    cdef vector[int] above_list, t_list
    cdef int x, y, z, new_index, index
    for index in indexes:
        t_list = index_list[index]
        x = t_list[0]
        y = t_list[1]
        z = t_list[2]
        new_index = grid_num(x, y, z + 1, unit_length)
        above_list.push_back(new_index)
    return above_list


cdef vector[int] find_lower_sites(vector[int] indexes, vector[vector[int]] index_list, int unit_length):
    cdef vector[int] below_list, t_list
    cdef int x, y, z, new_index, index
    for index in indexes:
        t_list = index_list[index]
        x = t_list[0]
        y = t_list[1]
        z = t_list[2]
        new_index = grid_num(x, y, z - 1, unit_length)
        below_list.push_back(new_index)
    return below_list

cdef vector[int] get_duplicate(vector[int] all_bonds):
    cdef vector[int] dup
    cdef int i, count, k
    for i in all_bonds:
        count = 0
        for k in all_bonds:
            if i == k:
                count += 1
        if count >= 2:
            dup.push_back(i)
    return remove_duplicate(dup)


cdef vector[int] find_shares(vector[int] sites_list, vector[vector[int]] bonds):
    cdef vector[int] all_bonds, shared_sites
    cdef int index, count, site
    for index in sites_list:
        all_bonds.insert(all_bonds.end(), bonds[index].begin(), bonds[index].end())
    shared_sites = get_duplicate(all_bonds)
    return shared_sites


cdef vector[int] find_candidates(vector[vector[int]] bonds, vector[vector[int]] index_list, int atom_x, int atom_y, int atom_z, int unit_length, int z_max):
    cdef vector[int] nnn_sites, nn_sites, above_layer, lower_layer, nn_vert, shared, shared_nnn_nighbor, candidates
    cdef int index
    index = grid_num(atom_x, atom_y, atom_z, unit_length)
    #
    nnn_sites.push_back(grid_num((atom_x - 1) % unit_length, atom_y, atom_z, unit_length))
    nnn_sites.push_back(grid_num(atom_x, (atom_y - 1) % unit_length, atom_z, unit_length))
    nnn_sites.push_back(grid_num((atom_x + 1) % unit_length, atom_y, atom_z, unit_length))
    nnn_sites.push_back(grid_num(atom_x, (atom_y + 1) % unit_length, atom_z, unit_length))
    nnn_sites.push_back(grid_num((atom_x + 1) % unit_length, (atom_y - 1) % unit_length, atom_z, unit_length))
    nnn_sites.push_back(grid_num((atom_x - 1) % unit_length, (atom_y + 1) % unit_length, atom_z, unit_length))
    #
    nn_sites = bonds[index]

    if atom_z == 0:
        above_layer = find_aboves(nn_sites, index_list, unit_length)
        shared = find_shares(above_layer, bonds)
    elif atom_z == 1:
        above_layer = find_aboves(nnn_sites, index_list, unit_length)
        nn_vert = bonds[grid_num(atom_x, atom_y, atom_z + 1, unit_length)]
    elif atom_z == z_max:
        lower_layer = find_lower_sites(nn_sites, index_list, unit_length)
        shared = find_shares(lower_layer, bonds)
    elif atom_z == z_max - 1:
        lower_layer = find_lower_sites(nnn_sites, index_list, unit_length)
        nn_vert = bonds[grid_num(atom_x, atom_y, atom_z - 1, unit_length)]
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
