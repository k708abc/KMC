# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from find_candidates cimport find_candidates
from Calc_grid_index cimport grid_num
from libcpp.vector cimport vector

cdef vector[vector[vector[double]]] form_first_3BL(double z_intra, double z_inter, int unit_length):
    cdef int x, y, index
    cdef vector[vector[vector[double]]] lattice_first = [[]]*unit_length**2
    for x in range(0, unit_length):
        for y in range(0, unit_length):
            index = grid_num(x, y, 0, unit_length)
            lattice_first[index] = [
                [x, y, 0],
                [x + 1 / 3.0, y + 1 / 3.0, z_intra],
                [x + 1 / 3.0, y + 1 / 3.0, z_intra + z_inter],
                [x + 2 / 3.0, y + 2 / 3.0, 2 * z_intra + z_inter],
                [x + 2 / 3.0, y + 2 / 3.0, 2 * (z_inter + z_intra)],
                [x, y, 2 * (z_inter + z_intra) + z_intra],
            ]
    return lattice_first



cdef vector[vector[double]] lattice_full_layers(double unit_height, int unit_length, int z_max, vector[vector[vector[double]]] lattice_first, int num_grids):
    cdef int x, y, z, index, index_xy
    cdef vector[vector[double]] lattice = [[]]*num_grids
    for x in range(0, unit_length):
        for y in range(0, unit_length):
            for z in range(0, z_max):
                index = grid_num(x, y, z, unit_length)
                index_xy = grid_num(x, y, 0, unit_length)
                lattice[index] = [
                    lattice_first[index_xy][z%6][0],
                    lattice_first[index_xy][z%6][1],
                    lattice_first[index_xy][z%6][2] + unit_height * (z//6),
                ]
    return lattice


cdef vector[int] neighbor_points(
    int x, int y, int z, int unit_length, int z_max
):
    cdef vector[int] neighbors
    cdef int z_judge
    z_judge = z % 6
    if z == z_max - 1:
        neighbors = [
            grid_num((x - 1) % unit_length, y, z-1, unit_length), 
            grid_num(x, (y - 1) % unit_length, z-1, unit_length), 
            grid_num((x - 1) % unit_length, (y - 1) % unit_length, z - 1, unit_length)
            ]
    elif z_judge in (0, 2):
        neighbors = [
            grid_num((x - 1) % unit_length, y, z + 1, unit_length),
            grid_num(x, (y - 1) % unit_length, z + 1, unit_length),
            grid_num(x, y, z + 1, unit_length),
        ]
        if z != 0:
            neighbors.push_back(grid_num(x, y, z - 1, unit_length))
    elif z_judge in (1, 3):
        neighbors = [
            grid_num((x + 1) % unit_length, y, z - 1, unit_length),
            grid_num(x, (y + 1) % unit_length, z - 1, unit_length),
            grid_num(x, y, z - 1, unit_length),
            grid_num(x, y, z + 1, unit_length),
        ]
    elif z_judge == 4:
        neighbors = [
            grid_num((x + 1) % unit_length, y, z + 1, unit_length),
            grid_num((x + 1) % unit_length, (y + 1) % unit_length, z + 1, unit_length),
            grid_num(x, (y + 1) % unit_length, z + 1, unit_length),
            grid_num(x, y, z - 1, unit_length),
        ]
    elif z_judge == 5:
        neighbors = [
            grid_num((x - 1) % unit_length, y, z - 1, unit_length),
            grid_num(x, (y - 1) % unit_length, z - 1, unit_length),
            grid_num((x - 1) % unit_length, (y - 1) % unit_length, z - 1, unit_length),
            grid_num(x, y, z + 1, unit_length),
        ]
    else:
        raise RuntimeError("Something wrong happens. check z_judge value")
    return neighbors


cdef vector[vector[int]] search_bond(int unit_length, int z_max, int num_grids):
    # Search for bonding atoms for all the atoms
    cdef int x, y, z, index
    cdef vector[vector[int]] bonds = [[]]*num_grids
    for x in range(0, unit_length):
        for y in range(0, unit_length):
            for z in range(0, z_max):
                index = grid_num(x, y, z, unit_length)
                bonds[index] = neighbor_points(x, y, z, unit_length, z_max)
    return bonds

cdef vector[vector[double]] lattice_form_lattice(
    int unit_length, 
    int z_units, 
    double z_intra, 
    double z_inter, 
    double unit_height, 
    int z_max,
    int num_one_layer,
    int num_grids
    ):
    cdef vector[vector[vector[double]]] lattice_first
    cdef vector[vector[double]] lattice
    #
    lattice_first = form_first_3BL(z_intra, z_inter, unit_length)
    lattice = lattice_full_layers(unit_height, unit_length, z_max, lattice_first, num_grids)

    return lattice

cdef vector[vector[int]] lattice_form_bonds(int unit_length, int num_grids, int z_max):
    cdef vector[vector[int]] bonds
    bonds = search_bond(unit_length, z_max, num_grids)
    return bonds

cdef vector[int] lattice_form_atom_set(int num_grids):
    cdef vector[int] atom_set
    atom_set = [0]*num_grids
    return atom_set

cdef vector[vector[int]] lattice_form_event(int num_grids):
    cdef vector[vector[int]] event
    event = [[]]*num_grids
    return event


cdef vector[vector[double]] lattice_form_event_time(int num_grids):
    cdef vector[vector[double]] event_time
    event_time = [[]]*num_grids
    return event_time

cdef vector[double] lattice_form_event_time_tot(int num_grids):
    cdef vector[double] event_time_tot
    event_time_tot = [0]*num_grids
    return event_time_tot


cdef vector[vector[int]] lattice_form_index_list(int unit_length, int num_grids, int z_max):
    cdef int x, y, z, index
    cdef vector[vector[int]] index_list
    index_list = [[]]*num_grids
    for x in range(unit_length):
        for y in range(unit_length):
            for z in range(z_max):
                index = grid_num(x, y, z, unit_length)
                index_list[index] = [x, y, z]
    return index_list


cdef vector[vector[int]] lattice_form_diffuse_candidates(
    int unit_length,
    int num_grids,
    int z_max, 
    vector[vector[int]] bonds, 
    vector[vector[int]] index_list
    ):
    cdef vector[vector[int]] diffuse_candidates
    cdef int x, y, z, index
    diffuse_candidates = [[]]*num_grids
    for x in range(unit_length):
        for y in range(unit_length):
            for z in range(z_max):
                index = grid_num(x, y, z, unit_length)
                diffuse_candidates[index] = find_candidates(bonds, index_list, x, y, z, unit_length, z_max)
    return diffuse_candidates


cdef vector[int] lattice_form_highest_atom(int num_one_layer):
    cdef vector[int] highest_atom
    highest_atom = [0]*num_one_layer
    return highest_atom

