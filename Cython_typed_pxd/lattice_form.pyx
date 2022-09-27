# cython: language_level=3, boundscheck=False, wraparound=False


from Modules.InputParameter import Params
from Modules.find_candidates import find_candidates
from Modules.Calc_grid_index import grid_num


cpdef list form_first_3BL(double z_intra, double z_inter, int unit_length):
    cdef int x, y, index
    cdef list lattice_first = [[]]*unit_length**2
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


cpdef list lattice_full_layers(double unit_height, int unit_length, int z_max, list lattice_first, int num_grids):
    cdef int x, y, z, index, index_xy
    cdef list lattice = [[]]*num_grids
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



cpdef list neighbor_points(
    int x, int y, int z, int unit_length, int z_max
):
    cdef list neighbors
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
            neighbors.append(grid_num(x, y, z - 1, unit_length))
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


cpdef list search_bond(int unit_length, int z_max, int num_grids):
    # Search for bonding atoms for all the atoms
    cdef int x, y, z, index
    cdef list bonds = [[]]*num_grids
    for x in range(0, unit_length):
        for y in range(0, unit_length):
            for z in range(0, z_max):
                index = grid_num(x, y, z, unit_length)
                bonds[index] = neighbor_points(x, y, z, unit_length, z_max)
    return bonds


cpdef tuple lattice_form(input_params):
    cdef int unit_length, z_units, z_max, num_one_layer, num_grids, i, j, k, index
    cdef double z_inter, z_intra, unit_height
    cdef list lattice_first
    cdef list lattice
    cdef list atom_set
    cdef list bonds
    cdef list event
    cdef list event_time
    cdef list event_time_tot
    cdef list site_list_correspondance
    cdef list list_site_correspondance
    cdef list diffuse_candidates
    cdef list highest_atom
    cdef list index_list
    #
    unit_length = input_params.cell_size_xy
    z_units = input_params.cell_size_z
    z_intra = float(input_params.distance_intra)
    z_inter = float(input_params.distance_inter)
    unit_height = 3 * (z_intra + z_inter)
    z_max = z_units * 6 -1
    #
    num_one_layer = unit_length**2
    num_grids = num_one_layer * z_max
    #
    highest_atom = [0]*num_one_layer
    atom_set = [0]*num_grids
    event = [[]]*num_grids
    event_time = [[]]*num_grids
    event_time_tot = [0]*num_grids
    diffuse_candidates = [[]]*num_grids
    index_list = [(int, int, int)]*num_grids
    #
    lattice_first = form_first_3BL(z_intra, z_inter, unit_length)
    lattice = lattice_full_layers(unit_height, unit_length, z_max, lattice_first, num_grids)
    bonds = search_bond(unit_length, z_max, num_grids)
    #
    for x in range(unit_length):
        for y in range(unit_length):
            for z in range(z_max):
                index = grid_num(x, y, z, unit_length)
                index_list[index] = (x, y, z)
    for x in range(unit_length):
        for y in range(unit_length):
            for z in range(z_max):
                index = grid_num(x, y, z, unit_length)
                diffuse_candidates[index] = find_candidates(bonds, index_list, x, y, z, unit_length, z_max)
    return (
        lattice,
        bonds,
        atom_set,
        event,
        event_time,
        event_time_tot,
        diffuse_candidates,
        highest_atom,
        index_list,
    )
