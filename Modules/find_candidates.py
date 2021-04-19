def find_aboves(indexes):
    return [(i[0], i[1], i[2] + 1) for i in indexes]


def find_lower_sites(indexes):
    return [(i[0], i[1], i[2] - 1) for i in indexes]


def find_candidates(bonds, target, unit_length, z_max):
    atom_x, atom_y, atom_z = target
    nnn_sites = [
        ((atom_x - 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y - 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y + 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, (atom_y - 1) % unit_length, atom_z),
        ((atom_x - 1) % unit_length, (atom_y + 1) % unit_length, atom_z),
    ]
    nn_sites = bonds[target]

    if atom_z == 0:
        above_layer = find_aboves(nn_sites)
        lower_layer = []
        nn_vert = []
    elif atom_z == 1:
        above_layer = find_aboves(nnn_sites)
        lower_layer = []
        nn_vert = bonds[atom_x, atom_y, atom_z + 1]

    elif atom_z == z_max:
        above_layer = []
        lower_layer = find_lower_sites(nn_sites)
        nn_vert = []
    elif atom_z == z_max - 1:
        above_layer = []
        lower_layer = find_lower_sites(nnn_sites)
        nn_vert = bonds[atom_x, atom_y, atom_z - 1]

    elif atom_z % 2 == 1:
        above_layer = find_aboves(nnn_sites)
        lower_layer = find_lower_sites(nn_sites)
        nn_vert = bonds[atom_x, atom_y, atom_z + 1]
    elif atom_z % 2 == 0:
        above_layer = find_aboves(nn_sites)
        lower_layer = find_lower_sites(nnn_sites)
        nn_vert = bonds[atom_x, atom_y, atom_z - 1]

    candidates = nnn_sites + nn_sites + above_layer + lower_layer + nn_vert

    return candidates
