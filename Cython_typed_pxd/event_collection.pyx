# cython: language_level=3, boundscheck=False, wraparound=False

from cal_rates cimport rate
from Calc_grid_index cimport grid_num
from InputParameter cimport Params

cdef double total_energy(
    list atom_set, list bonds, int target, list energy_bonding, list energy_diffuse, list highest_atom, list index_list, int unit_length
):
    cdef double bond_energy = 0
    cdef int x ,y, z, bond,b_x,b_y, b_z
    x, y, z = index_list[target]
    if z == 0:
        bond_energy += energy_bonding[0]
    # highest_atom >= 1 if sites above the transformation layer is occupied
    if highest_atom[grid_num(x, y, 0, unit_length)] == 0:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                b_x, b_y, b_z = index_list[bond]
                bond_energy += energy_bonding[max(z, b_z)]
        return energy_diffuse[z] + bond_energy
    elif highest_atom[grid_num(x, y, 0, unit_length)] >= 1:  # transformtion is activated
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[len(energy_bonding)-1]  # bulk bonding energy
        return energy_diffuse[len(energy_diffuse)-1] + bond_energy  # bulk diffusion energy


cdef list find_filled_sites(list atom_set, list indexes):
    cdef int atom_nn
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 1]


cdef list find_empty_sites(list atom_set, list indexes):
    cdef int atom_nn
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 0]


cdef tuple check_cluster(int c_target, list atom_set, list bonds, list cluster_start, list index_list):
    cdef list cluster, next_atoms, check_atoms
    cdef int bond, atom, c_x, c_y, c_z, b_x, b_y, b_z
    cluster = cluster_start.copy()
    cluster += [c_target]
    next_atoms = [c_target]
    c_x, c_y, c_z = index_list[c_target]
    if c_z == 0:
        return [], False
    while next_atoms != []:
        check_atoms = next_atoms.copy()
        next_atoms = []
        for atom in check_atoms:
            for bond in bonds[atom]:
                if (atom_set[bond] == 1) and (bond not in cluster):
                    cluster += [bond]
                    next_atoms += [bond]
                    b_x, b_y, b_z = index_list[bond]
                    if (b_z == 0) or (len(cluster) >= 5):
                        return cluster, False
    return cluster, True


# Remove events which cause isolated atoms
cdef list judge_isolation(
    list atom_set,
    list bonds,
    int target,
    list nn_atoms,
    list events,
    list index_list
):
    cdef list clusters, checked, remove, cluster, nn_nn_list, common, nonused, empty
    cdef int nn_atom, event, nn_nn, i ,e_x, e_y, e_z
    cdef bint z_judge
    atom_set[target] = 0
    clusters = []
    checked = []
    remove = []
    empty = []
    for nn_atom in nn_atoms:
        if nn_atom not in checked:
            cluster, z_judge = check_cluster(nn_atom, atom_set, bonds, empty, index_list)
            checked += cluster
            if z_judge is True:
                clusters.append(cluster)
        else:
            pass
    #

    if clusters == []:
        for event in events:
            e_x,e_y, e_z = index_list[event]
            nn_nn_list = [nn_nn for nn_nn in bonds[event] if atom_set[nn_nn] == 1]
            if nn_nn_list == [] and e_z != 0:
                remove += [event]
    else:
        for event in events:
            atom_set[event] = 1
            for cluster in clusters:
                common = [bond for bond in bonds[event] if bond in cluster]
                if common == []:
                    remove += [event]
                    break
                else:
                    nonused, z_judge = check_cluster(event, atom_set, bonds, cluster, index_list)
                    if z_judge:
                        remove += [event]
            atom_set[event] = 0
    remove = list(set(remove))
    for i in remove:
        events.remove(i)
    atom_set[target] = 1
    return events


cdef tuple possible_events(
    list atom_set,
    list bonds,
    int target,
    Params params,
    double energy,
    list diffuse_candidates,
    list index_list
):
    cdef double pre, kbt,rearange_rate
    cdef list eve_rate, nn_atom,event_f, rates_f
    cdef int cand, nonused
    pre = float(params.prefactor)
    kbt = params.temperature_eV()
    eve_rate = []
    rearange_rate = rate(pre, kbt, energy)
    nn_atom = find_filled_sites(atom_set, bonds[target])
    #
    eve_rate += [cand for cand in diffuse_candidates[target] if atom_set[cand] == 0]
    event_f = judge_isolation(atom_set, bonds, target, nn_atom, eve_rate, index_list)
    rates_f = [rearange_rate for nonused in event_f]

    return event_f, rates_f


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
):
    cdef list event_list = []
    cdef list rate_list = []
    cdef double energy
    energy = total_energy(
        atom_set,
        bonds,
        target,
        energy_bonding,
        energy_diffuse,
        highest_atom,
        index_list,
        params.cell_size_xy
    )
    # calculate possible events
    event_list, rate_list = possible_events(
        atom_set, bonds, target, params, energy, diffuse_candidates, index_list
    )
    return event_list, rate_list
