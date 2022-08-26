from Modules.cal_rates import rate


cpdef double total_energy(
    dict atom_set, dict bonds, tuple target, list energy_bonding, list energy_diffuse, dict highest_atom
):
    cdef double bond_energy
    cdef tuple bond
    bond_energy = 0
    if target[2] == 0:
        bond_energy += energy_bonding[0]
    # highest_atom >= 1 if sites above the transformation layer is occupied
    if highest_atom[(target[0], target[1])] == 0:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[max(target[2], bond[2])]

        return energy_diffuse[target[2]] + bond_energy
    elif highest_atom[(target[0], target[1])] >= 1:  # transformtion is activated
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[-1]  # bulk bonding energy
        return energy_diffuse[-1] + bond_energy  # bulk diffusion energy


cdef list find_filled_sites(dict atom_set, list indexes):
    cdef tuple atom_nn
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 1]


cdef list find_empty_sites(dict atom_set, list indexes):
    cdef tuple atom_nn
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 0]


cdef tuple check_cluster(tuple c_target, dict atom_set, dict bonds, list cluster_start):
    cdef list cluster, next_atoms, check_atoms
    cdef tuple bond, atom
    cluster = cluster_start.copy()
    cluster += [c_target]
    next_atoms = [c_target]
    if c_target[2] == 0:
        return [], False
    while next_atoms != []:
        check_atoms = next_atoms.copy()
        next_atoms = []
        for atom in check_atoms:
            for bond in bonds[atom]:
                if (atom_set[bond] == 1) and (bond not in cluster):
                    cluster += [bond]
                    next_atoms += [bond]
                    if (bond[2] == 0) or (len(cluster) >= 5):
                        return cluster, False
    return cluster, True


# Remove events which cause isolated atoms
cpdef list judge_isolation(
    dict atom_set,
    dict bonds,
    tuple target,
    list nn_atoms,
    list events,
):
    cdef list clusters, checked, remove, cluster, nn_nn_list, common, nonused, empty
    cdef tuple nn_atom, event, nn_nn, i 
    cdef bint z_judge
    atom_set[target] = 0
    clusters = []
    checked = []
    remove = []
    empty = []
    for nn_atom in nn_atoms:
        if nn_atom not in checked:
            cluster, z_judge = check_cluster(nn_atom, atom_set, bonds, empty)
            checked += cluster
            if z_judge is True:
                clusters.append(cluster)
        else:
            pass
    #

    if clusters == []:
        for event in events:
            nn_nn_list = [nn_nn for nn_nn in bonds[event] if atom_set[nn_nn] == 1]
            if nn_nn_list == [] and event[2] != 0:
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
                    nonused, z_judge = check_cluster(event, atom_set, bonds, cluster)
                    if z_judge:
                        remove += [event]
            atom_set[event] = 0
    remove = list(set(remove))
    for i in remove:
        events.remove(i)
    atom_set[target] = 1
    return events


cpdef tuple possible_events(
    dict atom_set,
    dict bonds,
    tuple target,
    params,
    double energy,
    dict diffuse_candidates,
):
    cdef double pre, kbt,rearange_rate
    cdef list eve_rate, nn_atom,event_f, rates_f
    cdef tuple cand, nonused
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    eve_rate = []
    # rearange_rate = Decimal(rate(pre, kbt, energy)).quantize(Decimal("0.00000001"))
    rearange_rate = rate(pre, kbt, energy)
    nn_atom = find_filled_sites(atom_set, bonds[target])
    #
    eve_rate += [cand for cand in diffuse_candidates[target] if atom_set[cand] == 0]
    event_f = judge_isolation(atom_set, bonds, target, nn_atom, eve_rate)
    rates_f = [rearange_rate for nonused in event_f]

    return event_f, rates_f


cpdef tuple site_events(
    dict atom_set,
    dict bonds,
    tuple target,
    params,
    list energy_bonding,
    list energy_diffuse,
    dict diffuse_candidates,
    dict highest_atom,
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
    )
    # calculate possible events
    event_list, rate_list = possible_events(
        atom_set, bonds, target, params, energy, diffuse_candidates
    )
    return event_list, rate_list
