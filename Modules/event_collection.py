from typing import List, Dict, Tuple
from Modules.cal_rates import rate
from decimal import Decimal

def total_energy(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    energy_bonding: List[float],
    energy_diffuse: List[float],
    highest_atom: Dict[Tuple, int],
):
    bond_energy = 0
    if target[2] == 0:
        bond_energy += energy_bonding[0]
    # highest_atom >= 1 if sites above the transformation layer is occupied
    if highest_atom[(target[0], target[1])] == 0:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[max(target[2], bond[2])]
        return energy_diffuse[target[2]] + bond_energy
    elif highest_atom[(target[0], target[1])] >= 1: # transformtion is activated
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[-1] #bulk bonding energy
        return energy_diffuse[-1] + bond_energy #bulk diffusion energy


def find_filled_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 1]


def find_empty_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 0]


def check_cluster(c_target, atom_set, bonds, cluster_start):
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
def judge_isolation(
    atom_set,
    bonds,
    target: Tuple[int, int, int],
    nn_atoms: List[Tuple[int, int, int]],
    events: List[Tuple[int, int, int]],
):
    atom_set[target] = 0
    clusters = []
    checked = []
    remove = []
    for nn_atom in nn_atoms:
        if nn_atom not in checked:
            cluster, z_judge = check_cluster(nn_atom, atom_set, bonds, [])
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
                    _, z_judge = check_cluster(event, atom_set, bonds, cluster)
                    if z_judge:
                        remove += [event]
            atom_set[event] = 0
    remove = list(set(remove))
    for i in remove:
        events.remove(i)
    atom_set[target] = 1
    return events


def possible_events(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    params,
    energy: float,
    diffuse_candidates: Dict,
):
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    eve_rate: List = []
    rearange_rate = Decimal(rate(pre, kbt, energy)).quantize(Decimal("0.00000001"))
    nn_atom = find_filled_sites(atom_set, bonds[target])
    #
    eve_rate += [cand for cand in diffuse_candidates[target] if atom_set[cand] == 0]
    event_f = judge_isolation(atom_set, bonds, target, nn_atom, eve_rate)
    rates_f: List = [rearange_rate for _ in event_f]
    return event_f, rates_f


def site_events(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    params,
    energy_bonding: List[float],
    energy_diffuse: List[float],
    diffuse_candidates: Dict,
    highest_atom: Dict,
):
    event_list: List[Tuple] = []
    rate_list: List[float] = []
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
