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
    # Ag-Si interaction
    if target[2] == 0:
        bond_energy += energy_bonding[0]
    if highest_atom[(target[0], target[1])] == 0:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[max(target[2], bond[2])]
        return energy_diffuse[target[2]] + bond_energy
    elif highest_atom[(target[0], target[1])] >= 1:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[-1]
        return energy_diffuse[-1] + bond_energy


def find_filled_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 1]


def find_empty_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 0]


# Remove events which cause isolated atoms
def judge_isolation(
    atom_set,
    bonds,
    target: Tuple[int, int, int],
    nn_atom: List[Tuple[int, int, int]],
    events: List[Tuple[int, int, int]],
):
    remove = []
    # for sites after move
    for check in events:
        bonds_sites = find_filled_sites(atom_set, bonds[check])
        num_bonds = len(bonds_sites)
        if (num_bonds >= 2) or (check[2] == 0):
            pass
        elif num_bonds == 0:
            remove.append(check)
        else:
            if bonds_sites[0] == target:
                remove.append(check)
    for i in remove:
        events.remove(i)
    # for neighbor atoms before move
    for check in nn_atom:
        nn_nn_atom = find_filled_sites(atom_set, bonds[check])
        if len(nn_nn_atom) >= 2:
            pass
        elif check[2] == 0:
            pass
        else:
            remove = []
            for post_move in events:
                # keep bonding after move
                if post_move in bonds[check]:
                    # nearest neighbor after move
                    post_nn_atom = find_filled_sites(atom_set, bonds[post_move])
                    # dimr isolation →remove
                    if len(post_nn_atom) <= 1:
                        remove.append(post_move)
                # No bond after move
                else:
                    remove.append(post_move)
            for i in remove:
                events.remove(i)
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


def rate_limit(rates, upper_limit) -> List:
    new_rate: List = []
    for rate_val in rates:
        if rate_val > upper_limit:
            new_rate.append(upper_limit)
        else:
            new_rate.append(rate_val)

    return new_rate


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
    if params.limit_check:
        rate_list = rate_limit(rate_list, Decimal(float(params.limit_num)))
    return event_list, rate_list
