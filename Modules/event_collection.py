from typing import List, Dict, Tuple
from Modules.cal_rates import rate
from InputParameter import Params
from decimal import Decimal


def total_energy(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    energy_bonding: List[float],
    energy_diffuse: List[float],
    trans_val: int,
    highest_atom: Dict[Tuple, int],
):
    bond_energy = 0
    # Ag-Si interaction
    if target[2] == 0:
        bond_energy += energy_bonding[0]
    if trans_val > highest_atom[(target[0], target[1])]:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[max(target[2], bond[2])]
        return energy_diffuse[target[2]] + bond_energy
    else:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[-1]
        return energy_diffuse[-1] + bond_energy


def find_filled_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 1]


def find_empty_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 0]


def find_aboves(indexes):
    return [(i[0], i[1], i[2] + 1) for i in indexes]


def find_lower_sites(indexes):
    return [(i[0], i[1], i[2] - 1) for i in indexes]


# target が移動することで、結合している周囲原子が孤立するような移動は省く
def judge_isolation(
    atom_set,
    bonds,
    target: Tuple[int, int, int],
    nn_atom: List[Tuple[int, int, int]],
    events: List[Tuple[int, int, int]],
):
    remove = []
    for check in events:
        bonds_sites = find_filled_sites(atom_set, bonds[check])
        num_bonds = len(bonds_sites)
        if num_bonds >= 2 or check[2] == 0:
            pass
        elif num_bonds == 0:
            remove.append(check)
        else:
            if bonds_sites[0] == target:
                remove.append(check)
    for i in remove:
        events.remove(i)

    for check in nn_atom:
        # 隣接原子の周囲の原子数
        nn_nn_atom = find_filled_sites(atom_set, bonds[check])
        # 二個以上なら移動後も孤立しない
        if len(nn_nn_atom) >= 2:
            pass
        # 隣接原子が基板上の時も孤立とは扱わない
        elif check[2] == 0:
            pass
        # 移動対象としか結合がない
        else:
            # 移動後の位置について
            remove = []
            for post_move in events:
                # 移動後も結合し続ける
                if post_move in bonds[check]:
                    # 移動後のサイトの隣接原子
                    post_nn_atom = find_filled_sites(atom_set, bonds[post_move])
                    # 2原子で孤立している→remove
                    if len(post_nn_atom) <= 1:
                        remove.append(post_move)
                # 結合がなくなる
                else:
                    remove.append(post_move)
            for i in remove:
                events.remove(i)

    return events


def judge_defect(target: Tuple[int, int, int], events):
    return [
        event
        for event in events
        if ((target[2] not in (0, 1)) and (event[2] in (0, 1)))
    ]


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
    #
    event_f = judge_isolation(atom_set, bonds, target, nn_atom, eve_rate)
    #
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
    trans_val: int,
):
    event_list: List[Tuple] = []
    rate_list: List[float] = []
    # states: List[int] = []
    # trans = params.trans_check
    # calculate total energy
    energy = total_energy(
        atom_set,
        bonds,
        target,
        energy_bonding,
        energy_diffuse,
        trans_val,
        highest_atom,
    )

    # calculate possible events
    event_list, rate_list = possible_events(
        atom_set, bonds, target, params, energy, diffuse_candidates
    )
    if params.limit_check:
        rate_list = rate_limit(rate_list, Decimal(float(params.limit_val)))
    return event_list, rate_list


def highest_z(atom_set):
    maxz = 1
    for index, state in atom_set.items():
        if (state != 0) and (index[2] + 1 > maxz):
            maxz = index[2] + 1
    return maxz


if __name__ == "__main__":
    from Modules.Test_modules.event_collection_check import (
        existing_atoms,
        event_check_poscar,
    )

    lattice = read_lattice()
    atom_set = read_atom_set()
    bonds = read_bonds()
    parameters = Params()
    unit_length = parameters.n_cell_init
    maxz = highest_z(atom_set)
    defect = True
    empty_first = 10
    #
    energy_bonding = [0.1 for i in range(0, 30)]
    energy_diffuse = [1 for i in range(0, 30)]

    dir_name = "Event_check/"
    os.makedirs(dir_name, exist_ok=True)
    target_cand = existing_atoms(atom_set)

    for target in target_cand:
        event_list, rate_list = site_events(
            atom_set, bonds, target, parameters, energy_bonding, energy_diffuse
        )
        #

        event_check_poscar(atom_set, event_list, lattice, unit_length, maxz, target)
    print("Poscars for event check are formed.")
