from typing import List
from cal_rates import rate


def get_energy(atom_set, bonds, target, params):
    z_judge = target[2]
    energy = 0
    if z_judge == 0:
        energy += params.binding_energies["AgSi"]
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                energy += params.binding_energies["Si12"]
    elif z_judge == 1:
        energy += params.binding_energies["AgSi"]
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                if bond[2] == 0:
                    energy += params.binding_energies["Si12"]
                elif bond[2] == 2:
                    energy += params.binding_energies["Si23"]
    elif z_judge == 2:
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                if bond[2] == 1:
                    energy += params.binding_energies["Si23"]
                elif bond[2] == 3:
                    energy += params.binding_energies["Si34"]

    elif z_judge == 3:
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                if bond[2] == 2:
                    energy += params.binding_energies["Si34"]
                elif bond[2] == 4:
                    energy += params.binding_energies["Si45"]

    elif z_judge == 4:
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                if bond[2] == 3:
                    energy += params.binding_energies["Si45"]
                elif bond[2] == 5:
                    energy += params.binding_energies["Si56"]

    elif z_judge == 5:
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                if bond[2] == 4:
                    energy += params.binding_energies["Si56"]
                elif bond[2] == 6:
                    energy += params.binding_energies["Si_intra"]

    elif z_judge % 2 == 0:
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                if bond[2] == z_judge - 1:
                    energy += params.binding_energies["Si_intra"]
                elif bond[2] == z_judge + 1:
                    energy += params.binding_energies["Si_inter"]

    elif z_judge % 2 == 1:
        for bond in bonds[target]:
            if atom_set[bond] != 0:
                if bond[2] == z_judge - 1:
                    energy += params.binding_energies["Si_inter"]
                elif bond[2] == z_judge + 1:
                    energy += params.binding_energies["Si_intra"]
    return energy


def possible_events(atom_set, bonds, target, params, energy, unit_length):
    events: List[tuple] = []
    rates: List[float] = []
    #
    atom_x = target[0]
    atom_y = target[1]
    atom_z = target[2]
    # nnn: next nearest neighbor
    site_nnn = [
        ((atom_x - 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y - 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y + 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, (atom_y - 1) % unit_length, atom_z),
        ((atom_x - 1) % unit_length, (atom_y + 1) % unit_length, atom_z),
        ]

    if z_judge == 0:    # 0層目
        nn_atoms = []
        nn_nn_atoms = []

        for bond_nn in bonds[target]:
            

        for bond_nn in bonds[target]:
            if atom_set[bond_nn] == 0:  # empty neighbor site
                events.append(bond_nn)
                rates.append(rate(energy))
            else:  # filled neighbor
                # 条件を満たすと上の層に登れる
                if atom_set[(bond_nn[0], bond_nn[1], bond_nn[2] + 1)] == 0:
                    # 隣接原子の上に原子がない
                    n_neighbor = 0
                    for bond_nnn in bonds[bond_nn]:
                        if atom_set[bond_nnn] != 0:
                            n_neighbor += 1
                    if n_neighbor >= 2:
                        # 隣接原子が二つ以上の結合原子を持つ
                        # 隣接原子の真上が候補
                        events.append((bond_nn[0], bond_nn[1], bond_nn[2] + 1))
                        rates.append(rate(energy))











    elif z_judge%2 == 1:
        for bond_nn in bonds[target]:
            if atom_set[bond_nn] == 0 and bond_nn[2] == 0:
                # empty neighbor site in lower layer
                events.append(bond_nn)
                rates.append(rate(energy))
            elif atom_set[bond_nn] == 0 and bond_nn[2] == 0:
                # empty neighbor in upperlayer
                n_neighbor = 0
                for bond_nnn in bonds[bond_nn]:
                    if atom_set[bond_nnn] != 0:
                        n_neighbor += 1
                if n_neighbor >= 2:
                    # 真上のサイトが二つ以上の隣接原子を持つ
                    # 真上が候補
                    events.append((bond_nn[0], bond_nn[1], bond_nn[2] + 1))
                    rates.append(rate(energy))






    elif z_judge == 2:


def get_events(atom_set, bonds, target, params):
    events: List[tuple] = []
    rates: List[float] = []
    unit_length = param.n_cell_init
    # calculate total energy
    energy = get_energy(atom_set, bonds, target, params)
    # calculate possible events
    event, rates = possible_events(atom_set, bonds, target, params, energy, unit_length)

    return events, rates
