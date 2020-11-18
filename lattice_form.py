from typing import List, Dict
import math
from InputParameter import Params
from lattice_form_check import check

lattice: Dict[str, list] = {}
atom_set: Dict[str, int] = {}
bonds: Dict[str, list] = {}
event: Dict[str, list] = {}
event_time: Dict[str, list] = {}
event_time_tot: Dict[str, float] = {}

def form_first_3BL(unit_length, zd1, zd2):
    lattice_first: List[float] = [
        [
            [
                [i, k, 0],
                [i + 1 / 3.0, k + 1 / 3.0, zd1],
                [i + 1 / 3.0, k + 1 / 3.0, zd1 + zd2],
                [i + 2 / 3.0, k + 2 / 3.0, 2 * zd1 + zd2],
                [i + 2 / 3.0, k + 2 / 3.0, 2 * (zd1 + zd2)],
                [i, k, 2 * (zd1 + zd2) + zd1],    # for accuracy
            ]
            for k in range(unit_length)
        ]
        for i in range(unit_length)
    ]
    return lattice_first

def lattice_full_layers(unit_length, z_units, lattice_first):
    for i in range(unit_length):
        for j in range(unit_length):
            for k in range(z_units):
                for l, first in enumerate(lattice_first[i][j]):
                    atom_index = (i, j, k*6+l)  # it's better to use tuple than str
                    lattice[atom_index] = [
                        round(first[0], 5),
                        round(first[1], 5),
                        round(first[2] + k * 2.448, 5),  # what does this 2.448 mean?
                    ]
                    atom_set[atom_index] = 0
                    event[atom_index] = []
                    event_time[atom_index] = []
                    event_time_tot[atom_index] = 0

def search_bond(unit_length, maxz):
    # Search for bonding atoms for all the atoms
    for i in range(unit_length):    #ここのfor loop意味が分からない。何してるの？
        for j in range(unit_length):
            for k in range(maxz):
                bond_with = []  ## what does "bond_with" mean?
                z_judge = k % 6
                ul_m = unit_length - 1
                if z_judge == 0:
                    if i == 0 and j == 0:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == 0:
                        bond_with = [[ul_m, j, k + 1], [0, j - 1, k + 1], [0, j, k + 1]]
                    elif j == 0:
                        bond_with = [[i - 1, 0, k + 1], [i, ul_m, k + 1], [i, 0, k + 1]]
                    else:
                        bond_with = [
                            [i - 1, j, k + 1],
                            [i, j - 1, k + 1],
                            [i, j, k + 1],
                        ]
                    if k != 0:
                        bond_with.append([i, j, k - 1])
                    else:
                        pass
                elif z_judge == 1:
                    if i == j == ul_m:
                        bond_with = [
                            [ul_m, 0, k - 1],
                            [0, ul_m, k - 1],
                            [ul_m, ul_m, k - 1],
                        ]
                    elif i == ul_m:
                        bond_with = [
                            [ul_m, j + 1, k - 1],
                            [0, j, k - 1],
                            [ul_m, j, k - 1],
                        ]
                    elif j == ul_m:
                        bond_with = [
                            [i + 1, ul_m, k - 1],
                            [i, 0, k - 1],
                            [i, ul_m, k - 1],
                        ]
                    else:
                        bond_with = [
                            [i + 1, j, k - 1],
                            [i, j + 1, k - 1],
                            [i, j, k - 1],
                        ]
                    bond_with.append([i, j, k + 1])
                elif z_judge == 2:
                    if i == j == 0:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == 0:
                        bond_with = [[ul_m, j, k + 1], [0, j - 1, k + 1], [0, j, k + 1]]
                    elif j == 0:
                        bond_with = [[i - 1, 0, k + 1], [i, ul_m, k + 1], [i, 0, k + 1]]
                    else:
                        bond_with = [
                            [i - 1, j, k + 1],
                            [i, j - 1, k + 1],
                            [i, j, k + 1],
                        ]
                    bond_with.append([i, j, k - 1])
                elif z_judge == 3:
                    if i == j == ul_m:
                        bond_with = [
                            [ul_m, 0, k - 1],
                            [0, ul_m, k - 1],
                            [ul_m, ul_m, k - 1],
                        ]
                    elif i == ul_m:
                        bond_with = [
                            [ul_m, j + 1, k - 1],
                            [0, j, k - 1],
                            [ul_m, j, k - 1],
                        ]
                    elif j == ul_m:
                        bond_with = [
                            [i + 1, ul_m, k - 1],
                            [i, 0, k - 1],
                            [i, ul_m, k - 1],
                        ]
                    else:
                        bond_with = [
                            [i + 1, j, k - 1],
                            [i, j + 1, k - 1],
                            [i, j, k - 1],
                        ]
                    bond_with.append([i, j, k + 1])
                elif z_judge == 4:
                    if i == j == ul_m:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == ul_m:
                        bond_with = [
                            [0, j, k + 1],
                            [0, j + 1, k + 1],
                            [ul_m, j + 1, k + 1],
                        ]
                    elif j == ul_m:
                        bond_with = [
                            [i + 1, 0, k + 1],
                            [i, 0, k + 1],
                            [i + 1, ul_m, k + 1],
                        ]
                    else:
                        bond_with = [
                            [i + 1, j, k + 1],
                            [i + 1, j + 1, k + 1],
                            [i, j + 1, k + 1],
                        ]
                    bond_with.append([i, j, k - 1])
                elif z_judge == 5:
                    if i == j == 0:
                        bond_with = [
                            [ul_m, ul_m, k - 1],
                            [0, ul_m, k - 1],
                            [ul_m, 0, k - 1],
                        ]
                    elif i == 0:
                        bond_with = [
                            [0, j - 1, k - 1],
                            [ul_m, j, k - 1],
                            [ul_m, j - 1, k - 1],
                        ]
                    elif j == 0:
                        bond_with = [
                            [i, ul_m, k - 1],
                            [i - 1, ul_m, k - 1],
                            [i - 1, 0, k - 1],
                        ]
                    else:
                        bond_with = [
                            [i - 1, j, k - 1],
                            [i, j - 1, k - 1],
                            [i - 1, j - 1, k - 1],
                        ]
                    bond_with.append([i, j, k + 1])
                #
                atom_index = (i, j, k)
                bonds[atom_index] = bond_with


def lattice_form(input_params):  # 　ここ長すぎるので、少なくとも3つか四つの関数に分ける。
    # 荒船だったら，多分5個ぐらい。
    unit_length = input_params.n_cell_init
    z_units = input_params.z_unit_init
    maxz = z_units * 6 - 1
    zd1 = float(input_params.intra_distance)
    zd2 = float(input_params.inter_distance)
    #
    lattice_first = form_first_3BL(unit_length, zd1, zd2)
    #
    lattice_full_layers(unit_length, z_units, lattice_first)
    #
    search_bond(unit_length, maxz)
    return lattice, bonds, atom_set, event, event_time, event_time_tot



if __name__ == "__main__":
    init_values= Params()
    lattice_formed = lattice_form(init_values)
    check(init_values, lattice_formed[0], lattice_formed[1])
