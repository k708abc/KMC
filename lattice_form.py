from typing import List, Dict
import math
from InputParameter import Params
from lattice_form_check import check
import itertools

lattice: Dict[tuple, list] = {}
atom_set: Dict[tuple, int] = {}
bonds: Dict[tuple, list] = {}
event: Dict[tuple, list] = {}
event_time: Dict[tuple, list] = {}
event_time_tot: Dict[tuple, float] = {}


def form_first_3BL(unit_length, zd1, zd2):
    lattice_first: List[float] = [
        [
            [
                [i, k, 0],
                [i + 1 / 3.0, k + 1 / 3.0, zd1],
                [i + 1 / 3.0, k + 1 / 3.0, zd1 + zd2],
                [i + 2 / 3.0, k + 2 / 3.0, 2 * zd1 + zd2],
                [i + 2 / 3.0, k + 2 / 3.0, 2 * (zd1 + zd2)],
                [i, k, 2 * (zd1 + zd2) + zd1],  # for accuracy
            ]
            for k in range(unit_length)
        ]
        for i in range(unit_length)
    ]
    return lattice_first


def lattice_full_layers(lattice_first):
    for key in lattice:
        lattice[key] = [
            lattice_first[key[0]][key[1]][key[2] % 6][0],
            lattice_first[key[0]][key[1]][key[2] % 6][1],
            lattice_first[key[0]][key[1]][key[2] % 6][2] + 2.448 * key[2] // 6,
        ]
    """六方晶系でSiを記述したとき、3BL分のハニカム構造がz方向の単位構造になります。
    この3BL分の単位格子の高さが2.448(nm)です。"""  # ←そういった細かい情報の結果ならば、コードにその旨書いておかないと絶対ダメ。


def search_bond(unit_length, maxz):
    # Search for bonding atoms for all the atoms
    for key in bonds:
        bond_with = []
        z_judge = key[2] % 6
        ul_m = unit_length - 1
        i = key[0]
        j = key[1]
        k = key[2]

        if z_judge == 0:
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
        bonds[key] = bond_with


def lattice_form(input_params):  # 　ここ長すぎるので、少なくとも3つか四つの関数に分ける。
    # 荒船だったら，多分5個ぐらい。
    unit_length = input_params.n_cell_init
    z_units = input_params.z_unit_init
    maxz = z_units * 6 - 1
    zd1 = float(input_params.intra_distance)
    zd2 = float(input_params.inter_distance)
    #
    for i in range(unit_length):
        for j in range(unit_length):
            for k in range(z_units * 6):
                lattice[(i, j, k)] = []
                atom_set[(i, j, k)] = 0
                bonds[(i, j, k)] = []
                event[(i, j, k)] = []
                event_time[(i, j, k)] = []
                event_time_tot[(i, j, k)] = 0.0
    #
    lattice_first = form_first_3BL(unit_length, zd1, zd2)
    #
    lattice_full_layers(lattice_first)
    #
    search_bond(unit_length, maxz)
    return lattice, bonds, atom_set, event, event_time, event_time_tot


if __name__ == "__main__":
    init_values = Params()
    lattice_formed = lattice_form(init_values)
    check(init_values, lattice_formed[0], lattice_formed[1])
