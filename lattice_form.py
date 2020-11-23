from typing import List, Dict, Tuple
import math
from InputParameter import Params
from lattice_form_check import check

lattice_first: Dict[Tuple, List] = {}
lattice: Dict[Tuple, List] = {}
atom_set: Dict[Tuple, int] = {}
bonds: Dict[Tuple, List] = {}
event: Dict[Tuple, List] = {}
event_time: Dict[Tuple, List] = {}
event_time_tot: Dict[Tuple, float] = {}


def reset_dicts():
    lattice_first.clear()
    lattice.clear()
    atom_set.clear()
    bonds.clear()
    event.clear()
    event_time.clear()
    event_time_tot.clear()


def form_first_3BL(unit_length, zd1, zd2):
    for key in lattice_first:
        i = key[0]
        k = key[1]
        lattice_first[key] = [
            [i, k, 0],
            [i + 1 / 3.0, k + 1 / 3.0, zd1],
            [i + 1 / 3.0, k + 1 / 3.0, zd1 + zd2],
            [i + 2 / 3.0, k + 2 / 3.0, 2 * zd1 + zd2],
            [i + 2 / 3.0, k + 2 / 3.0, 2 * (zd1 + zd2)],
            [i, k, 2 * (zd1 + zd2) + zd1],  # for accuracy
        ]


def lattice_full_layers(unit_height: int):
    for key in lattice:
        key_first = (key[0], key[1])
        lattice[key] = [
            lattice_first[key_first][key[2] % 6][0],
            lattice_first[key_first][key[2] % 6][1],
            lattice_first[key_first][key[2] % 6][2] + unit_height * key[2] // 6,
        ]
    """六方晶系でSiを記述したとき、3BL分のハニカム構造がz方向の単位構造になります。
    この3BL分の単位格子の高さが2.448(nm)です。"""  # ←そういった細かい情報の結果ならば、コードにその旨書いておかないと絶対ダメ。


def neighbor_points(
    i: int, j: int, k: int, z_judge: int, unit_length: int
) -> List[List[int]]:
    neighbors: List[List[int]]
    if z_judge in (0, 2):
        neighbors = [
            [(i - 1) % unit_length, j, k + 1],
            [i, (j - 1) % unit_length, k + 1],
            [i, j, k + 1],
        ]
        if k != 0:
            neighbors.append([i, j, k - 1])
    elif z_judge in (1, 3):
        neighbors = [
            [(i + 1) % unit_length, j, k - 1],
            [i, (j + 1) % unit_length, k - 1],
            [i, j, k - 1],
            [i, j, k + 1],
        ]
    elif z_judge == 4:
        neighbors = [
            [(i + 1) % unit_length, j, k + 1],
            [(i + 1) % unit_length, (j + 1) % unit_length, k + 1],
            [i, (j + 1) % unit_length, k + 1],
            [i, j, k - 1],
        ]
    elif z_judge == 5:
        neighbors = [
            [(i - 1) % unit_length, j, k - 1],
            [i, (j - 1) % unit_length, k - 1],
            [(i - 1) % unit_length, (j - 1) % unit_length, k - 1],
            [i, j, k + 1],
        ]
    else:
        raise RuntimeError("Something wrong happens. check z_judge value")
    return neighbors


def search_bond(unit_length: int):
    # Search for bonding atoms for all the atoms
    for key in bonds:
        z_judge = key[2] % 6
        i = key[0]
        j = key[1]
        k = key[2]
        bonds[key] = neighbor_points(i, j, k, z_judge, unit_length)


def lattice_form(input_params):  # 　ここ長すぎるので、少なくとも3つか四つの関数に分ける。
    # 荒船だったら，多分5個ぐらい。
    unit_length: int = input_params.n_cell_init
    z_units: int = input_params.z_unit_init
    zd1: int = float(input_params.intra_distance)
    zd2 = float(input_params.inter_distance)
    unit_height = 3 * (zd1 + zd2)
    reset_dicts()
    #
    for i in range(unit_length):
        for j in range(unit_length):
            lattice_first[(i, j)] = []
            for k in range(z_units * 6):
                lattice[(i, j, k)] = []
                atom_set[(i, j, k)] = 0
                bonds[(i, j, k)] = []
                event[(i, j, k)] = []
                event_time[(i, j, k)] = []
                event_time_tot[(i, j, k)] = 0.0
    #
    form_first_3BL(unit_length, zd1, zd2)
    #
    lattice_full_layers(unit_height)
    #
    search_bond(unit_length)
    return lattice, bonds, atom_set, event, event_time, event_time_tot


if __name__ == "__main__":
    init_values = Params()
    lattice_formed = lattice_form(init_values)
    check(init_values, lattice_formed[0], lattice_formed[1])
