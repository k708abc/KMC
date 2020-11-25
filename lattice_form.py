from typing import List, Dict, Tuple
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


def form_first_3BL(unit_length: int, z_intra: float, z_inter: float):
    for index in lattice_first:
        x, y = index
        lattice_first[index] = [
            [x, y, 0],
            [x + 1 / 3.0, y + 1 / 3.0, z_intra],
            [x + 1 / 3.0, y + 1 / 3.0, z_intra + z_inter],
            [x + 2 / 3.0, y + 2 / 3.0, 2 * z_intra + z_inter],
            [x + 2 / 3.0, y + 2 / 3.0, 2 * (z_inter + z_intra)],
            [x, y, 2 * (z_inter + z_intra) + z_intra],  # for accuracy
        ]


def lattice_full_layers(unit_height: int):
    for index in lattice:
        index_xy = (index[0], index[1])
        lattice[index] = [
            lattice_first[index_xy][index[2] % 6][0],
            lattice_first[index_xy][index[2] % 6][1],
            lattice_first[index_xy][index[2] % 6][2] + unit_height * (index[2] // 6),
        ]
    """六方晶系でSiを記述したとき、3BL分のハニカム構造がz方向の単位構造になります。
    この3BL分の単位格子の高さが2.448(nm)です。"""  # ←そういった細かい情報の結果ならば、コードにその旨書いておかないと絶対ダメ。


def neighbor_points(
    index: tuple, z_judge: int, unit_length: int, z_max: int
) -> List[List[int]]:
    neighbors: List[List[int]]
    x, y, z = index
    if z == z_max:
        neighbors = [
            ((x - 1) % unit_length, y, z - 1),
            (x, (y - 1) % unit_length, z - 1),
            ((x - 1) % unit_length, (y - 1) % unit_length, z - 1),
        ]
    elif z_judge in (0, 2):
        neighbors = [
            ((x - 1) % unit_length, y, z + 1),
            (x, (y - 1) % unit_length, z + 1),
            (x, y, z + 1),
        ]
        if z != 0:
            neighbors.append((x, y, z - 1))
    elif z_judge in (1, 3):
        neighbors = [
            ((x + 1) % unit_length, y, z - 1),
            (x, (y + 1) % unit_length, z - 1),
            (x, y, z - 1),
            (x, y, z + 1),
        ]
    elif z_judge == 4:
        neighbors = [
            ((x + 1) % unit_length, y, z + 1),
            ((x + 1) % unit_length, (y + 1) % unit_length, z + 1),
            (x, (y + 1) % unit_length, z + 1),
            (x, y, z - 1),
        ]
    elif z_judge == 5:
        neighbors = [
            ((x - 1) % unit_length, y, z - 1),
            (x, (y - 1) % unit_length, z - 1),
            ((x - 1) % unit_length, (y - 1) % unit_length, z - 1),
            (x, y, z + 1),
        ]
    else:
        raise RuntimeError("Something wrong happens. check z_judge value")
    return neighbors


def search_bond(unit_length: int, z_max: int):
    # Search for bonding atoms for all the atoms
    for index in bonds:
        z_judge = index[2] % 6
        bonds[index] = neighbor_points(index, z_judge, unit_length, z_max)


def lattice_form(input_params):  # 　ここ長すぎるので、少なくとも3つか四つの関数に分ける。
    # 荒船だったら，多分5個ぐらい。
    unit_length: int = input_params.n_cell_init
    z_units: int = input_params.z_unit_init
    z_intra: int = float(input_params.intra_distance)
    z_inter = float(input_params.inter_distance)
    unit_height = 3 * (z_intra + z_inter)
    reset_dicts()
    z_max = z_units * 6 - 1
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
    form_first_3BL(unit_length, z_intra, z_inter)
    #
    lattice_full_layers(unit_height)
    #
    search_bond(unit_length, z_max)
    return lattice, bonds, atom_set, event, event_time, event_time_tot


if __name__ == "__main__":
    init_values = Params()
    lattice_formed = lattice_form(init_values)
    check(init_values, lattice_formed[0], lattice_formed[1])
