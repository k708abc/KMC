from typing import List, Dict, Tuple
from Modules.cal_rates import rate
import random
from Test_modules.read_examples import read_lattice, read_bonds, read_atom_set
from InputParameter import Params
from Test_modules.event_collection_check import random_target, event_check_poscar
import os
import time
import decimal

"""
# バイレイヤー内のみを考慮
def total_energy(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    energy_bonding: List[float],
    energy_diffuse: List[float],
):
    num_bond = 0
    for bond in bonds[target]:
        if (atom_set[bond] == 1) and (target[2] // 2 == bond[2] // 2):
            num_bond += 1
    return (
        energy_diffuse[target[2] // 2] + num_bond * energy_bonding[int(target[2] // 2)]
    )
"""


# バイレイヤー間も含める
def total_energy(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    energy_bonding: List[float],
    energy_diffuse: List[float],
):
    num_bond = 0
    for bond in bonds[target]:
        if atom_set[bond] == 1:
            num_bond += 1
    # Ag-Siも結合数に含める
    if target[0] == 0:
        num_bond += 1
    return (
        energy_diffuse[target[2] // 2] + num_bond * energy_bonding[int(target[2] // 2)]
    )


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
    events: List[Tuple[Tuple[int, int, int], float]],
):
    for check in nn_atom:
        # 隣接原子の周囲の原子数
        nn_nn_atom = find_filled_sites(atom_set, bonds[check])
        # 二個以上なら移動後も孤立しない
        if len(nn_nn_atom) >= 2:
            pass
        # 移動対象としか結合がない
        else:
            # 移動後の位置について
            for post_move in events:
                # 移動後も結合し続ける
                if post_move[0] in bonds[check]:
                    # 移動後のサイトの隣接原子
                    post_nn_atom = find_filled_sites(atom_set, bonds[post_move[0]])
                    # 2原子で孤立している→remove
                    if len(post_nn_atom) <= 1:
                        events.remove(post_move)
                # 結合がなくなる
                else:
                    events.remove(post_move)
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
):
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    eve_rate: List = []
    rearange_rate = decimal.Decimal(rate(pre, kbt, energy))
    unit_length = params.n_cell_init
    #
    atom_x, atom_y, atom_z = target
    # nnn: next nearest neighbor
    nnn_sites = [
        ((atom_x - 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y - 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y + 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, (atom_y - 1) % unit_length, atom_z),
        ((atom_x - 1) % unit_length, (atom_y + 1) % unit_length, atom_z),
    ]
    # 最近接の空きサイト
    nn_empty = find_empty_sites(atom_set, bonds[target])
    # 次近接の空きサイト
    nnn_empty = find_empty_sites(atom_set, nnn_sites)
    # 最近接の存在原子
    nn_atom = find_filled_sites(atom_set, bonds[target])
    # BL内での移動
    # 最近接空きサイトへの移動
    eve_rate += [
        (empty, rearange_rate)
        for empty in nn_empty
        if (len(find_filled_sites(atom_set, bonds[empty])) >= 2 or empty[2] == 0)
    ]
    # 同高さの次近接への移動
    eve_rate += [
        (empty, rearange_rate)
        for empty in nnn_empty
        if (len(find_filled_sites(atom_set, bonds[empty])) >= 1 or empty[2] == 0)
    ]
    #
    # BLの上り下り
    # BLの下層原子
    if atom_z % 2 == 0:
        # BLを上る判定
        eve_rate += [
            ((filled[0], filled[1], filled[2] + 1), rearange_rate)
            for filled in nn_atom
            if (
                atom_set[(filled[0], filled[1], filled[2] + 1)] == 0
                and len(find_filled_sites(atom_set, bonds[filled])) >= 2
            )
        ]
        # BLを下る判定
        # Ag直上原子は下れない
        if atom_z == 0:
            pass
        else:
            # 直下原子のインデックス
            direct_below = (atom_x, atom_y, atom_z - 1)
            # 直下に原子がないとき
            if (
                atom_set[direct_below] == 0
                and len(find_filled_sites(atom_set, bonds[direct_below])) >= 2
            ):
                pass
            # 直下に原子があるとき
            else:
                eve_rate += [
                    (empty, rearange_rate)
                    for empty in find_empty_sites(atom_set, bonds[direct_below])
                ]

            # 次近接の下へも移動可能
            # 次近接空きサイト
            nnn_empty = find_empty_sites(atom_set, nnn_sites)
            # 次近接空きサイト下サイト
            nnn_lower_site = find_lower_sites(nnn_empty)
            # 次近接空きサイト下空きサイト
            nnn_lower_empty = find_empty_sites(atom_set, nnn_lower_site)
            #
            eve_rate += [
                (cand, rearange_rate)
                for cand in nnn_lower_empty
                if len(find_filled_sites(atom_set, bonds[cand])) >= 1
            ]
    # BLの上層原子
    # 次近接の真上か、最近接真上と次近接真上の共通サイト
    elif atom_z % 2 == 1:
        # BLを上る判定
        # 直上原子のインデックス
        direct_above = (atom_x, atom_y, atom_z + 1)
        # 直上原子があるときpass
        if atom_set[direct_above] == 1:
            pass
        else:
            # 次近接の原子
            nnn_atoms = find_filled_sites(atom_set, nnn_sites)
            # 次近接の直上サイト
            nnn_above_site = find_aboves(nnn_atoms)
            # 次近接の直上の空きサイト
            nnn_above_empty = find_empty_sites(atom_set, nnn_above_site)
            # 次近接直上の原子
            nnn_above_atoms = find_filled_sites(atom_set, nnn_above_site)
            # 次近接原子直上の空きサイト→候補
            eve_rate += [
                (above_empty, rearange_rate) for above_empty in nnn_above_empty
            ]

            # 次近接直上に原子がある場合
            for above_filled in nnn_above_atoms:
                # 次近接直上の結合サイト
                above_filled_nn = bonds[above_filled]
                # 対象原子直上の結合サイト
                dir_above_nn = bonds[direct_above]
                # 共通項のサイトの抜出
                common_site = list(set(above_filled_nn) & set(dir_above_nn))
                # 共通サイトに原子がない場合→候補
                if atom_set[common_site[0]] == 0:
                    eve_rate.append((common_site[0], rearange_rate))

        #
        # BLを下る判定
        if atom_z == 1:
            pass
        else:
            # 最近接の空きサイト
            nn_empty = find_empty_sites(atom_set, bonds[target])
            # 最近接空きサイトの下サイト
            nn_lower = find_lower_sites(nn_empty)
            # 最近接空きサイトの下空きサイト
            nn_lower_enpty = find_empty_sites(atom_set, nn_lower)
            #
            eve_rate += [
                (cand, rearange_rate)
                for cand in nn_lower_enpty
                if len(find_filled_sites(atom_set, bonds[cand])) >= 1
            ]

    # イベントリストから、孤立原子を生じるイベントを抽出
    eve_rate = judge_isolation(atom_set, bonds, target, nn_atom, eve_rate)
    #
    event_f: List = [eves[0] for eves in eve_rate]
    rates_f: List = [eves[1] for eves in eve_rate]
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
):
    event_list: List[Tuple] = []
    rate_list: List[float] = []
    # states: List[int] = []
    # trans = params.trans_check
    # calculate total energy
    energy = total_energy(atom_set, bonds, target, energy_bonding, energy_diffuse)

    # calculate possible events
    event_list, rate_list = possible_events(
        atom_set,
        bonds,
        target,
        params,
        energy,
    )
    if params.limit_check:
        rate_list = rate_limit(rate_list, decimal.Decimal(float(params.limit_val)))
    return event_list, rate_list


def highest_z(atom_set):
    maxz = 1
    for index, state in atom_set.items():
        if (state != 0) and (index[2] + 1 > maxz):
            maxz = index[2] + 1
    return maxz


if __name__ == "__main__":
    lattice = read_lattice()
    atom_set = read_atom_set()
    bonds = read_bonds()
    parameters = Params()
    unit_length = parameters.n_cell_init
    maxz = highest_z(atom_set)
    trans = True
    defect = True
    empty_first = 10
    #
    start = time.time()

    dir_name = "Event_check/"
    os.makedirs(dir_name, exist_ok=True)
    target_cand = random_target(atom_set)
    for target in target_cand:
        event_list, rate_list = site_events(
            atom_set,
            bonds,
            target,
            parameters,
        )
        #
        event_check_poscar(atom_set, event_list, lattice, unit_length, maxz, target)
    print("Poscars for event check are formed.")
    end_time = time.time() - start
    print(end_time)
