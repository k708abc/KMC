from typing import List, Dict, Tuple
from Modules.cal_rates import rate
import random
from Test_modules.read_examples import read_lattice, read_bonds, read_atom_set
from InputParameter import Params
from Test_modules.event_collection_check import random_target, event_check_poscar
import os
import time
import decimal


def bond_energy_same_state(
    target: Tuple[int, int, int],
    atom_state: int,
    energy2D,
    energy3D,
):
    if atom_state == 2:
        return energy2D[int(target[2] // 2)]
    else:
        return energy3D[int(target[2] // 2)]


def bond_energy_diff_state(target: Tuple[int, int, int], energy2D3D: List[float]):
    return energy2D3D[int(target[2] // 2)]


def total_energy_trans(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    energy2D: List[float],
    energy2D3D: List[float],
    energy3D: List[float],
    diffuse_energy: List[float],
):
    target_state = atom_set[target]
    energy = diffuse_energy[target[2] // 2]  # base energy for diffusion
    for bond in bonds[target]:
        bond_state = atom_set[bond]
        if bond_state == 0:
            pass
        elif target[2] // 2 != bond[2] // 2:
            pass
        elif bond_state != target_state:
            energy += bond_energy_diff_state(target, energy2D3D)
        elif bond_state == target_state:
            energy += bond_energy_same_state(target, target_state, energy2D, energy3D)
    return energy


def total_energy_wo_trans(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    energy2D: List[float],
    diffuse_energy: List[float],
):
    energy = diffuse_energy[target[2] // 2]
    num_bond = 0
    for bond in bonds[target]:
        if (atom_set[bond] != 0) and (target[2] // 2 == bond[2] // 2):
            num_bond += 1
    energy += num_bond * bond_energy_same_state(target, 2, energy2D, 0)
    return energy


def find_filled_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] != 0]


def find_empty_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 0]


def find_aboves(indexes):
    return [(i[0], i[1], i[2] + 1) for i in indexes]


def find_lower_sites(indexes):
    return [(i[0], i[1], i[2] - 1) for i in indexes]


# 移動後に原子が孤立するような移動は省く
def judge_isolation(atom_set, bonds, target: Tuple[int, int, int], nn_atom, events):
    for check in nn_atom:
        # 隣接原子の周囲の原子数
        nn_nn_atom = find_lower_sites(bonds[check])
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
                    post_nn_atom = find_filled_sites(atom_set, post_move)
                    # 2原子で孤立している→remove
                    if len(post_nn_atom) <= 1:
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
    # BLの上り下り
    rearange_rate = decimal.Decimal(rate(pre, kbt, energy))
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
            # 直下に原子がないときパス
            if atom_set[direct_below] == 0:
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
        if atom_set[direct_above] != 0:
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


def state_after_move(atom_set, bonds, event, params, energy2D, energy3D):
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    E2 = 0
    E3 = 0
    # 移動後の隣接原子の状態を確認
    for bond in bonds[event]:
        if atom_set[bond] == 2:
            E2 += bond_energy_same_state(event, atom_set[bond], energy2D, energy3D)
        elif atom_set[bond] == 3:
            E3 += bond_energy_same_state(event, atom_set[bond], energy2D, energy3D)

    if E2 == E3 == 0:
        return 2
    elif E2 == 0:
        return 3
    elif E3 == 0:
        return 2
    else:
        # 結合原子の状態ごとの速度定数逆数を計算
        rates = [rate(pre, kbt, -E2), rate(pre, kbt, -E3)]
        states = [2, 3]
        # 速度定数の大きい状態を取りやすい
        det_state = random.choices(states, weights=rates)
        return det_state[0]


def state_determinate(
    atom_set, bonds, event_list, params, energy2D, energy2D3D, energy3D
):
    # 原子が動いた後の構造を設定
    return [
        state_after_move(atom_set, bonds, event, params, energy2D, energy3D)
        for event in event_list
    ]


def state_change_to_neighbor(
    atom_set, bonds, target: Tuple[int, int, int], params, energy2D, energy3D
):
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    t_state = atom_set[target]
    E_change = 0
    diff_i = 0
    num_bond = 0
    for bond in bonds[target]:
        if atom_set[bond] == t_state:
            E_change += bond_energy_same_state(
                target, bond, t_state, energy2D, energy3D
            )
            num_bond += 1
        elif atom_set[bond] != 0:
            diff_i += 1
            num_bond += 1
    change_rate = decimal.Decimal(rate(pre, kbt, E_change))
    if num_bond == 2 and diff_i == 1:
        return 3, 0
    elif diff_i == 0:
        return t_state, 0

    elif t_state == 2:
        return 3, change_rate

    elif t_state == 3:
        return 2, change_rate

    else:
        print("target state = " + str(t_state))
        raise RuntimeError("Some error in state change to neighbor")


def state_change_new(atom_set, bonds, target: Tuple[int, int, int], params):
    target_state = atom_set[target]
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    neighbor_i = 0
    neighbor_2 = 0
    for bond in bonds[target]:
        if atom_set[bond] != 0:
            neighbor_i += 1
        if atom_set[bond] == 2:
            neighbor_2 += 1
    #
    if target_state == 2 and neighbor_i == 4:
        return 4, decimal.Decimal(rate(pre, kbt, params.transformation))
    elif target_state == 3 and target[2] != 0 and neighbor_i != 4 and neighbor_2 != 0:
        return 5, decimal.Decimal(rate(pre, kbt, params.transformation))
    else:
        return target_state, 0

    """
    neighbor_i = 0
    target_state = atom_set[target]
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    for bond in bonds[target]:
        if atom_set[bond] != 0:
            neighbor_i += 1
    # エネルギーの設定については要再考
    # 2次元から3次元：周辺原子が多いと起きやすい
    if target_state == 2 and target[2] >= 3:
        E_trans = (5 - neighbor_i) * params.transformation
        trans_rate = decimal.Decimal(rate(pre, kbt, E_trans))
        return 4, trans_rate
    # 3次元から2次元：周辺原子が少ないと起きやすい
    elif target_state == 3:
        if neighbor_i == 0:
            neighbor_i == 1
        E_trans = neighbor_i * params.transformation
        trans_rate = decimal.Decimal(rate(pre, kbt, E_trans))
        return 2, trans_rate
    elif target_state == 2:
        return 2, decimal.Decimal(0)
    else:
        raise RuntimeError("Some error in state change to neighbor")
    """


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
    energy2D: List[float],
    energy2D3D: List[float],
    energy3D: List[float],
    diffuse_energy: List[float],
):
    event_list: List[Tuple] = []
    rate_list: List[float] = []
    states: List[int] = []
    trans = params.trans_check
    # calculate total energy
    if trans is False:
        energy = total_energy_wo_trans(
            atom_set, bonds, target, energy2D, diffuse_energy
        )

    else:
        energy = total_energy_trans(
            atom_set, bonds, target, energy2D, energy2D3D, energy3D, diffuse_energy
        )
    # calculate possible events
    event_list, rate_list = possible_events(
        atom_set,
        bonds,
        target,
        params,
        energy,
    )

    # event_list: List[tuple[int, int, int, int], rate_list: List[float]
    if trans is True:
        # 各イベントの構造判定
        states = state_determinate(
            atom_set, bonds, event_list, params, energy2D, energy2D3D, energy3D
        )
        # サイトの変わらない構造変化の設定
        """
        # 隣接原子の状態による変化
        state, rate = state_change_to_neighbor(atom_set, bonds, target, params, energy2D, energy3D)
        event_list.append(target)
        rate_list.append(rate)
        states.append(state)
        """

        # 隣接原子の数による変化
        state, rate = state_change_new(atom_set, bonds, target, params)
        event_list.append(target)
        rate_list.append(rate)
        states.append(state)

    else:
        states = [2 for _ in range(len(event_list))]

    if params.limit_check is True:
        rate_list = rate_limit(rate_list, decimal.Decimal(float(params.limit_val)))

    return event_list, rate_list, states


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
        event_list, rate_list, states = site_events(
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
