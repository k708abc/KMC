from typing import List, Dict, Tuple
from cal_rates import rate
import random
from read_examples import read_lattice, read_bonds, read_atom_set
from InputParameter import Params
from event_collection_check import random_target, event_check_poscar
import os


def bond_energy_same_state(
    target: Tuple[int, int, int], bond: Tuple, params, atom_state: int
):
    z_target = target[2]
    z_bond = bond[2]
    if atom_state == 2:
        if z_target in (1, 0) and z_bond in (1, 0):
            return params.binding_energies["Si12"]
        elif z_target in (1, 2) and z_bond in (1, 2):
            return params.binding_energies["Si23"]
        elif z_target in (2, 3) and z_bond in (2, 3):
            return params.binding_energies["Si34"]
        elif z_target in (3, 4) and z_bond in (3, 4):
            return params.binding_energies["Si45"]
        elif z_target in (4, 5) and z_bond in (4, 5):
            return params.binding_energies["Si56"]
        elif z_target in (5, 6) and z_bond in (5, 6):
            return params.binding_energies["Si_intra"]
        elif z_target % 2 == 0 and z_bond == z_target + 1:
            return params.binding_energies["Si_intra"]
        elif z_target % 2 == 0 and z_bond == z_target - 1:
            return params.binding_energies["Si_inter"]
        elif z_target % 2 == 1 and z_bond == z_target + 1:
            return params.binding_energies["Si_inter"]
        elif z_target % 2 == 1 and z_bond == z_target - 1:
            return params.binding_energies["Si_intra"]
        else:
            raise RuntimeError("Something wrong in 2D energy")
    else:
        if z_target % 2 == 0 and z_bond == z_target + 1:
            return params.binding_energies["Si_intra"]
        elif z_target % 2 == 0 and z_bond == z_target - 1:
            return params.binding_energies["Si_inter"]
        elif z_target % 2 == 1 and z_bond == z_target + 1:
            return params.binding_energies["Si_inter"]
        elif z_target % 2 == 1 and z_bond == z_target - 1:
            return params.binding_energies["Si_intra"]
        else:
            raise RuntimeError("Something wrong in 3D energy")
        return 0


def bond_energy_diff_state(target: Tuple[int, int, int], bond: Tuple, params):
    z_target = target[2]
    z_bond = bond[2]
    if z_target in (1, 0) and z_bond in (1, 0):
        return (
            params.binding_energies["Si12"] + params.binding_energies["Si_inter"]
        ) / 2
    elif z_target in (1, 2) and z_bond in (1, 2):
        return (
            params.binding_energies["Si23"] + params.binding_energies["Si_intra"]
        ) / 2
    elif z_target in (2, 3) and z_bond in (2, 3):
        return (
            params.binding_energies["Si34"] + params.binding_energies["Si_inter"]
        ) / 2
    elif z_target in (3, 4) and z_bond in (3, 4):
        return (
            params.binding_energies["Si45"] + params.binding_energies["Si_intra"]
        ) / 2
    elif z_target in (4, 5) and z_bond in (4, 5):
        return (
            params.binding_energies["Si56"] + params.binding_energies["Si_inter"]
        ) / 2
    elif z_target in (5, 6) and z_bond in (5, 6):
        return params.binding_energies["Si_intra"]
    elif z_target % 2 == 0 and z_bond == z_target + 1:
        return params.binding_energies["Si_intra"]
    elif z_target % 2 == 0 and z_bond == z_target - 1:
        return params.binding_energies["Si_inter"]
    elif z_target % 2 == 1 and z_bond == z_target + 1:
        return params.binding_energies["Si_inter"]
    elif z_target % 2 == 1 and z_bond == z_target - 1:
        return params.binding_energies["Si_intra"]
    else:
        raise RuntimeError("Something wrong in 2D energy")


def total_energy_trans(
    atom_set: Dict, bonds: Dict, target: Tuple[int, int, int], params
):
    target_z = target[2]
    target_state = atom_set[target]
    energy = 0
    if target_z == 0:
        energy += params.binding_energies["AgSi"]
    for bond in bonds[target]:
        bond_state = atom_set[bond]
        if bond_state not in (0, target_state):
            energy += bond_energy_diff_state(target, bond, params)
        elif bond_state == target_state:
            energy += bond_energy_same_state(target, bond, params, target_state)
    return energy


def total_energy_wo_trans(atom_set: Dict, bonds: Dict, target: Tuple, params):
    target_z = target[2]
    energy = 0
    if target_z == 0:
        energy += params.binding_energies["AgSi"]
    for bond in bonds[target]:
        if atom_set[bond] != 0:
            energy += bond_energy_same_state(target, bond, params, 2)
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
    rem_eve = []
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
                if post_move in bonds[check]:
                    # 移動後のサイトの隣接原子
                    post_nn_atom = find_filled_sites(atom_set, post_move)
                    # ダイマーを形成（2原子で孤立）→孤立
                    if len(post_nn_atom) <= 1:
                        rem_eve.append(post_move)
                # 移動後に孤立
                else:
                    rem_eve.append(post_move)
    return rem_eve


def judge_defect(target: Tuple[int, int, int], events):
    remove = []
    for event in events:
        if target[2] not in (0, 1) and event[2] in (0, 1):
            remove.append(event)
    return remove


def possible_events(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    params,
    energy: float,
    unit_length: int,
    defect: bool,
    empty_first: int,
    trans: bool,
):
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    events: List[Tuple] = []
    rates: List[float] = []
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
    #
    # BL内での移動
    # 最近接空きサイトへの移動
    for empty in nn_empty:
        # 最近接空きサイトの周辺原子
        nn_nn_site = find_filled_sites(atom_set, bonds[empty])
        # 移動先でも結合原子があるか、移動先がAg直上→候補
        if len(nn_nn_site) >= 2 or empty[2] == 0:
            events.append(empty)
            rates.append(rate(pre, kbt, energy))
    # 同高さの次近接への移動
    for empty in nnn_empty:
        # 次近接空きサイトの周辺原子
        nn_nnn_site = find_filled_sites(atom_set, bonds[empty])
        # 移動後に隣接原子があるか、Ag直上のサイト→候補
        if len(nn_nnn_site) >= 1 or empty[2] == 0:
            events.append(empty)
            rates.append(rate(pre, kbt, energy))
    #
    # BLの上り下り
    # BLの下層原子
    if atom_z % 2 == 0:
        # BLを上る判定
        for filled in nn_atom:
            # 最近接原子の周囲原子
            nn_nn_atom = find_filled_sites(atom_set, bonds[filled])
            # 隣接原子の上が空いていて、かつ隣接原子が別原子でも支えられている→候補
            above_site = (filled[0], filled[1], filled[2] + 1)
            if atom_set[above_site] == 0 and len(nn_nn_atom) >= 2:
                events.append(above_site)
                rates.append(rate(pre, kbt, energy))
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
            # 直下に原子があるおき
            else:
                # 直下原子周辺の空きサイト
                d_step_empty = find_empty_sites(atom_set, bonds[direct_below])
                # 直下原子の周辺空きサイト→候補
                for empty in d_step_empty:
                    events.append(empty)
                    rates.append(rate(pre, kbt, energy))
            # 次近接の下へも移動可能
            # 次近接空きサイト
            nnn_empty = find_empty_sites(atom_set, nnn_sites)
            # 次近接空きサイト下サイト
            nnn_lower_site = find_lower_sites(nnn_empty)
            # 次近接空きサイト下空きサイト
            nnn_lower_empty = find_empty_sites(atom_set, nnn_lower_site)
            for cand in nnn_lower_empty:
                # 次近接空きサイト下空きサイトの周辺原子
                cand_nn_atom = find_filled_sites(atom_set, bonds[cand])
                # 移動後に隣接原子がある→候補
                if len(cand_nn_atom) >= 1:
                    events.append(cand)
                    rates.append(rate(pre, kbt, energy))
    #
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
            for above_empty in nnn_above_empty:
                events.append(above_empty)
                rates.append(rate(pre, kbt, energy))
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
                    events.append(common_site[0])
                    rates.append(rate(pre, kbt, energy))
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
            for cand in nn_lower_enpty:
                # 候補サイトの周辺原子
                cand_nn = find_filled_sites(atom_set, bonds[cand])
                # 候補サイト移動後何らかの結合がある→候補
                if len(cand_nn) >= 1:
                    events.append(cand)
                    rates.append(rate(pre, kbt, energy))
    # イベントリストから、孤立原子を生じるイベントを抽出
    remove = judge_isolation(atom_set, bonds, target, nn_atom, events)
    # first BLにデフェクトを残す場合
    if (defect is True) and (empty_first == 1):
        remove.extend(judge_defect(target, events))
    # 削除
    event_f: List = []
    rates_f: List = []
    for eve, rat in zip(events, rates):
        if eve in remove:
            pass
        else:
            event_f.append(eve)
            rates_f.append(rat)
    return event_f, rates_f


def state_after_move(atom_set, bonds, event, params):
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    E2 = 0
    E3 = 0
    # 移動後の隣接原子の状態を確認
    for bond in bonds[event]:
        if atom_set[bond] == 2:
            E2 += bond_energy_same_state(event, bond, params, atom_set[bond])
        elif atom_set[bond] == 3:
            E3 += bond_energy_same_state(event, bond, params, atom_set[bond])

    if E2 == E3 == 0:
        E2 = -1  # 蒸着時、周辺原子がない場合2次元にする
    # 結合原子の状態ごとの速度定数逆数を計算
    rates = [rate(pre, kbt, -E2), rate(pre, kbt, -E3)]
    states = [2, 3]
    # 速度定数の大きい状態を取りやすい
    det_state = random.choices(states, weights=rates)
    return det_state[0]


def state_determinate(atom_set, bonds, event_list, params):
    # 原子が動いた後の構造を設定
    states: List[int] = []
    for event in event_list:
        state = state_after_move(atom_set, bonds, event, params)
        states.append(state)
    return states


def state_change_to_neighbor(atom_set, bonds, target: Tuple[int, int, int], params):
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    t_state = atom_set[target]
    E_change = 0
    diff_i = 0
    for bond in bonds[target]:
        if atom_set[bond] == t_state:
            E_change += bond_energy_same_state(target, bond, params, t_state)
        elif atom_set[bond] != 0:
            diff_i += 1
    change_rate = rate(pre, kbt, E_change)
    if diff_i == 0:
        return t_state, 0
    elif t_state == 2:
        return 3, change_rate
    elif t_state == 3:
        return 2, change_rate
    else:
        print("target state = " + str(t_state))
        raise RuntimeError("Some error in state change to neighbor")


def state_change_new(atom_set, bonds, target: Tuple[int, int, int], params):
    neighbor_i = 0
    target_state = atom_set[target]
    pre = float(params.prefactor)
    kbt = params.temperature_eV
    for bond in bonds[target]:
        if atom_set[bond] != 0:
            neighbor_i += 1
    # エネルギーの設定については要再考
    # 2次元から3次元：周辺原子が多いと起きやすい
    if target_state == 2:
        E_trans = (5 - neighbor_i) * params.transformation
        trans_rate = rate(pre, kbt, E_trans)
        return 3, trans_rate
    # 3次元から2次元：周辺原子が少ないと起きやすい
    elif target_state == 3:
        if neighbor_i == 0:
            neighbor_i == 1
        E_trans = neighbor_i * params.transformation
        trans_rate = rate(pre, kbt, E_trans)
        return 2, trans_rate
    else:
        raise RuntimeError("Some error in state change to neighbor")


def site_events(
    atom_set: Dict,
    bonds: Dict,
    target: Tuple[int, int, int],
    params,
    defect: bool,
    empty_first: int,
    trans: bool,
):
    event_list: List[Tuple] = []
    rate_list: List[float] = []
    unit_length = params.n_cell_init
    states: List[int] = []
    # calculate total energy
    if params.trans_check is False:
        energy = total_energy_wo_trans(atom_set, bonds, target, params)
    elif params.trans_check is True:
        energy = total_energy_trans(atom_set, bonds, target, params)

    # calculate possible events
    event_list, rate_list = possible_events(
        atom_set, bonds, target, params, energy, unit_length, defect, empty_first, trans
    )
    # event_list: List[tuple[int, int, int, int], rate_list: List[float]
    if trans is True:
        # 各イベントの構造判定
        states = state_determinate(atom_set, bonds, event_list, params)
        # サイトの変わらない構造変化の設定
        # 隣接原子の状態による変化
        state, rate = state_change_to_neighbor(atom_set, bonds, target, params)
        event_list.append(target)
        rate_list.append(rate)
        states.append(state)
        # 隣接原子の数による変化
        state, rate = state_change_new(atom_set, bonds, target, params)
        event_list.append(target)
        rate_list.append(rate)
        states.append(state)

    else:
        states = [2 for _ in range(len(event_list))]
    """
    for rate in rate_list:
        if rate > 10:
            print("too large rate constant")
            print(target)
            print(event_list)
            print(rate_list)
    """

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
    dir_name = "Event_check/"
    os.makedirs(dir_name, exist_ok=True)
    target_cand = random_target(atom_set)
    for target in target_cand:
        event_list, rate_list, states = site_events(
            atom_set, bonds, target, parameters, defect, empty_first, trans
        )
        #
        event_check_poscar(atom_set, event_list, lattice, unit_length, maxz, target)
    print("Poscars for event check are formed.")