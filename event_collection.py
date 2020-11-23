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


def get_filled_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] != 0]


def get_empty_sites(atom_set, indexes):
    return [atom_nn for atom_nn in indexes if atom_set[atom_nn] == 0]


def get_aboves(indexes):
    return [(i[0], i[1], i[2] + 1) for i in indexes]


def get_lowers_sites(indexes):
    return [(i[0], i[1], i[2] - 1) for i in indexes]


# 移動後に原子が孤立するような移動は省く
def judge_isolation(atom_set, bonds, target, nn_atom, events):
    rem_eve = []
    for check in nn_atom:
        # 隣接原子の周囲の原子数
        nn_nn_atom = get_lowers_sites(atom_set, bonds[check])
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
                    post_nn_atom = get_filled_sites(atom_set, post_move)
                    # ダイマーを形成（2原子で孤立）→孤立
                    if len(post_nn_atom) <= 1:
                        rem_eve.append(post_move)
                # 移動後に孤立
                else:
                    rem_eve.append(post_move)
    return rem_eve


def possible_events(atom_set, bonds, target, params, energy, unit_length):
    events: List[tuple] = []
    rates: List[float] = []
    #
    atom_x = target[0]
    atom_y = target[1]
    atom_z = target[2]
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
    nn_empty = get_empty_sites(atom_set, bonds[target])
    # 次近接の空きサイト
    nnn_empty = get_empty_sites(atom_set, nnn_sites)
    # 最近接の存在原子
    nn_atom = get_filled_sites(atom_set, bonds[target])
    #
    # BL内での移動
    # 最近接空きサイトへの移動
    for empty in nn_empty:
        # 最近接空きサイトの周辺原子
        nn_nn_site = get_filled_sites(atom_set, bonds[empty])
        # 移動先でも結合原子があるか、移動先がAg直上→候補
        if len(nn_nn_site) >= 2 or empty[2] == 0:
            events.append(empty)
            rates.append(rate(energy))
    # 同高さの次近接への移動
    for empty in nnn_empty:
        # 次近接空きサイトの周辺原子
        nn_nnn_site = get_filled_sites(atom_set, bonds[empty])
        # 移動後に隣接原子があるか、Ag直上のサイト→候補
        if len(nn_nnn_site) >= 1 or empty[2] == 0:
            events.append(empty)
            rates.append(rate(energy))
    #
    # BLの上り下り
    # BLの下層原子
    if atom_z % 2 == 0:
        # BLを上る判定
        for filled in nn_atom:
            # 最近接原子の周囲原子
            nn_nn_atom = get_filled_sites(atom_set, bonds[filled])
            # 隣接原子の上が空いていて、かつ隣接原子が別原子でも支えられている→候補
            if atom_set(bonds[filled][3]) == 0 and len(nn_nn_atom) >= 2:
                events.append(filled)
                rates.append(rate(energy))
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
                d_step_empty = get_empty_sites(atom_set, bonds[direct_below])
                # 直下原子の周辺空きサイト→候補
                for empty in d_step_site:
                    events.append(empty)
                    rates.append(rate(energy))
            # 次近接の下へも移動可能
            # 次近接空きサイト
            nnn_empty = get_empty_sites(atom_set, nnn_sites)
            # 次近接空きサイト下サイト
            nnn_lower_site = get_lowers_sites(nnn_empty)
            # 次近接空きサイト下空きサイト
            nnn_lower_empty = get_empty_sites(atom_set, nnn_lower_site)
            for cand in nnn_lower_empty:
                # 次近接空きサイト下空きサイトの周辺原子
                cand_nn_atom = get_filled_sites(atom_set, bonds[cand])
                # 移動後に隣接原子がある→候補
                if len(cand_nn_atom) >= 1:
                    events.append(cand)
                    rates.append(rate(energy))
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
            nnn_atoms = get_filled_sites(atom_set, nnn_sites)
            # 次近接の直上サイト
            nnn_above_site = get_aboves(nnn_atoms)
            # 次近接の直上の空きサイト
            nnn_above_empty = get_empty_sites(atom_set, nnn_above_site)
            # 次近接直上の原子
            nnn_above_atoms = get_filled_sites(atom_set, nnn_above_site)
            # 次近接原子直上の空きサイト→候補
            for above_empty in nnn_above_empty:
                events.append(above_empty)
                rates.append(rate(energy))
            # 次近接直上に原子がある場合
            for above_filled in nnn_above_atoms:
                # 次近接直上の結合サイト
                above_filled_nn = bonds[above_filled]
                # 対象原子直上の結合サイト
                dir_above_nn = bonds[direct_above]
                # 共通項のサイトの抜出
                common_site = list(set(above_filled) & set(dir_above_nn))
                # 共通サイトに原子がない場合→候補
                if atom_set[common_site[0]] == 0:
                    events.append(common_site[0])
                    rates.append(rate(energy))
        #
        # BLを下る判定
        if atom_z == 1:
            pass
        else:
            # 最近接の空きサイト
            nn_empty = get_empty_sites(atom_set, bonds[target])
            # 最近接空きサイトの下サイト
            nn_lower = get_lowers_sites(nn_empty)
            # 最近接空きサイトの下空きサイト
            nn_lower_enpty = get_empty_sites(atom_set, nn_lower)
            for cand in nn_lower_enpty:
                # 候補サイトの周辺原子
                cand_nn = get_filled_sites(atom_set, bonds[cand])
                # 候補サイト移動後何らかの結合がある→候補
                if len(cand_nn) >= 1:
                    events.append(cand)
                    rates.append(rate(energy))
    # イベントリストから、孤立原子を生じるイベントを抽出
    remove = judge_isolation(atom_set, bonds, target, nn_atom, events)
    # 削除
    event_f = []
    rates_f = []
    for eve, rat in zip(events, rates):
        if eve in remove:
            pass
        else:
            event_f.append(eve)
            rates_f.append(rat)
    return event_f, rates_f


def get_events(atom_set, bonds, target, params):
    event_list: List[tuple] = []
    rate_list: List[float] = []
    unit_length = param.n_cell_init
    # calculate total energy
    energy = get_energy(atom_set, bonds, target, params)
    # calculate possible events
    event_list, rate_list = possible_events(
        atom_set, bonds, target, params, energy, unit_length
    )

    return event_list, rate_list
