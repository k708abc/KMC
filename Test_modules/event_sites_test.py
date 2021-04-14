from typing import Dict, Tuple


def find_lower_sites(indexes):
    return [(i[0], i[1], i[2] - 1) for i in indexes]


def find_aboves(indexes):
    return [(i[0], i[1], i[2] + 1) for i in indexes]


def possible_events(
    bonds: Dict,
    target: Tuple[int, int, int],
    params,
):
    candidate_sites = []
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
    #
    #
    # 最近接の空きサイト
    # nn_empty = find_empty_sites(atom_set, bonds[target])
    # 次近接の空きサイト
    # nnn_empty = find_empty_sites(atom_set, nnn_sites)
    # 最近接の存在原子
    # nn_atom = find_filled_sites(atom_set, bonds[target])
    #
    #
    # BL内での移動
    # 最近接空きサイトへの移動
    candidate_sites += bonds[target]
    # 同高さの次近接への移動
    candidate_sites += nnn_sites

    #
    # BLの上り下り
    # BLの下層原子
    if atom_z % 2 == 0:
        # BLを上る判定
        for nn in bonds[target]:
            candidate_sites += [(nn[0], nn[1], nn[2] + 1)]

        """
        eve_rate += [
            ((filled[0], filled[1], filled[2] + 1), rearange_rate)
            for filled in nn_atom
            if (
                atom_set[(filled[0], filled[1], filled[2] + 1)] == 0
                and len(find_filled_sites(atom_set, bonds[filled])) >= 2
            )
        ]
        """
        # BLを下る判定
        # Ag直上原子は下れない
        if atom_z == 0:
            pass
        else:
            # 直下原子のインデックス
            direct_below = (atom_x, atom_y, atom_z - 1)
            candidate_sites += bonds[direct_below]

            """
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
            """
            # 次近接の下へも移動可能
            nnn_lower_site = find_lower_sites(nnn_sites)
            candidate_sites += nnn_lower_site

            """
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
            """
    # BLの上層原子
    # 次近接の真上か、最近接真上と次近接真上の共通サイト
    elif atom_z % 2 == 1:
        # BLを上る判定
        # 直上原子のインデックス
        direct_above = (atom_x, atom_y, atom_z + 1)
        #
        # 次近接原子直上の空きサイト→候補
        nnn_above_atoms = find_aboves(nnn_sites)
        candidate_sites += nnn_above_atoms
        # 共通サイトに原子がない場合→候補
        for above_filled in nnn_above_atoms:
            # 次近接直上の結合サイト
            above_filled_nn = bonds[above_filled]
            # 対象原子直上の結合サイト
            dir_above_nn = bonds[direct_above]
            # 共通項のサイトの抜出
            common_site = list(set(above_filled_nn) & set(dir_above_nn))
            candidate_sites += common_site

        """
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
        """

        #
        # BLを下る判定
        if atom_z == 1:
            pass
        else:
            candidate_sites += find_lower_sites(bonds[target])

            """
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
            """
    """
    # イベントリストから、孤立原子を生じるイベントを抽出
    eve_rate = judge_isolation(atom_set, bonds, target, nn_atom, eve_rate)
    #
    event_f: List = [eves[0] for eves in eve_rate]
    rates_f: List = [eves[1] for eves in eve_rate]
    return event_f, rates_f
    """
    return candidate_sites
