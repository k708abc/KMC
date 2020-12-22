from InputParameter import Params
from typing import List, Dict, Tuple


def recalculate(
    target: Tuple[int, int, int],
    bonds: Dict[Tuple[int, int, int], List[Tuple[int, int, int]]],
    atom_set: Dict[Tuple[int, int, int], int],
    params: Params,
) -> List[Tuple[int, int, int]]:
    recal_list: List[Tuple[int, int, int]] = [target]
    atom_x, atom_y, atom_z = target
    unit_length = params.n_cell_init
    # 最近接
    nn_sites = bonds[target]
    #
    recal_list += [nn for nn in nn_sites if atom_set[nn] != 0]
    """
    for nn in nn_sites:
        if atom_set[nn] != 0:
            recal_list.append(nn)
    """
    # 次近接
    nnn_sites = [
        ((atom_x - 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y - 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, atom_y, atom_z),
        (atom_x, (atom_y + 1) % unit_length, atom_z),
        ((atom_x + 1) % unit_length, (atom_y - 1) % unit_length, atom_z),
        ((atom_x - 1) % unit_length, (atom_y + 1) % unit_length, atom_z),
    ]
    #
    recal_list += [nnn for nnn in nnn_sites if atom_set[nnn] != 0]
    """
    for nnn in nnn_sites:
        if atom_set[nnn] != 0:
            recal_list.append(nnn)
    """
    # BL下層に対する移動
    if atom_z % 2 == 0:
        # BL上から下ってくる
        # 最近接周囲原子の上
        recal_list += [
            (nn[0], nn[1], nn[2] + 1)
            for nn in nn_sites
            if atom_set[(nn[0], nn[1], nn[2] + 1)] != 0
        ]
        """
        for nn in nn_sites:
            above_site = (nn[0], nn[1], nn[2] + 1)
            if atom_set[above_site] != 0:
                recal_list.append(above_site)
        """
        # BL下から上ってくる
        if atom_z == 0:
            pass
        else:
            recal_list += [
                nn_d_b
                for nn_d_b in bonds[(atom_x, atom_y, atom_z - 1)]
                if atom_set[nn_d_b] != 0
            ]
            """
            direct_below = (atom_x, atom_y, atom_z - 1)
            nn_direct_below = bonds[direct_below]
            for nn_d_b in nn_direct_below:
                if atom_set[nn_d_b] != 0:
                    recal_list.append(nn_d_b)
            """
            recal_list += [
                (nnn[0], nnn[1], nnn[2] - 1)
                for nnn in nnn_sites
                if atom_set[(nnn[0], nnn[1], nnn[2] - 1)] != 0
            ]
            """
            # 次近接下から登ってくる
            for nnn in nnn_sites:
                nnn_below = (nnn[0], nnn[1], nnn[2] - 1)
                if atom_set[nnn_below] != 0:
                    recal_list.append(nnn_below)
            """
    # BL上層に対する移動
    elif atom_z % 2 == 1:
        # 次近接真上から下ってくる
        recal_list += [
            (nnn[0], nnn[1], nnn[2] + 1)
            for nnn in nnn_sites
            if atom_set[(nnn[0], nnn[1], nnn[2] + 1)] != 0
        ]

        """
        for nnn in nnn_sites:
            nnn_above = (nnn[0], nnn[1], nnn[2] + 1)
            if atom_set[nnn_above] != 0:
                recal_list.append(nnn_above)
        """
        # 直上近接から下ってくる
        recal_list += [
            nn_d_a
            for nn_d_a in bonds[(atom_x, atom_y, atom_z + 1)]
            if atom_set[nn_d_a] != 0
        ]
        """
        direct_above = (atom_x, atom_y, atom_z + 1)
        nn_direct_above = bonds[direct_above]
        for nn_d_a in nn_direct_above:
            if atom_set[nn_d_a] != 0:
                recal_list.append(nn_d_a)
        """
        # 上ってくる
        if atom_z == 1:
            pass
        else:
            # 隣接直下から登る
            recal_list += [
                (nn[0], nn[1], nn[2] - 1)
                for nn in nn_sites
                if atom_set[(nn[0], nn[1], nn[2] - 1)] != 0
            ]

            """
            for nn in nn_sites:
                nn_below = (nn[0], nn[1], nn[2] - 1)
                if atom_set[nn_below] != 0:
                    recal_list.append(nn_below)
            """
    return recal_list
