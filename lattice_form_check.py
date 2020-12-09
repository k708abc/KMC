from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from InputParameter import Params

unit_x: List[float] = [1, 0, 0]
unit_y: List[float] = [0.5, 0.866, 0]


def reconstuct_representation(lattice, maxz):
    x_list: List[float] = [[] for _ in range(maxz + 1)]
    y_list: List[float] = [[] for _ in range(maxz + 1)]
    for site, atom_pos in lattice.items():
        xpos = atom_pos[0] * unit_x[0] + atom_pos[1] * unit_y[0]
        ypos = atom_pos[0] * unit_x[1] + atom_pos[1] * unit_y[1]
        x_list[site[2]].append(xpos)
        y_list[site[2]].append(ypos)
    return x_list, y_list


def represent_bond_for_figure(lattice, bonds, maxz):
    bond_list: List[float] = [[] for _ in range(maxz + 1)]
    for site, atom_pos in lattice.items():
        for bond in bonds[site]:
            if bond[2] >= maxz:
                pass
            else:
                xpos = atom_pos[0] * unit_x[0] + atom_pos[1] * unit_y[0]
                ypos = atom_pos[0] * unit_x[1] + atom_pos[1] * unit_y[1]
                atom_pos_bond = lattice[bond]
                xpos_b = atom_pos_bond[0] * unit_x[0] + atom_pos_bond[1] * unit_y[0]
                ypos_b = atom_pos_bond[0] * unit_x[1] + atom_pos_bond[1] * unit_y[1]
                bond_list[site[2]].append([[xpos, ypos], [xpos_b, ypos_b]])
    return bond_list


def figure_formation(unit_length, x_list, y_list, bond_list, BL_num):
    # drawing range
    min_x = -1
    max_x = (unit_x[0] + unit_y[0]) * unit_length + 1
    min_y = -1
    max_y = (unit_x[1] + unit_y[1]) * unit_length + 1
    #
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(xlim=(min_x, max_x), ylim=(min_y, max_y))
    # Draw repeatition unit
    ax.add_line(mlines.Line2D([0, unit_x[0] * unit_length], [0, 0], c="black"))
    ax.add_line(
        mlines.Line2D([unit_x[0] * unit_length, max_x - 1], [0, max_y - 1], c="black")
    )
    ax.add_line(
        mlines.Line2D(
            [max_x - 1, unit_y[0] * unit_length], [max_y - 1, max_y - 1], c="black"
        )
    )
    ax.add_line(mlines.Line2D([unit_y[0] * unit_length, 0], [max_y - 1, 0], c="black"))

    # draw bonding
    for bond_i in bond_list[BL_num * 2]:
        line = mlines.Line2D([bond_i[0][0], bond_i[1][0]], [bond_i[0][1], bond_i[1][1]])
        ax.add_line(line)
    for bond_i in bond_list[BL_num * 2 + 1]:
        line = mlines.Line2D([bond_i[0][0], bond_i[1][0]], [bond_i[0][1], bond_i[1][1]])
        ax.add_line(line)
    # draw atoms
    ax.scatter(x_list[BL_num * 2], y_list[BL_num * 2], c="b", s=30)
    ax.scatter(x_list[BL_num * 2 + 1], y_list[BL_num * 2 + 1], c="r", s=30)
    return fig


def check_lattice(
    input_params: Params, lattice: Dict[Tuple, List], bonds: Dict[Tuple, List]
):
    unit_length = input_params.n_cell_init
    z_units = input_params.z_unit_init
    maxz = z_units * 6 - 1
    #
    x_list, y_list = reconstuct_representation(lattice, maxz)
    bond_list = represent_bond_for_figure(lattice, bonds, maxz)
    # input BL number to be checked
    BL_num = int(input("Input layer number to be checked: "))
    #
    fig = figure_formation(unit_length, x_list, y_list, bond_list, BL_num)
    plt.show()
