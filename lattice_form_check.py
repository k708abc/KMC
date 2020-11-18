from typing import List, Dict
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def check(input_params, lattice_v, bonds_v):  # 多分別モジュールの方がいい。
    unit_x: List[float] = [1, 0, 0]
    unit_y: List[float] = [0.5, 0.866, 0]
    unit_length = input_params.n_cell_init
    z_units = input_params.z_unit_init
    maxz = z_units * 6 - 1
    # drawing range
    min_x = -1
    max_x = (unit_x[0] + unit_y[0]) * unit_length + 1
    min_y = -1
    max_y = (unit_x[1] + unit_y[1]) * unit_length + 1
    x_list: List[float] = [[] for _ in range(maxz + 1)]
    y_list: List[float] = [[] for _ in range(maxz + 1)]
    bond_list: List[float] = [[] for _ in range(maxz + 1)]
    #
    for key, atom_pos in lattice_v.items():
        xpos = atom_pos[0] * unit_x[0] + atom_pos[1] * unit_y[0]
        ypos = atom_pos[0] * unit_x[1] + atom_pos[1] * unit_y[1]
        x_list[key[2]].append(xpos)
        y_list[key[2]].append(ypos)
        bond_with = bonds_v[key]

        for u in bond_with:
            if u[2] >= maxz:
                pass
            else:
                atom_index_bond = (u[0], u[1], u[2])
                atom_pos_bond = lattice_v[atom_index_bond]
                xpos_b = (
                    atom_pos_bond[0] * unit_x[0] + atom_pos_bond[1] * unit_y[0]
                )
                ypos_b = (
                    atom_pos_bond[0] * unit_x[1] + atom_pos_bond[1] * unit_y[1]
                )

                bond_list[key[2]].append([[xpos, ypos], [xpos_b, ypos_b]])


    # input BL number to be checked
    BL_num = 0
    #
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(xlim=(min_x, max_x), ylim=(min_y, max_y))
    # Draw repeatition unit
    l1 = mlines.Line2D([0, unit_x[0] * unit_length], [0, 0], c="black")
    l2 = mlines.Line2D([unit_x[0] * unit_length, max_x - 1], [0, max_y - 1], c="black")
    l3 = mlines.Line2D(
        [max_x - 1, unit_y[0] * unit_length], [max_y - 1, max_y - 1], c="black"
    )
    l4 = mlines.Line2D([unit_y[0] * unit_length, 0], [max_y - 1, 0], c="black")
    ax.add_line(l1)
    ax.add_line(l2)
    ax.add_line(l3)
    ax.add_line(l4)
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
    plt.show()