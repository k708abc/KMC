from typing import List, Dict
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

lattice: Dict[str, list] = {}
atom_set: Dict[str, int] = {}
bonds: Dict[str, list] = {}
event: Dict[str, list] = {}
event_time: Dict[str, list] = {}
event_time_tot: Dict[str, float] = {}

def lattice_form (input_params):
    unit_length = int(input_params["n_cell_init"])
    z_units = int(input_params["z_unit_init"])
    maxz = z_units * 6 - 1
    zd1 = float(input_params["intra_distance"])
    zd2 = float(input_params["inter_distance"])
    #
    lattice_first: List[float] = [
        [[
            [i, k, 0],
            [i+0.333, k+0.333, zd1],
            [i+0.333, k+0.333, zd1+zd2],
            [i+0.666, k+0.666, 2*zd1+zd2],
            [i+0.666, k+0.666, 2*(zd1+zd2)],
            [i, k, 2*(zd1+zd2)+zd1],            
        ] for k in range(unit_length)
        ] for i in range(unit_length)
    ]
    #
    for i in range(unit_length):
        for j in range(unit_length):
            for k in range(z_units):
                for l, first in enumerate(lattice_first[i][j]):
                    atom_index = str(i) + str(j) + str(k*6 + l)
                    lattice[atom_index] = [
                        round(first[0], 5),
                        round(first[1], 5),
                        round(first[2] + k * 2.448, 5),
                    ]
                    atom_set[atom_index] = 0
                    event[atom_index] = []
                    event_time[atom_index] = []
                    event_time_tot[atom_index] = 0
    # Search for bonding atoms for all the atoms
    for i in range(unit_length):
        for j in range(unit_length):
            for k in range(maxz):
                bond_with = []
                z_judge = int(k%6)
                ul_m = unit_length - 1
                if z_judge == 0:
                    if i == 0 and j == 0:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == 0:
                        bond_with = [[ul_m, j, k + 1], [0, j - 1, k + 1], [0, j, k + 1]]
                    elif j == 0:
                        bond_with = [[i - 1, 0, k + 1], [i, ul_m, k + 1], [i, 0, k + 1]]
                    else:
                        bond_with = [
                            [i - 1, j, k + 1],
                            [i, j - 1, k + 1],
                            [i, j, k + 1],
                        ]
                    if k != 0:
                        bond_with.append([i, j, k - 1])
                    else:
                        pass
                elif z_judge == 1:
                    if i == ul_m and j == ul_m:
                        bond_with = [
                            [ul_m, 0, k - 1],
                            [0, ul_m, k - 1],
                            [ul_m, ul_m, k - 1],
                        ]
                    elif i == ul_m:
                        bond_with = [
                            [ul_m, j + 1, k - 1],
                            [0, j, k - 1],
                            [ul_m, j, k - 1],
                        ]
                    elif j == ul_m:
                        bond_with = [
                            [i + 1, ul_m, k - 1],
                            [i, 0, k - 1],
                            [i, ul_m, k - 1],
                        ]
                    else:
                        bond_with = [
                            [i + 1, j, k - 1],
                            [i, j + 1, k - 1],
                            [i, j, k - 1],
                        ]
                    bond_with.append([i, j, k + 1])
                elif z_judge == 2:
                    if i == 0 and j == 0:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == 0:
                        bond_with = [[ul_m, j, k + 1], [0, j - 1, k + 1], [0, j, k + 1]]
                    elif j == 0:
                        bond_with = [[i - 1, 0, k + 1], [i, ul_m, k + 1], [i, 0, k + 1]]
                    else:
                        bond_with = [
                            [i - 1, j, k + 1],
                            [i, j - 1, k + 1],
                            [i, j, k + 1],
                        ]
                    bond_with.append([i, j, k - 1])
                elif z_judge == 3:
                    if i == ul_m and j == ul_m:
                        bond_with = [
                            [ul_m, 0, k - 1],
                            [0, ul_m, k - 1],
                            [ul_m, ul_m, k - 1],
                        ]
                    elif i == ul_m:
                        bond_with = [
                            [ul_m, j + 1, k - 1],
                            [0, j, k - 1],
                            [ul_m, j, k - 1],
                        ]
                    elif j == ul_m:
                        bond_with = [
                            [i + 1, ul_m, k - 1],
                            [i, 0, k - 1],
                            [i, ul_m, k - 1],
                        ]
                    else:
                        bond_with = [
                            [i + 1, j, k - 1],
                            [i, j + 1, k - 1],
                            [i, j, k - 1],
                        ]
                    bond_with.append([i, j, k + 1])
                elif z_judge == 4:
                    if i == ul_m and j == ul_m:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == ul_m:
                        bond_with = [
                            [0, j, k + 1],
                            [0, j + 1, k + 1],
                            [ul_m, j + 1, k + 1],
                        ]
                    elif j == ul_m:
                        bond_with = [
                            [i + 1, 0, k + 1],
                            [i, 0, k + 1],
                            [i + 1, ul_m, k + 1],
                        ]
                    else:
                        bond_with = [
                            [i + 1, j, k + 1],
                            [i + 1, j + 1, k + 1],
                            [i, j + 1, k + 1],
                        ]
                    bond_with.append([i, j, k - 1])
                elif z_judge == 5:
                    if i == 0 and j == 0:
                        bond_with = [
                            [ul_m, ul_m, k - 1],
                            [0, ul_m, k - 1],
                            [ul_m, 0, k - 1],
                        ]
                    elif i == 0:
                        bond_with = [
                            [0, j - 1, k - 1],
                            [ul_m, j, k - 1],
                            [ul_m, j - 1, k - 1],
                        ]
                    elif j == 0:
                        bond_with = [
                            [i, ul_m, k - 1],
                            [i - 1, ul_m, k - 1],
                            [i - 1, 0, k - 1],
                        ]
                    else:
                        bond_with = [
                            [i - 1, j, k - 1],
                            [i, j - 1, k - 1],
                            [i - 1, j - 1, k - 1],
                        ]
                    bond_with.append([i, j, k + 1])
                #
                atom_index = str(i) + str(j) + str(k)
                bonds[atom_index] = bond_with
    return lattice, bonds, atom_set, event, event_time, event_time_tot

def lattice_visual (input_params, lattice_v, bonds_v):
    unit_x: List[float] = [1, 0, 0]
    unit_y: List[float] = [0.5, 0.866, 0]
    unit_length = int(input_params["n_cell_init"])
    z_units = int(input_params["z_unit_init"])
    maxz = z_units * 6 - 1
    # drawing range
    min_x = -1
    max_x = (unit_x[0] + unit_y[0]) * unit_length + 1
    min_y = -1
    max_y = (unit_x[1] + unit_y[1]) * unit_length + 1
    x_list: List[float] = [[] for _ in range(maxz)]
    y_list: List[float] = [[] for _ in range(maxz)]
    bond_list: List[float] = [[] for _ in range(maxz)]
    #
    for k in range(maxz):
        for i in range(unit_length):
            for j in range(unit_length):
                atom_index = str(i) + str(j) + str(k)
                atom_pos = lattice_v[atom_index]
                #
                xpos = atom_pos[0] * unit_x[0] + atom_pos[1] * unit_y[0]
                ypos = atom_pos[0] * unit_x[1] + atom_pos[1] * unit_y[1]
                #
                x_list[k].append(xpos)
                y_list[k].append(ypos)
                bond_with = bonds_v[atom_index]
                for u in bond_with:
                    if u[2] >= maxz:
                        pass
                    else:
                        atom_index_bond = str(u[0]) + str(u[1]) + str(u[2])
                        atom_pos_bond = lattice_v[atom_index_bond]
                        xpos_b = (
                            atom_pos_bond[0] * unit_x[0] + atom_pos_bond[1] * unit_y[0]
                        )
                        ypos_b = (
                            atom_pos_bond[0] * unit_x[1] + atom_pos_bond[1] * unit_y[1]
                        )

                        bond_list[k].append([[xpos, ypos], [xpos_b, ypos_b]])
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

"""
if __name__ == "__main__":
    init_values: Dict = dict(
        n_cell_init = 5, 
        z_unit_init = 5,
        temperature = 550,
        dep_rate = 0.4,
        dep_time = 5,
        post_anneal = 0,
        prefactor = "1E+13",
        AgSi = -1.4,
        Si12 = -1.3,
        Si23 = -1.4,
        Si34 = -1.4,
        Si45 = -1.2,
        Si56 = -1.4,
        Si_intra = -1.4,
        Si_inter = -1.4,
        Ag_top = -1.4,
        transformation = -0.3,
        record_name = "KMC_Si_rec",
        img_per = 10,
        comments = "No comments",
        intra_distance = 0.25,
        inter_distance = 0.22
    )
    lattice_formed = lattice_form(init_values)
    lattice_visual(init_values, lattice_formed[0], lattice_formed[1])
"""
