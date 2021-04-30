from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys

sys.path.append("../../KMC/")

from Modules.InputParameter import Params
from Modules.lattice_form import lattice_form
import os

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


def rec_image(fig, BL_num):
    file_name = "Lattice_check/" + str(BL_num + 1) + "BL.png"
    fig.savefig(file_name)


def lattice_poscar(lattice, input_params: Params):
    xp: List[float] = []
    yp: List[float] = []
    zp: List[float] = []
    atom_i = 0
    unit_length = input_params.cell_size_xy
    maxz = input_params.cell_size_z * 6 - 1

    for index, atom_state in lattice.items():
        xp.append(lattice[index][0] / unit_length)
        yp.append(lattice[index][1] / unit_length)
        zp.append(lattice[index][2] / maxz / 2.448)
        atom_i += 1

    file_data = open("Lattice_check/lattice_poscar.vasp", "w")
    file_data.write("lattice check" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Si" + "\t" + "O" + "\n")
    file_data.write(str(atom_i) + "\n")
    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()


def check_lattice(
    input_params: Params, lattice: Dict[Tuple, List], bonds: Dict[Tuple, List]
):
    unit_length = input_params.cell_size_xy
    z_units = input_params.cell_size_z
    maxz = z_units * 6 - 1
    #
    x_list, y_list = reconstuct_representation(lattice, maxz)
    bond_list = represent_bond_for_figure(lattice, bonds, maxz)
    #
    for BL_num in range(0, 6):
        fig = figure_formation(unit_length, x_list, y_list, bond_list, BL_num)
        rec_image(fig, BL_num)


def diffuse_candidates(target, lattice, candidate, unit_length, maxz):
    xp: List[float] = [lattice[target][0] / unit_length]
    yp: List[float] = [lattice[target][1] / unit_length]
    zp: List[float] = [lattice[target][2] / maxz / 2.448]
    atom_cand = 0
    atom_lattice = 0
    for index in candidate[target]:
        if index != target:
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_cand += 1
    for index in lattice:
        if (index not in candidate[target]) and (index != target):
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_lattice += 1

    file_name = "Possible_diffusion/pos_" + str(target) + ".vasp"
    file_data = open(file_name, "w")
    file_data.write("check possible diffusion" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Bi" + "\t" + "Si" + "\t" + "O" + "\n")

    file_data.write(str(1) + "\t" + str(atom_cand) + "\t" + str(atom_lattice) + "\n")

    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()


if __name__ == "__main__":
    init_values = Params("../kmc_input.yml")
    if os.path.exists("Lattice_check") is False:
        os.mkdir("Lattice_check")
    lattice_formed = lattice_form(init_values)
    check_lattice(init_values, lattice_formed[0], lattice_formed[1])
    lattice_poscar(
        lattice_formed[0],
        init_values,
    )
    #
    if os.path.exists("Possible_diffusion") is False:
        os.mkdir("Possible_diffusion")
    unit_length = init_values.cell_size_xy
    maxz = init_values.cell_size_z
    for key in lattice_formed[0]:
        if key[2] > maxz:
            pass
        else:
            diffuse_candidates(
                key, lattice_formed[0], lattice_formed[7], unit_length, maxz
            )
    print("Lattice formation checked")