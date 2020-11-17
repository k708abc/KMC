#!/usr/bin/env python3

import numpy as np
import copy
import math
import time
from typing import List, Dict
import tkinter as tk
from tkinter import ttk
import os
import random
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.lines as mlines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageGrab
from pptx import Presentation
from pptx.util import Inches, Pt
from preference_window import Window

input_params: Dict[str, float] = {}
lattice: Dict[str, list] = {}
atom_set: Dict[str, int] = {}
bonds: Dict[str, list] = {}
event: Dict[str, list] = {}
event_time: Dict[str, list] = {}
event_time_tot: Dict[str, float] = {}


def update_values():
    input_params["Unit_length"] = float(entry_unit_length.get())
    input_params["Number_of_z_units"] = float(entry_zunit.get())
    input_params["Temperature"] = float(entry_temperature.get())
    input_params["kbt"] = float(entry_temperature.get())*8.617e-5
    input_params["deposition_rate"] = float(entry_dep_rate.get())
    input_params["deposition_time"] = float(entry_dep_time.get())
    input_params["post_anneal"] = float(entry_post_anneal.get())
    input_params["prefactor"] = float(entry_prefactor.get())
    # calculate deposition speed in atoms/sec
    unit_length = int(entry_unit_length.get())
    deposition_rate = float(entry_dep_rate.get())  # ML/min
    deposition_rate = deposition_rate/60  # ML/s
    atoms_in_BL = 2 *unit_length**2  # atoms/ML
    deposition_rate = deposition_rate*atoms_in_BL  # atoms/s
    input_params["Atoms_in_BL"] = float(atoms_in_BL)
    input_params["dep_rate_(atom/s)"] = float(deposition_rate)
    text_dep_atom_persec["text"] = str("{:.5f}".format(deposition_rate))
    #update rates
    for i in range(len(rates_list)):
        rates_list[i]["text"] = str("{:.5f}".format(cal_rate(float(bonding_entry_list[i].get()))))
    #
    root.update()

def update(event):
    update_values()



def show_current():
    global atom_set, lattice, c_num, max_layer
    print("current saving")
    # image
    fig = plt.figure()
    ax = fig.add_subplot(111)
    length = nl
    p = pat.Polygon(
        xy=[
            (0, 0),
            (unit_x[0] * length, unit_x[1] * length),
            ((unit_x[0] + unit_y[0]) * length, (unit_x[1] + unit_y[1]) * length),
            (unit_y[0] * length, unit_y[1] * length),
        ],
        fc=(0.5, 0.5, 0.5),
        ec=(0, 0, 0),
    )
    ax.add_patch(p)
    atom_set_n = atom_set
    for z in range(0, max_layer):
        for i in range(len(atom_set_n)):
            for k in range(len(atom_set_n[i])):
                if atom_set_n[i][k][z] == 0:
                    pass
                else:
                    xp = lattice[i][k][z][0]
                    yp = lattice[i][k][z][1]
                    zp = lattice[i][k][z][2]

                    if z % 2 == 0:
                        t1x = (
                            xp * unit_x[0]
                            + yp * unit_y[0]
                            - (unit_x[0] + unit_y[0]) / 3
                        )
                        t1y = (
                            xp * unit_x[1]
                            + yp * unit_y[1]
                            - (unit_x[1] + unit_y[1]) / 3
                        )
                        t2x = t1x + unit_x[0]
                        t2y = t1y + unit_x[1]
                        t3x = t1x + unit_y[0]
                        t3y = t1y + unit_y[1]
                    else:
                        t1x = (
                            xp * unit_x[0]
                            + yp * unit_y[0]
                            + (unit_x[0] + unit_y[0]) / 3
                        )
                        t1y = (
                            xp * unit_x[1]
                            + yp * unit_y[1]
                            + (unit_x[1] + unit_y[1]) / 3
                        )

                        t2x = t1x - unit_x[0]
                        t2y = t1y - unit_x[1]

                        t3x = t1x - unit_y[0]
                        t3y = t1y - unit_y[1]

                    color_num = math.floor(z / 2) * 2
                    if z == 1 or z == 0:
                        color = [0, 1, 0]
                    else:
                        color = [color_num / max_layer, 0, 1 - color_num / max_layer]

                    p = pat.Polygon(
                        xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                    )

                    ax.add_patch(p)

                    if i == 0:
                        if z % 2 == 0:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                - (unit_x[0] + unit_y[0]) / 3
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                - (unit_x[1] + unit_y[1]) / 3
                                + unit_x[1] * length
                            )

                            t2x = t1x + unit_x[0]
                            t2y = t1y + unit_x[1]

                            t3x = t1x + unit_y[0]
                            t3y = t1y + unit_y[1]
                        else:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                + (unit_x[0] + unit_y[0]) / 3
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                + (unit_x[1] + unit_y[1]) / 3
                                + unit_x[1] * length
                            )

                            t2x = t1x - unit_x[0]
                            t2y = t1y - unit_x[1]

                            t3x = t1x - unit_y[0]
                            t3y = t1y - unit_y[1]

                        p = pat.Polygon(
                            xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                        )
                        ax.add_patch(p)

                    if k == 0:
                        if z % 2 == 0:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                - (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                - (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                            )

                            t2x = t1x + unit_x[0]
                            t2y = t1y + unit_x[1]

                            t3x = t1x + unit_y[0]
                            t3y = t1y + unit_y[1]
                        else:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                + (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                + (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                            )

                            t2x = t1x - unit_x[0]
                            t2y = t1y - unit_x[1]

                            t3x = t1x - unit_y[0]
                            t3y = t1y - unit_y[1]

                        p = pat.Polygon(
                            xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                        )
                        ax.add_patch(p)

                    if (i == 0) and (k == 0):

                        if z % 2 == 0:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                - (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                - (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                                + unit_x[1] * length
                            )

                            t2x = t1x + unit_x[0]
                            t2y = t1y + unit_x[1]

                            t3x = t1x + unit_y[0]
                            t3y = t1y + unit_y[1]
                        else:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                + (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                + (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                                + unit_x[1] * length
                            )

                            t2x = t1x - unit_x[0]
                            t2y = t1y - unit_x[1]

                            t3x = t1x - unit_y[0]
                            t3y = t1y - unit_y[1]

                        p = pat.Polygon(
                            xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                        )
                        ax.add_patch(p)
    p = pat.Polygon(
        xy=[(0, 0), (unit_y[0] * length, unit_y[1] * length), (0, unit_y[1] * length)],
        fc=(1, 1, 1),
        ec=(0, 0, 0),
    )
    ax.add_patch(p)
    p = pat.Polygon(
        xy=[
            (unit_x[0] * length, 0),
            ((unit_x[0] + unit_y[0]) * length, 0),
            ((unit_x[0] + unit_y[0]) * length, (unit_x[1] + unit_y[1]) * length),
        ],
        fc=(1, 1, 1),
        ec=(0, 0, 0),
    )
    ax.add_patch(p)
    ax.set_xlim(0, (unit_x[0] + unit_y[0]) * length)
    ax.set_ylim(0, (unit_x[1] + unit_y[1]) * length)
    ax.set_aspect("equal")
    file_name = entry_rec.get() + "_current_" + str(c_num) + ".png"
    fig.savefig(file_name)

    # poscar

    xp = []
    yp = []
    zp = []
    num = 0
    s = atom_set

    for i in range(len(s)):
        for k in range(len(s[i])):
            for z in range(len(s[i][k])):
                if s[i][k][z] != 0:
                    xp.append(lattice[i][k][z][0] / nl)
                    yp.append(lattice[i][k][z][1] / nl)
                    zp.append(lattice[i][k][z][2] / zl / 2.448)
                    num = num + 1

    file_name = entry_rec.get() + "_" + "current" + str(c_num) + ".vasp"
    file_data = open(file_name, "w")
    file_data.write(entry_rec.get() + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(nl) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(str(nl / 2) + "\t" + str(nl / 2 * 1.732) + "\t" + "0" + "\n")
    file_data.write("0" + "\t" + "0" + "\t" + str(zl * 2.448) + "\n")
    file_data.write("Si" + "\n")
    file_data.write(str(num) + "\n")
    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()
    c_num = c_num + 1


def lattice_form():
    #global lattice, atom_set, bonds, event, event_time, event_tot, maxz, max_atom
    unit_length = int(input_params["Unit_length"])
    z_units = int(input_params["Number_of_z_units"])
    zd1 = input_params["intra_distance"]
    zd2 = input_params["inter_distance"]
    lattice_first = []
    # Calculating the position of atoms in the first unit in z (first 3 BL)
    for i in range(0, unit_length):
        lattice_first.append([])
        for k in range(0, unit_length):
            lattice_first[-1].append([])
            unit_pos = [i, k, 0]  # coodination of the unit cell
            # 6 atoms in a unit in z direction
            # first BL
            atom1 = [unit_pos[0], unit_pos[1], unit_pos[2]]
            atom2 = [atom1[0] + 0.333, atom1[1] + 0.333, atom1[2] + zd1]
            # second BL
            atom3 = [atom2[0], atom2[1], atom2[2] + zd2]
            atom4 = [atom3[0] + 0.334, atom3[1] + 0.334, atom3[2] + zd1]
            # Third BL
            atom5 = [atom4[0], atom4[1], atom4[2] + zd2]
            atom6 = [atom1[0], atom1[1], atom5[2] + zd1]
            # latice_first: coodination of all atoms
            lattice_first[-1][-1].extend([atom1, atom2, atom3, atom4, atom5, atom6])
    # Calculating the position of all atoms
    for i in range(0, unit_length):
        for j in range(0, unit_length):
            for k in range(0, z_units):
                for l in range(len(lattice_first[i][k])):
                    x_index = i
                    y_index = j
                    z_index = k*len(lattice_first[i][k]) + l
                    atom_index = str(x_index) + str(y_index) + str(z_index)
                    lattice[atom_index] = [
                            round(lattice_first[i][j][l][0], 5),
                            round(lattice_first[i][j][l][1], 5),
                            round(lattice_first[i][j][l][2] + k * 2.448, 5),
                        ]
                    atom_set[atom_index] = 0
                    event[atom_index] = []
                    event_time[atom_index] = []
                    event_time_tot[atom_index] = 0
    maxz = z_units*6 - 1
    input_params["maxz"] = maxz
    # Search for bonding atoms for all the atoms
    for i in range(0, unit_length):
        for j in range(0, unit_length):
            for k in range(0, maxz):
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
                        bond_with = [[i - 1, j, k + 1], [i, j - 1, k + 1], [i, j, k + 1]]
                    if k != 0:
                        bond_with.append([i, j, k - 1])
                    else:
                        pass
                elif z_judge == 1:
                    if i == ul_m and j == ul_m:
                        bond_with = [[ul_m, 0, k - 1], [0, ul_m, k - 1], [ul_m, ul_m, k - 1]]
                    elif i == ul_m:
                        bond_with = [[ul_m, j + 1, k - 1], [0, j, k - 1], [ul_m, j, k - 1]]
                    elif j == ul_m:
                        bond_with = [[i + 1, ul_m, k - 1], [i, 0, k - 1], [i, ul_m, k - 1]]
                    else:
                        bond_with = [[i + 1, j, k - 1], [i, j + 1, k - 1], [i, j, k - 1]]
                    bond_with.append([i, j, k + 1])
                elif z_judge == 2:
                    if i == 0 and j == 0:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == 0:
                        bond_with = [[ul_m, j, k + 1], [0, j - 1, k + 1], [0, j, k + 1]]
                    elif j == 0:
                        bond_with = [[i - 1, 0, k + 1], [i, ul_m, k + 1], [i, 0, k + 1]]
                    else:
                        bond_with = [[i - 1, j, k + 1], [i, j - 1, k + 1], [i, j, k + 1]]
                    bond_with.append([i, j, k - 1])
                elif z_judge == 3:
                    if i == ul_m and j == ul_m:
                        bond_with = [[ul_m, 0, k - 1], [0, ul_m, k - 1], [ul_m, ul_m, k - 1]]
                    elif i == ul_m:
                        bond_with = [[ul_m, j + 1, k - 1], [0, j, k - 1], [ul_m, j, k - 1]]
                    elif j == ul_m:
                        bond_with = [[i + 1, ul_m, k - 1], [i, 0, k - 1], [i, ul_m, k - 1]]
                    else:
                        bond_with = [[i + 1, j, k - 1], [i, j + 1, k - 1], [i, j, k - 1]]
                    bond_with.append([i, j, k + 1])
                elif z_judge == 4:
                    if i == ul_m and j == ul_m:
                        bond_with = [[ul_m, 0, k + 1], [0, ul_m, k + 1], [0, 0, k + 1]]
                    elif i == ul_m:
                        bond_with = [[0, j, k + 1], [0, j + 1, k + 1], [ul_m, j + 1, k + 1]]
                    elif j == ul_m:
                        bond_with = [[i + 1, 0, k + 1], [i, 0, k + 1], [i + 1, ul_m, k + 1]]
                    else:
                        bond_with = [
                            [i + 1, j, k + 1],
                            [i + 1, j + 1, k + 1],
                            [i, j + 1, k + 1],
                        ]
                    bond_with.append([i, j, k - 1])
                elif z_judge == 5:
                    if i == 0 and j == 0:
                        bond_with = [[ul_m, ul_m, k - 1], [0, ul_m, k - 1], [ul_m, 0, k - 1]]
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

def lattice_form_check():
    #For tsting lattice_form
    lattice_form()
    unit_length = int(input_params["Unit_length"])
    maxz = int(input_params["maxz"])
    #drawing range
    min_x = -1
    max_x = (unit_x[0] + unit_y[0]) * unit_length + 1
    min_y = -1
    max_y = (unit_x[1] + unit_y[1]) * unit_length + 1
    x_list = []
    y_list = []
    bond_list = []
    #
    for k in range(maxz):
        x_list.append([])
        y_list.append([])
        bond_list.append([])
        for i in range(unit_length):
            for j in range(unit_length):
                atom_index = str(i) + str(j) + str(k)
                atom_pos = lattice[atom_index]
                #
                xpos = atom_pos[0] * unit_x[0] + atom_pos[1] * unit_y[0]
                ypos = atom_pos[0] * unit_x[1] + atom_pos[1] * unit_y[1]
                #
                x_list[-1].append(xpos)
                y_list[-1].append(ypos)
                bond_with = bonds[atom_index]
                for u in bond_with:
                    if u[2] >= maxz:
                        pass
                    else:
                        atom_index_bond = str(u[0]) + str(u[1]) + str(u[2])
                        atom_pos_bond = lattice[atom_index_bond]
                        xpos_b = (
                            atom_pos_bond[0] * unit_x[0]
                            + atom_pos_bond[1] * unit_y[0]
                        )
                        ypos_b = (
                            atom_pos_bond[0] * unit_x[1]
                            + atom_pos_bond[1] * unit_y[1]
                        )

                        bond_list[-1].append([[xpos, ypos], [xpos_b, ypos_b]])
    # input BL number to be checked
    BL_num = 3
    #
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(xlim=(min_x, max_x), ylim=(min_y, max_y))
    # Draw repeatition unit
    l1 = mlines.Line2D([0, unit_x[0] * unit_length], [0, 0], c="black")
    l2 = mlines.Line2D([unit_x[0] * unit_length, max_x - 1], [0, max_y - 1], c="black")
    l3 = mlines.Line2D([max_x - 1, unit_y[0] * unit_length], [max_y - 1, max_y - 1], c="black")
    l4 = mlines.Line2D([unit_y[0] * unit_length, 0], [max_y - 1, 0], c="black")
    ax.add_line(l1)
    ax.add_line(l2)
    ax.add_line(l3)
    ax.add_line(l4)
    # draw bonding
    for d in bond_list[BL_num * 2]:
        l = mlines.Line2D([d[0][0], d[1][0]], [d[0][1], d[1][1]])
        ax.add_line(l)
    for d in bond_list[BL_num * 2 + 1]:
        l = mlines.Line2D([d[0][0], d[1][0]], [d[0][1], d[1][1]])
        ax.add_line(l)
    # draw atoms
    ax.scatter(x_list[BL_num * 2], y_list[BL_num * 2], c="b", s=30)
    ax.scatter(x_list[BL_num * 2 + 1], y_list[BL_num * 2 + 1], c="r", s=30)
    plt.show()


def deposition():
    # choose a position to deposit an atom
    global mod, atom_set, n_atoms, cov, max_layer

    # check nuber of atoms

    num_atom = 0
    for i in range(len(atom_set)):
        for k in range(len(atom_set[i])):
            for z in range(len(atom_set[i][k])):
                if atom_set[i][k][z] != 0:
                    num_atom += 1
    if num_atom != n_atoms:
        print("Mass error")
        print("deposited atom = " + str(n_atoms))
        print("caculated atom = " + str(num_atom))
    candidate = []
    mod = []
    n_atoms = n_atoms + 1
    cov = n_atoms / NBL
    if n_atoms in (5, 20, 300):
        show_current()

    if n_atoms >= max_atom:
        print("max_reached")

    # first bilayer (directlly attached to substrate) is always the candidate
    for i in range(len(lattice)):
        for k in range(len(lattice[i])):
            if atom_set[i][k][0] == 0:
                candidate.append([i, k, 0])
            if atom_set[i][k][1] == 0:
                candidate.append([i, k, 1])

    # the neighbor sites of existing atoms are candidates
    for i in range(len(atom_set)):
        for k in range(len(atom_set[i])):
            for z in range(len(atom_set[i][k])):
                if atom_set[i][k][z] == 0:
                    pass
                else:
                    for u in bonds[i][k][z]:
                        if atom_set[u[0]][u[1]][u[2]] == 0:
                            candidate.append([u[0], u[1], u[2]])
                        else:
                            pass

    # Delete duplication
    candidate2 = []
    for i in candidate:
        if (i in candidate2) == False:
            candidate2.append(i)

    # Choose a site
    r = random.random()
    num = math.floor(r * len(candidate2))
    dep = candidate2[num]

    # mod: Sites that should change the probability of events
    mod.append(dep)

    # Check the structure of surrounding atoms
    st = 2
    if dep[2] == maxz:
        print("atom at maxz")
    if dep[2] > max_layer:
        max_layer = dep[2]

    for u in bonds[dep[0]][dep[1]][dep[2]]:
        # if at leat one atom around the deposition site is 3D, deposit atom is also 3D.
        if atom_set[u[0]][u[1]][u[2]] == 3:
            st = 3
            mod.append(u)
        elif atom_set[u[0]][u[1]][u[2]] == 2:
            mod.append(u)

    atom_set[dep[0]][dep[1]][dep[2]] = st


def params():
    global pre, Eagsi, E12, E23, E34, E45, E56, Eintra, Einter, Eagtp, Etr, Vene, rec, rec1
    global d_rate, NBL
    pre = float(entry_pre.get())
    Eagsi = float(entry_AgSi.get())
    E12 = float(entry_Si12.get())
    E23 = float(entry_Si23.get())
    E34 = float(entry_Si34.get())
    E45 = float(entry_Si45.get())
    E56 = float(entry_Si56.get())
    Eintra = float(entry_Siel.get())
    Einter = float(entry_Sielin.get())
    Eagtp = float(entry_Agtp.get())
    Etr = float(entry_tr.get())

    Vene = [E12, E23, E34, E45, E56, Eintra, Einter]
    rec = float(entry_img.get())
    rec1 = float(entry_img.get())

    # calculate deposition speed
    d_rate = float(entry_rate.get())  # ML/min
    d_rate = d_rate / 60  # ML/s
    NBL = 2 * nl * nl  # atoms/ML
    d_rate = d_rate * NBL  # atoms/s

    global ttime, dep_time
    dep_time = float(entry_time.get()) * 60
    post_time = float(entry_post.get()) * 60
    ttime = dep_time + post_time


def update_progress():
    text_time_c["text"] = str("{:.1f}".format(rel_time)) + " s"
    pbval = int(rel_time / ttime * 100)
    pb.configure(value=pbval)
    text_count["text"] = str(pbval) + " %"
    text_event["text"] = str(event_num) + " events"
    text_coverage["text"] = str("{:.2f}".format(cov)) + " ML"
    text_atoms["text"] = str("{:.2f}".format(n_atoms)) + " atoms"
    root.update()


def time_progress():
    global rel_time, tot, pbval, rec, rec1
    r2 = 0
    while r2 == 0:
        r2 = random.random()
    ratio = int(rel_time / ttime * 100)
    dt = -math.log(r2) / tot
    rel_time = rel_time + dt
    # recording
    if ratio >= rec:
        rec_atoms()
        rec = rec + rec1


def rec_atoms():
    global a_set_rec, cov_rec, t_rec, cov
    a_set_rec.append([])

    a_set_rec[-1] = copy.deepcopy(atom_set)

    cov_rec.append(cov)
    t_rec.append(rel_time)


def cal_rate(E):
    pre = input_params["prefactor"]
    kbt = input_params["kbt"]
    rate = pre * math.exp(E/kbt)
    return rate


def update_events():
    global mod, event, event_time, event_tot, tot
    # Modify the event at the site in mod
    for e in mod:
        i = e[0]
        k = e[1]
        z = e[2]
        tot = tot - event_tot[i][k][z]
        event[i][k][z] = []
        event_time[i][k][z] = []
        event_tot[i][k][z] = 0
        tot_rec = 0
        # 最初はシンプルに組んで後で考え直す
        if atom_set[i][k][z] == 0:
            pass
        elif z == 0:
            # calculate total energy
            E = Eagsi
            for u in bonds[i][k][z]:
                s = atom_set[u[0]][u[1]][u[2]]
                if s != 0:
                    E += E12
            # get rate constant
            kr = cal_rate(E)
            # Record possible events and rate
            for u in bonds[i][k][z]:
                if atom_set[u[0]][u[1]][u[2]] == 0:
                    # empty neighbor site
                    event[i][k][z].append(["m", u[0], u[1], u[2]])
                    event_time[i][k][z].append(kr)
                    tot_rec += kr
                else:
                    # filled neighbor site
                    # Jucge upclimb event
                    sx = u[0]
                    sy = u[1]
                    sz = u[2] + 1
                    if atom_set[sx][sy][sz] == 0:
                        # upclimb possible
                        event[i][k][z].append(["u", sx, sy, sz])
                        event_time[i][k][z].append(kr)
                        tot_rec += kr
                    else:
                        # upclimb impossible
                        pass
        elif z == 1:
            # calculate total energy
            E = Eagsi
            for u in bonds[i][k][z]:
                s = atom_set[u[0]][u[1]][u[2]]
                if (s != 0) and (u[2] == 1):
                    # bonding with 0th atoms
                    E += E12
                elif (s != 0) and (u[2] == 3):
                    # bonding with next bilayer
                    E += E23
            # get rate constant
            kr = cal_rate(E)
            # Record possible events and rate
            for u in bonds[i][k][z]:
                if (atom_set[u[0]][u[1]][u[2]] == 0) and (u[2] == 0):
                    # empty neighbor site
                    event[i][k][z].append(["m", u[0], u[1], u[2]])
                    event_time[i][k][z].append(kr)
                    tot_rec += kr
                else:
                    pass
        elif (z == 2) or (z == 4):
            # calculate total energy
            E = 0
            for u in bonds[i][k][z]:
                s = atom_set[u[0]][u[1]][u[2]]
                if (s != 0) and (u[2] == z - 1):
                    # bonding with under BL
                    E += Vene[u[2]]
                elif (s != 0) and (u[2] == z + 1):
                    # bonding in same BL
                    E += Vene[u[2] + 1]
            # get rate constant
            kr = cal_rate(E)
            # Record possible events and rates
            for u in bonds[i][k][z]:
                if (atom_set[u[0]][u[1]][u[2]] == 0) and (u[2] == z + 1):
                    # empty neighbor site in same BL
                    event[i][k][z].append(["m", u[0], u[1], u[2]])
                    event_time[i][k][z].append(kr)
                    tot_rec += kr
                elif (atom_set[u[0]][u[1]][u[2]] == 0) and (u[2] == z - 1):
                    # empty neighbor site in lower BL
                    # Judge neighbor in lower site
                    judge = 0
                    for j in bonds[u[0]][u[1]][u[2]]:
                        if atom_set[j[0]][j[1]][j[2]] != 0:
                            judge += 1

                    if judge >= 2:
                        # the site has neighbor, which can make a bond. Candidate.
                        event[i][k][z].append(["m", u[0], u[1], u[2]])
                        event_time[i][k][z].append(kr)
                        tot_rec += kr
                elif (atom_set[u[0]][u[1]][u[2]] == 1) and (u[2] == z - 1):
                    # Downstair
                    for j in bonds[u[0]][u[1]][u[2]]:
                        if atom_set[j[0]][j[1]][j[2]] == 0:
                            # down stair possible
                            event[i][k][z].append(["d", j[0], j[1], j[2]])
                            event_time[i][k][z].append(kr)
                            tot_rec += kr
                        else:
                            pass
                elif (atom_set[u[0]][u[1]][u[2]] == 1) and (u[2] == z + 1):
                    # upstair
                    judge = 0
                    for j in bonds[u[0]][u[1]][[2]]:
                        if atom_set[j[0]][j[1]][j[2]] != 0:
                            judge += 1
                        else:
                            pass
                    if (judge >= 2) and (atom_set[u[0]][u[1]][u[2] + 1] == 0):
                        # upstair possible
                        event[i][k][z].append(["u", u[0], u[1], u[2] + 1])
                        event_time[i][k][z].append(kr)
                        tot_rec += kr
        elif z % 2 == 0:
            # calculate total energy
            E = 0
            for u in bonds[i][k][z]:
                s = atom_set[u[0]][u[1]][u[2]]
                if (s != 0) and (u[2] == z - 1):
                    # bonding with under BL
                    E += Vene[5]
                elif (s != 0) and (u[2] == z + 1):
                    # bonding in same BL
                    E += Vene[6]
            # get rate constant
            kr = cal_rate(E)
            # Record possible events and rates
            for u in bonds[i][k][z]:
                if (atom_set[u[0]][u[1]][u[2]] == 0) and (u[2] == z + 1):
                    # empty neighbor site in same BL
                    event[i][k][z].append(["m", u[0], u[1], u[2]])
                    event_time[i][k][z].append(kr)
                    tot_rec += kr
                elif (atom_set[u[0]][u[1]][u[2]] == 0) and (u[2] == z - 1):
                    # empty neighbor site in lower BL
                    # Judge neighbor in lower site
                    judge = 0
                    for j in bonds[u[0]][u[1]][u[2]]:
                        if atom_set[j[0]][j[1]][j[2]] != 0:
                            judge += 1

                    if judge >= 2:
                        # the site has neighbor, which can make a bond. Candidate.
                        event[i][k][z].append(["m", u[0], u[1], u[2]])
                        event_time[i][k][z].append(kr)
                        tot_rec += kr
                elif (atom_set[u[0]][u[1]][u[2]] == 1) and (u[2] == z - 1):
                    # Downstair
                    for j in bonds[u[0]][u[1]][u[2]]:
                        if atom_set[j[0]][j[1]][j[2]] == 0:
                            # down stair possible
                            event[i][k][z].append(["d", j[0], j[1], j[2]])
                            event_time[i][k][z].append(kr)
                            tot_rec += kr
                        else:
                            pass
                elif (atom_set[u[0]][u[1]][u[2]] == 1) and (u[2] == z + 1):
                    # upstair
                    judge = 0
                    for j in bonds[u[0]][u[1]][[2]]:
                        if atom_set[j[0]][j[1]][j[2]] != 0:
                            judge += 1
                        else:
                            pass

                    if (judge >= 2) and (atom_set[u[0]][u[1]][u[2] + 1] == 0):
                        # upstair possible
                        if u[2] + 1 >= maxz:
                            print("maxz reached")
                        event[i][k][z].append(["u", u[0], u[1], u[2] + 1])
                        event_time[i][k][z].append(kr)
                        tot_rec += kr

        elif (z == 3) or (z == 5):
            # calculate total energy
            E = 0
            for u in bonds[i][k][z]:
                s = atom_set[u[0]][u[1]][u[2]]
                if (s != 0) and (u[2] == z - 1):
                    # bonding in same BL
                    E += Vene[z]
                elif (s != 0) and (u[2] == z + 1):
                    # bonding with next bilayer
                    E += Vene[z + 1]
            # get rate constant
            kr = cal_rate(E)

            # Record possible events and rate
            for u in bonds[i][k][z]:
                if atom_set[u[0]][u[1]][u[2]] == 0:
                    # empty neighbor site
                    judge = 0
                    for j in bonds[u[0]][u[1]][u[2]]:
                        if atom_set[j[0]][j[1]][j[2]] != 0:
                            judge += 1

                    if judge >= 2:
                        # candidate
                        event[i][k][z].append(["m", u[0], u[1], u[2]])
                        event_time[i][k][z].append(kr)
                        tot_rec += kr

                else:
                    pass

        elif z % 2 == 1:
            # calculate total energy
            E = 0
            for u in bonds[i][k][z]:
                s = atom_set[u[0]][u[1]][u[2]]
                if (s != 0) and (u[2] == z - 1):
                    # bonding in same BL
                    E += Vene[6]
                elif (s != 0) and (u[2] == z + 1):
                    # bonding with next bilayer
                    E += Vene[5]
            # get rate constant
            kr = cal_rate(E)

            # Record possible events and rate
            for u in bonds[i][k][z]:
                if u[2] >= maxz:
                    print("maxz reached")
                elif atom_set[u[0]][u[1]][u[2]] == 0:
                    # empty neighbor site
                    judge = 0
                    for j in bonds[u[0]][u[1]][u[2]]:
                        if atom_set[j[0]][j[1]][j[2]] != 0:
                            judge += 1

                    if judge >= 2:
                        # candidate
                        event[i][k][z].append(["m", u[0], u[1], u[2]])
                        event_time[i][k][z].append(kr)
                        tot_rec += kr
                else:
                    pass

        else:
            print("No updates")
        event_tot[i][k][z] = tot_rec
        tot = tot + tot_rec

    mod = []


def event_progress(i, k, z, j):
    global atom_set, mod, max_layer

    if (
        (i >= len(event))
        or (k >= len(event[i]))
        or (z >= len(event[i][k]))
        or (j >= len(event[i][k][z]))
    ):
        print("Over")
        print("i = " + str(i))
        print("len(i) = " + str(len(event)))
        print("k = " + str(k))
        print("len(k) = " + str(len(event[i])))
        print("z = " + str(z))
        print("len(z) = " + str(len(event[i][k])))
        print("j = " + str(j))
        print("len(j) = " + str(len(event[i][k][z])))
        print("atomset = " + str(atom_set[i][k][z]))
        print("time = " + str(event_tot[i][k][z]))

    s = event[i][k][z][j]
    mod = []
    if atom_set[i][k][z] == 0:
        print("some mistake")

    if s[0] == "m":
        atom_set[i][k][z] = 0
        atom_set[s[1]][s[2]][s[3]] = 2
    elif s[0] == "u":
        atom_set[i][k][z] = 0
        atom_set[s[1]][s[2]][s[3]] = 2
    elif s[0] == "d":
        atom_set[i][k][z] = 0
        atom_set[s[1]][s[2]][s[3]] = 2
    elif s[0] == "t":
        pass
    if s[3] > max_layer:
        max_layer = s[3]

    mod.append([i, k, z])
    mod.append([s[1], s[2], s[3]])

    if bonds[i][k][z] == []:
        pass
    else:
        for a in bonds[i][k][z]:
            if atom_set[a[0]][a[1]][a[2]] != 0:
                mod.append([a[0], a[1], a[2]])

    if bonds[s[1]][s[2]][s[3]] == []:
        pass
    else:
        for a in bonds[s[1]][s[2]][s[3]]:
            if atom_set[a[0]][a[1]][a[2]] != 0:
                mod.append([a[0], a[1], a[2]])


def kMC():
    global tot, d_rate, rel_time, event_num
    event_num = 1
    while rel_time <= dep_time:
        update_progress()

        event_num = event_num + 1

        r1 = random.random()
        rt = r1 * tot
        # judge deposition
        if rt <= d_rate:
            # deposition happens
            deposition()
            # time progess
            time_progress()
            # update events and rates
            update_events()

        else:
            ie = 100
            rt = rt - d_rate
            for i in range(len(event_tot)):
                for k in range(len(event_tot[i])):
                    for z in range(len(event_tot[i][k])):
                        if rt <= 0:
                            pass
                        elif event_tot[i][k][z] == 0:
                            pass
                        elif rt <= event_tot[i][k][z]:
                            ie = i
                            ke = k
                            ze = z
                            trec = rt
                            rt = rt - event_tot[i][k][z]

                        else:
                            rt = rt - event_tot[i][k][z]

            if rt > 0:
                print("strange")
                print("tot = " + str(tot))
                print("remain = " + str(rt))
                print("ie =" + str(ie))
            else:
                for j in range(len(event_time[ie][ke][ze])):
                    # print(event_time[ie][ke][ze][j])
                    if trec <= 0:
                        pass
                    elif event_time[ie][ke][ze][j] == 0:
                        pass
                    elif trec <= event_time[ie][ke][ze][j]:
                        je = j
                        trec = trec - event_time[ie][ke][ze][j]
                    else:
                        trec = trec - event_time[ie][ke][ze][j]

                if trec > 0:
                    print("strange")
                    print("tot = " + str(tot))
                    print("remain = " + str(trec))

                else:

                    # event[ie][ke][ze][je] happens
                    event_progress(ie, ke, ze, je)
                    # time progess
                    time_progress()
                    # update events and rates
                    update_events()

    # deposition finish
    # post annealing
    tot = tot - d_rate

    # repetition untill post annealing finish
    while rel_time <= ttime:
        r1 = random.random()
        rt = r1 * tot

        for i in range(len(event_tot)):
            for k in range(len(event_tot[i])):
                for z in range(len(event_tot[i][k])):
                    if rt <= 0:
                        pass
                    elif event_tot[i][k][z] == 0:
                        pass
                    elif rt <= event_tot[i][k][z]:
                        ie = i
                        ke = k
                        ze = z
                        trec = rt
                        rt = rt - event_tot[i][k][z]
                    else:
                        rt = rt - event_tot[i][k][z]

        if rt > 0:
            print("strange")
            print("tot = " + str(tot))
            print("remain = " + str(rt))
            print("ie =" + str(ie))

        else:

            for j in range(len(event_time[ie][ke][ze])):
                if trec <= 0:
                    pass
                elif event_time[ie][ke][ze] == 0:
                    pass
                elif trec <= event_time[ie][ke][ze][j]:
                    je = j
                    trec = trec - event_time[ie][ke][ze][j]

            if trec > 0:
                pass
            else:
                # event[ie][ke][ze][je] happens
                event_progress(ie, ke, ze, je)
                # time progess
                time_progress()
                # update events and rates
                update_events()


def rec_pos():
    # record the position of atoms in poscar form
    for n in range(len(a_set_rec)):
        xp = []
        yp = []
        zp = []
        num = 0
        s = a_set_rec[n]

        for i in range(len(s)):
            for k in range(len(s[i])):
                for z in range(len(s[i][k])):
                    if s[i][k][z] != 0:
                        xp.append(lattice[i][k][z][0] / nl)
                        yp.append(lattice[i][k][z][1] / nl)
                        zp.append(lattice[i][k][z][2] / zl / 2.448)
                        num = num + 1

        file_name = (
            entry_rec.get() + "_" + str(int(t_rec[n])) + "s" + "_" + str(n) + ".vasp"
        )
        file_data = open(file_name, "w")
        file_data.write(entry_rec.get() + "\n")
        file_data.write("10" + "\n")

        file_data.write(str(nl) + "\t" + "0" + "\t" + "0" + "\n")
        file_data.write(str(nl / 2) + "\t" + str(nl / 2 * 1.732) + "\t" + "0" + "\n")
        file_data.write("0" + "\t" + "0" + "\t" + str(zl * 2.448) + "\n")
        file_data.write("Si" + "\n")
        file_data.write(str(num) + "\n")
        file_data.write("direct" + "\n")
        for i in range(len(xp)):
            file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
        file_data.close()


def figure_draw_rec():
    for n in range(len(a_set_rec)):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        kx = figure_formation(ax, n)
        file_name = (
            entry_rec.get()
            + "_middle_"
            + str(int(t_rec[n]))
            + "s"
            + "_"
            + str(n)
            + ".png"
        )
        fig.savefig(file_name)
        images.append(file_name)


def coverage_color():
    global images, gray, green
    number = 0
    distribution = []
    for i in images:
        distribution.append([0, 0])
        number = number + 1
        img = Image.open(i)
        img.convert("RGB")
        psize = img.size
        for s in range(psize[0]):
            for t in range(psize[1]):
                r, g, b, no = img.getpixel((s, t))
                col = [r, g, b]

                if col == [128, 128, 128]:
                    distribution[-1][0] += 1
                elif col == [0, 255, 0]:
                    distribution[-1][1] += 1
    gray = []
    green = []
    for i in distribution:
        gray.append(i[0])
        green.append(i[1])


def hist_rec():
    global images_h
    images_h = []
    for k in range(len(layers_each)):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        kx = hist_formation(ax, k)

        file_name = (
            entry_rec.get()
            + "_hist_"
            + str(int(t_rec[k]))
            + "s"
            + "_"
            + str(k)
            + ".png"
        )
        fig.savefig(file_name)
        images_h.append(file_name)


def ppt_form():
    global images, images_h
    k = os.path.exists("Layer_analysis_results.pptx")
    if k == True:
        prs = Presentation("Layer_analysis_results.pptx")
    else:
        prs = Presentation()
        prs.slide_height = Inches(7.5)
        prs.slide_width = Inches(13.33)
        title_slide_layout = prs.slide_layouts[0]
        slide = prs.slides.add_slide(title_slide_layout)
        title = slide.shapes.title
        title.text = "Layer analysis calculation results"
    # First page
    blank_slide_layout = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank_slide_layout)
    width = height = Inches(1)
    top = Inches(0)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    p = tf.add_paragraph()
    p.text = "Parameters"
    p.font.size = Pt(28)
    shapes = slide.shapes
    left = Inches(0.5)
    top = Inches(1.5)
    rows = 2
    cols = 11
    width = Inches(12)
    height = Inches(1)
    # record parameters
    table0 = shapes.add_table(rows, cols, left, top, width, height).table
    table0.cell(0, 0).text = "Num. cells"
    table0.cell(0, 1).text = "Z units"
    table0.cell(0, 2).text = "T (K)"
    table0.cell(0, 3).text = "kbT"
    table0.cell(0, 4).text = "dep. rate (ML/min)"
    table0.cell(0, 5).text = "dep. time (min)"
    table0.cell(0, 6).text = "post anneal(min)"
    table0.cell(0, 7).text = "prefactor"
    table0.cell(0, 8).text = "transform"
    table0.cell(0, 9).text = "keep defect"
    table0.cell(0, 10).text = "Cal. time (s)"
    table0.cell(1, 0).text = str(nl)
    table0.cell(1, 1).text = str(zl)
    table0.cell(1, 2).text = entry_kbT.get()
    table0.cell(1, 3).text = str("{:.3g}".format(kbt))
    table0.cell(1, 4).text = entry_rate.get()
    table0.cell(1, 5).text = entry_time.get()
    table0.cell(1, 6).text = entry_post.get()
    table0.cell(1, 7).text = entry_pre.get()
    if bln_tr.get() == True:
        table0.cell(1, 8).text = "on"
    else:
        table0.cell(1, 8).text = "off"
    if bln_def.get() == True:
        table0.cell(1, 9).text = "on"
    else:
        table0.cell(1, 9).text = "off"
    table0.cell(1, 10).text = str(minute) + " m " + str(second) + " s"
    left = Inches(0.5)
    top = Inches(3.5)
    rows = 2
    cols = 11
    width = Inches(12)
    height = Inches(1)
    table = shapes.add_table(rows, cols, left, top, width, height).table
    table.cell(0, 1).text = "Ag-Si"
    table.cell(0, 2).text = "Si(1-2)"
    table.cell(0, 3).text = "Si(2-3)"
    table.cell(0, 4).text = "Si(3-4)"
    table.cell(0, 5).text = "Si(4-5)"
    table.cell(0, 6).text = "Si(5-6)"
    table.cell(0, 7).text = "Si(intra)"
    table.cell(0, 8).text = "Si(inter)"
    table.cell(0, 9).text = "Ag(top)"
    table.cell(0, 10).text = "Trans."
    table.cell(1, 0).text = "Energy"
    table.cell(1, 1).text = entry_AgSi.get()
    table.cell(1, 2).text = entry_Si12.get()
    table.cell(1, 3).text = entry_Si23.get()
    table.cell(1, 4).text = entry_Si34.get()
    table.cell(1, 5).text = entry_Si45.get()
    table.cell(1, 6).text = entry_Si56.get()
    table.cell(1, 7).text = entry_Siel.get()
    table.cell(1, 8).text = entry_Sielin.get()
    table.cell(1, 9).text = entry_Agtp.get()
    table.cell(1, 10).text = entry_tr.get()
    """
    left = Inches(0.5)
    top = Inches(6.5)
    rows = 2
    cols = 3
    width = Inches(8)
    height = Inches(0.5)

    table = shapes.add_table(rows, cols, left, top, width, height).table
    table.cell(0, 0).text = "Total energy"
    table.cell(0, 1).text = "1 ML"
    table.cell(0, 2).text = "Final"

    table.cell(1, 0).text = str('{:.3g}'.format(E))
    table.cell(1, 1).text = str('{:.3g}'.format(ML_check))
    table.cell(1, 2).text = str('{:.3g}'.format(final_check))
    """

    width = height = Inches(1)
    top = Inches(6)
    left = Inches(0.5)

    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame

    p = tf.add_paragraph()
    p.text = entry_text.get()
    p.font.size = Pt(20)

    # second slide
    slide = prs.slides.add_slide(blank_slide_layout)

    width = height = Inches(1)
    top = Inches(-0.1)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame

    p = tf.add_paragraph()
    p.text = "Results"
    p.font.size = Pt(28)

    # put images

    left0 = 0.2
    top0 = 0.7
    height = Inches(2.5)
    height_w = Inches(1)

    num_ims = len(images)
    imn = 0
    slides = 0

    fin = 0

    while fin == 0:
        slides = slides + 1
        if slides != 1:
            slide = prs.slides.add_slide(blank_slide_layout)

            width = height = Inches(1)
            top = Inches(-0.1)
            left = Inches(0.5)
            txBox = slide.shapes.add_textbox(left, top, width, height)
            tf = txBox.text_frame

            p = tf.add_paragraph()
            p.text = "Results"
            p.font.size = Pt(28)

            left0 = 0.2
            top0 = 0.7
            height = Inches(2.5)
            height_w = Inches(1)

        for i in range(0, 3):
            for k in range(0, 4):
                if imn >= num_ims:
                    fin = 1
                else:
                    left = Inches(left0 + k * 3.2)
                    top = Inches(top0 + i * 2.1)
                    slide.shapes.add_picture(images[imn], left, top, height=height)

                    txBox = slide.shapes.add_textbox(
                        left, top - Inches(0.2), width, height_w
                    )
                    tf = txBox.text_frame

                    p = tf.add_paragraph()
                    p.text = (
                        str(int(t_rec[imn]))
                        + " s, "
                        + str("{:.2f}".format(cov_rec[imn]))
                        + " ML"
                    )
                    p.font.size = Pt(20)

                    imn = imn + 1
        if imn == num_ims:
            fin = 1

    # Third slide
    slide = prs.slides.add_slide(blank_slide_layout)

    width = height = Inches(1)
    top = Inches(-0.1)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame

    p = tf.add_paragraph()
    p.text = "Results: layer analysis"
    p.font.size = Pt(28)

    # put images

    left0 = 0.2
    top0 = 0.7
    height = Inches(2)
    height_w = Inches(1)

    num_ims = len(images_h)
    imn = 0
    slides = 0

    fin = 0

    while fin == 0:
        slides = slides + 1
        if slides != 1:
            slide = prs.slides.add_slide(blank_slide_layout)
            width = height = Inches(1)
            top = Inches(-0.1)
            left = Inches(0.5)
            txBox = slide.shapes.add_textbox(left, top, width, height)
            tf = txBox.text_frame
            p = tf.add_paragraph()
            p.text = "Results: layer analysis"
            p.font.size = Pt(28)
            left0 = 0.2
            top0 = 0.7
            height = Inches(2)
            height_w = Inches(1)
        for i in range(0, 3):
            for k in range(0, 4):
                if imn >= num_ims:
                    fin = 1
                else:
                    left = Inches(left0 + k * 3.2)
                    top = Inches(top0 + i * 2.1)
                    slide.shapes.add_picture(
                        images_h[imn], left, top + Inches(0.2), height=height
                    )

                    txBox = slide.shapes.add_textbox(
                        left, top - Inches(0.2), width, height_w
                    )
                    tf = txBox.text_frame

                    p = tf.add_paragraph()
                    p.text = (
                        str(int(t_rec[imn]))
                        + " s, "
                        + str("{:.2f}".format(cov_rec[imn]))
                        + " ML"
                    )
                    p.font.size = Pt(20)

                    imn = imn + 1
        if imn == num_ims:
            fin = 1
    """
    slide = prs.slides.add_slide(blank_slide_layout)

    width = height = Inches(1)
    top = Inches(-0.1)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame

    p = tf.add_paragraph()
    p.text = "Results: Coverage dependence"
    p.font.size = Pt(28)

    height = Inches(5.5)
    top = Inches(1)
    left = Inches(0.7)

    file_name = entry_rec.get() + "_coverage" + ".png"

    slide.shapes.add_picture(file_name, left, top, height=height)
    """
    prs.save("Layer_analysis_results.pptx")


def layer_sum():
    global layers_each, max_layer
    layers_each = []
    max_layer = 0
    for s in range(len(a_set_rec)):
        layers_each.append([])
        co = a_set_rec[s]
        for z in range(len(co[0][0])):
            c_each = 0
            if z % 2 == 0:
                count = 0
            for i in range(len(co)):
                for k in range(len(co[i])):
                    if co[i][k][z] != 0:
                        count += 1
                        c_each += 1
            if z % 2 == 1:
                layers_each[-1].append(count)
            if (c_each != 0) and (z >= max_layer):
                max_layer = z


def figure_formation(ax, n):
    length = nl
    p = pat.Polygon(
        xy=[
            (0, 0),
            (unit_x[0] * length, unit_x[1] * length),
            ((unit_x[0] + unit_y[0]) * length, (unit_x[1] + unit_y[1]) * length),
            (unit_y[0] * length, unit_y[1] * length),
        ],
        fc=(0.5, 0.5, 0.5),
        ec=(0, 0, 0),
    )
    ax.add_patch(p)
    atom_set_n = a_set_rec[n]
    for z in range(0, max_layer):
        for i in range(len(atom_set_n)):
            for k in range(len(atom_set_n[i])):
                if atom_set_n[i][k][z] == 0:
                    pass
                else:
                    xp = lattice[i][k][z][0]
                    yp = lattice[i][k][z][1]
                    zp = lattice[i][k][z][2]
                    if z % 2 == 0:
                        t1x = (
                            xp * unit_x[0]
                            + yp * unit_y[0]
                            - (unit_x[0] + unit_y[0]) / 3
                        )
                        t1y = (
                            xp * unit_x[1]
                            + yp * unit_y[1]
                            - (unit_x[1] + unit_y[1]) / 3
                        )

                        t2x = t1x + unit_x[0]
                        t2y = t1y + unit_x[1]

                        t3x = t1x + unit_y[0]
                        t3y = t1y + unit_y[1]

                    else:
                        t1x = (
                            xp * unit_x[0]
                            + yp * unit_y[0]
                            + (unit_x[0] + unit_y[0]) / 3
                        )
                        t1y = (
                            xp * unit_x[1]
                            + yp * unit_y[1]
                            + (unit_x[1] + unit_y[1]) / 3
                        )

                        t2x = t1x - unit_x[0]
                        t2y = t1y - unit_x[1]

                        t3x = t1x - unit_y[0]
                        t3y = t1y - unit_y[1]

                    color_num = math.floor(z / 2) * 2
                    if z == 1 or z == 0:
                        color = [0, 1, 0]
                    else:
                        color = [color_num / max_layer, 0, 1 - color_num / max_layer]

                    p = pat.Polygon(
                        xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                    )

                    ax.add_patch(p)

                    if i == 0:
                        if z % 2 == 0:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                - (unit_x[0] + unit_y[0]) / 3
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                - (unit_x[1] + unit_y[1]) / 3
                                + unit_x[1] * length
                            )

                            t2x = t1x + unit_x[0]
                            t2y = t1y + unit_x[1]

                            t3x = t1x + unit_y[0]
                            t3y = t1y + unit_y[1]
                        else:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                + (unit_x[0] + unit_y[0]) / 3
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                + (unit_x[1] + unit_y[1]) / 3
                                + unit_x[1] * length
                            )

                            t2x = t1x - unit_x[0]
                            t2y = t1y - unit_x[1]

                            t3x = t1x - unit_y[0]
                            t3y = t1y - unit_y[1]

                        p = pat.Polygon(
                            xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                        )
                        ax.add_patch(p)

                    if k == 0:
                        if z % 2 == 0:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                - (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                - (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                            )

                            t2x = t1x + unit_x[0]
                            t2y = t1y + unit_x[1]

                            t3x = t1x + unit_y[0]
                            t3y = t1y + unit_y[1]
                        else:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                + (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                + (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                            )

                            t2x = t1x - unit_x[0]
                            t2y = t1y - unit_x[1]

                            t3x = t1x - unit_y[0]
                            t3y = t1y - unit_y[1]

                        p = pat.Polygon(
                            xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                        )
                        ax.add_patch(p)

                    if (i == 0) and (k == 0):

                        if z % 2 == 0:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                - (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                - (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                                + unit_x[1] * length
                            )

                            t2x = t1x + unit_x[0]
                            t2y = t1y + unit_x[1]

                            t3x = t1x + unit_y[0]
                            t3y = t1y + unit_y[1]
                        else:
                            t1x = (
                                xp * unit_x[0]
                                + yp * unit_y[0]
                                + (unit_x[0] + unit_y[0]) / 3
                                + unit_y[0] * length
                                + unit_x[0] * length
                            )
                            t1y = (
                                xp * unit_x[1]
                                + yp * unit_y[1]
                                + (unit_x[1] + unit_y[1]) / 3
                                + unit_y[1] * length
                                + unit_x[1] * length
                            )

                            t2x = t1x - unit_x[0]
                            t2y = t1y - unit_x[1]

                            t3x = t1x - unit_y[0]
                            t3y = t1y - unit_y[1]

                        p = pat.Polygon(
                            xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                        )
                        ax.add_patch(p)

    p = pat.Polygon(
        xy=[(0, 0), (unit_y[0] * length, unit_y[1] * length), (0, unit_y[1] * length)],
        fc=(1, 1, 1),
        ec=(0, 0, 0),
    )
    ax.add_patch(p)

    p = pat.Polygon(
        xy=[
            (unit_x[0] * length, 0),
            ((unit_x[0] + unit_y[0]) * length, 0),
            ((unit_x[0] + unit_y[0]) * length, (unit_x[1] + unit_y[1]) * length),
        ],
        fc=(1, 1, 1),
        ec=(0, 0, 0),
    )
    ax.add_patch(p)

    ax.set_xlim(0, (unit_x[0] + unit_y[0]) * length)
    ax.set_ylim(0, (unit_x[1] + unit_y[1]) * length)
    ax.set_aspect("equal")

    return ax


def hist_formation(bx, n):
    numbers = layers_each[n]
    lay = []
    for i in range(len(numbers)):
        lay.append(i + 1)

    left = np.arange(len(numbers))

    bx.barh(left, numbers)
    bx.set_yticks(left)
    bx.set_yticklabels(lay)
    bx.set_xlabel("Coverage")
    bx.set_ylabel("layer number")

    bx.set_xlim(0, NBL)
    return bx


def figure_draw(n):
    ax = fig.add_subplot(211)
    ax.cla()
    ax = fig.add_subplot(211)
    kx = figure_formation(ax, n)
    # hist
    bx = fig.add_subplot(212)
    bx.cla()
    bx = fig.add_subplot(212)
    jx = hist_formation(bx, n)


def show_pictures():
    global canvas, imnum, fig
    root_p = tk.Toplevel()
    root_p.geometry("650x630")

    def p1_clicked():
        global imnum, cov_rec, canvas, fig
        imnum = imnum - 1
        if imnum < 0:
            imnum = len(cov_rec) - 1
        figure_draw(imnum)
        canvas.draw()
        canvas.get_tk_widget().pack()
        text_p8["text"] = str(int(t_rec[imnum])) + " s"
        text_p4["text"] = str(cov_rec[imnum])
        text_p6["text"] = str(imnum + 1) + "/" + str(len(cov_rec))

    def p2_clicked():
        global imnum, cov_rec, canvas, fig
        imnum = imnum + 1
        if imnum >= len(cov_rec):
            imnum = 0
        figure_draw(imnum)
        canvas.draw()
        canvas.get_tk_widget().pack()
        text_p8["text"] = str(int(t_rec[imnum])) + " s"
        text_p4["text"] = str(cov_rec[imnum])
        text_p6["text"] = str(imnum + 1) + "/" + str(len(cov_rec))

    button_p1 = tk.Button(
        root_p, text="Previous", command=p1_clicked, height=2, width=25
    )
    button_p1.place(x=30, y=560)
    button_p2 = tk.Button(
        root_p, text="Next", command=p2_clicked, height=2, width=25
    )
    button_p2.place(x=230, y=560)
    text_p7 = tk.Label(root_p, text="Time: ", font=("", 16))
    text_p7.place(x=30, y=500)
    text_p8 = tk.Label(root_p, text=str(cov_rec[0]), font=("", 16))
    text_p8.place(x=130, y=500)
    text_p3 = tk.Label(root_p, text="Coverage: ", font=("", 16))
    text_p3.place(x=230, y=500)
    text_p4 = tk.Label(root_p, text=str(cov_rec[0]), font=("", 16))
    text_p4.place(x=350, y=500)
    text_p6 = tk.Label(root_p, text="1/" + str(len(cov_rec)), font=("", 16))
    text_p6.place(x=450, y=500)
    imnum = len(cov_rec) - 1
    fig = plt.Figure()
    text_p8["text"] = str(int(t_rec[imnum])) + " s"
    text_p4["text"] = str(cov_rec[imnum])
    text_p6["text"] = str(imnum + 1) + "/" + str(len(cov_rec))
    figure_draw(imnum)
    canvas = FigureCanvasTkAgg(fig, master=root_p)
    canvas.draw()
    canvas.get_tk_widget().pack()

    def button_closep_clicked():
        plt.close("all")
        root_p.destroy()

    button_closep = tk.Button(
        root_p, text="Close", command=button_closep_clicked, height=2, width=25
    )
    button_closep.place(x=430, y=560)
    root_p.mainloop()


def cal_start():
    global lattice, atom_set, bonds, event, mod, rel_time, tot
    global a_set_rec, cov_rec, t_rec, images, d_rate
    start = time.time()
    pbval = 0
    pb.configure(value=pbval)
    text_count["text"] = "Started"
    text_coverage["text"] = "0 ML"
    root.update()
    # lists for recording
    a_set_rec = []
    cov_rec = []
    t_rec = []
    images = []
    # form lattice
    lattice_form()
    params()
    # tot: sum of rate constants. First, only deposition
    tot = d_rate
    # evaporation start
    rel_time = 0
    # first deposition
    global n_atoms
    n_atoms = 0
    deposition()
    text_count["text"] = "0 %"
    root.update()
    time_progress()
    # record first deposition
    rec_atoms()
    # update events and rates
    update_events()
    # second atomdepsoition
    deposition()
    time_progress()
    update_events()
    # repetition
    kMC()
    # record final structure
    text_count["text"] = "Saving"
    rec_atoms()
    # make poscar
    rec_pos()
    # record figures
    layer_sum()
    figure_draw_rec()
    coverage_color()
    hist_rec()
    global minute, second
    elapsed_time = time.time() - start
    minute = math.floor(elapsed_time / 60)
    second = int(elapsed_time % 60)
    text_count["text"] = str(minute) + " min " + str(second) + " sec"
    text_count["text"] = "PPT forming"
    ppt_form()
    text_count["text"] = "Finished"
    root.update()
    # calculation end
    # Show results
    show_pictures()


def button_start_clicked():
    lattice_form_check()
    #cal_start()

def button_close_clicked():
    plt.close("all")
    root.destroy()

class App(Window):
    def __init__(self, master):
        super().__init__(master)

if __name__ == "__main__":
    unit_x: List[float] = [1, 0, 0]
    unit_y: List[float] = [0.5, 0.866, 0]
    unit_z: List[float] = [0, 0, 1]

    application = tk.Tk()
    app = App(application)
    app.run()

    """
    application = tk.Tk()
    application.title("kMC_Si")
    Window(application)
    application.mainloop()
    """

    """
    input_params["intra_distance"] = 0.204
    input_params["inter_distance"] = 0.612
    input_params["max_layer"] = 0
    #max_layer = 0
    input_params["Saved_images"] = 0    #Use instead of "c_num"
    #Initial values
    input_params["Unit_length"] = 5
    input_params["Number_of_z_units"] = 5
    input_params["Temperature"]  = 550
    input_params["kbt"] = input_params["Temperature"] * 8.617e-5
    input_params["deposition_rate"] = 0.4
    input_params["deposition_time"] = 5
    input_params["post_anneal"] = 0
    input_params["prefactor"] = 1e+13
    #bonding parameters
    input_params["bond_Ag_Si"] = -1.2
    input_params["bond_Si1_Si2"] = -1.2
    input_params["bond_Si2_Si3"] = -1.2
    input_params["bond_Si3_Si4"] = -1.2
    input_params["bond_Si4_Si5"] = -1.2
    input_params["bond_Si5_Si6"] = -1.2
    input_params["bond_Si_intra"] = -1.2
    input_params["bond_Si_inter"] = -1.2
    input_params["bond_Ag_top"] = -1.2
    input_params["transformation_energy"] = -1.0
    #other inouts
    input_params["record_name"] = "KMC"
    input_params["record_image_every"] = 10
    input_params["comments"] = "No comment"
    #For adjusting the position of entries
    arx: List[int] = [0, 20, 140, 200, 300, 380, 500, 560, 660]
    ary: List[int] = [0, 20, 50, 100, 130, 160, 0, 190, 220, 250, 280, 310]
    #
    root = tk.Tk()
    root.title("kMC_Si")
    root.geometry("800x550")
    #
    text_unit_length = tk.Label(root, text="Number of cells")
    text_unit_length.place(x=arx[1], y=ary[1])
    entry_unit_length = tk.Entry(root, text="Number of cells", width=7)
    entry_unit_length.place(x=arx[2], y=ary[1])
    entry_unit_length.delete(0, tk.END)
    entry_unit_length.insert(tk.END, input_params["Unit_length"] )
    entry_unit_length.bind("<Return>", update)
    #
    text_zunit = tk.Label(root, text="Z units")
    text_zunit.place(x=arx[3], y=ary[1])
    entry_zunit = tk.Entry(root, text="Z units", width=7)
    entry_zunit.place(x=arx[4], y=ary[1])
    entry_zunit.delete(0, tk.END)
    entry_zunit.insert(tk.END, input_params["Number_of_z_units"])
    entry_zunit.bind("<Return>", update)
    #
    text_temperature = tk.Label(root, text="T (K)")
    text_temperature.place(x=arx[5], y=ary[1])
    entry_temperature = tk.Entry(root, text="T (K)", width=7)
    entry_temperature.place(x=arx[6], y=ary[1])
    entry_temperature.delete(0, tk.END)
    entry_temperature.insert(tk.END, input_params["Temperature"])
    entry_temperature.bind("<Return>", update)
    #
    text_kbt = tk.Label(root, text="kbT")
    text_kbt.place(x=arx[7], y=ary[1])
    text_kbt = tk.Label(root, text=str("{:.3g}".format(input_params["kbt"])))
    text_kbt.place(x=arx[8], y=ary[1])
    #
    text_dep_rate = tk.Label(root, text="dep_rate (ML/min)")
    text_dep_rate.place(x=arx[1], y=ary[2])
    entry_dep_rate = tk.Entry(root, text="dep_rate", width=7)
    entry_dep_rate.place(x=arx[2], y=ary[2])
    entry_dep_rate.delete(0, tk.END)
    entry_dep_rate.insert(tk.END, input_params["deposition_rate"])
    entry_dep_rate.bind("<Return>", update)
    #
    text_dep_atom_persec = tk.Label(root, text="0")
    text_dep_atom_persec.place(x=arx[2], y=ary[2] + 30)
    #
    text_dep_time = tk.Label(root, text="Dep.time (min)")
    text_dep_time.place(x=arx[3], y=ary[2])
    entry_dep_time = tk.Entry(root, text="Dep.time", width=7)
    entry_dep_time.place(x=arx[4], y=ary[2])
    entry_dep_time.delete(0, tk.END)
    entry_dep_time.insert(tk.END, input_params["deposition_time"])
    entry_dep_time.bind("<Return>", update)
    #
    text_post_anneal = tk.Label(root, text="Post anneal (min)")
    text_post_anneal.place(x=arx[5], y=ary[2])
    entry_post_anneal = tk.Entry(root, text="Post annealing", width=7)
    entry_post_anneal.place(x=arx[6], y=ary[2])
    entry_post_anneal.delete(0, tk.END)
    entry_post_anneal.insert(tk.END, input_params["post_anneal"] )
    entry_post_anneal.bind("<Return>", update)
    #
    text_prefactor = tk.Label(root, text="prefactor (1/s)")
    text_prefactor.place(x=arx[7], y=ary[2])
    entry_prefactor = tk.Entry(root, text="Prefactor (1/s)", width=7)
    entry_prefactor.place(x=arx[8], y=ary[2])
    entry_prefactor.delete(0, tk.END)
    entry_prefactor.insert(tk.END, '{:.1e}'.format(input_params["prefactor"]))
    entry_prefactor.bind("<Return>", update)    
    #set bonding parameters
    text_bonding = tk.Label(root, text="Energy")
    text_bonding.place(x=20, y=ary[4])
    # Ag-Si
    text_AgSi = tk.Label(root, text="Ag-Si")
    text_AgSi.place(x=100, y=ary[3])
    entry_AgSi = tk.Entry(root, text="Ag-Si", width=7)
    entry_AgSi.place(x=100, y=ary[4])
    entry_AgSi.delete(0, tk.END)
    entry_AgSi.insert(tk.END, input_params["bond_Ag_Si"] )
    entry_AgSi.bind("<Return>", update)
    # Si1-2
    text_Si12 = tk.Label(root, text="Si(1-2)")
    text_Si12.place(x=180, y=ary[3])
    entry_Si12 = tk.Entry(root, text="Si(1-2)", width=7)
    entry_Si12.place(x=180, y=ary[4])
    entry_Si12.delete(0, tk.END)
    entry_Si12.insert(tk.END, input_params["bond_Si1_Si2"])
    entry_Si12.bind("<Return>", update)
    # Si2-3
    text_Si23 = tk.Label(root, text="Si(2-3)")
    text_Si23.place(x=260, y=ary[3])
    entry_Si23 = tk.Entry(root, text="Si(2-3)", width=7)
    entry_Si23.place(x=260, y=ary[4])
    entry_Si23.delete(0, tk.END)
    entry_Si23.insert(tk.END, input_params["bond_Si2_Si3"])
    entry_Si23.bind("<Return>", update)
    # Si3-4
    text_Si34 = tk.Label(root, text="Si(3-4)")
    text_Si34.place(x=340, y=ary[3])
    entry_Si34 = tk.Entry(root, text="Si(3-4)", width=7)
    entry_Si34.place(x=340, y=ary[4])
    entry_Si34.delete(0, tk.END)
    entry_Si34.insert(tk.END, input_params["bond_Si3_Si4"])
    entry_Si34.bind("<Return>", update)
    # Si4-5
    text_Si45 = tk.Label(root, text="Si(4-5)")
    text_Si45.place(x=420, y=ary[3])
    entry_Si45 = tk.Entry(root, text="Si(4-5)", width=7)
    entry_Si45.place(x=420, y=ary[4])
    entry_Si45.delete(0, tk.END)
    entry_Si45.insert(tk.END, input_params["bond_Si4_Si5"])
    entry_Si45.bind("<Return>", update)
    # Si5-6
    text_Si56 = tk.Label(root, text="Si(5-6)")
    text_Si56.place(x=500, y=ary[3])
    entry_Si56 = tk.Entry(root, text="Si(5-6)", width=7)
    entry_Si56.place(x=500, y=ary[4])
    entry_Si56.delete(0, tk.END)
    entry_Si56.insert(tk.END, input_params["bond_Si5_Si6"])
    entry_Si56.bind("<Return>", update)
    # else intra layers
    text_Si_intra = tk.Label(root, text="Si(intra)")
    text_Si_intra.place(x=580, y=ary[3])
    entry_Si_intra = tk.Entry(root, text="Si(intra)", width=7)
    entry_Si_intra.place(x=580, y=ary[4])
    entry_Si_intra.delete(0, tk.END)
    entry_Si_intra.insert(tk.END, input_params["bond_Si_intra"])
    entry_Si_intra.bind("<Return>", update)
    # else inter layers
    text_Si_inter = tk.Label(root, text="Si(inter)")
    text_Si_inter.place(x=660, y=ary[3])
    entry_Si_inter = tk.Entry(root, text="Si(inter)", width=7)
    entry_Si_inter.place(x=660, y=ary[4])
    entry_Si_inter.delete(0, tk.END)
    entry_Si_inter.insert(tk.END, input_params["bond_Si_inter"])
    entry_Si_inter.bind("<Return>", update)
    # Agtop
    text_Ag_top = tk.Label(root, text="Ag(top)")
    text_Ag_top.place(x=740, y=ary[3])
    entry_Ag_top = tk.Entry(root, text="Ag(top)", width=7)
    entry_Ag_top.place(x=740, y=ary[4])
    entry_Ag_top.delete(0, tk.END)
    entry_Ag_top.insert(tk.END, input_params["bond_Ag_top"])
    entry_Ag_top.bind("<Return>", update)
    #bond entry lists
    bonding_entry_list = [
        entry_AgSi, entry_Si12, entry_Si23, entry_Si34, entry_Si45,
        entry_Si56, entry_Si_intra, entry_Si_inter
        ]
    # Reference rates
    text_rates = tk.Label(root, text="Rates/bond")
    text_rates.place(x=20, y=ary[5])
    # Ag-Si
    text_AgSi_rates = tk.Label(root, text="0")
    text_AgSi_rates.place(x=100, y=ary[5])
    # Si1-2
    text_Si12_rates = tk.Label(root, text="0")
    text_Si12_rates.place(x=180, y=ary[5])
    # Si2-3
    text_Si23_rates = tk.Label(root, text="0")
    text_Si23_rates.place(x=260, y=ary[5])
    # Si3-4
    text_Si34_rates = tk.Label(root, text="0")
    text_Si34_rates.place(x=340, y=ary[5])
    # Si4-5
    text_Si45_rates = tk.Label(root, text="0")
    text_Si45_rates.place(x=420, y=ary[5])
    # Si5-6
    text_Si56_rates = tk.Label(root, text="0")
    text_Si56_rates.place(x=500, y=ary[5])
    # else between layers
    text_Si_intra_rates = tk.Label(root, text="0")
    text_Si_intra_rates.place(x=580, y=ary[5])
    # else inter layers
    text_Si_inter_rates = tk.Label(root, text="0")
    text_Si_inter_rates.place(x=660, y=ary[5])
    # Agtop
    text_Agtp_rates = tk.Label(root, text="0")
    text_Agtp_rates.place(x=740, y=ary[5])
    #rates list
    rates_list = [
        text_AgSi_rates, text_Si12_rates, text_Si23_rates, text_Si34_rates,
        text_Si45_rates, text_Si56_rates, text_Si_intra_rates, text_Si_inter_rates
    ]
    #set transformation
    bln_transformation = tk.BooleanVar()
    bln_transformation.set(True)
    chk_transformation = tk.Checkbutton(root, variable=bln_transformation, text="Transformation")
    chk_transformation.place(x=20, y=ary[7])
    entry_transformation = tk.Entry(root, text="transformation", width=7)
    entry_transformation.place(x=200, y=ary[7] + 5)
    entry_transformation.delete(0, tk.END)
    entry_transformation.insert(tk.END, input_params["transformation_energy"])
    entry_transformation.bind("<Return>", update)
    #set keep defects
    bln_keep_defect = tk.BooleanVar()
    bln_keep_defect.set(True)
    chk_keep_defect = tk.Checkbutton(
        root, variable=bln_keep_defect, text="Keep a defect in first layer"
    )
    chk_keep_defect.place(x=20, y=ary[8])
    #Set recording settings
    text_record = tk.Label(root, text="Record")
    text_record.place(x=20, y=ary[9])
    entry_record = tk.Entry(root, text="Name", width=50)
    entry_record.place(x=100, y=ary[9])
    entry_record.delete(0, tk.END)
    entry_record.insert(tk.END, input_params["record_name"])
    #Image recording setting
    text_img_every = tk.Label(root, text="Image rec. (%) :  ")
    text_img_every.place(x=450, y=ary[9])
    entry_img_every = tk.Entry(root, text="img", width=10)
    entry_img_every.place(x=570, y=ary[9])
    entry_img_every.delete(0, tk.END)
    entry_img_every.insert(tk.END, input_params["record_image_every"])
    #Set comments
    text_comment = tk.Label(root, text="Comments")
    text_comment.place(x=20, y=ary[10])
    entry_comment = tk.Entry(root, text="Comments", width=110)
    entry_comment.place(x=100, y=ary[10])
    entry_comment.delete(0, tk.END)
    entry_comment.insert(tk.END, input_params["comments"])
    #Set progress bar
    input_params["progress_value"] = 0
    progress_bar = ttk.Progressbar(root, orient=tk.HORIZONTAL, length=350, mode="determinate")
    progress_bar.configure(maximum=100, value=input_params["progress_value"])
    progress_bar.place(x=350, y=ary[11] + 5)
    #Other statuses
    text_count_progress = tk.Label(root, text="Waiting")
    text_count_progress.place(x=720, y=ary[11] + 5)
    #
    text_time_progress = tk.Label(root, text="time (s)")
    text_time_progress.place(x=370, y=ary[11] + 35)
    #
    text_event = tk.Label(root, text="events")
    text_event.place(x=470, y=ary[11] + 35)
    #
    text_number_atoms = tk.Label(root, text="Num. atoms")
    text_number_atoms.place(x=570, y=ary[11] + 35)
    #
    text_coverage = tk.Label(root, text="Coverage")
    text_coverage.place(x=670, y=ary[11] + 35)
    #
    button_start = tk.Button(
        root, text="Start", command=button_start_clicked, height=1, width=20
    )
    button_start.place(x=20, y=ary[11])
    #
    button_close = tk.Button(
        root, text="Close", command=button_close_clicked, height=1, width=20
    )
    button_close.place(x=180, y=ary[11])
    update_values()
    root.mainloop()
    """
    
