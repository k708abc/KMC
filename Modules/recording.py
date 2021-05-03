from typing import List, Dict
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import math
from Modules.record_ppt import rec_ppt
import os
import glob
from multiprocessing import Pool

unit_x: List[float] = [1, 0, 0]
unit_y: List[float] = [0.5, 0.866, 0]


def highest_z(pos_all: List[Dict]) -> int:
    maxz = 0
    for positions in pos_all:
        for index, state in positions.items():
            if (state != 0) and (index[2] > maxz):
                maxz = index[2]
    return maxz


def triangle(xp, yp, z):
    if z % 2 == 0:
        t1x = xp * unit_x[0] + yp * unit_y[0] - (unit_x[0] + unit_y[0]) / 3
        t1y = xp * unit_x[1] + yp * unit_y[1] - (unit_x[1] + unit_y[1]) / 3
        t2x = t1x + unit_x[0]
        t2y = t1y + unit_x[1]
        t3x = t1x + unit_y[0]
        t3y = t1y + unit_y[1]
    else:
        t1x = xp * unit_x[0] + yp * unit_y[0] + (unit_x[0] + unit_y[0]) / 3
        t1y = xp * unit_x[1] + yp * unit_y[1] + (unit_x[1] + unit_y[1]) / 3
        t2x = t1x - unit_x[0]
        t2y = t1y - unit_x[1]
        t3x = t1x - unit_y[0]
        t3y = t1y - unit_y[1]
    return [(t1x, t1y), (t2x, t2y), (t3x, t3y)]


def color_determinate(z, maxz):

    color_num = math.floor(z / 2)
    maxz_BL = math.floor(maxz / 2) - 0.999
    if color_num == 0:
        color = [0, 1, 0]
    else:
        color = [(color_num - 1) / maxz_BL, 0, 1 - (color_num - 1) / maxz_BL]
    return color

def image_formaiton(args):
    pos, lattice, length, maxz, img_name = args

    unit_x: List[float] = [1, 0, 0]
    unit_y: List[float] = [0.5, 0.866, 0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #
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
    #
    for z in range(maxz + 1):
        for x in range(length):
            for y in range(length):
                if pos[(x, y, z)] == 0:
                    pass
                else:
                    xp = lattice[(x, y, z)][0]
                    yp = lattice[(x, y, z)][1]
                    tri_pos = triangle(xp, yp, z)
                    color = color_determinate(z, maxz)
                    p = pat.Polygon(xy=tri_pos, fc=color, ec=color)
                    ax.add_patch(p)

                    if x == 0:
                        xpl = xp + length
                        ypl = yp
                        tri_pos = triangle(xpl, ypl, z)

                        p = pat.Polygon(xy=tri_pos, fc=color, ec=color)
                        ax.add_patch(p)

                    if y == 0:
                        xpl = xp
                        ypl = yp + length
                        tri_pos = triangle(xpl, ypl, z)
                        p = pat.Polygon(xy=tri_pos, fc=color, ec=color)
                        ax.add_patch(p)

                    if (x == 0) and (y == 0):
                        xpl = xp + length
                        ypl = yp + length
                        tri_pos = triangle(xpl, ypl, z)

                        p = pat.Polygon(xy=tri_pos, fc=color, ec=color)
                        ax.add_patch(p)

                    if x == length - 1:
                        xpl = xp - length
                        ypl = yp
                        tri_pos = triangle(xpl, ypl, z)

                        p = pat.Polygon(xy=tri_pos, fc=color, ec=color)
                        ax.add_patch(p)

                    if y == length - 1:
                        xpl = xp
                        ypl = yp - length

                        tri_pos = triangle(xpl, ypl, z)

                        p = pat.Polygon(xy=tri_pos, fc=color, ec=color)
                        ax.add_patch(p)

                    if (x == length - 1) and (y == length - 1):
                        xpl = xp - length
                        ypl = yp - length
                        tri_pos = triangle(xpl, ypl, z)

                        p = pat.Polygon(xy=tri_pos, fc=color, ec=color)
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
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.set_aspect("equal")
    fig.savefig(img_name)


def hist_formation(args):
    pos, maxz, n_BL, img_name = args
    hist: List[float] = [0 for _ in range(math.floor(maxz / 2))]
    left: List[float] = [i + 1 for i in range(math.floor(maxz / 2))]
    for pos_index, atom_state in pos.items():
        if atom_state == 1:
            hist[pos_index[2] // 2] += 1 / n_BL * 100

    fig = plt.figure()
    bx = fig.add_subplot(111)
    bx.barh(left, hist, label="Number")
    bx.set_yticks(left)
    # bx.legend(fontsize=16)
    # bx.set_yticklabels(lay)
    bx.set_xlim(0, 100)
    plt.xlabel("Coverage", fontsize=20)
    plt.ylabel("layer number", fontsize=20)
    plt.tick_params(labelsize=18)
    fig.subplots_adjust(bottom=0.2)
    #
    fig.savefig(img_name)


def rec_poscar(pos: Dict, lattice: Dict, unit_length: int, maxz: int, rec_name: str):
    xp: list[float] = [[]]
    yp: list[float] = [[]]
    zp: list[float] = [[]]
    zmax = 1
    atom_list = [
        "B",
        "C",
        "N",
        "O",
        "F",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
    ]
    #
    x_Ag: list[float] = []
    y_Ag: list[float] = []
    Si_unit = 0.384
    Ag_unit = 0.289
    unit_size = unit_length * Si_unit
    Ag_num = round(unit_size / Ag_unit)
    x_Ag = [i / Ag_num for i in range(Ag_num) for _ in range(Ag_num)]
    y_Ag = [i / Ag_num for _ in range(Ag_num) for i in range(Ag_num)]
    #

    for index, atom_state in pos.items():
        if atom_state != 0:
            z_val = index[2]
            if z_val > zmax:
                while zmax < z_val:
                    xp.append([])
                    yp.append([])
                    zp.append([])
                    zmax += 2
            xp[z_val // 2].append(lattice[index][0] / unit_length)
            yp[z_val // 2].append(lattice[index][1] / unit_length)
            zp[z_val // 2].append((lattice[index][2] + 0.7) / maxz / 2.448)
    #
    file_data = open(rec_name, "w")
    file_data.write(rec_name + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    #
    file_data.write("Ag" + "\t")
    for i in range(len(zp)):
        file_data.write(atom_list[i] + "\t")
    file_data.write("\n")
    file_data.write(str(len(x_Ag)) + "\t")
    for i in range(len(zp)):
        file_data.write(str(len(zp[i])) + "\t")
    file_data.write("\n")
    file_data.write("direct" + "\n")
    #
    for i in range(len(x_Ag)):
        file_data.write(str(x_Ag[i]) + "\t" + str(y_Ag[i]) + "\t" + str(0) + "\n")
    for (i, j, k) in zip(xp, yp, zp):
        for (xval, yval, zval) in zip(i, j, k):
            file_data.write(str(xval) + "\t" + str(yval) + "\t" + str(zval) + "\n")
    file_data.close()


def dir_formarion(name: str):
    os.makedirs(name, exist_ok=True)


def rec_events_per_dep_fig(args):
    atoms, num_events, fig_name = args
    fig = plt.figure()
    plt.plot(atoms, num_events)
    plt.xlabel("Num. Atoms")
    plt.ylabel("Num. Events")
    plt.rcParams["font.size"] = 18
    fig.subplots_adjust(bottom=0.2, left=0.2)
    plt.savefig(fig_name)

def rec_events_per_dep(num_events: List[int], atoms: List[int], params):
    dir_name = "Record/" + params.record_name + "/"
    if os.path.exists(dir_name) is False:
        os.mkdir(dir_name)
    fig_name = dir_name + "Num_events_per_dep.png"
    p = Pool(1)
    p.map(rec_events_per_dep_fig, [[atoms, num_events, fig_name]])
    p.close()
    #
    file_data = open(dir_name + "Num_events_per_dep.txt", "w")
    file_data.write("atoms" + "\t" + "events" + "\n")

    for i in range(len(atoms)):
        file_data.write(str(atoms[i]) + "\t" + str(num_events[i]) + "\n")
    file_data.close()


def rec_growth_mode(args):
    growth_mode, coverage, params = args
    ag = []
    first = []
    multi = []
    for i in growth_mode:
        ag.append(i[0])
        first.append(i[1])
        multi.append(i[2])
    #
    dir_name = "Record/" + params.record_name + "/"
    if os.path.exists(dir_name) is False:
        os.mkdir(dir_name)
    fig = plt.figure()
    plt.plot(coverage, ag, label="Ag")
    plt.plot(coverage, first, label="First ML")
    plt.plot(coverage, multi, label="Multi layer")
    plt.legend(bbox_to_anchor=(1, 1), loc="upper right", borderaxespad=0, fontsize=14)
    plt.ylim(0, 100)
    plt.xlabel("Total amount (ML)")
    plt.ylabel("Coverage (%)")
    plt.rcParams["font.size"] = 18
    fig.subplots_adjust(bottom=0.2, left=0.2)
    plt.savefig(dir_name + "Coverage_change.png")


def height_check(pos_x, pos_y, pos, maxz):
    max_pos = -1
    for z in range(maxz + 1):
        if pos[(pos_x, pos_y, z)] != 0:
            max_pos = z
    return max_pos


def growth_check(pos, unit_length, maxz, atom_BL):
    num_ag = 0
    num_1st = 0
    num_2nd = 0
    num_multi = 0
    for i in range(unit_length):
        for k in range(unit_length):
            height = height_check(i, k, pos, maxz)
            if height == -1:
                num_ag += 2
            elif height == 0:
                num_ag += 1
                num_1st += 1
            elif height == 1:
                num_1st += 2
            elif height == 2:
                num_1st += 1
                num_2nd += 1
            elif height == 3:
                num_2nd += 2
            elif height == 4:
                num_2nd += 1
                num_multi += 1
            else:
                num_multi += 2
    return (
        round(num_ag / atom_BL * 100, 2),
        round(num_1st / atom_BL * 100, 2),
        round(num_2nd / atom_BL * 100, 2),
        round(num_multi / atom_BL * 100, 2),
    )


def num_check(coverage, num):
    value = 100
    number = 0
    for i in range(len(coverage)):
        cf_val = abs(num - coverage[i])
        if value > cf_val:
            number = i
            value = cf_val
    return number


def growth_val(growth_mode, coverage, num_1ML, num_2ML):
    First_at_1ML = growth_mode[num_1ML][1]
    Second_at_2ML = growth_mode[num_2ML][2]
    if growth_mode[num_2ML][3] == 100:
        First_at_second = 0
    else:
        First_at_second = (
            growth_mode[num_2ML][1] / (100 - growth_mode[num_2ML][3]) * 100
        )
    if First_at_1ML < 50:
        if Second_at_2ML < 50:
            mode = "VW"
        else:
            mode = "BL"
    else:
        if Second_at_2ML >= 50:
            mode = "FM"
        else:
            if First_at_second >= 50:
                mode = "SK"
            else:
                mode = "DW"

    return [
        First_at_1ML,
        Second_at_2ML,
        First_at_second,
        mode,
        coverage[num_1ML],
        coverage[num_2ML],
    ]


def mode_check(growth_mode, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    return growth_val(growth_mode, coverage, num_1ML, num_2ML)


def delete_images(dir_name, recursive=True):
    dir_name = dir_name + "*.png"
    for p in glob.glob(dir_name, recursive=recursive):
        if os.path.isfile(p):
            os.remove(p)


def record_data(
    pos_all: List[dict],
    time: List,
    coverage: List,
    lattice: Dict,
    params,
    minute: int,
    second: float,
    time_per_eve,
): 

    maxz = highest_z(pos_all)
    maxz_unit = params.cell_size_z * 6
    unit_length = params.cell_size_xy
    rec_name_body = params.record_name
    n_BL = params.atoms_in_BL
    img_names: List[str] = []
    hist_names: List[str] = []
    dir_name = "Record/" + params.record_name + "/"

    dir_formarion(dir_name)
    growth_mode = []

    for rec_num, (pos_i, time_i, cov_i) in enumerate(zip(pos_all, time, coverage)):
        rec_name = (
            str(rec_name_body)
            + "_"
            + str(rec_num)
            + "_"
            + str(int(time_i))
            + "s_"
            + str(round(cov_i, 2))
            + "_ML"
        )
        #
        img_name = dir_name + rec_name + ".png"
        #
        p = Pool(1)
        p.map(image_formaiton, [[pos_i, lattice, unit_length, maxz, img_name]])
        p.close()
        img_names.append(img_name)
        #
        hist_name = dir_name + rec_name + "_hist.png"
        p = Pool(1)
        p.map(hist_formation, [[pos_i, maxz_unit, n_BL, hist_name]])
        p.close()
        hist_names.append(hist_name)
        #
        poscar_name = dir_name + rec_name + "_poscar.vasp"
        rec_poscar(pos_i, lattice, unit_length, maxz, poscar_name)
        #
        growth_mode.append(growth_check(pos_i, unit_length, maxz, params.atoms_in_BL))
    #
    p = Pool(1)
    p.map(rec_growth_mode, [[growth_mode, coverage, params]])
    p.close()
    mode_val = mode_check(growth_mode, coverage)
    #
    rec_ppt(
        params,
        minute,
        second,
        img_names,
        hist_names,
        time,
        coverage,
        dir_name,
        time_per_eve,
        growth_mode,
        mode_val,
    )
    delete_images(dir_name)
    return mode_val
    #
