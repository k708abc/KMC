# cython: language_level=3, boundscheck=False, wraparound=False

import matplotlib.pyplot as plt
import matplotlib.patches as pat
import math
from record_ppt cimport rec_ppt
from Calc_grid_index cimport grid_num
import os
import glob
import yaml
from growth_mode_determination import (
    growth_check,
    mode_check_prev,
    mode_check_1,
    mode_check_2,
    mode_check_3,
    mode_check_4,
    mode_check_5,
    mode_check_6,
)
from InputParameter cimport Params


cdef int highest_z(list pos_all, list index_list):
    cdef int maxz
    cdef list positions
    cdef int index, state, _, z
    maxz = 0
    for positions in pos_all:
        for index, state in enumerate(positions):
            _, _, z = index_list[index]
            if (state != 0) and (z > maxz):
                maxz = z
    return maxz


cdef dir_formarion(str name):
    os.makedirs(name, exist_ok=True)


cdef list[tuple] triangle(double xp, double yp, int z):
    cdef double t1x, t1y, t2x, t2y, t3x, t3y
    cdef list unit_x = [1.0, 0, 0]
    cdef list unit_y = [0.5, 0.866, 0]
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


cdef list color_determinate(int z, int maxz):
    cdef double color_num, maxz_BL
    cdef list color
    color_num = math.floor(z / 2)
    maxz_BL = math.floor(maxz / 2) - 0.999
    if color_num == 0:
        color = [0, 1, 0]
    else:
        color = [(color_num - 1) / maxz_BL, 0, 1 - (color_num - 1) / maxz_BL]
    return color


cdef image_formaiton(list pos, list lattice, int length, int maxz, str img_name):

    cdef int x, y, z, max_z, index
    cdef double xp, yp
    cdef list tri_pos
    cdef list color
    cdef list unit_x = [1.0, 0, 0]
    cdef list unit_y = [0.5, 0.866, 0]
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
                index = grid_num(x, y, z, length)
                if pos[index] == 0:
                    pass
                else:
                    xp = lattice[index][0]
                    yp = lattice[index][1]
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


cdef list occupation_of_layers(int maxz, pos, int n_BL, list index_list):
    cdef int _, 
    cdef list hist = [0 for _ in range(math.floor(maxz / 2))]
    for pos_index, atom_state in enumerate(pos):
        if atom_state == 1:
            _, _, z = index_list[pos_index]
            hist[z // 2] += 1 / n_BL * 100
    return hist


cdef hist_formation(list pos, int maxz, int n_BL, str img_name, list index_list):
    cdef list hist
    cdef list left
    left= [i + 1 for i in range(math.floor(maxz / 2))]
    #
    hist = occupation_of_layers(maxz, pos, n_BL, index_list)
    #
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


cdef rec_poscar(list pos, list lattice, int unit_length, int maxz, str rec_name, list index_list):
    cdef list xp = [[]]
    cdef list yp = [[]]
    cdef list zp = [[]]
    cdef int zmax
    cdef list atom_list
    cdef list x_Ag = []
    cdef list y_Ag = []
    cdef double Si_unit = 0.384
    cdef double Ag_unit = 0.289
    cdef double unit_size

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
    unit_size = unit_length * Si_unit
    Ag_num = round(unit_size / Ag_unit)
    x_Ag = [i / Ag_num for i in range(Ag_num) for _ in range(Ag_num)]
    y_Ag = [i / Ag_num for _ in range(Ag_num) for i in range(Ag_num)]
    #

    for index, atom_state in enumerate(pos):
        x_val, y_val, z_val = index_list[index]
        if atom_state != 0:
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


cdef rec_events_per_dep_fig(list atoms, list num_events, str fig_name):
    fig = plt.figure()
    plt.plot(atoms, num_events)
    plt.xlabel("Num. Atoms")
    plt.ylabel("Num. Events")
    plt.rcParams["font.size"] = 18
    fig.subplots_adjust(bottom=0.2, left=0.2)
    plt.savefig(fig_name)


cdef rec_events_per_dep(list num_events, list atoms, str rec_name):
    cdef str dir_name, fig_name
    cdef int i
    dir_name = "Record/" + rec_name + "/"
    if os.path.exists(dir_name) is False:
        os.mkdir(dir_name)
    fig_name = dir_name + "Num_events_per_dep.png"
    #
    """
    if params.start_from_middle is False:
        p_eve = Pool(1)
        p_eve.map(rec_events_per_dep_fig, [[atoms, num_events, fig_name]])
        p_eve.close()
    else:
        rec_events_per_dep_fig([atoms, num_events, fig_name])
    """
    rec_events_per_dep_fig(atoms, num_events, fig_name)
    #
    file_data = open(dir_name + "Num_events_per_dep.txt", "w")
    file_data.write("atoms" + "\t" + "events" + "\n")

    for i in range(len(atoms)):
        file_data.write(str(atoms[i]) + "\t" + str(num_events[i]) + "\n")
    file_data.close()


cdef rec_growth_mode(list growth_mode, list coverage, Params params):
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


cdef delete_images(str dir_name):
    cdef bint recursive=True
    cdef dir_name2 = dir_name + "*.png"
    for p in glob.glob(dir_name2, recursive=recursive):
        if os.path.isfile(p):
            os.remove(p)


cdef rec_yaml(Params init_value, str dir_name):
    cdef str file_name
    file_name = dir_name + "kmc_input.yml"
    yml_write = {key: val for key, val in init_value.__dict__.items()}
    with open(file_name, "w") as file:
        yaml.dump(yml_write, file)


cdef record_data(
    list pos_all,
    list time,
    list coverage,
    list lattice,
    Params params,
    int minute,
    int second,
    time_per_eve,
    Params init_value,
    list index_list,
):
    cdef list pos_i
    cdef double time_i
    cdef double cov_i
    cdef int rec_num
    maxz = highest_z(pos_all, index_list)
    maxz_unit = params.cell_size_z * 6
    unit_length = params.cell_size_xy
    rec_name_body = params.record_name
    n_BL = params.atoms_in_BL()
    img_names: List[str] = []
    hist_names: List[str] = []
    dir_name = "Record/" + params.record_name + "/"

    dir_formarion(dir_name)
    growth_mode = []
    occupation = []

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
        """
        if params.start_from_middle is False:
            p_image = Pool(1)
            p_image.map(
                image_formaiton, [[pos_i, lattice, unit_length, maxz, img_name]]
            )
            p_image.close()
        else:
            image_formaiton([pos_i, lattice, unit_length, maxz, img_name])
            plt.clf()
            plt.close()
        # p.join()
        """
        image_formaiton(pos_i, lattice, unit_length, maxz, img_name)
        plt.clf()
        plt.close()

        img_names.append(img_name)
        #
        hist_name = dir_name + rec_name + "_hist.png"
        """
        if params.start_from_middle is False:
            p_hist = Pool(1)
            p_hist.map(hist_formation, [[pos_i, maxz_unit, n_BL, hist_name, index_list]])
            p_hist.close()
            # p.join()
        else:
            hist_formation([pos_i, maxz_unit, n_BL, hist_name, index_list])
            plt.clf()
            plt.close()
        """
        hist_formation(pos_i, maxz_unit, n_BL, hist_name, index_list)
        plt.clf()
        plt.close()

        hist_names.append(hist_name)
        #
        poscar_name = dir_name + rec_name + "_poscar.vasp"
        rec_poscar(pos_i, lattice, unit_length, maxz, poscar_name, index_list)
        #
        growth_mode.append(growth_check(pos_i, unit_length, maxz, params.atoms_in_BL(), index_list))
        occupation.append(occupation_of_layers(maxz_unit, pos_i, n_BL, index_list))
    #
    """
    if params.start_from_middle is False:
        p_mode = Pool(1)
        p_mode.map(rec_growth_mode, [[growth_mode, coverage, params]])
        p_mode.close()
        # p.join()
    else:
        rec_growth_mode([growth_mode, coverage, params])
        plt.clf()
        plt.close()
    """
    rec_growth_mode(growth_mode, coverage, params)
    plt.clf()
    plt.close()
    mode_val = mode_check_prev(growth_mode, coverage)
    #
    # 2022/1/7 Other judgement
    mode_val_1 = mode_check_1(growth_mode, coverage)
    mode_val_2 = mode_check_2(growth_mode, coverage)
    mode_val_3 = mode_check_3(occupation, coverage)
    mode_val_4 = mode_check_4(occupation, coverage)
    mode_val_5 = mode_check_5(occupation, coverage)
    mode_val_6 = mode_check_6(occupation, coverage)
    other_mode_vals = [
        mode_val_1,
        mode_val_2,
        mode_val_3,
        mode_val_4,
        mode_val_5,
        mode_val_6,
    ]
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
        other_mode_vals,
    )
    delete_images(dir_name)
    #
    #rec_yaml(init_value, dir_name)
    return mode_val, other_mode_vals
