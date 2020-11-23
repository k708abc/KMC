from typing import List
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import math
from record_ppt import rec_ppt


def highest_z(pos):
    maxz = 0
    for positions in pos:
        for key in positions:
            if key[2] > maxz:
                maxz = key[2]
    return maxz


def image_formaiton(pos, lattice, length, maxz):
    unit_x: List[float] = [1, 0, 0]
    unit_y: List[float] = [0.5, 0.866, 0]
    unit_z: List[float] = [0, 0, 1]
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
    for z in range(maxz):
        for x in range(length):
            for y in range(length):
                if pos[(x, y, z)] == 0:
                    pass
                else:
                    xp = lattice[(x, y, z)][0]
                    yp = lattice[(x, y, z)][1]
                    zp = lattice[(x, y, z)][2]
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
                        color = [color_num / maxz, 0, 1 - color_num / maxz]

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

    return fig


def rec_img(img, name):
    img.savefig(name)


def hist_formation(pos, maxz, n_BL):
    hist: List[float] = [0 for _ in range(math.ceil(maxz / 2))]
    left: List[float] = [i for i in range(math.ceil(maxz / 2))]
    for pos_index, atom_state in pos.items():
        if atom_state != 0:
            hist[pos_index[3] // 2] += 1 / n_BL * 100
    fig = plt.figure()
    bx = fig.add_subplot(111)
    bx.barh(left, hist)
    bx.set_yticks(left)
    # bx.set_yticklabels(lay)
    bx.set_xlabel("Coverage")
    bx.set_ylabel("layer number")
    bx.set_xlim(0, 100)
    return fig


def rec_poscar(pos, unit_length, maxz, rec_name):
    xp: list[float] = []
    yp: list[float] = []
    zp: list[float] = []
    atom_i = 0
    for index, atom_state in pos.items():
        if atom_state != 0:
            xp.append(index[0] / unit_length)
            yp.append(index[1] / unit_length)
            zp.append(index[2] / maxz / 2.448)
            atom_i += 1
    file_data = open(rec_name, "w")
    file_data.write(rec_name + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Si" + "\n")
    file_data.write(str(atom_i) + "\n")
    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()


def record_data(position, time, coverage, lattice, params, minute, second):
    maxz = highest_z(position)
    rec_num = 0
    unit_length = params.n_cell_init
    rec_name = params.record_name
    n_BL = params.atoms_in_BL
    img_names: List[str] = []
    hist_names: List[str] = []
    for pos_i, time_i, cov_i in zip(position, time, coverage):
        rec_name = (
            str(rec_name)
            + "_"
            + str(rec_num)
            + "_"
            + str(int(time_i))
            + "s_"
            + str(round(cov_i, 2))
            + "_ML"
        )
        #
        img = image_formaiton(pos_i, lattice, unit_length, maxz)
        rec_img(img, rec_name)
        img_names.append(rec_name + ".png")
        #
        hist = hist_formation(pos_i, maxz, n_BL)
        hist_name = rec_name + "_hist"
        rec_img(hist, hist_name)
        hist_names.append(hist_name + ".png")
        #
        poscar_name = rec_name + "_poscar.vasp"
        rec_poscar(position, unit_length, maxz, poscar_name)
    #
    rec_ppt(params, minute, second, img_names, hist_names, time, coverage)
