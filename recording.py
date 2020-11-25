from typing import List, Dict
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import math
from record_ppt import rec_ppt
import os


def highest_z(pos_all: List[dict]) -> int:
    maxz = 1
    for positions in pos_all:
        for index, state in positions.items():
            if (state != 0) and (index[2] + 1 > maxz):
                maxz = index[2] + 1
    return maxz


def image_formaiton(pos: Dict, lattice: Dict, length: int, maxz: int):
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
    for z in range(maxz):
        for x in range(length):
            for y in range(length):
                if pos[(x, y, z)] == 0:
                    pass
                else:
                    xp = lattice[(x, y, z)][0]
                    yp = lattice[(x, y, z)][1]
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

                    color_num = math.floor(z / 2) - 1
                    maxz_BL = math.floor(maxz / 2) - 1

                    if color_num == -1:
                        color = [0, 1, 0]
                    else:
                        color = [color_num / maxz_BL, 0, 1 - color_num / maxz_BL]

                    p = pat.Polygon(
                        xy=[(t1x, t1y), (t2x, t2y), (t3x, t3y)], fc=color, ec=color
                    )
                    ax.add_patch(p)

                    if x == 0:
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

                    if y == 0:
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

                    if (x == 0) and (y == 0):

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


def rec_img(img, name: str):
    img.savefig(name)


def hist_formation(pos: Dict, maxz: int, n_BL: int):
    hist: List[float] = [0 for _ in range(math.ceil(maxz / 2))]
    left: List[float] = [i for i in range(math.ceil(maxz / 2))]
    for pos_index, atom_state in pos.items():
        if atom_state != 0:
            hist[pos_index[2] // 2] += 1 / n_BL * 100
    fig = plt.figure()
    bx = fig.add_subplot(111)
    bx.barh(left, hist)
    bx.set_yticks(left)
    # bx.set_yticklabels(lay)
    bx.set_xlabel("Coverage")
    bx.set_ylabel("layer number")
    bx.set_xlim(0, 100)
    return fig


def rec_poscar(pos: Dict, lattice: Dict, unit_length: int, maxz: int, rec_name: str):
    xp: list[float] = []
    yp: list[float] = []
    zp: list[float] = []
    atom_i = 0
    for index, atom_state in pos.items():
        if atom_state != 0:
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
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


def dir_formarion(name: str):
    os.makedirs(name, exist_ok=True)


def record_data(
    pos_all: List[dict],
    time: List,
    coverage: List,
    lattice: Dict,
    params,
    minute: int,
    second: float,
    defect,
):
    maxz = highest_z(pos_all)
    maxz_unit = params.z_unit_init * 6
    unit_length = params.n_cell_init
    rec_name_body = params.record_name
    n_BL = params.atoms_in_BL
    img_names: List[str] = []
    hist_names: List[str] = []
    dir_name = params.record_name + "/"
    dir_formarion(dir_name)

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
        img = image_formaiton(pos_i, lattice, unit_length, maxz)
        img_name = dir_name + rec_name + ".png"
        rec_img(img, img_name)
        img_names.append(img_name)
        #
        hist = hist_formation(pos_i, maxz_unit, n_BL)
        hist_name = dir_name + rec_name + "_hist.png"
        rec_img(hist, hist_name)
        hist_names.append(hist_name)
        #
        poscar_name = dir_name + rec_name + "_poscar.vasp"
        rec_poscar(pos_i, lattice, unit_length, maxz, poscar_name)
    #
    rec_ppt(
        params, minute, second, img_names, hist_names, time, coverage, dir_name, defect
    )
