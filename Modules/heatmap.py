import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import mathtext

# from pptx import Presentation
# from pptx.util import Inches, Pt
# import os


def heatmap_image(dx, dy, xlist, ylist, modelist):
    mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
    plt.rcParams.update(
        {
            "mathtext.default": "default",
            "mathtext.fontset": "stix",
            "font.family": "Times New Roman",
            "font.size": 12,
            "figure.figsize": (6, 6),
        }
    )

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    for x, y, mode in zip(ylist, xlist, modelist):
        RGB = [0, 0, 0]
        if mode == "DW":
            RGB = [1, 0, 0]

        elif mode == "BL":
            RGB = [0.1, 0.1, 0.1]

        elif mode == "SK":
            RGB = [0, 1, 0]

        elif mode == "FM":
            RGB = [0, 0, 1]

        elif mode == "VW":
            RGB = [1, 1, 0]

        r = patches.Rectangle(
            xy=(x - dx / 2, y - dy / 2),
            width=dx,
            height=dy,
            facecolor=RGB,
            edgecolor="white",
            fill=True,
        )

        ax.add_patch(r)

    ax.set_xlim(min(ylist) - dy / 2, max(ylist) + dy / 2)
    ax.set_ylim(min(xlist) - dx / 2, max(xlist) + dx / 2)
    ax.set_xticks(ylist)
    ax.set_yticks(xlist)
    ax.set_xlabel("$E_{a,bulk}$ (eV)", fontsize=24)
    ax.set_ylabel("$E_{a,2}$ (eV)", fontsize=24)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=45, fontsize=22)
    plt.tick_params(labelsize=22)
    plt.tight_layout()
    #
    heat_name = "Record/heatmap.png"
    plt.savefig(heat_name)


def rec_modes(xlist, ylist, modelist):
    file_name = "Record/growth_modes.txt"
    with open(file_name, "a") as f:
        for x, y, mode in zip(xlist, ylist, modelist):
            print(str(x) + "\t" + str(y) + "\t" + mode, file=f)


def form_heatmap(growth_mode, rec_name, E1, E2, diff1, diff2):
    xlist = []
    ylist = []
    modelist = []
    rec_name = "Record/" + rec_name
    for i in range(len(E1)):
        for k in range(len(E2)):
            xlist.append(E1[i])
            ylist.append(E2[k])
            modelist.append(growth_mode[i][k][3])
    # heatmap_image(diff1, diff2, xlist, ylist, modelist)
    rec_modes(xlist, ylist, modelist)
