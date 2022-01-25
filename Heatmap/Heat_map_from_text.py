import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import mathtext
import glob
import os


def form_examples(RGB_set, patterns):
    for i in range(len(RGB_set)):
        RGB = RGB_set[i]
        pattern = patterns[i]
        fig = plt.figure(figsize=(1, 1))
        ax = fig.add_subplot(111)
        r = patches.Rectangle(
            xy=(0, 0),
            width=1,
            height=1,
            facecolor=RGB,
            edgecolor="white",
            hatch=pattern,
            fill=True,
        )
        ax.add_patch(r)
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        heat_name = "example_" + str(i) + ".png"
        plt.savefig(heat_name)


def heatmap_image(dx, dy, xlist, ylist, modelist, RGB_set, patterns, data_name):
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
            RGB = RGB_set[0]
            pattern = patterns[0]

        elif mode == "BL":
            RGB = RGB_set[1]
            pattern = patterns[1]

        elif mode == "SK":
            RGB = RGB_set[2]
            pattern = patterns[2]

        elif mode == "FM":
            RGB = RGB_set[3]
            pattern = patterns[3]

        elif mode == "VW":
            RGB = RGB_set[4]
            pattern = patterns[4]

        r = patches.Rectangle(
            xy=(x - dx / 2, y - dy / 2),
            width=dx,
            height=dy,
            facecolor=RGB,
            edgecolor="white",
            hatch=pattern,
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

    heat_name = os.path.splitext(os.path.basename(data_name))[0] + ".pdf"

    plt.savefig(heat_name)


if __name__ == "__main__":

    xlist = []
    ylist = []
    modelist = []
    RGB_set = [
        [1, 0, 0],
        [0.1, 0.1, 0.1],
        [0.21, 0.75, 0.4],
        [0.7, 0.7, 0.7],
        [0, 0, 1],
    ]
    patterns = ["/", None, "\\", ".", "x"]

    """
        def read_image_list(self, dir_name):
            image_list = glob.glob("*")
            image_list = [os.path.basename(pathname) for pathname in image_list]
            image_list2 = []
            for file in image_list:
                data_path = dir_name + file
                imtype = get_datatypes(data_path)
                if imtype[0] == "folder":
                    if (data_path + "/") not in self.folder_list:
                        self.folder_list.append(data_path + "/")
                else:
                    image_list2.append(file)
            image_list2 = self.folder_list + image_list2
            return image_list2



    """
    data_list = glob.glob("*.txt")

    for data in data_list:
        with open(data, "r") as f:
            for line in f:
                line_s = line.split("\t")
                if line_s[0] == "dx":
                    dx = float(line_s[1])
                elif line_s[0] == "dy":
                    dy = float(line_s[1])
                else:
                    xlist.append(float(line_s[0]))
                    ylist.append(float(line_s[1]))
                    modelist.append(line_s[2].replace("\n", ""))

        heatmap_image(dx, dy, xlist, ylist, modelist, RGB_set, patterns, data)
    form_examples(RGB_set, patterns)
    print("heat map is formed")
    input()
