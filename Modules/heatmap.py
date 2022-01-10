import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import mathtext


def rec_modes(xlist, ylist, modelist, diff1, diff2, num):
    file_name = "Record/growth_modes" + str(num) + ".txt"
    with open(file_name, "a") as f:
        f.write("dx" + "\t" + str(diff1) + "\n")
        f.write("dy" + "\t" + str(diff2) + "\n")
        for x, y, mode in zip(xlist, ylist, modelist):
            print(str(x) + "\t" + str(y) + "\t" + mode, file=f)


def form_heatmap(growth_mode, rec_name, E1, E2, diff1, diff2, num):
    xlist = []
    ylist = []
    modelist = []
    rec_name = "Record/" + rec_name
    for i in range(len(E1)):
        for k in range(len(E2)):
            xlist.append(E1[i])
            ylist.append(E2[k])
            modelist.append(growth_mode[i][k][3])
    rec_modes(xlist, ylist, modelist, diff1, diff2, num)
