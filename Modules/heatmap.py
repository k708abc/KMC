import matplotlib.patches as patches
import matplotlib.pyplot as plt


def num_check(coverage, num):
    value = 100
    number = 0
    for i in range(len(coverage)):
        cf_val = abs(num - coverage[i])
        if value > cf_val:
            number = i
            value = cf_val
    return number


def reorder_growth(growth_mode, num_1ML, num_2ML):
    new_growth = []
    for i in range(len(growth_mode)):
        new_growth.append([])
        for k in range(len(growth_mode[i])):
            val1 = growth_mode[i][k][num_1ML][1] / (100 - growth_mode[i][k][num_1ML][2])
            if 100 - growth_mode[i][k][num_2ML][2] == 0:
                val2 = 0
            else:
                val2 = growth_mode[i][k][num_2ML][1] / (
                    100 - growth_mode[i][k][num_2ML][2]
                )
            new_growth[-1].append([val1, val2])
    return new_growth


def form_heatmap(growth_mode, coverage, E1, E2, diff1, diff2):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    #
    growth_judge = reorder_growth(growth_mode, num_1ML, num_2ML)
    #
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    for i in range(len(E1)):
        for k in range(len(E2)):
            RGB = [0, 0, 0]
            s1 = growth_judge[i][k][0]
            s2 = growth_judge[i][k][1]

            sn1 = s1 + s2 - 100
            sn2 = -s1 + s2

            if sn1 >= 0:
                RGB[1] = sn1 / 100
            else:
                RGB[2] = abs(sn1 / 100)

            if sn2 <= 0:
                RGB[0] = abs(sn2 / 100)
            """
            if i == k == 0:
                RGB = [0.2, 0.2, 0.2]
            elif i == k == 1:
                RGB = [0.9, 0.9, 0.9]
            elif i == k == 2:
                RGB = [0.5, 0.5, 0.5]
            """
            r = patches.Rectangle(
                xy=(E1[i] - diff1 / 2, E2[k] - diff2 / 2),
                width=diff1,
                height=diff2,
                color=RGB,
                fill=True,
            )

            ax.add_patch(r)

    ax.set_xlim(min(E1) + diff1 / 2, max(E1) - diff1 / 2)
    ax.set_ylim(min(E2) + diff2 / 2, max(E2) - diff2 / 2)

    E1_2 = []
    E2_2 = []

    if len(E1) >= 10:
        for i in range(len(E1)):
            if i % 2 == 0:
                E1_2.append(E1[i])
    else:
        for i in range(len(E1)):
            E1_2.append(E1[i])

    if len(E2) >= 10:
        for i in range(len(E2)):
            if i % 2 == 0:
                E2_2.append(E2[i])
    else:
        for i in range(len(E2)):
            E2_2.append(E2[i])

    ax.set_xticks(E1_2)
    ax.set_yticks(E2_2)
    ax.set_xlabel("First layer")
    ax.set_ylabel("Multi layer")

    ax.invert_xaxis()
    ax.invert_yaxis()
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=45, fontsize=16)
    plt.tick_params(labelsize=16)
    plt.savefig("heatmap.png")
    #
    file_name = "heat_values.txt"
    file_data = open(file_name, "w")

    file_data.write("E1" + "\t" + "E2" + "\t" + "ML" + "\t" + "BL" + "\n")

    for i in range(len(E1)):
        for k in range(len(E2)):
            file_data.write(
                str(E1[i])
                + "\t"
                + str(E2[k])
                + "\t"
                + str(growth_judge[i][k][0])
                + "\t"
                + str(growth_judge[i][k][1])
                + "\n"
            )

    file_data.close()
