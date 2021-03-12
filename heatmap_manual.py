from Modules.heatmap import form_heatmap

if __name__ == "__main__":
    file_name = "heat_values.txt"
    num = 0
    growth_sum = []
    E1_vals = []
    E2_vals = []
    with open(file_name, "r", encoding="utf-8_sig") as f:
        for row in f:
            row_sp = row.split()
            if num == 0:
                if "diff1" in row_sp:
                    diff1 = float(row_sp[1])
                elif "diff2" in row_sp:
                    diff2 = float(row_sp[1])
                elif "values" in row_sp:
                    num = 1
            elif num == 1:
                num = 2
            elif num == 2:
                E1 = float(row_sp[0])
                E2 = float(row_sp[1])
                value1 = float(row_sp[2])
                value2 = float(row_sp[3])
                if E1 not in E1_vals:
                    E1_vals.append(E1)
                if E2 not in E2_vals:
                    E2_vals.append(E2)
                growth_sum.append([E1, E2, value1, value2])
    #
    E1_vals.sort()
    E2_vals.sort()
    growth_mode = [[[] for i in range(len(E2_vals))] for j in range(len(E1_vals))]
    for val in growth_sum:
        pos1 = E1_vals.index(val[0])
        pos2 = E2_vals.index(val[1])
        growth_mode[pos1][pos2] = [val[2], val[3]]
    """
    print(E1_vals)
    print(E2_vals)
    print(diff1)
    print(diff2)
    print(growth_mode)
    """

    form_heatmap(growth_mode, "heatmap_summary", E1_vals, E2_vals, diff1, diff2)
    print("heatmap formed")