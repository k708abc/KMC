from Rejection_free_kmc import rejection_free
from InputParameter import Params
from Modules.heatmap import form_heatmap

if __name__ == "__main__":
    first_input = Params()
    start_E1 = first_input.start_E1
    end_E1 = first_input.end_E1
    diff_E1 = first_input.diff_E1
    start_E2 = first_input.start_E2
    end_E2 = first_input.end_E2
    diff_E2 = first_input.diff_E2
    or_rec_name = first_input.record_name
    #
    if diff_E1 != 0:
        E1_list = [
            round(start_E1 + diff_E1 * i, 3)
            for i in range(int(round((end_E1 - start_E1) / diff_E1)) + 1)
        ]
    else:
        E1_list = [start_E1]
    if diff_E2 != 0:
        E2_list = [
            round(start_E2 + diff_E2 * i, 3)
            for i in range(int(round((end_E2 - start_E2) / diff_E2)) + 1)
        ]
    else:
        E2_list = [start_E2]
    #
    print(E1_list)
    print(E2_list)
    print("repetetion start")

    growth_rec = []
    total_cal = len(E1_list) * len(E2_list)
    num_cal = 0

    #
    for i in E1_list:
        growth_rec.append([])
        for k in E2_list:
            num_cal += 1
            print("starting:" + str(num_cal) + "/" + str(total_cal))
            print("val 1 = " + str(i))
            print("val 2 = " + str(k))
            rec_name = or_rec_name + "_" + str(i) + "_" + str(k)
            rf_class = rejection_free(1)
            rf_class.init_value.binding_energies["Si_second"] = i
            rf_class.init_value.binding_energies["Si_third"] = k
            rf_class.init_value.binding_energies["Si_else"] = k
            rf_class.init_value.record_name = rec_name
            rf_class.start()
            growth_rec[-1].append(rf_class.mode_val)

    form_heatmap(
        growth_rec, rf_class.init_value.record_name, E1_list, E2_list, diff_E1, diff_E2
    )
    print("heatmap formed")
    input()
