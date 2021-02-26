from Rejection_free_kmc import rejection_free
from Modules.heatmap import form_heatmap

if __name__ == "__main__":
    start_E1 = -0.2
    end_E1 = -0.6
    diff_E1 = -0.1
    start_E2 = -0.2
    end_E2 = -0.6
    diff_E2 = -0.1
    #
    E1_list = [
        start_E1 + diff_E1 * i for i in range(int((end_E1 - start_E1) / diff_E1))
    ]
    E2_list = [
        start_E2 + diff_E2 * i for i in range(int((end_E2 - start_E2) / diff_E2))
    ]

    growth_rec = []
    #
    for i in E1_list:
        growth_rec.append([])
        for k in E2_list:
            rf_class = rejection_free(1)
            rf_class.init_value.binding_energies["Si_2nd"] = i
            rf_class.init_value.binding_energies["Si_3rd"] = k
            rf_class.init_value.binding_energies["Si_else"] = k
            rf_class.start()
            growth_rec[-1].append(rf_class.growth_mode)

    form_heatmap(growth_rec, rf_class.coverage, E1_list, E2_list, diff_E1, diff_E2)
    print("heatmap formed")
