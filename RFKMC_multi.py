from concurrent import futures
from Rejection_free_kmc import rejection_free
from Modules.InputParameter import Params
from Modules.heatmap import form_heatmap
import os


def function(first_input, val1, val2, i, total_cals):
    print(
        "("
        + str(val1)
        + ", "
        + str(val2)
        + ")"
        + " started: "
        + str(i + 1)
        + "/"
        + str(total_cals)
    )
    rec_name = first_input.record_name + "_" + str(val1) + "_" + str(val2)
    rf_class = rejection_free(2)

    if first_input.repeat_val == 0:
        rf_class.init_value.energies_binding["second"] = val1
        rf_class.init_value.energies_binding["third"] = val2
        rf_class.init_value.energies_binding["upper"] = val2
        val0 = rf_class.init_value.energies_binding["first"]

    elif first_input.repeat_val == 1:
        rf_class.init_value.energies_diffusion["second"] = val1
        rf_class.init_value.energies_diffusion["third"] = val2
        rf_class.init_value.energies_diffusion["upper"] = val2
        val0 = rf_class.init_value.energies_diffusion["first"]
    #
    """
    trans_energy_gain = (val2 - val0) + (val2 - val1)
    if trans_energy_gain < 0:
        rf_class.init_value.trans_check = False
    """

    rf_class.init_value.record_name = rec_name
    rf_class.start()
    print("(" + str(val1) + ", " + str(val2) + ")" + " finished")
    return rf_class.mode_val, rf_class.other_modes


def run_multi():
    first_input = Params("kmc_input.yml")
    start_E1 = first_input.repeat_E1_start
    end_E1 = first_input.repeat_E1_end
    diff_E1 = first_input.repeat_E1_diff
    start_E2 = first_input.repeat_E2_start
    end_E2 = first_input.repeat_E2_end
    diff_E2 = first_input.repeat_E2_diff
    para_val = first_input.max_workers_val
    if os.path.exists("Record") is False:
        os.mkdir("Record")
    E1_range = abs(end_E1 - start_E1)
    E2_range = abs(end_E2 - start_E2)
    if diff_E1 != 0:
        E1_list = [
            round(start_E1 + diff_E1 * i, 3)
            for i in range(int(round(E1_range / diff_E1)) + 1)
        ]
    else:
        E1_list = [start_E1]
    if diff_E2 != 0:
        E2_list = [
            round(start_E2 + diff_E2 * i, 3)
            for i in range(int(round(E2_range / diff_E2)) + 1)
        ]
    else:
        E2_list = [start_E2]
    #
    energy_list = [(val1, val2) for val1 in E1_list for val2 in E2_list]
    growth_list = [[0 for k in E2_list] for i in E1_list]
    growth_list_1 = [[0 for k in E2_list] for i in E1_list]
    growth_list_2 = [[0 for k in E2_list] for i in E1_list]
    growth_list_3 = [[0 for k in E2_list] for i in E1_list]
    growth_list_4 = [[0 for k in E2_list] for i in E1_list]
    growth_list_5 = [[0 for k in E2_list] for i in E1_list]
    growth_list_6 = [[0 for k in E2_list] for i in E1_list]
    #
    print("RFKMC_multi started")
    print("E1: " + str(E1_list))
    print("E2: " + str(E2_list))
    #
    with futures.ProcessPoolExecutor(max_workers=para_val) as executor:
        future_dict = {
            executor.submit(
                function, first_input, var[0], var[1], i, len(energy_list)
            ): var
            for i, var in enumerate(energy_list)
        }

        for future in futures.as_completed(future_dict):
            values = future_dict[future]
            val1_index = E1_list.index(values[0])
            val2_index = E2_list.index(values[1])
            growth_list[val1_index][val2_index] = future.result()[0]
            growth_list_1[val1_index][val2_index] = future.result()[1][0]
            growth_list_2[val1_index][val2_index] = future.result()[1][1]
            growth_list_3[val1_index][val2_index] = future.result()[1][2]
            growth_list_4[val1_index][val2_index] = future.result()[1][3]
            growth_list_5[val1_index][val2_index] = future.result()[1][4]
            growth_list_6[val1_index][val2_index] = future.result()[1][5]

    form_heatmap(
        growth_list, first_input.record_name, E1_list, E2_list, diff_E1, diff_E2, 0
    )
    form_heatmap(
        growth_list_1, first_input.record_name, E1_list, E2_list, diff_E1, diff_E2, 1
    )
    form_heatmap(
        growth_list_2, first_input.record_name, E1_list, E2_list, diff_E1, diff_E2, 2
    )
    form_heatmap(
        growth_list_3, first_input.record_name, E1_list, E2_list, diff_E1, diff_E2, 3
    )
    form_heatmap(
        growth_list_4, first_input.record_name, E1_list, E2_list, diff_E1, diff_E2, 4
    )
    form_heatmap(
        growth_list_5, first_input.record_name, E1_list, E2_list, diff_E1, diff_E2, 5
    )
    form_heatmap(
        growth_list_6, first_input.record_name, E1_list, E2_list, diff_E1, diff_E2, 6
    )
    print("complete")


if __name__ == "__main__":
    run_multi()
