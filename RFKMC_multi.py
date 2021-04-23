from concurrent import futures
from Rejection_free_kmc import rejection_free
from InputParameter import Params
from Modules.heatmap import form_heatmap


def function(first_input, val1, val2):
    # print("starting:" + str(val1) + "/" + str(total_cal))
    print("val 1 = " + str(val1) + ", " + "val 2 = " + str(val2))
    rec_name = first_input.record_name + "_" + str(val1) + "_" + str(val2)
    rf_class = rejection_free(2)

    if first_input.repeat_val == 0:
        rf_class.init_value.binding_energies["Si_second"] = val1
        rf_class.init_value.binding_energies["Si_third"] = val2
        rf_class.init_value.binding_energies["Si_else"] = val2

    elif first_input.repeat_val == 1:
        rf_class.init_value.diffusion_barriers["Si_second"] = val1
        rf_class.init_value.diffusion_barriers["Si_third"] = val2
        rf_class.init_value.diffusion_barriers["Si_else"] = val2

    rf_class.init_value.record_name = rec_name
    rf_class.start()
    print("(" + str(val1) + "," + str(val2) + ")" + " finished")
    # growth_rec[(val1, val2)] = rf_class.mode_val


def run_multi():
    first_input = Params()
    start_E1 = first_input.start_E1
    end_E1 = first_input.end_E1
    diff_E1 = first_input.diff_E1
    start_E2 = first_input.start_E2
    end_E2 = first_input.end_E2
    diff_E2 = first_input.diff_E2
    # or_rec_name = first_input.record_name

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
    energy_list = [(val1, val2) for val1 in E1_list for val2 in E2_list]

    future_list = []
    with futures.ProcessPoolExecutor(max_workers=2) as executor:
        for var in energy_list:
            future = executor.submit(function, first_input, var[0], var[1])
            future_list.append(future)
        for i in futures.as_completed(future_list):
            print(i.result())
    print("complete")

    """

    growth_rec = []
    total_cal = len(E1_list) * len(E2_list)
    num_cal = 0
    for i in E1_list:
        growth_rec.append([])
        for k in E2_list:
            num_cal += 1
            print("starting:" + str(num_cal) + "/" + str(total_cal))
            print("val 1 = " + str(i))
            print("val 2 = " + str(k))
            rec_name = or_rec_name + "_" + str(i) + "_" + str(k)
            rf_class = rejection_free(1)

            if first_input.repeat_val == 0:
                rf_class.init_value.binding_energies["Si_second"] = i
                rf_class.init_value.binding_energies["Si_third"] = k
                rf_class.init_value.binding_energies["Si_else"] = k
            elif first_input.repeat_val == 1:
                rf_class.init_value.diffusion_barriers["Si_second"] = i
                rf_class.init_value.diffusion_barriers["Si_third"] = k
                rf_class.init_value.diffusion_barriers["Si_else"] = k

            rf_class.init_value.record_name = rec_name
            rf_class.start()
            growth_rec[-1].append(rf_class.mode_val)
    """


if __name__ == "__main__":
    run_multi()