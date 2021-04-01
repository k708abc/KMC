import os
from InputParameter import Params
from Modules.cal_rates import rate
from Modules.choose_site import choose_atom
from Modules.deposition import *
from Modules.event_collection import site_events
from Test_modules.read_examples import *
from Test_modules.deposition_check import dep_check_poscar
from Test_modules.event_collection_check import *
from Modules.rejection_free_choose import rejection_free_choise


def check_cal_rates():
    print("Chek cal_rates.py")
    init_values = Params()
    pre = float(init_values.prefactor)
    kbt = init_values.temperature_eV
    energy = float(input("Input energy: "))
    print("prefactor: " + str(pre))
    print("kbt: " + str(kbt))
    print("energy: " + str(energy))
    print("rate = " + str(rate(pre, kbt, energy)))


def check_choose_site():
    print("Chek choose_site.py")
    num_cand = 10
    repetition = 1000
    candidates = [(0, 0, i) for i in range(num_cand)]
    results = [0 for _ in range(num_cand)]
    for _ in range(repetition):
        result = choose_atom(candidates)
        results[result[2]] += 1
    print("repetition : " + str(repetition))
    print("candidates : " + str(candidates))
    print("results : " + str(results))


def check_deposition():
    print("Chek deposition.py")
    init_values = Params()
    defect = init_values.keep_defect_check
    empty_first = init_values.num_defect
    unit_length = init_values.n_cell_init
    #
    atom_set = read_atom_set()
    bonds = read_bonds()
    lattice = read_lattice()
    maxz = highest_z(atom_set)

    candidate = find_candidates(atom_set, bonds)
    if (defect is True) and (empty_first == 1):
        candidate = remove_first(candidate)

    selected = [0 for _ in range(len(candidate))]
    for i in range(1000):
        deposited_pos = dep_position(candidate)
        index = candidate.index(deposited_pos)
        selected[index] += 1
    dep_check_poscar(atom_set, candidate, lattice, unit_length, maxz)
    print("candidates :" + str(candidate))
    print("deposited position:" + str(selected))
    print("A poscar is formed to check the candidate of deposition")


def check_event_collection():
    lattice = read_lattice()
    atom_set = read_atom_set()
    bonds = read_bonds()
    parameters = Params()
    unit_length = parameters.n_cell_init
    maxz = highest_z(atom_set)
    defect = parameters.keep_defect_check
    empty_first = parameters.num_defect
    #
    energy_bonding = [0.1 for i in range(0, 30)]
    energy_diffuse = [1 for i in range(0, 30)]

    dir_name = "./Test_modules/Event_check/"
    os.makedirs(dir_name, exist_ok=True)
    target_cand = existing_atoms(atom_set)

    if os.path.exists("Test_modules/Event_check/Check_rates.txt"):
        os.remove("Test_modules/Event_check/Check_rates.txt")

    for target in target_cand:
        event_list, rate_list = site_events(
            atom_set, bonds, target, parameters, energy_bonding, energy_diffuse
        )
        #
        event_check_poscar(atom_set, event_list, lattice, unit_length, maxz, target)
        rate_check(atom_set, bonds, rate_list, target)
    print("Poscars for event check are formed.")


"""
def check_rf_choise():
    event_time = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9, 1.0]]
    event_time_tot = []
    total_event_time = 0.05
    selected_num = []
    for i in event_time:
        site_tot_time = 0
        selected_num.append([])
        for k in i:
            selected_num[-1].append(0)
            total_event_time += k
            site_tot_time += k
        event_time_tot.append(site_tot_time)

    site, num = rejection_free_choise(total_event_time, event_time, event_time_tot)
"""

if __name__ == "__main__":
    cal_rates = False
    choose_site = False
    deposition = False
    event_collection = False
    rf_choise = True
    # cal_rates

    if cal_rates:
        check_cal_rates()
    elif choose_site:
        check_choose_site()
    elif deposition:
        check_deposition()
    elif event_collection:
        check_event_collection()
    elif rf_choise:
        pass
        # check_rf_choise()
