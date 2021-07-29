from typing import List, Dict, Tuple
import sys
import os

sys.path.append("../../KMC/")
from Modules.InputParameter import Params
from Modules.event_collection import site_events
from Modules.find_candidates import find_candidates
from Test_conditions.read_examples import (
    read_atom_set,
    read_bonds,
    read_lattice,
    read_candidate,
    lattice_size,
)



def existing_atoms(atom_set: Dict):
    candidate = []
    for index, state in atom_set.items():
        if state != 0:
            candidate.append(index)
    return candidate


def candidate_check_poscar(lattice, unit_length, maxz, target, candidates):
    #

    candidate_layers = []
    min_layer = target[2] - 3
    max_layer = target[2] + 3
    lattice_layers = [[] for i in range(7)]

    for index, coodinate in lattice.items():
        if index[2] >= min_layer and index[2] <= max_layer:
            if index not in candidates and index != target:
                lattice_layers[index[2] - min_layer].append(coodinate)
    #
    candidate_layer_num = []
    for index in candidates:
        if index[2] in candidate_layer_num:
            pass
        else:
            candidate_layer_num.append(index[2])
            candidate_layers.append([])
    min_val = min(candidate_layer_num)
    for index in candidates:
        if index != target:
            candidate_layers[index[2] - min_val].append(lattice[index])

    #
    file_name = "Event_check_2/Check_candidates_" + str(target) + ".vasp"
    file_data = open(file_name, "w")
    file_data.write("check events" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    atom_list = ["Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Ar", "Si", "P", "S", "Cr", "Ar", "K", "Ka"]
    atom_list2 = ["Ar", "Kr", "Bi", "Sn", "Fe", "Ni"]
    for i in range(len(lattice_layers)):
        file_data.write(atom_list[i] + "\t")
    for i in range(len(candidate_layers)):
        file_data.write(atom_list2[i] + "\t")
    file_data.write("Ti" + "\n")
    #
    for i in lattice_layers:
        file_data.write(str(len(i)) + "\t")
    for i in candidate_layers:
        file_data.write(str(len(i)) + "\t")
    file_data.write("1" + "\n")
    #
    file_data.write("direct" + "\n")
    for i in range(len(lattice_layers)):
        for k in lattice_layers[i]:
            file_data.write(
                str(k[0] / unit_length)
                + "\t"
                + str(k[1] / unit_length)
                + "\t"
                + str(k[2] / maxz / 2.448)
                + "\n"
            )
    for i in range(len(candidate_layers)):
        for k in candidate_layers[i]:
            file_data.write(
                str(k[0] / unit_length)
                + "\t"
                + str(k[1] / unit_length)
                + "\t"
                + str(k[2] / maxz / 2.448)
                + "\n"
            )
    file_data.write(
        str(lattice[target][0] / unit_length)
        + "\t"
        + str(lattice[target][1] / unit_length)
        + "\t"
        + str(lattice[target][2] / maxz / 2.448)
        + "\n"
    )
    file_data.close()


def highest_z(atom_set):
    maxz = 1
    for index, state in atom_set.items():
        if (state != 0) and (index[2] + 1 > maxz):
            maxz = index[2] + 1
    return maxz


def trans_formation(atom_set, unit_xy, trans_num):
    heighest_atom: Dict = {
        (i, k): 0 for i in range(0, unit_xy) for k in range(0, unit_xy)
    }
    for index, val in atom_set.items():
        if index[2] >= trans_num and val == 1:
            heighest_atom[(index[0], index[1])] += 1

    return heighest_atom


if __name__ == "__main__":
    lattice = read_lattice()
    atom_set = read_atom_set()
    bonds = read_bonds()
    candidates = read_candidate()
    parameters = Params("../kmc_input.yml")
    unit_length = lattice_size(lattice)
    maxz = highest_z(atom_set)
    defect = parameters.keep_defect_check
    empty_first = parameters.keep_defect_num
    trans_num = parameters.trans_num
    highest_vals = trans_formation(atom_set, unit_length, trans_num)
    #
    energy_bonding = [0.1 for i in range(0, 30)]
    energy_diffuse = [1 for i in range(0, 30)]
    #
    if os.path.exists("Event_check_2") is False:
        os.mkdir("Event_check_2")

    target = (5, 5, 5)
    #
    candidates = find_candidates(bonds, target, unit_length, maxz)
    candidate_check_poscar(lattice, unit_length, maxz, target, candidates)

    print("Poscars for event check are formed.")
