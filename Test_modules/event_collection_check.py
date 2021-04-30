from typing import List, Dict, Tuple
import sys
import os

sys.path.append("../../KMC/")
from Modules.InputParameter import Params
from Modules.event_collection import site_events
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


def rate_check(atom_set, bonds, rate_list, target):
    bond = []
    for i in bonds[target]:
        if atom_set[i] == 1:
            bond.append(1)
        else:
            bond.append(0)

    file_name = "Event_check/Check_rates.txt"
    file_data = open(file_name, "a")
    if len(rate_list) > 0:
        file_data.write(
            str(target)
            + "\t"
            + "rate: "
            + str(rate_list[0])
            + "\t"
            + "bonds: "
            + str(bond)
            + "\n"
        )
    else:
        file_data.write(
            str(target) + "\t" + "rate: 0" + "\t" + "bonds: " + str(bond) + "\n"
        )
    file_data.close()


def event_check_poscar(
    atom_set: Dict[Tuple[int, int, int], int],
    events: List,
    lattice: Dict,
    unit_length: int,
    maxz: int,
    target: Tuple[int, int, int],
):
    xp: List[float] = [lattice[target][0] / unit_length]
    yp: List[float] = [lattice[target][1] / unit_length]
    zp: List[float] = [lattice[target][2] / maxz / 2.448]
    atom_i = 0
    eve_i = 0
    for index, atom_state in atom_set.items():
        if atom_state != 0:
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_i += 1
    for cand in events:
        xp.append(lattice[cand][0] / unit_length)
        yp.append(lattice[cand][1] / unit_length)
        zp.append(lattice[cand][2] / maxz / 2.448)
        eve_i += 1
    #
    file_name = "Event_check/Check_events_" + str(target) + ".vasp"
    file_data = open(file_name, "w")
    file_data.write("check events" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Bi" + "\t" + "Si" + "\t" + "O" + "\n")
    file_data.write(str(1) + "\t" + str(atom_i) + "\t" + str(eve_i) + "\n")
    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
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
    if os.path.exists("Event_check") is False:
        os.mkdir("Event_check")

    target_cand = existing_atoms(atom_set)
    #
    for target in target_cand:
        event_list, rate_list = site_events(
            atom_set,
            bonds,
            target,
            parameters,
            energy_bonding,
            energy_diffuse,
            candidates,
            highest_vals,
        )
        #
        event_check_poscar(atom_set, event_list, lattice, unit_length, maxz, target)
        rate_check(atom_set, bonds, rate_list, target)
    print("Poscars for event check are formed.")