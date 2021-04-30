from typing import List, Dict, Tuple
import sys
import os

sys.path.append("../../KMC/")
from Modules.InputParameter import Params
from Modules.deposition import find_candidates, remove_first, dep_position
from Test_conditions.read_examples import (
    read_atom_set,
    read_bonds,
    read_lattice,
    lattice_size,
)


def dep_check_poscar(
    atom_set: Dict[Tuple[int, int, int], int],
    candidate: List[Tuple[int, int, int]],
    lattice: Dict,
    unit_length: int,
    maxz: int,
):
    xp: List[float] = []
    yp: List[float] = []
    zp: List[float] = []
    atom_i = 0
    cand_i = 0
    for index, atom_state in atom_set.items():
        if atom_state != 0:
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_i += 1
    for cand in candidate:
        xp.append(lattice[cand][0] / unit_length)
        yp.append(lattice[cand][1] / unit_length)
        zp.append(lattice[cand][2] / maxz / 2.448)
        cand_i += 1
    #
    file_name = "Deposition_check\poscar_for_check_depsotion.vasp"
    file_data = open(file_name, "w")
    file_data.write("check depsotion" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Si" + "\t" + "O" + "\n")
    file_data.write(str(atom_i) + "\t" + str(cand_i) + "\n")
    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()


def highest_z(atom_set: Dict):
    maxz = 1
    for index, state in atom_set.items():
        if (state != 0) and (index[2] + 1 > maxz):
            maxz = index[2] + 1
    return maxz


if __name__ == "__main__":
    print("Chek deposition.py")
    if os.path.exists("Deposition_check") is False:
        os.mkdir("Deposition_check")
    #
    init_values = Params("../kmc_input.yml")
    defect = init_values.keep_defect_check
    empty_first = init_values.keep_defect_num
    #
    atom_set = read_atom_set()
    bonds = read_bonds()
    lattice = read_lattice()
    maxz = highest_z(atom_set)
    unit_length = lattice_size(lattice)

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