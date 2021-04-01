from typing import List, Tuple, Dict
import random
import math


# 20210401 checked


def find_candidates(atom_set: Dict, bonds: Dict) -> List[Tuple[int, int, int]]:
    candidate: List[Tuple[int, int, int]] = []
    for atom_index, condition in atom_set.items():
        if condition == 1:
            pass
        elif atom_index[2] == 0:
            candidate.append(atom_index)
        else:
            for bond in bonds[atom_index]:
                if atom_set[bond] in (1, 2):
                    candidate.append(atom_index)
                    break
    return list(set(candidate))


def dep_position(candidate: List) -> Tuple[int, int, int]:
    return candidate[math.floor(random.random() * len(candidate))]


def remove_first(candidate) -> List[Tuple]:
    return [site for site in candidate if site[2] not in (0, 1)]


def deposit_an_atom(atom_set: Dict, bonds: Dict, params, empty_first: int) -> Tuple:
    candidate = find_candidates(atom_set, bonds)
    if (params.keep_defect_check is True) and (empty_first == int(params.num_defect)):
        candidate = remove_first(candidate)
    return dep_position(candidate)


def highest_z(atom_set: Dict):
    maxz = 1
    for index, state in atom_set.items():
        if (state != 0) and (index[2] + 1 > maxz):
            maxz = index[2] + 1
    return maxz


"""
if __name__ == "__main__":
    defect = True
    empty_first = 2
    parameter = Params()
    unit_length = parameter.n_cell_init
    #
    atom_set = read_atom_set()
    bonds = read_bonds()
    lattice = read_lattice()
    maxz = highest_z(atom_set)
    candidate = find_candidates(atom_set, bonds)
    if (defect is True) and (empty_first == 1):
        candidate = remove_first(candidate)
    dep_check_poscar(atom_set, candidate, lattice, unit_length, maxz)
    print("Test candidates for depositing an atom")
    print("candidates :" + str(candidate))
    print("A poscar is formed to check the candidate of deposition")
"""